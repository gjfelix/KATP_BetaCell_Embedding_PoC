"""
The human beta-cell model of Pedersen (2010), re-implemented in Python from the XPPAUT file human_beta_model.ode,
with the option of replacing its constant KATP conductance by the channel scheme of katp_scheme.py.

Units are those of the .ode file: V in mV, t in ms, g in nS/pF, I in pA/pF. The membrane capacitance (9.9 pF)
converts between the conductance density of the file and a number of open channels:
N_open * gamma [pS] = g_KATP * Cm * 1000.
The original model has no Ca2+ or nucleotide dynamics, so the nucleotide ratio a = ADP/ATP is an input that stays
fixed in each run. To compare with the original model we assume g_KATP(a) = 0.015 nS/pF * (a / 0.5), which is the
value in the file at a = 0.5.

This only illustrates how a channel scheme can be embedded in a cell model. The channel parameters are fitted or
assumed (see the README); nothing here comes from protein structures or molecular simulations.
"""
import numpy as np
from numba import njit
from katp_scheme import po_of_a

CM = 9.9            # pF (the .ode: num cm=9.9; conductances are per pF)
GK_REF = 0.015      # nS/pF, the .ode value (gkatp) taken at a = 0.5
A_REF = 0.5

# parameters of human_beta_model.ode (order used by the kernel)
P = dict(vK=-75.0, vCa=65.0, vNa=70.0, gl=0.015, vleak=-30.0, gk=1.0, vnK=0.0, nnK=-10.0, taun0=2.0,
         gKCa=0.02, taum=2.0, VmK=0.0, nmK=-10.0, basal=20.0, gERG=0.2, VxERG=-30.0, nxERG=-10.0, tauxERG=100.0,
         tauyERG=50.0, VyERG=-42.0, nyERG=17.5, gNa=0.4, VhNa=-42.0, nhNa=6.0, tauh=2.0, VkNa=-18.0, nkNa=-5.0,
         gCaL=0.14, VxL=-25.0, nxL=-6.0, tauL=20.0, gCaPQ=0.17, VxPQ=-10.0, nxPQ=-10.0,
         gCaT=0.05, VxT=-40.0, nxT=-4.0, VyT=-64.0, nyT=8.0, tauT=7.0)
_KEYS = list(P.keys())


def pvec(**over):
    d = dict(P); d.update(over)
    return np.array([d[k] for k in _KEYS], dtype=np.float64)


_I = {k: i for i, k in enumerate(_KEYS)}


@njit(cache=True)
def _sig(x):
    return 1.0 / (1.0 + np.exp(x))


@njit(cache=True)
def rhs_h(v, n, m, h, yL, yT, xE, yE, gkatp, p):
    vK = p[0]; vCa = p[1]; vNa = p[2]; gl = p[3]; vleak = p[4]; gk = p[5]; vnK = p[6]; nnK = p[7]; taun0 = p[8]
    gKCa = p[9]; taum = p[10]; VmK = p[11]; nmK = p[12]; basal = p[13]; gERG = p[14]; VxERG = p[15]; nxERG = p[16]
    tauxERG = p[17]; tauyERG = p[18]; VyERG = p[19]; nyERG = p[20]; gNa = p[21]; VhNa = p[22]; nhNa = p[23]
    tauh = p[24]; VkNa = p[25]; nkNa = p[26]; gCaL = p[27]; VxL = p[28]; nxL = p[29]; tauL = p[30]
    gCaPQ = p[31]; VxPQ = p[32]; nxPQ = p[33]; gCaT = p[34]; VxT = p[35]; nxT = p[36]; VyT = p[37]; nyT = p[38]
    tauT = p[39]
    ninf = _sig((v - vnK) / nnK)
    taun = taun0 + 10.0 * np.exp(min(np.log(3.0), (-20.0 - v) / 6.0))
    minf = _sig((v - VmK) / nmK)
    xLinf = _sig((v - VxL) / nxL)
    yLinf = max(0.0, min(1.0, 1.0 + xLinf * (v - vCa) / 57.0))
    xPQinf = _sig((v - VxPQ) / nxPQ)
    xTinf = _sig((v - VxT) / nxT)
    yTinf = _sig((v - VyT) / nyT)
    hinf = _sig((v - VhNa) / nhNa)
    kinf = _sig((v - VkNa) / nkNa)
    xEinf = _sig((v - VxERG) / nxERG)
    yEinf = _sig((v - VyERG) / nyERG)
    ICaL = gCaL * xLinf * yL * (v - vCa)
    ICaPQ = gCaPQ * xPQinf * (v - vCa)
    ICaT = gCaT * xTinf * yT * (v - vCa)
    ICa = ICaL + ICaPQ + ICaT
    IKv = gk * n * (v - vK)
    IKCa = gKCa * (-ICa + basal) * m * (v - vK)
    IERG = gERG * xE * yE * (v - vK)
    INa = gNa * h * kinf * (v - vNa)
    IKATP = gkatp * (v - vK)
    Ileak = gl * (v - vleak)
    dv = -(IKv + IKCa + IERG + INa + ICaL + ICaPQ + ICaT + IKATP + Ileak)
    dn = (ninf - n) / taun
    dm = (minf - m) / taum
    dh = (hinf - h) / tauh
    dyL = (yLinf - yL) / tauL
    dyT = (yTinf - yT) / tauT
    dxE = (xEinf - xE) / tauxERG
    dyE = (yEinf - yE) / tauyERG
    return dv, dn, dm, dh, dyL, dyT, dxE, dyE


@njit(cache=True)
def run_h(mode, T_ms, dt, stride, seed, p, N, gamma, tau_ch, L, KC, KO, KD, Tot, avail, a_in, sigma,
          gk_fixed, maxspk):
    """
    mode 0: original model, g_KATP = GK_REF * a_in / A_REF (or gk_fixed if gk_fixed >= 0)
    mode 1: deterministic embedding (mean-field open fraction P of the two-state channel)
    mode 2: stochastic embedding (binomial update of the open-channel number)
    sigma:  white-noise amplitude of the .ode (v' += sigma*q); 0 in all KATP runs
    Heun step for the 8 voltage-dependent variables; Ito noise term added after the step.
    Returns out (5, nout): t (s), V, max V in the interval, spikes in the interval, mean open channels;
            spike_times (ms) and their number.
    """
    np.random.seed(seed)
    nsteps = int(T_ms / dt)
    nout = nsteps // stride + 1
    out = np.zeros((5, nout))
    spk = np.zeros(maxspk)
    nspk = 0
    v = -49.0; n = 0.0; m = 0.0; h = 1.0; yL = 1.0; yT = 1.0; xE = 0.0; yE = 0.0
    Nav = int(round(N * avail))
    Pinf = po_of_a(a_in, L, KC, KO, KD, Tot)
    Pc = Pinf
    k = int(round(Pinf * Nav))
    kc = 1.0 / tau_ch
    ko = kc * Pinf / (1.0 - Pinf)
    armed = True
    sumopen = 0.0; cnt = 0; vmax = -1e9; spikes_int = 0; j = 0
    sq = np.sqrt(dt)
    for i in range(nsteps + 1):
        if mode == 0:
            if gk_fixed >= 0.0:
                g = gk_fixed
            else:
                g = GK_REF * a_in / A_REF
            nopen = g * CM * 1000.0 / gamma
        elif mode == 1:
            nopen = Nav * Pc
            g = nopen * gamma / (CM * 1000.0)
        else:
            nopen = float(k)
            g = nopen * gamma / (CM * 1000.0)
        sumopen += nopen; cnt += 1
        if v > vmax:
            vmax = v
        if armed and v > -20.0:
            armed = False
            spikes_int += 1
            if nspk < maxspk:
                spk[nspk] = i * dt
                nspk += 1
        elif (not armed) and v < -30.0:
            armed = True
        if i % stride == 0:
            out[0, j] = i * dt * 1e-3; out[1, j] = v; out[2, j] = vmax; out[3, j] = spikes_int
            out[4, j] = sumopen / cnt
            sumopen = 0.0; cnt = 0; vmax = -1e9; spikes_int = 0; j += 1
        if i == nsteps:
            break
        d1 = rhs_h(v, n, m, h, yL, yT, xE, yE, g, p)
        v2 = v + dt * d1[0]; n2 = n + dt * d1[1]; m2 = m + dt * d1[2]; h2 = h + dt * d1[3]
        yL2 = yL + dt * d1[4]; yT2 = yT + dt * d1[5]; xE2 = xE + dt * d1[6]; yE2 = yE + dt * d1[7]
        d2 = rhs_h(v2, n2, m2, h2, yL2, yT2, xE2, yE2, g, p)
        v += 0.5 * dt * (d1[0] + d2[0]); n += 0.5 * dt * (d1[1] + d2[1]); m += 0.5 * dt * (d1[2] + d2[2])
        h += 0.5 * dt * (d1[3] + d2[3]); yL += 0.5 * dt * (d1[4] + d2[4]); yT += 0.5 * dt * (d1[5] + d2[5])
        xE += 0.5 * dt * (d1[6] + d2[6]); yE += 0.5 * dt * (d1[7] + d2[7])
        if sigma > 0.0:
            v += sigma * sq * np.random.randn()
        if mode == 1:
            Pc += dt * (ko * (1.0 - Pc) - kc * Pc)
        elif mode == 2:
            kclose = np.random.binomial(k, min(1.0, kc * dt))
            kopen = np.random.binomial(Nav - k, min(1.0, ko * dt))
            k = k - kclose + kopen
    return out, spk[:nspk]


def simulate_h(mode, T_s=60.0, dt=0.02, stride=250, seed=1, N=1000, gamma=70.0, tau_ch=5.0, avail=1.0, a_in=0.5,
               sigma=0.0, scheme=None, gk_fixed=-1.0, cell=None, maxspk=400000):
    sc = dict(L=0.89, KC=18.1, KO=163.0, KD=100.0, Tot=5000.0)
    sc.update(scheme or {})
    out, spk = run_h(mode, T_s * 1e3, dt, stride, seed, pvec(**(cell or {})), N, gamma, tau_ch, sc['L'], sc['KC'],
                     sc['KO'], sc['KD'], sc['Tot'], avail, a_in, sigma, gk_fixed, maxspk)
    return dict(t=out[0], V=out[1], Vmax=out[2], nspk=out[3], nopen=out[4], spikes_ms=spk)
