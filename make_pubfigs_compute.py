"""Runs all the simulations behind Figure 5 and the supplementary figures and writes results_pub.json.
Variant A (K_D = 100 uM, assumed) is the case used in the paper. Variant B, in which K_D is fitted as well, is only
used to show what that would require (K_D = 3.8 uM). Run: python make_pubfigs_compute.py (takes about 10 minutes).
"""
import json, time, sys
import numpy as np
import katp_scheme as kp
import human_poc as hp
import analysis_h as ah

T0 = time.time()
N = 10000
GBAR = hp.GK_REF * hp.CM * 1000.0 / hp.A_REF       # 297 pS per unit a
A_LO, A_HI = 0.40, 0.90
SKIP = 10.0
_cal = {}


def cal(gamma, variant='A'):
    key = (gamma, variant)
    if key not in _cal:
        _cal[key] = kp.calibrate(N, gamma, KD=100.0, free_KD=(variant == 'B'), gbar=GBAR, a_lo=A_LO, a_hi=A_HI)
    return _cal[key]


def sim(mode, a=0.5, gamma=70.0, tau_ch=5.0, avail=1.0, seed=1, scheme_over=None, variant='A', sigma=0.0, T_s=60.0, gk_fixed=-1.0):
    c = cal(gamma, variant)
    scheme = dict(L=c['L'], KC=c['KC'], KO=c['KO'], KD=c['KD'], Tot=c['Tot'])
    scheme.update(scheme_over or {})
    return hp.simulate_h(mode, T_s=T_s, seed=seed, N=N, gamma=gamma, tau_ch=tau_ch, avail=avail, a_in=a, sigma=sigma,
                         scheme=scheme, gk_fixed=gk_fixed)


def st(r):
    s = ah.spike_stats(r, t_skip=SKIP)
    return {k: (None if (isinstance(v, float) and np.isnan(v)) else v) for k, v in s.items()}


def npo_of(gamma, a, variant='A'):
    c = cal(gamma, variant)
    return N * kp.po_of_a(a, c['L'], c['KC'], c['KO'], c['KD'], c['Tot'])


def log(msg):
    print(f'[{time.time()-T0:6.0f}s] {msg}', flush=True)


def main():
    R = {}
    # ---- S1: implementation check, original model, rate versus g_KATP
    gk = [0.005, 0.0075, 0.010, 0.0125, 0.015, 0.0165, 0.018, 0.0185, 0.019, 0.0195, 0.020, 0.025, 0.030, 0.040]
    R['S1'] = [dict(g=g, **st(sim(0, gk_fixed=g, T_s=40.0))) for g in gk]
    log('S1 done')

    # ---- scheme parameters and patch titration (panel A, S2)
    cA, cB = cal(70.0, 'A'), cal(70.0, 'B')
    A = np.logspace(-1, 5, 300)
    f = 298.0 / 7.0
    R['scheme'] = dict(
        A=dict(cA, elasticity=kp.elasticity(cA)), B=dict(cB, elasticity=kp.elasticity(cB)),
        atp=A.tolist(),
        po_wt=[kp.po_scheme(x, 0.0, cA['L'], cA['KC'], cA['KO'], 1e12) for x in A],
        po_mut=[kp.po_scheme(x, 0.0, cA['L'], cA['KC'] * f, cA['KO'] * f, 1e12) for x in A],
        a=np.linspace(0.3, 1.0, 100).tolist(),
        target=(GBAR * np.linspace(0.3, 1.0, 100) / (N * 70.0)).tolist(),
        poA=[kp.po_of_a(x, cA['L'], cA['KC'], cA['KO'], cA['KD'], cA['Tot']) for x in np.linspace(0.3, 1.0, 100)],
        poB=[kp.po_of_a(x, cB['L'], cB['KC'], cB['KO'], cB['KD'], cB['Tot']) for x in np.linspace(0.3, 1.0, 100)],
        ic50_wt=kp.patch_ic50(cA['L'], cA['KC'], cA['KO']), ic50_mut=kp.patch_ic50(cA['L'], cA['KC'] * f, cA['KO'] * f))
    log('scheme done')

    # ---- deterministic sweep (panel B, S3)
    a_grid = sorted(set(np.round(np.concatenate([np.arange(0.30, 0.60, 0.05), np.arange(0.60, 0.701, 0.01), np.arange(0.75, 1.001, 0.05)]), 2)))
    R['sweep'] = dict(a=a_grid, ref=[], A=[], B=[])
    for a in a_grid:
        R['sweep']['ref'].append(st(sim(0, a=a, T_s=40.0)))
        R['sweep']['A'].append(st(sim(1, a=a, variant='A', T_s=40.0)))
        R['sweep']['B'].append(st(sim(1, a=a, variant='B', T_s=40.0)))
    log('sweep done')

    # ---- stochastic gating versus the number of open channels (panels C, D, E)
    GAM = [140.0, 100.0, 70.0, 50.0, 35.0, 25.0, 20.0, 14.0, 10.0, 7.0, 5.0, 3.5, 2.5, 2.0, 1.0, 0.5]
    R['noise'] = {}
    traces = {}
    for a in (0.5, 0.7):
        rows = []
        for g in GAM:
            ss = []
            for seed in (1, 2, 3):
                r = sim(2, a=a, gamma=g, seed=seed, T_s=120.0)
                ss.append(st(r))
                if seed == 1 and a == 0.7 and g in (70.0, 20.0, 7.0, 2.0):
                    m = (r['t'] >= 20) & (r['t'] < 40)
                    traces[str(g)] = dict(t=r['t'][m].tolist(), v=r['Vmax'][m].tolist(), npo=npo_of(g, a))
            rows.append(dict(gamma=g, npo=npo_of(g, a), runs=ss))
        R['noise'][str(a)] = rows
        log(f'noise a={a} done')
    R['traces_a07'] = traces
    R['det_ref'] = {str(a): st(sim(1, a=a, T_s=60.0)) for a in (0.5, 0.7)}

    # ---- channel kinetics (S5)
    R['tau'] = {}
    for a in (0.5, 0.7):
        rows = []
        for tau in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0):
            rows.append(dict(tau=tau, runs=[st(sim(2, a=a, gamma=70.0, tau_ch=tau, seed=s, T_s=120.0)) for s in (1, 2, 3)]))
        R['tau'][str(a)] = rows
    log('tau done')

    # ---- perturbations (panel F, S6)
    c = cal(70.0, 'A')
    conds = [('wild type', 0.5, {}), ('R201H-like', 0.5, dict(scheme_over=dict(KC=c['KC'] * 298 / 7, KO=c['KO'] * 298 / 7))),
             ('Q52R-like', 0.5, dict(scheme_over=dict(L=c['L'] / 4))), ('wild type', 0.7, {}),
             ('PIP2 depletion', 0.7, dict(scheme_over=dict(L=c['L'] * 4))), ('sulfonylurea', 0.7, dict(avail=0.1))]
    R['perturb'] = []
    for name, a, kw in conds:
        rd = sim(1, a=a, T_s=60.0, **kw)
        runs_s = []
        tr = None
        for seed in (1, 2, 3):
            rs = sim(2, a=a, gamma=70.0, seed=seed, T_s=120.0, **kw)
            runs_s.append(st(rs))
            if seed == 1:
                m = (rs['t'] >= 20) & (rs['t'] < 35)
                tr = dict(t=rs['t'][m].tolist(), v=rs['V'][m].tolist())
        m = (rd['t'] >= 20) & (rd['t'] < 35)
        R['perturb'].append(dict(name=name, a=a, det=st(rd), stoch=runs_s, trace_det=dict(t=rd['t'][m].tolist(), v=rd['V'][m].tolist()), trace_stoch=tr))
    log('perturbations done')

    # ---- the model's own white noise (S7)
    R['sigma'] = {}
    for a in (0.5, 0.7):
        R['sigma'][str(a)] = [dict(sigma=sg, **st(sim(0, a=a, sigma=sg, T_s=120.0))) for sg in (0.0, 0.1, 0.3, 1.0, 2.0)]
    log('sigma done')

    json.dump(R, open('results_pub.json', 'w'))
    log('saved results_pub.json')


if __name__ == '__main__':
    main()
