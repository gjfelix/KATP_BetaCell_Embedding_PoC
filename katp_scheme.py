"""
Reduced allosteric KATP scheme (concerted MWC form with four ATP sites; Enkvetchakul and Nichols, 2003):

    P_O = 1 / (1 + L_eff * [(1 + A/K_C) / (1 + A/K_O)]^4),   L_eff = L / (1 + D/K_D),   K_O > K_C

and its calibration (patch IC50 and the conductance of the original cell model).

This only illustrates the embedding. L and the patch IC50 come from electrophysiology (Proks et al., 2004),
K_C and K_O are fitted, and K_D, the adenine pool and the size of the perturbations are assumptions.
Nothing here is derived from structures or simulations.

Units: concentrations in uM (A, D, K_C, K_O, K_D, Tot), conductances in pS.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def po_scheme(A, D, L, KC, KO, KD):
    Leff = L / (1.0 + D / KD)
    return 1.0 / (1.0 + Leff * ((1.0 + A / KC) / (1.0 + A / KO)) ** 4)


@njit(cache=True)
def po_of_a(a, L, KC, KO, KD, Tot):
    A = Tot / (1.0 + a)
    D = Tot * a / (1.0 + a)
    return po_scheme(A, D, L, KC, KO, KD)


# ----------------------------------------------------------------------------------------
# Calibration helpers (plain Python, scipy)
# ----------------------------------------------------------------------------------------
def patch_ic50(L, KC, KO, KD=1e12):
    """ATP concentration (uM) giving half of Po,zero in an excised patch (ADP = 0)."""
    from scipy.optimize import brentq
    p0 = po_scheme(0.0, 0.0, L, KC, KO, KD)
    return brentq(lambda A: po_scheme(A, 0.0, L, KC, KO, KD) - p0 / 2, 1e-4, 1e7)


def calibrate(N, gamma, KD=100.0, L=0.89, Tot=5000.0, ic50_target=7.0, a_lo=0.40, a_hi=0.60, free_KD=False, gbar=500.0):
    """Fit KC and KO (log space) so that
       (i) the patch IC50 (ADP = 0) is ic50_target uM (Proks 2004, wild type), and
       (ii) N * gamma * Po(a) = gbar * a (reference line of the original model; 297 pS in make_pubfigs_compute.py) for a in [a_lo, a_hi] at the nucleotide
            levels ATP = Tot/(1+a), ADP = Tot*a/(1+a).
    Returns dict(KC, KO, KD, L, Tot, max_rel_err, ic50)."""
    from scipy.optimize import least_squares
    agrid = np.linspace(a_lo, a_hi, 21)
    target = gbar * agrid / (N * gamma)

    def resid(x):
        KC, KO = np.exp(x[:2])
        kd = np.exp(x[2]) if free_KD else KD
        po = np.array([po_of_a(a, L, KC, KO, kd, Tot) for a in agrid])
        return np.concatenate([np.log(po) - np.log(target), [2.0 * np.log(patch_ic50(L, KC, KO) / ic50_target)]])

    best = None
    for KC0 in (5.0, 20.0, 80.0):
        for KO0 in (50.0, 300.0, 3000.0):
            for KD0 in ((0.1, 10.0, 1000.0) if free_KD else (None,)):
                x0 = [np.log(KC0), np.log(KO0)] + ([np.log(KD0)] if free_KD else [])
                try:
                    s_ = least_squares(resid, x0)
                except Exception:
                    continue
                if s_.x[1] > s_.x[0] and (best is None or s_.cost < best.cost):
                    best = s_
    KC, KO = np.exp(best.x[:2])
    if free_KD:
        KD = float(np.exp(best.x[2]))
    po = np.array([po_of_a(a, L, KC, KO, KD, Tot) for a in agrid])
    return dict(KC=KC, KO=KO, KD=KD, L=L, Tot=Tot, max_rel_err=float(np.max(np.abs(po / target - 1))),
                ic50=patch_ic50(L, KC, KO))


def elasticity(c, a=0.5, h=1e-4):
    """d ln(Po) / d ln(a) of the reduced scheme at nucleotide ratio a (the original model has exactly 1)."""
    f = lambda x: np.log(po_of_a(x, c['L'], c['KC'], c['KO'], c['KD'], c['Tot']))
    return (f(a + h) - f(a - h)) / (np.log(a + h) - np.log(a - h))
