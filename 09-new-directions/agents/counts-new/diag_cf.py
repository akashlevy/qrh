"""Closed form of the DIAG optimum (slack-audit diag_closed.py) for an ARBITRARY x = 0 count R(delta):
theta(R) = min_ell max(1 - (1+3ell)/12, max_delta (1+delta)(1+ell)/4 - ell/2 + ((1+3ell)/2)(R(delta) + delta/2 - 1/6)).
Valid while the DIAG geometry (b = 0, M + ell = 1, binding bin at x = 0) persists (checked against the LP for kA).
Count families (all at x = 0; alpha = 5/6 copies slope):
  R0(d)                 paper / long-inverse closed form
  R_kappa(d, k)         M_kappa: long Moebius mean square diagonal for H >= D^k (with copies): e(r) = max(1,(1+5kr)/6)
  R_r1(d, r1)           diagonal sum_{u<=U}|M_r|^2 << U up to r <= r1, then copies slope 5/6
  R_joint(d, w)         joint sum_u |M_r S_m|^2 << U for r+m <= 1+w, MVT slope 1 beyond (no copies)
  R_joint_slope(d, l)   joint with copies-type slope l beyond s = 1
  R_cube(d, w)          plain sixth moment sum_u |S_m|^6 << U for 3m <= 1+w (4th moment beyond)
  R_p6full(d)           plain sixth moment for all m <= 1/2
"""
import numpy as np
from scipy.optimize import brentq

A = 5 / 6
def R0(d, a=A): return 1 - 2 * d * a / (3 * a - d)
def tstar(d, a=A): return a / (a - d / 3)
def R_kappa(d, k): return np.minimum(R0(d), 1 - 10 * d / (3 * (5 * k - 2 * d)))
def R_r1(d, r1): return np.minimum(R0(d), 1 - 2 * d * r1 * tstar(d) / 3)
def R_joint(d, w): return np.minimum(R0(d), 1 - 2 * d * (1 + w) / (3 - d))
def R_joint_slope(d, l): return R0(d, min(A, l))
def R_cube(d, w): return np.minimum(R0(d), np.maximum(1 - 2 * d * (1 + w) / 3, 1 - 3 * d * A / (4 * A - d)))
def R_p6full(d): return 1 - 3 * d * A / (4 * A - d)

DS = np.linspace(1e-3, 0.75, 7501)
def high0(ell, Rfun):
    h = (1 + 3 * ell) / 2
    v = (1 + DS) * (1 + ell) / 4 - ell / 2 + h * (Rfun(DS) + DS / 2 - 1 / 6)
    k = v.argmax(); return v[k], DS[k]
def theta(Rfun):
    f = lambda ell: (1 - (1 + 3 * ell) / 12) - high0(ell, Rfun)[0]
    ell = brentq(f, 0.05, 0.4)
    return 1 - (1 + 3 * ell) / 12, ell, high0(ell, Rfun)[1]
def need(fam, lo, hi):
    """parameter value of a one-parameter family reaching 27/31"""
    return brentq(lambda p: theta(lambda d: fam(d, p))[0] - 27 / 31, lo, hi)

if __name__ == "__main__":
    print("DIAG (R0):", theta(R0))
    print("check M_kappa k=0.9627 (LP: 0.870928):", theta(lambda d: R_kappa(d, 0.9627)))
    print("needed for 27/31 with DIAG:")
    print("  M_kappa  kappa =", round(need(R_kappa, 0.85, 0.999), 4))
    print("  diag Moebius up to r1 =", round(need(R_r1, 1.0001, 1.2), 4))
    print("  joint M.S diagonal up to length U^(1+w), w =", round(need(R_joint, 0.0, 0.3), 4))
    print("  joint M.S slope beyond s=1 (copies-type), lambda =", round(need(R_joint_slope, 0.3, 0.8333), 4))
    print("  plain sixth moment diagonal up to 3m <= 1+w, w =", round(need(R_cube, 0.0, 0.5), 4))
    print("full beyond-edge plain sixth moment (m <= 1/2):", theta(R_p6full))
    print("DH counts:", theta(lambda d: 1 - d))
