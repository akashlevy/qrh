"""Plain-cube beyond-edge slope lam6: sum_u |S_m|^6 << U^{1 + lam6 (3m-1)} for 3m > 1 (count 1 + lam6(3m-1) - 3 d m).
Brute force with the detector's witness freedom (r + m >= t).  Claim: R(delta,0) = R0(delta; min(5/6, lam6))."""
import numpy as np, counts_lab as CL
def R_pc(d, lam6, nt=101, nr=801, nm=401, alpha=5/6):
    best = np.inf
    for t in np.linspace(1, 1.5, nt):
        r = np.linspace(max(t - .5, 0), t, nr)[:, None]; m = np.linspace(0, .5, nm)[None, :]
        ok = r + m >= t - 1e-12
        fA = np.maximum(1, 1 + alpha * (r - 1)) - d * r
        fP = 1 - 2 * d * m
        f6 = np.where(3 * m > 1, 1 + lam6 * (3 * m - 1) - 3 * d * m, 1 - 3 * d * m)
        F = np.where(ok, np.minimum(np.minimum(fA, fP), f6), -np.inf)
        best = min(best, F.max())
    return best
for d in (0.27, 0.33, 0.39):
    for lam6 in (1.0, 5/6, 0.80, 0.70):
        print(f"delta={d} lam6={lam6:.3f}: R={R_pc(d, lam6):.5f}  R0(min)={CL.R0(d, min(5/6, lam6)):.5f}")
