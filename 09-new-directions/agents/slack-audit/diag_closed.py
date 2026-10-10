"""Closed-form reduction of the DIAG optimum: b = 0, wall M + ell = 1, binding high bin at x = 0.
theta(ell) = max(1 - h/6, max_delta high0(delta, ell)),  h = (1+3 ell)/2,  l_x = l_y = (1-ell)/2,
high0 = (1+delta)(1+ell)/4 - ell/2 + h (R0(delta) + delta/2 - 1/6),  R0 = 1 - 2 delta alpha/(3 alpha - delta)
(R0 = exact P78 count at x = 0, long-inverse closed form). alpha = amplification slope (5/6 paper)."""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
def high0(ell, alpha=5/6, d0=1e-3):
    h = (1 + 3 * ell) / 2
    ds = np.linspace(d0, 0.75, 20001)
    R0 = 1 - 2 * ds * alpha / (3 * alpha - ds)
    v = (1 + ds) * (1 + ell) / 4 - ell / 2 + h * (R0 + ds / 2 - 1 / 6)
    k = v.argmax(); return v[k], ds[k]
def theta(alpha=5/6):
    f = lambda ell: (1 - (1 + 3 * ell) / 12) - high0(ell, alpha)[0]
    ell = brentq(f, 0.05, 0.4)
    return 1 - (1 + 3 * ell) / 12, ell, high0(ell, alpha)[1]
if __name__ == "__main__":
    th, ell, dl = theta()
    print(f"DIAG closed form: theta = {th:.7f}, ell = {ell:.6f}, h = {(1+3*ell)/2:.6f}, binding delta = {dl:.4f}")
    # alpha needed for 27/31 on this line
    a = brentq(lambda al: theta(al)[0] - 27/31, 0.5, 5/6)
    print(f"alpha needed for 27/31 (x=0 long-count slope; = (5/6) kA): {a:.5f}  -> kA = {a/(5/6):.4f}")
    # needed h for 27/31: h = 6(1 - 27/31) = 24/31
    ell_t = (2 * 24 / 31 - 1) / 3
    v, dl = high0(ell_t)
    print(f"target ell = {ell_t:.5f}, high0 at paper counts = {v:.6f} (needs <= {27/31:.6f}), at delta = {dl:.3f}")
