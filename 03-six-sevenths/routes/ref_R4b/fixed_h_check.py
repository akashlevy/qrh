# Referee: best theta when the row range is FORCED to h = 6/7, for W(w)+refined count with the paper's counts.
import sys, numpy as np
sys.path.insert(0, '/private/tmp/claude-501/-Users-akashlevy-Documents-riemann/eab6979c-f095-4107-a475-7e51d88a9fe7/scratchpad/wf67/r4')
from model import R_of, DEFAULT
from scipy.optimize import minimize
def theta_geom(lx, ly, ell, w, caseA="refined", cfg=None):
    cfg = dict(DEFAULT) if cfg is None else cfg
    h = 1 - lx + ell; b = ly - lx; M = lx + ly
    if not (lx >= ell and ly - ell >= 0.0 and M - 2*ell >= 0.0 and ell >= 0): return 9.0
    lows = []
    for d in (0.0, ell):
        g = max(0.0, b/6, 2*b - ly + d) if b > 0 else 0.0
        Mp, lp = M - 2*d, ell - d
        ca = (1 + 2*lp - 2*Mp)/3 if caseA == "refined" else (1 + 3*lp - 2*Mp)/4
        e = max(0.0, ca, Mp + lp - 1 - w, (Mp + 4*lp - 2)/3 - w)
        lows.append(5/6 + lx/6 - ell/6 - d + g/2 + e/2)
    S = ell / h
    dl = np.linspace(1/50, 5/6, 81); xs = np.linspace(0, 0.5, 11)
    D, X = np.meshgrid(dl, xs, indexing="ij")
    R = R_of(D, X, S, cfg, nt=31)
    high = ((1 + D)*(1 - ly)/2 - ell/2 + X*D*ell + h*(R + D/2 - 1/6)).max()
    d0 = 1/50
    floor = (1 + d0)*(1 - ly)/2 - ell/2 + d0*ell/2 + h*(1 + d0/2 - 1/6)
    return max(max(lows), high, floor)
for w in (2/7, 0.28493):
    best = (9, None)
    for ell in np.linspace(0.15, 0.40, 26):
        lx = 1/7 + ell
        for ly in np.linspace(0.30, 0.80, 51):
            t = theta_geom(lx, ly, ell, w)
            if t < best[0]: best = (t, (lx, ly, ell))
    # local polish with h fixed
    f = lambda v: theta_geom(1/7 + v[0], v[1], v[0], w)
    r = minimize(f, x0=[best[1][2], best[1][1]], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-8, maxiter=400))
    print(f"w={w:.5f}: best theta with h forced to 6/7 = {min(best[0], r.fun):.6f} at ell,ly = {r.x}  (6/7 = {6/7:.6f})")
