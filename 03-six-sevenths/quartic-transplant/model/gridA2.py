import numpy as np, sys, json
sys.path.insert(0, '.')
import qmodel as Q
kappa = 1.0
cfg = Q.cfg_quartic("A", kappa=kappa, dmax=kappa)
m = Q.Model(cfg, kappa=kappa, dgrid=61, xgrid=6, Sgrid=np.linspace(0, 0.6, 31), nt=61, refine=10)
rng = np.random.default_rng(1)
best = (9, None)
cnt = 0
for _ in range(60000):
    lx, ly = rng.uniform(0, 1.1, 2); ell = rng.uniform(0, min(lx, ly))
    if ly - ell < 0.05 or lx + ly - 2 * ell < 0.05: continue
    low, hi, fl, sm = m.parts((lx, ly, ell), nd=5)
    v = max(hi, fl)
    if low < 0.999:
        cnt += 1
        th = max(low, v)
        if th < best[0]: best = (th, (lx, ly, ell), (low, hi, fl))
print("points with low<0.999:", cnt, " min max(low,high,floor) there:", best)
# polish from best
from scipy.optimize import minimize
def f(p):
    lx, ly, ell = p
    if min(lx, ly, ell) < 0 or ell > lx or ly - ell < 0.05 or lx + ly - 2 * ell < 0.05: return 9
    low, hi, fl, sm = m.parts(p, nd=5)
    return max(low, hi, fl) + 10 * max(0, low - 0.999)
r = minimize(f, best[1], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-8, maxiter=3000))
print("polished (restricted to non-trivial low side):", r.fun, r.x, m.parts(r.x, nd=21, exact=True))
