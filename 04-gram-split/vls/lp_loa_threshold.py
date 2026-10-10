"""LP consequence of the Mellin-shift reformulation of (V'):
   if  sum_{N(mu)~L} |Pi(5/12+it; mu)|^2 << L^{th} (1+|t|)^A,  Pi(z;mu) = sum_m pi(m) (mu/m)_2 N(m)^{-z}
   (conductor N(mu)^3 [H], centre 5/12 [R-]), then  -2 s(M') = max(M'-1, (2 th - 1) M' - th)
   (Q = Z^{M'}, L = Q^2/Z).  th = 3/2 is the quadratic large sieve (no saving for M' > 2/3... see report),
   th = 1 is Lindelof-on-average (= the 'random' row).  P^2/Y Gram piece dropped (LoA6), as in lp_const_saving.py."""
import sys, inspect
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import numpy as np
import theta_split_model as T
src = inspect.getsource(T.lp_for_S)
src = src.replace('def lp_for_S(S, cfg, mode, return_geom=False, Rcache=None):',
                  'def lp3(S, cfg, mode, return_geom=False, Rcache=None, th=1.5, s0=0.0):')
src = src.replace('        gpieces.append(add(B, LY, d, coef=(2, -1, 1)))', '        pass')
src = src.replace('        if mode == "none":',
                  '        if mode == "loa":\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 1, -1)))\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 2 * th - 1, -th)))\n'
                  '        elif mode == "loamarked":\n'
                  '            gpieces.append(add(B, Mp, lp, ONE, coef=(1 / 6, 1, 1, -1)))\n'
                  '            gpieces.append(add(B, Mp, lp, ONE, coef=(1 / 6, 2 * th - 1, th, -th)))\n'
                  '            gpieces.append(add(B, Mp, coef=(1 / 6, -1 / 6)))\n'
                  '        elif mode == "const":\n'
                  '            gpieces.append(add(B, ONE, coef=(1 / 6, -2 * s0)))\n'
                  '        elif mode == "none":')
assert 'mode == "loamarked"' in src and 'mode == "loa"' in src
ns = dict(T.__dict__); exec(src, ns); lp3 = ns["lp3"]
import model_r4 as MR
from scipy.optimize import minimize_scalar
def th_of(mode, **kw):
    cfg = dict(MR.DEFAULT); cache = {}
    Sg = np.linspace(0.0005, 0.6, 60)
    vals = [lp3(S, cfg, mode, Rcache=cache, **kw) for S in Sg]
    k = int(np.argmin(vals)); a = Sg[max(k-1, 0)]; c = Sg[min(k+1, len(Sg)-1)]
    r = minimize_scalar(lambda S: lp3(S, cfg, mode, **kw), bounds=(a, c), method="bounded", options=dict(xatol=1e-7))
    Sb = r.x if r.fun < vals[k] else Sg[k]
    return lp3(Sb, cfg, mode, return_geom=True, **kw)
sg, g = th_of("const", s0=0.0)
print("check: constant saving s0=0 -> theta=%.6f" % sg, {k: round(v, 4) for k, v in g.items()})
import sys as _s
modes = _s.argv[1:] or ["loa", "loamarked"]
for mode in modes:
    for th in [1.0, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.5]:
        sg, g = th_of(mode, th=th)
        print("%s th=%.2f: theta=%.6f  geometry" % (mode, th, sg), {k: round(v, 4) for k, v in g.items()}, flush=True)
