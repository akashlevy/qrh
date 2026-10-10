"""Referee check: per-block route mixing for the zeta-opening ('zexp') route of B_bilinear.md.
modes:
  zexp0_zero : zexp theta pieces at d=0 (lam=0), NO theta piece at d=ell (lam=1)   [lower bound for any (zexp, *) mix]
  minroute   : at each block, g = min(CS, zexp) handled by enumerating the 4 (route0, route1) combos -> min
"""
import sys, inspect
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import numpy as np
import theta_split_model as T
import model_r4 as M
from scipy.optimize import minimize_scalar
src = inspect.getsource(T.lp_for_S)
src = src.replace('def lp_for_S(S, cfg, mode, return_geom=False, Rcache=None):',
                  'def lpm(S, cfg, mode, return_geom=False, Rcache=None, dropq=True):')
src = src.replace('        gpieces.append(add(B, LY, d, coef=(2, -1, 1)))',
                  '        if not dropq: gpieces.append(add(B, LY, d, coef=(2, -1, 1)))')
src = src.replace('        if mode == "none":',
                  '        r = mode[0] if lam == 0.0 else mode[1]\n'
                  '        if r == "cs":\n'
                  '            gpieces.append(add(B, coef=(1 / 6,)))\n'
                  '        elif r == "zero":\n'
                  '            pass\n'
                  '        elif r == "zexp":\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 1.5, -1)))\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 3, -2)))\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 10 / 3, -2)))\n'
                  '        if False:')
ns = dict(T.__dict__); exec(src, ns); lpm = ns["lpm"]
def th(mode):
    cfg = dict(M.DEFAULT); cache = {}
    Sg = np.concatenate([np.linspace(1e-5, 0.0005, 6), np.linspace(0.001, 0.6, 50)])
    vals = [lpm(S, cfg, mode, Rcache=cache) for S in Sg]
    k = int(np.argmin(vals)); a = Sg[max(k - 1, 0)]; c = Sg[min(k + 1, len(Sg) - 1)]
    r = minimize_scalar(lambda S: lpm(S, cfg, mode), bounds=(a, c), method="bounded", options=dict(xatol=1e-8))
    Sb = r.x if r.fun < vals[k] else Sg[k]
    return lpm(Sb, cfg, mode, return_geom=True)
if __name__ == "__main__":
    for mode in [("zexp", "zero"), ("zexp", "zexp"), ("cs", "zero")]:
        sg, g = th(mode)
        print(mode, "theta=%.6f" % sg, {k: round(float(v), 4) for k, v in g.items()}, flush=True)
