"""Per-block route mixing in the theta piece: the LP low side has two blocks d = 0 (M' = M) and d = ell (M' = M - 2 ell).
Earlier LP runs (lp_dual_route.py) used one route for both blocks.  Here: Cauchy-Schwarz (b/6) at d = 0 and
route X at d = ell, where X in
   'zero'   : theta piece 0 at d = ell (infinite saving there; diagnostic: is the d = ell block binding at all?)
   'dualKL' : A2 dual route with trivial sublattice factor, -2s = 8M'-5
   'dualK1' : A2 dual route with K = 1, -2s = 6M'-4
   'cs'     : Cauchy-Schwarz (reproduces Liu)."""
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
                  '        if lam == 0.0 or mode == "cs":\n'
                  '            gpieces.append(add(B, coef=(1 / 6,)))\n'
                  '        elif mode == "zero":\n'
                  '            pass\n'
                  '        elif mode == "dualKL":\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 8, -5)))\n'
                  '        elif mode == "dualK1":\n'
                  '            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 6, -4)))\n'
                  '        elif mode == "none":')
ns = dict(T.__dict__); exec(src, ns); lpm = ns["lpm"]

def th(mode):
    cfg = dict(M.DEFAULT); cache = {}
    Sg = np.linspace(0.0005, 0.6, 60)
    vals = [lpm(S, cfg, mode, Rcache=cache) for S in Sg]
    k = int(np.argmin(vals)); a = Sg[max(k - 1, 0)]; c = Sg[min(k + 1, len(Sg) - 1)]
    r = minimize_scalar(lambda S: lpm(S, cfg, mode), bounds=(a, c), method="bounded", options=dict(xatol=1e-7))
    Sb = r.x if r.fun < vals[k] else Sg[k]
    return lpm(Sb, cfg, mode, return_geom=True)

if __name__ == "__main__":
    for mode in sys.argv[1:] or ["cs", "zero"]:
        sg, g = th(mode)
        print("%-7s theta=%.6f" % (mode, sg), {k: round(v, 4) for k, v in g.items()}, flush=True)
