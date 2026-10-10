"""LP consequence of the A2-residue functional-equation route for (V').
The residual series has conductor N(mu)^3, so F'(mu) has a dual sum of length X = L^3/Q.
Bounding the dual trivially (Cauchy-Schwarz in s, sublattice mean-square factor K):
    T' << Q^{-1/12} L^2 K^{1/2},   versus Cauchy-Schwarz scale Q^{11/12}.
Ratio exponent (in Z): 3M-2 (K=1, uniform sublattice mean square) or 4M-5/2 (K=L, trivial).
Theta piece in the LP: g = max(0, b/6 - 2s) with -2s = 6M'-4 (K=1) or 8M'-5 (K=L).
Also 'dual_LS': hypothetical cubic-type large sieve on the dual (not available; for orientation only).
theta_true = min(LP_CS, LP_dual) by distributivity (IDEA1_NOTES section 6)."""
import sys, inspect
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import numpy as np
import theta_split_model as T
src = inspect.getsource(T.lp_for_S)
src = src.replace('def lp_for_S(S, cfg, mode, return_geom=False, Rcache=None):',
                  'def lp3(S, cfg, mode, return_geom=False, Rcache=None, dropq=True):')
src = src.replace('        gpieces.append(add(B, LY, d, coef=(2, -1, 1)))',
                  '        if not dropq: gpieces.append(add(B, LY, d, coef=(2, -1, 1)))')
src = src.replace('        if mode == "none":',
                  '        if mode == "dualK1":\n            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 6, -4)))\n'
                  '        elif mode == "dualKL":\n            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 8, -5)))\n'
                  '        elif mode == "none":')
ns = dict(T.__dict__); exec(src, ns); lp3 = ns["lp3"]
import model_r4 as M
from scipy.optimize import minimize_scalar
def th(mode, dropq=True):
    cfg = dict(M.DEFAULT); cache = {}
    Sg = np.linspace(0.0005, 0.6, 60)
    vals = [lp3(S, cfg, mode, Rcache=cache, dropq=dropq) for S in Sg]
    k = int(np.argmin(vals)); a = Sg[max(k-1,0)]; c = Sg[min(k+1,len(Sg)-1)]
    r = minimize_scalar(lambda S: lp3(S, cfg, mode, dropq=dropq), bounds=(a,c), method="bounded", options=dict(xatol=1e-7))
    Sb = r.x if r.fun < vals[k] else Sg[k]
    return lp3(Sb, cfg, mode, return_geom=True, dropq=dropq)
for mode in ["none", "dualK1", "dualKL"]:
    sg, g = th(mode)
    print(mode, "theta=%.6f" % sg, {k: round(v, 4) for k, v in g.items()})
