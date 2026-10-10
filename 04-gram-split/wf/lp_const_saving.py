"""theta as a function of a hypothetical uniform power saving s0 (in Z-exponent) for the theta pairing
<rho,B> over its Cauchy-Schwarz bound, with the P^2/Y Gram piece dropped (LoA6).  Also the LS-shaped saving
scaled by a factor kappa.  Copied from theta_split_model.lp_for_S with a modified theta piece."""
import sys, inspect
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import numpy as np
import theta_split_model as T
src = inspect.getsource(T.lp_for_S)
src = src.replace('def lp_for_S(S, cfg, mode, return_geom=False, Rcache=None):',
                  'def lp2(S, cfg, mode, return_geom=False, Rcache=None, s0=0.0, dropq=True):')
src = src.replace('        gpieces.append(add(B, LY, d, coef=(2, -1, 1)))',
                  '        if not dropq: gpieces.append(add(B, LY, d, coef=(2, -1, 1)))')
src = src.replace('        if mode == "none":',
                  '        if mode == "const":\n            gpieces.append(add(B, ONE, coef=(1 / 6, -2 * s0)))\n        elif mode == "none":')
ns = dict(T.__dict__); exec(src, ns); lp2 = ns["lp2"]
import model_r4 as M
from scipy.optimize import minimize_scalar
def th(mode, s0=0.0, dropq=True):
    cfg = dict(M.DEFAULT); cache = {}
    Sg = np.linspace(0.0005, 0.6, 60)
    vals = [lp2(S, cfg, mode, Rcache=cache, s0=s0, dropq=dropq) for S in Sg]
    k = int(np.argmin(vals)); a = Sg[max(k-1,0)]; c = Sg[min(k+1,len(Sg)-1)]
    r = minimize_scalar(lambda S: lp2(S, cfg, mode, s0=s0, dropq=dropq), bounds=(a,c), method="bounded", options=dict(xatol=1e-7))
    Sb = r.x if r.fun < vals[k] else Sg[k]
    return lp2(Sb, cfg, mode, return_geom=True, s0=s0, dropq=dropq)
for s0 in [0.0, 0.002, 0.005, 0.01, 0.02, 0.04]:
    sg, g = th("const", s0)
    print("uniform saving s0=%.3f: theta=%.6f" % (s0, sg), {k: round(v, 4) for k, v in g.items()})
sg, g = th("LSmarked"); print("LSmarked (dropq):", round(sg, 6), {k: round(v, 4) for k, v in g.items()})
