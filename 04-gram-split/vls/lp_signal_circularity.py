"""C.3 check: bounding the theta pairing T by evaluating its signal vertex S and using an a-priori
zero-free half-plane Re s <= th_k gives the saving
    s0(M') = 5/6 + M'/12 - th_k      (Z-units, relative to the CS scale Q^{11/12}),
i.e. (since sigma_theta = base + b/12 with base = 5/6 + l_x/6 - l/6 - d the LP low-side base) a theta Gram
 piece g = b/6 - 2 s0 = 2 (th_k - base); the LP constraint becomes sigma >= base + max(0,g,...)/2 + e/2 >= th_k.
We feed this into the same LP as wf/lp_const_saving.py (P^2/Y piece dropped by LoA6) and check that the
LP returns theta >= th_k, i.e. the a-priori region is returned unchanged (circularity).
For comparison we also print the constant-s0 rows of lp_const_saving.py."""
import sys, inspect
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf")
import numpy as np
import theta_split_model as T
src = inspect.getsource(T.lp_for_S)
src = src.replace('def lp_for_S(S, cfg, mode, return_geom=False, Rcache=None):',
                  'def lp3(S, cfg, mode, return_geom=False, Rcache=None, thk=0.875, s0=0.0):')
src = src.replace('        gpieces.append(add(B, LY, d, coef=(2, -1, 1)))', '        pass')
src = src.replace('        if mode == "none":',
                  '        if mode == "zf":\n'
                  '            gpieces.append(add(ONE, LX, ell, d, coef=(2 * thk - 5 / 3, -1 / 3, 1 / 3, 2)))\n'
                  '        elif mode == "const":\n'
                  '            gpieces.append(add(B, ONE, coef=(1 / 6, -2 * s0)))\n'
                  '        elif mode == "none":')
ns = dict(T.__dict__); exec(src, ns); lp3 = ns["lp3"]
import model_r4 as M
from scipy.optimize import minimize_scalar
def th(mode, **kw):
    cfg = dict(M.DEFAULT); cache = {}
    Sg = np.linspace(0.0005, 0.6, 60)
    vals = [lp3(S, cfg, mode, Rcache=cache, **kw) for S in Sg]
    k = int(np.argmin(vals)); a = Sg[max(k-1,0)]; c = Sg[min(k+1,len(Sg)-1)]
    r = minimize_scalar(lambda S: lp3(S, cfg, mode, **kw), bounds=(a,c), method="bounded", options=dict(xatol=1e-7))
    Sb = r.x if r.fun < vals[k] else Sg[k]
    return lp3(Sb, cfg, mode, return_geom=True, **kw)
for s0 in [0.027]:
    sg, g = th("const", s0=s0); print("constant s0=%.3f      : theta=%.6f" % (s0, sg), {k: round(v, 4) for k, v in g.items()})
for thk in [0.874957, 0.86]:
    sg, g = th("zf", thk=thk); print("signal-vertex bound, a-priori th_k=%.6f : theta=%.6f" % (thk, sg), {k: round(v, 4) for k, v in g.items()})
