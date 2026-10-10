import sys, inspect
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf")
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
which = sys.argv[1]
if which == "drop16":
    sg, g = th("const", 5.0, dropq=False); print("theta piece removed, P2/Y kept:", round(sg,6), {k: round(v,4) for k,v in g.items()})
elif which == "drop16q":
    sg, g = th("const", 5.0, dropq=True); print("theta piece removed, P2/Y dropped:", round(sg,6), {k: round(v,4) for k,v in g.items()})
elif which == "random":
    sg, g = th("random", 0.0, dropq=False); print("random, P2/Y kept:", round(sg,6), {k: round(v,4) for k,v in g.items()})
elif which == "rand_dropq":
    sg, g = th("random", 0.0, dropq=True); print("random, P2/Y dropped:", round(sg,6), {k: round(v,4) for k,v in g.items()})
