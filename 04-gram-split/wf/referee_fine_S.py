import types, sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import os; os.chdir("/Users/akashlevy/Documents/riemann/gram_split")
import model_r4 as M
from scipy.optimize import minimize_scalar
base = open("theta_split_model.py").read().split("if __name__")[0]
def mk(dropq):
    src = base.replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass") if dropq else base
    mod = types.ModuleType("t"); exec(compile(src, "t", "exec"), mod.__dict__); return mod
cfg = dict(M.DEFAULT)
for dropq in (False, True):
    mod = mk(dropq)
    cache = {}
    Sg = np.linspace(0.10, 0.30, 41)
    vals = [mod.lp_for_S(S, cfg, "none", Rcache=cache) for S in Sg]
    k = int(np.argmin(vals))
    res = minimize_scalar(lambda S: mod.lp_for_S(S, cfg, "none"), bounds=(Sg[max(k-1,0)], Sg[min(k+1,len(Sg)-1)]), method="bounded", options=dict(xatol=1e-7))
    sg, g = mod.lp_for_S(res.x, cfg, "none", return_geom=True)
    print("dropP2Y=%s grid-min %.6f at S=%.4f ; refined %.6f at S=%.5f geom lx=%.4f ly=%.4f ell=%.4f" % (dropq, min(vals), Sg[k], sg, res.x, g['lx'], g['ly'], g['ell']))
    # also check whether P^2/Y piece is slack at optimum
    lx, ly, ell = g['lx'], g['ly'], g['ell']
    b = ly - lx
    for d in (0.0, ell):
        print("   d=%.4f  2b-ly+d=%.5f  b/6=%.5f" % (d, 2*b - ly + d, b/6))
