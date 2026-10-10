import numpy as np, theta_split_model as T, model_r4 as M
from scipy.optimize import minimize_scalar
cfg = dict(M.DEFAULT)
for mode in ["LSmarked"]:
    f = lambda S: T.lp_for_S(S, cfg, mode)
    Sg = np.linspace(0.02, 0.3, 57); v = [f(S) for S in Sg]; k = int(np.argmin(v))
    r = minimize_scalar(f, bounds=(Sg[max(k-1,0)], Sg[min(k+1,len(Sg)-1)]), method="bounded", options=dict(xatol=1e-9))
    sg, g = T.lp_for_S(r.x, cfg, mode, return_geom=True)
    print(mode, repr(sg), {k_: float(v_) for k_, v_ in g.items()})
