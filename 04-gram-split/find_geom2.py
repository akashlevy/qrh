import types, numpy as np
from scipy.optimize import minimize_scalar
import model_r4 as M
src = open("theta_split_model.py").read().split("if __name__")[0].replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass")
mod = types.ModuleType("t"); exec(compile(src, "t", "exec"), mod.__dict__)
cfg = dict(M.DEFAULT); f = lambda S: mod.lp_for_S(S, cfg, "LSmarked")
Sg = np.linspace(0.02, 0.3, 57); v = [f(S) for S in Sg]; k = int(np.argmin(v))
r = minimize_scalar(f, bounds=(Sg[max(k-1,0)], Sg[min(k+1,56)]), method="bounded", options=dict(xatol=1e-9))
print(mod.lp_for_S(r.x, cfg, "LSmarked", return_geom=True))
