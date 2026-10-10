import types, numpy as np
import model_r4 as M
base = open("theta_split_model.py").read().split("if __name__")[0]
for dropq in (False, True):
    src = base.replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass") if dropq else base
    mod = types.ModuleType("t"); exec(compile(src, "t", "exec"), mod.__dict__)
    cfg = dict(M.DEFAULT)
    best = None
    for S in [0.0005, 0.001, 0.002, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.055, 0.07]:
        sg, g = mod.lp_for_S(S, cfg, "LSmarked", return_geom=True)
        print("dropq=%s S=%.4f theta=%.6f lx=%.4f ly=%.4f ell=%.4f" % (dropq, S, sg, g['lx'], g['ly'], g['ell']))
