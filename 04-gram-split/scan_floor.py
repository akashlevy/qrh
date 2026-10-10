import types, numpy as np
import model_r4 as M
base = open("theta_split_model.py").read().split("if __name__")[0]
src = base.replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass")
mod = types.ModuleType("t"); exec(compile(src, "t", "exec"), mod.__dict__)
for d0 in (1/50, 1/200, 1/1000):
    cfg = dict(M.DEFAULT); cfg["delta0"] = d0
    best = (9, None)
    for S in [0.0, 0.0005, 0.002, 0.01, 0.03, 0.055]:
        sg, g = mod.lp_for_S(max(S, 1e-6), cfg, "LSmarked", return_geom=True)
        if sg < best[0]: best = (sg, g, S)
    sg, g, S = best
    print("delta0=%.4f  best theta=%.6f at S=%.4f lx=%.4f ly=%.4f ell=%.4f" % (d0, sg, S, g['lx'], g['ly'], g['ell']))
