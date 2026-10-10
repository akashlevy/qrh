import types, sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import os; os.chdir("/Users/akashlevy/Documents/riemann/gram_split")
import model_r4 as M
base = open("theta_split_model.py").read().split("if __name__")[0]
def mk(dropq):
    src = base.replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass") if dropq else base
    mod = types.ModuleType("t"); exec(compile(src, "t", "exec"), mod.__dict__); return mod
cfg = dict(M.DEFAULT)
print("DEFAULT cfg:", cfg)
for dropq in (False, True):
    mod = mk(dropq)
    for mode in ("none", "LSmarked"):
        best = None
        for S in [1e-6, 1e-4, 1e-3, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]:
            sg, g = mod.lp_for_S(S, cfg, mode, return_geom=True)
            if best is None or sg < best[0]: best = (sg, S, g)
        sg0, g0 = mod.lp_for_S(1e-9, cfg, mode, return_geom=True)
        print("dropP2Y=%-5s mode=%-9s  S~0: theta=%.6f (lx=%.4f ly=%.4f) | best over S grid: theta=%.6f at S=%.4g lx=%.4f ly=%.4f ell=%.4f"
              % (dropq, mode, sg0, g0['lx'], g0['ly'], best[0], best[1], best[2]['lx'], best[2]['ly'], best[2]['ell']))
