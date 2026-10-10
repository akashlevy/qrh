import types, numpy as np
import model_r4 as M
base = open("theta_split_model.py").read().split("if __name__")[0]
src = base.replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass")
mod = types.ModuleType("t"); exec(compile(src, "t", "exec"), mod.__dict__)
cfg = dict(M.DEFAULT)
# mode 'none' with b/6 piece removed  => monkeypatch: use mode 'zero' (no theta piece at all)
src2 = src.replace('elif mode == "random":', 'elif mode == "zero":\n            pass\n        elif mode == "random":')
mod2 = types.ModuleType("t2"); exec(compile(src2, "t2", "exec"), mod2.__dict__)
for S in [0.0005, 0.01, 0.03, 0.06, 0.1, 0.15, 0.2, 0.25]:
    sg, g = mod2.lp_for_S(S, cfg, "zero", return_geom=True)
    print("theta-part negligible, no P2/Y: S=%.4f theta=%.6f lx=%.4f ly=%.4f ell=%.4f h=%.4f" % (S, sg, g['lx'], g['ly'], g['ell'], g['h']))
