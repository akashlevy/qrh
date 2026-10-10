import types, numpy as np
import model_r4 as M
base = open("theta_split_model.py").read().split("if __name__")[0]
def make(dropq, zero_slot_R=False):
    src = base.replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass") if dropq else base
    if zero_slot_R:
        src = src.replace("RR = M.R_of(DD, XX, S, cfg)", "RR = 1 - 2 * DD / 3")
    mod = types.ModuleType("t"); exec(compile(src, "t", "exec"), mod.__dict__); return mod
for dropq in (False, True):
    mod = make(dropq)
    for d0 in (1/50, 1/1000):
        cfg = dict(M.DEFAULT); cfg["delta0"] = d0
        res = []
        for S in np.linspace(7/37 + 0.002, 0.4, 40):
            sg, g = mod.lp_for_S(S, cfg, "LSmarked", return_geom=True); res.append((sg, S, g))
        sg, S, g = min(res, key=lambda r: r[0])
        print("S>=7/37  dropP2Y=%-5s delta0=%.3f  theta=%.6f S=%.4f lx=%.4f ly=%.4f ell=%.4f" % (dropq, d0, sg, S, g['lx'], g['ly'], g['ell']))
    for d0 in (1/50, 1/1000):
        mod0 = make(dropq, zero_slot_R=True)
        cfg = dict(M.DEFAULT); cfg["delta0"] = d0
        sg, g = mod0.lp_for_S(1e-9, cfg, "LSmarked", return_geom=True)
        print("S=0 (R=1-2d/3) dropP2Y=%-5s delta0=%.3f  theta=%.6f lx=%.4f ly=%.4f" % (dropq, d0, sg, g['lx'], g['ly']))
# paper baseline restricted to S>=7/37 for reference
mod = make(False)
src = open("theta_split_model.py").read().split("if __name__")[0]
cfg = dict(M.DEFAULT)
res = [mod.lp_for_S(S, cfg, "none") for S in np.linspace(7/37 + 0.002, 0.4, 40)]
print("paper lemmas, S>=7/37: theta=%.6f" % min(res))
