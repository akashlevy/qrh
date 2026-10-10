import numpy as np
import theta_split_model as T
import model_r4 as M
# variants: drop the P^2/Y' piece as well (gram_q=0) under each theta mode
def theta2(mode, dropq):
    orig = T.lp_for_S
    def lp(S, cfg, mode_, return_geom=False, Rcache=None):
        return orig(S, cfg, mode_, return_geom, Rcache)
    return T.theta(mode)
import types
src = open("theta_split_model.py").read()
for mode in ["none", "LS", "random"]:
    for dropq in [False, True]:
        s2 = src.replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass" if dropq else "gpieces.append(add(B, LY, d, coef=(2, -1, 1)))")
        mod = types.ModuleType("tsm"); exec(compile(s2.split('if __name__')[0], "tsm", "exec"), mod.__dict__)
        sg, g = mod.theta(mode)
        print("%-7s drop P^2/Y=%-5s theta=%.6f  geom=%s" % (mode, dropq, sg, {k: round(float(v), 4) for k, v in g.items()}))
