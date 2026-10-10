import types, numpy as np
import model_r4 as M
src = open("theta_split_model.py").read().split("if __name__")[0]
def run(mode, over=None, dropq=False, extra=""):
    s2 = src
    if dropq:
        s2 = s2.replace("gpieces.append(add(B, LY, d, coef=(2, -1, 1)))", "pass")
    mod = types.ModuleType("t"); exec(compile(s2, "t", "exec"), mod.__dict__)
    sg, g = mod.theta(mode, over or {})
    print("%-12s %-28s dropP2Y=%-5s theta=%.6f  lx=%.4f ly=%.4f ell=%.4f M=%.4f b=%.4f" % (mode, over, dropq, sg, g['lx'], g['ly'], g['ell'], g['M'], g['b']))
run("LSmarked")
run("LSmarked", dropq=True)
run("LSmarked", {"caseA": "refined"})
run("LSmarked", {"wE": 0.1})
run("LSmarked", {"wE": 0.1, "caseA": "refined"})
run("LSmarked", {"delta0": 0.001})
run("LSmarked", {"kA": 8/9, "delta0": 0.001})
run("LSmarked", {"perfect": True, "delta0": 0.001})
run("randommarked")
run("randommarked", dropq=True)
