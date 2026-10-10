"""Relative saving eta on the energy-wall excess: e = (1-eta)(M'+l'-1) instead of (M'+l'-1).
Also the variant 'abs': the excess is reduced by eta*(M'+l') (saving N^eta on the column term)."""
import sys, json, numpy as np
sys.path.insert(0, "/home/user/qrh/03-six-sevenths/routes/r4")
import model as M
src = open("/home/user/qrh/03-six-sevenths/routes/r4/model.py").read()
old = "            epieces.append(add(Mp, lp, ONE, coef=(1, 1, -1 - wE)))"
assert old in src
new = ("            eta = cfg.get('eta', 0.0)\n"
       "            epieces.append(add(Mp, lp, ONE, coef=(1 - eta, 1 - eta, -(1 - eta) - wE)))")
src = src.replace(old, new)
ns = {}
exec(compile(src, "model_eta", "exec"), ns)
for eta in [float(a) for a in sys.argv[1:]]:
    sg, geom = ns["theta"]({"eta": eta})
    print(json.dumps(dict(eta=eta, theta=sg, geom={k: float(v) for k, v in geom.items()})), flush=True)
