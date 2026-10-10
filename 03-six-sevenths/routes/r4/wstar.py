import json, sys
from model import theta
mode = sys.argv[1]
def th(w):
    over = {"wE": w}
    if mode == "refined": over["caseA"] = "refined"
    if mode == "gb0": over["gram_b"] = 0.0
    sg, g = theta(over)
    return sg, g
lo, hi = (0.27, 0.2857142857) if mode == "refined" else (0.03, 0.04)
for _ in range(9):
    mid = (lo + hi) / 2
    sg, g = th(mid)
    if sg > 6/7: lo = mid
    else: hi = mid
sg, g = th(hi)
print(json.dumps(dict(mode=mode, wstar_upper=hi, wstar_lower=lo, theta_at_upper=sg, geom={k: float(v) for k, v in g.items()})))
