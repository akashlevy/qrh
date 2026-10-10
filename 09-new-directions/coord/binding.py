"""Report binding high-side (delta,x) bins and their LP duals at the optimum of the r4 model."""
import sys, json, numpy as np
sys.path.insert(0, "/home/user/qrh/03-six-sevenths/routes/r4")
import model as M
from scipy.optimize import linprog

def analyze(over):
    cfg = dict(M.DEFAULT); cfg.update(over)
    sg, geom = M.theta(over)
    S = geom["S"]
    # rebuild the high-side constraint values at the optimal geometry
    lx, ly, ell, h = geom["lx"], geom["ly"], geom["ell"], geom["h"]
    deltas = np.linspace(cfg["delta0"], cfg["dmax"], cfg["dgrid"]); xs = np.linspace(0, 0.5, cfg["xgrid"])
    DD, XX = np.meshgrid(deltas, xs, indexing="ij")
    RR = M.R_of(DD, XX, S, cfg)
    val = (1 + DD) * (1 - ly) / 2 - ell / 2 + XX * DD * ell + h * (RR + DD / 2 - 1 / 6)
    out = []
    idx = np.argsort(-val.ravel())[:12]
    for k in idx:
        i, j = np.unravel_index(k, val.shape)
        out.append((round(float(DD[i, j]), 4), float(XX[i, j]), round(float(RR[i, j]), 5), round(float(val[i, j]), 6)))
    # low side value at geometry
    return sg, geom, out

if __name__ == "__main__":
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    sg, geom, out = analyze(over)
    print("theta", sg, geom)
    print("top high-side bins (delta, x, R, value):")
    for o in out: print(o)
