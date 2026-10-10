import numpy as np, sys, json, time
sys.path.insert(0, '.')
import qmodel as Q
from energy_closed import e_quartic
def build(spec, kappa):
    base = spec["base"]
    if base == "theta":
        cfg = Q.cfg_theta_route(kappa=kappa, dmax=kappa)
    else:
        cfg = Q.cfg_quartic(base, kappa=kappa, dmax=kappa)
    if "counts_from" in spec:
        c2 = Q.cfg_quartic(spec["counts_from"], kappa=kappa, dmax=kappa)
        for k in ("zI", "rI", "amp", "dirs_inv", "dirs_plain"): cfg[k] = c2[k]
    if "energy" in spec:
        en = spec["energy"]
        cfg["energy"] = (lambda M, l, s=en: 0.0) if en == "zero" else (lambda M, l, s=en: e_quartic(M, l, s, True))
    if "margin" in spec: cfg["margin"] = spec["margin"]
    if "ell_frac" in spec: cfg["ell_frac"] = spec["ell_frac"]
    if "dmax1" in spec:
        cfg["dmax"] = 1.0
    return cfg
def solve(spec, theta0=0.9, iters=5):
    th = theta0; hist = []
    for it in range(iters):
        kappa = min(2 * th - 1, 1.0)
        cfg = build(spec, kappa)
        m = Q.Model(cfg, kappa=kappa, dgrid=81, xgrid=11, nt=61, refine=12)
        res = m.optimize(); res["kappa"] = kappa; res["cP"] = m.cP
        hist.append(round(res["theta"], 7))
        if abs(res["theta"] - th) < 3e-6: break
        th = res["theta"]
    res["history"] = hist
    return res
if __name__ == "__main__":
    spec = json.loads(sys.argv[1]); t0 = time.time()
    r = solve(spec, float(sys.argv[2]) if len(sys.argv) > 2 else 0.9)
    r["spec"] = spec; r["sec"] = round(time.time() - t0, 1)
    print(json.dumps(r), flush=True)
