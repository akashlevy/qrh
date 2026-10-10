import numpy as np, sys, json, time
sys.path.insert(0, '.')
import qmodel as Q
def build(scen, kappa, tmin, tmax, dmax_mode):
    dmax = kappa if dmax_mode == "kappa" else 1.0
    if scen == "theta":
        return Q.cfg_theta_route(kappa=kappa, dmax=dmax, tmin=tmin, tmax=tmax)
    return Q.cfg_quartic(scen, kappa=kappa, dmax=dmax, tmin=tmin, tmax=tmax)
def solve(scen, theta0, tmin=0.3, tmax=3.0, dmax_mode="kappa", iters=5, dgrid=81, nt=61):
    th = theta0; hist = []
    for it in range(iters):
        kappa = min(2 * th - 1, 1.0)
        cfg = build(scen, kappa, tmin, tmax, dmax_mode)
        m = Q.Model(cfg, kappa=kappa, dgrid=dgrid, xgrid=11, nt=nt, refine=12)
        res = m.optimize(); res["kappa"] = kappa; res["cP"] = m.cP
        hist.append(round(res["theta"], 7))
        if abs(res["theta"] - th) < 3e-6: break
        th = res["theta"]
    res["history"] = hist
    return res
if __name__ == "__main__":
    scen = sys.argv[1]; theta0 = float(sys.argv[2]); kw = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    t0 = time.time(); r = solve(scen, theta0, **kw); r["scenario"] = scen; r["kw"] = kw; r["sec"] = round(time.time() - t0, 1)
    print(json.dumps(r), flush=True)
