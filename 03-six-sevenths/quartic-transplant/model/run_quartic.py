import numpy as np, sys, json, time
sys.path.insert(0, '.')
import qmodel as Q
def solve(sieve, theta0=0.95, dmax_mode="kappa", iters=4, **kw):
    th = theta0; out = None
    for it in range(iters):
        kappa = 2 * th - 1
        cfg = Q.cfg_quartic(sieve, kappa=kappa, dmax=(kappa if dmax_mode == "kappa" else 1.0), **kw)
        m = Q.Model(cfg, kappa=kappa, dgrid=81, xgrid=11)
        res = m.optimize()
        res["kappa"] = kappa; res["cP"] = m.cP
        if abs(res["theta"] - th) < 2e-6:
            out = res; break
        th = res["theta"]; out = res
    return out
if __name__ == "__main__":
    sieve = sys.argv[1]; extra = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    t0 = time.time()
    r = solve(sieve, **extra)
    r["sieve"] = sieve; r["extra"] = extra; r["sec"] = round(time.time() - t0, 1)
    print(json.dumps(r))
