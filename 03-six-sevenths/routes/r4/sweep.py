import json, sys, numpy as np, multiprocessing as mp
from scipy.optimize import minimize_scalar
from model import DEFAULT, lp_for_S

def sweep(args):
    count_cfg, low_cfgs = args
    Sgrid = np.round(np.arange(0.05, 0.62, 0.01), 4)
    cache = {}
    out = []
    for low in low_cfgs:
        cfg = dict(DEFAULT); cfg.update(count_cfg); cfg.update(low)
        vals = np.array([lp_for_S(S, cfg, Rcache=cache) for S in Sgrid])
        k = int(np.argmin(vals))
        a = Sgrid[max(k - 1, 0)]; c = Sgrid[min(k + 1, len(Sgrid) - 1)]
        res = minimize_scalar(lambda S: lp_for_S(S, cfg), bounds=(a, c), method="bounded", options=dict(xatol=1e-7))
        Sb = res.x if res.fun < vals[k] else Sgrid[k]
        sg, geom = lp_for_S(Sb, cfg, return_geom=True)
        out.append(dict(count=count_cfg, low=low, theta=float(sg), geom={kk: float(v) for kk, v in (geom or {}).items()}))
    return out

if __name__ == "__main__":
    spec = json.load(open(sys.argv[1]))
    jobs = [(c, spec["low"]) for c in spec["counts"]]
    with mp.Pool(min(len(jobs), 11)) as p:
        res = p.map(sweep, jobs)
    flat = [r for rr in res for r in rr]
    json.dump(flat, open(sys.argv[2], "w"), indent=1)
    for r in flat:
        print(f"{json.dumps(r['count']):45s} {json.dumps(r['low']):55s} theta={r['theta']:.7f}  6/7-gap={r['theta']-6/7:+.2e}  geom={ {k: round(v,5) for k,v in r['geom'].items()} }")
