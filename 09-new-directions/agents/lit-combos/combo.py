"""Batch driver for model.py (copied unchanged from 03-six-sevenths/routes/r4/model.py).

Exactness: identical to model.theta(), except that the row-count table R(delta,x;S) is cached
per (count-knob key, S) inside a worker, so configs that differ only in LOW-side knobs
(wE, gram_b, caseA) reuse the 45 grid-S tables. No approximation is introduced.

usage: python3 combo.py configs.json out.jsonl [nproc]
"""
import json, sys, os, time, multiprocessing as mp
import numpy as np
from scipy.optimize import minimize_scalar
import importlib
model = importlib.import_module(os.environ.get("MODEL", "model"))

COUNT_KEYS = ("kA", "w1", "w2", "cP", "mk", "w2scale", "perfect", "grh", "delta0", "dmax", "tmax", "xgrid", "dgrid")

def theta_cached(over, cache):
    cfg = dict(model.DEFAULT); cfg.update(over)
    Sgrid = np.linspace(0.05, 0.6, 45)
    f = lambda S: model.lp_for_S(S, cfg, Rcache=cache)
    vals = [f(S) for S in Sgrid]
    k = int(np.argmin(vals))
    a = Sgrid[max(k - 1, 0)]; c = Sgrid[min(k + 1, len(Sgrid) - 1)]
    res = minimize_scalar(f, bounds=(a, c), method="bounded", options=dict(xatol=1e-6))
    Sb = res.x if res.fun < vals[k] else Sgrid[k]
    sg, geom = model.lp_for_S(Sb, cfg, return_geom=True, Rcache=cache)
    return sg, geom

def countkey(over):
    cfg = dict(model.DEFAULT); cfg.update(over)
    return json.dumps({k: cfg[k] for k in COUNT_KEYS}, sort_keys=True)

def run_chunk(args):
    chunk, out = args
    cache = {}
    res = []
    for over in chunk:
        t0 = time.time()
        try:
            sg, g = theta_cached(over, cache)
            r = dict(cfg=over, theta=float(sg), geom={k: float(v) for k, v in g.items()}, sec=round(time.time() - t0, 1))
        except Exception as e:
            r = dict(cfg=over, error=repr(e))
        with open(out, "a") as fh:
            fh.write(json.dumps(r) + "\n")
        res.append(r)
    return res

if __name__ == "__main__":
    cfgs = json.load(open(sys.argv[1])); out = sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    chunksize = int(os.environ.get("CHUNK", "9"))
    groups = {}
    for c in cfgs:
        groups.setdefault(countkey(c), []).append(c)
    chunks = []
    for g in groups.values():
        for i in range(0, len(g), chunksize):
            chunks.append((g[i:i + chunksize], out))
    chunks.sort(key=lambda ch: -len(ch[0]))
    with mp.Pool(nproc) as p:
        for _ in p.imap_unordered(run_chunk, chunks):
            pass
