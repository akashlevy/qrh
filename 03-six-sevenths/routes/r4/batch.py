import json, sys, multiprocessing as mp
from model import theta
def run(over):
    try:
        sg, g = theta(over)
        return dict(cfg=over, theta=float(sg), geom={k: float(v) for k, v in g.items()})
    except Exception as e:
        return dict(cfg=over, error=repr(e))
if __name__ == "__main__":
    cfgs = json.load(open(sys.argv[1])); out = sys.argv[2]
    with mp.Pool(min(len(cfgs), 10)) as p:
        res = p.map(run, cfgs)
    json.dump(res, open(out, "w"), indent=1)
    for r in res:
        print(json.dumps(r["cfg"]), "->", r.get("theta"), r.get("geom"))
