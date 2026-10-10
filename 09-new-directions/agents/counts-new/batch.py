"""Batch driver: python3 batch.py configs.json out.jsonl [nproc].  Each config: model knobs + counts_lab knobs;
special keys _tag, _duals."""
import json, sys, time, multiprocessing as mp
import counts_lab as CL
import model_sa as M
CL.install()

def run(over):
    over = dict(over); tag = over.pop("_tag", None); duals = over.pop("_duals", False)
    t0 = time.time()
    try:
        if duals:
            sg, geom, du = M.theta(over, duals=True)
            r = dict(tag=tag, cfg=over, theta=float(sg), geom={k: float(v) for k, v in geom.items()},
                     duals=sorted([(d[0], d[1]) for d in du], key=lambda z: -z[1])[:6])
        else:
            sg, geom = M.theta(over)
            r = dict(tag=tag, cfg=over, theta=float(sg), geom={k: float(v) for k, v in geom.items()})
    except Exception as e:
        r = dict(tag=tag, cfg=over, error=repr(e))
    r["sec"] = round(time.time() - t0, 1)
    return r

if __name__ == "__main__":
    cfgs = json.load(open(sys.argv[1])); out = sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    with mp.Pool(nproc) as p:
        for r in p.imap_unordered(run, cfgs):
            with open(out, "a") as fh:
                fh.write(json.dumps(r) + "\n")
            print(r.get("tag"), r.get("theta"), r.get("sec"), flush=True)
