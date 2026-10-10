"""mark-design model scan (model_md.py): wall_l (conductor cost per unit mark length), dens (restricted prime density).
Outputs runs/<tag>.json; summary in runs/summary.txt; the 27/31 threshold runs (wall_l 0.69/0.685/0.68) are in runs/wl_threshold.log."""
import json, sys
from multiprocessing import Pool
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/mark-design")
import model_md

JOBS = [({}, "base"), ({"wE": 0.064}, "wE0.064"),
        ({"wall_l": 0.9}, "wl0.9"), ({"wall_l": 0.8}, "wl0.8"), ({"wall_l": 0.75}, "wl0.75"),
        ({"wall_l": 0.7}, "wl0.7"), ({"wall_l": 0.6}, "wl0.6"), ({"wall_l": 0.5}, "wl0.5"),
        ({"wall_l": 0.0}, "wl0.0"), ({"wall_l": 2.0}, "wl2.0"),
        ({"dens": 0.9}, "dens0.9"), ({"dens": 0.8}, "dens0.8"), ({"dens": 2 / 3}, "dens0.667")]

def one(job):
    cfg, tag = job
    th, geom = model_md.theta(cfg)
    rec = dict(tag=tag, cfg=cfg, theta=th, geom=geom)
    json.dump(rec, open(f"runs/{tag}.json", "w"))
    print(json.dumps(rec), flush=True)
    return rec

if __name__ == "__main__":
    jobs = JOBS if len(sys.argv) < 2 else [j for j in JOBS if j[1] in sys.argv[1:]]
    with Pool(2) as p:
        p.map(one, jobs)
