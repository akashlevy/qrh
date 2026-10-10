import json, sys
from frontier import slack
base = {"caseA": "refined", "delta0": 0.001}
for cap in ["paper", "none"]:
    for mk in [0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8]:
        for cP in [0.5, 0.55, 0.6, 0.65, 0.7, 0.75]:
            o = dict(base, mk=mk, cP=cP)
            if cap == "none": o["w2scale"] = 10.0
            s, d, x, R = slack(o)
            print(json.dumps(dict(cap=cap, mk=mk, cP=cP, minslack=s, delta=d, x=x, R=R))); sys.stdout.flush()
