"""Slack of every high-side bin at a FIXED geometry, for given count knobs.
slack(delta,x) = R_req(delta,x) - R(delta,x),  R_req from sigma0 = theta_target at that geometry.
Used to find the cheapest partial-DIAG (mk, cP, 2r+8z cap) that already reaches the DIAG optimum."""
import numpy as np, json, sys
import model_sa as M

GEOM_DIAG = dict(lx=0.41176619673983683, ly=0.4117663818739259, ell=0.17646742138623728, h=0.7647012246464004, S=0.23076649506849714)
TH_DIAG = 0.8725498113201074

def slack(over, geom=GEOM_DIAG, th=TH_DIAG, dgrid=321, xgrid=21):
    cfg = dict(M.DEFAULT); cfg.update(over)
    ds = np.linspace(cfg["delta0"], cfg["dmax"], dgrid); xs = np.linspace(0, 0.5, xgrid)
    DD, XX = np.meshgrid(ds, xs, indexing="ij")
    R = M.R_of(DD, XX, geom["S"], cfg)
    lx, ly, ell, h = geom["lx"], geom["ly"], geom["ell"], geom["h"]
    Rreq = (th - (1 + DD) * (1 - ly) / 2 + ell / 2 - XX * DD * ell) / h - DD / 2 + 1 / 6
    sl = Rreq - R
    k = np.unravel_index(sl.argmin(), sl.shape)
    return float(sl.min()), float(DD[k]), float(XX[k]), float(R[k])

if __name__ == "__main__":
    base = {"caseA": "refined", "delta0": 0.001}
    rows = []
    for cap in ["paper", "none"]:
        for mk in [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0]:
            for cP in [2/9, 0.25, 0.3, 1/3, 0.4, 0.5, 0.6, 0.75, 1.0]:
                o = dict(base, mk=mk, cP=cP)
                if cap == "none": o["w2scale"] = 10.0
                s, d, x, R = slack(o)
                rows.append(dict(cap=cap, mk=mk, cP=cP, minslack=s, delta=d, x=x, R=R))
                print(json.dumps(rows[-1])); sys.stdout.flush()
