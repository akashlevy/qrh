"""Slack map R_req - R over (delta, x) at a fixed geometry (default: DIAG optimum)."""
import numpy as np, json, sys
import model_sa as M
from frontier import GEOM_DIAG, TH_DIAG
def slackmap(over, geom=GEOM_DIAG, th=TH_DIAG, ds=None, xs=None):
    cfg = dict(M.DEFAULT); cfg.update(over)
    ds = np.linspace(cfg["delta0"], cfg["dmax"], 161) if ds is None else ds
    xs = np.linspace(0, 0.5, 11) if xs is None else xs
    DD, XX = np.meshgrid(ds, xs, indexing="ij")
    R = M.R_of(DD, XX, geom["S"], cfg)
    lx, ly, ell, h = geom["lx"], geom["ly"], geom["ell"], geom["h"]
    Rreq = (th - (1 + DD) * (1 - ly) / 2 + ell / 2 - XX * DD * ell) / h - DD / 2 + 1 / 6
    return DD, XX, R, Rreq
if __name__ == "__main__":
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    DD, XX, R, Rreq = slackmap(over)
    sl = Rreq - R
    print("min slack per x (delta at min, R, Rreq):")
    for j in range(DD.shape[1]):
        i = sl[:, j].argmin()
        print(f"  x={XX[0,j]:.2f}: min slack {sl[i,j]:+.5f} at delta={DD[i,j]:.3f}  R={R[i,j]:.4f} Rreq={Rreq[i,j]:.4f}  DH={1-DD[i,j]:.4f}")
