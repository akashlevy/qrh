"""Anatomy of the binding row-count bins: which witness lengths (t, r, m) and slot lengths z are used.
usage: python3 diag_anatomy.py '<cfg json>'   (computes theta, then reports top bins)"""
import sys, json, numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/diag")
import model as M

def pieces(delta, q, S, t, r, cfg):
    mk, w2s, w1, w2, cP, kA = cfg["mk"], cfg["w2scale"], cfg["w1"], cfg["w2"], cfg["cP"], cfg["kA"]
    wz = cfg.get("wz", 0.0)
    fA = max(1.0, (1 + 5 * kA * r) / 6) - delta * r
    zI = max(0.0, min(mk * (1 + w1 + wz - r), w2s * (3 + w2 - 2 * r), S))
    fI = 1 - delta * r - 2 * q * zI if r < 1 + w1 else np.inf
    m = max(t - r, 0.0)
    zP = min(cP * max(1 - 2 * m, 0.0), S)
    fP = 1 - 2 * delta * m - 2 * q * zP if m < 0.5 else np.inf
    fP0 = 1 - 2 * delta * m
    return dict(fA=fA, fI=fI, zI=zI, fP=fP, zP=zP, fP0=fP0, m=m)

def anatomy(over, geom=None, top=6):
    cfg = dict(M.DEFAULT); cfg.update(over)
    if geom is None:
        sg, geom = M.theta(over)
    else:
        sg = None
    S = geom["S"]; lx, ly, ell, h = geom["lx"], geom["ly"], geom["ell"], geom["h"]
    deltas = np.linspace(cfg["delta0"], cfg["dmax"], cfg["dgrid"]); xs = np.linspace(0, 0.5, cfg["xgrid"])
    DD, XX = np.meshgrid(deltas, xs, indexing="ij")
    RR = M.R_of(DD, XX, S, cfg)
    val = (1 + DD) * (1 - ly) / 2 - ell / 2 + XX * DD * ell + h * (RR + DD / 2 - 1 / 6)
    idx = np.argsort(-val.ravel())[:top]
    rows = []
    for k in idx:
        i, j = np.unravel_index(k, val.shape)
        d, x = float(DD[i, j]), float(XX[i, j]); q = x * d
        # find t* (argmin over t of max over r), and r*
        best = None
        for t in np.linspace(1, cfg["tmax"], 501):
            rs = np.linspace(max(t - 0.5, 0), t, 1001)
            vals = []
            for r in rs:
                p = pieces(d, q, S, t, r, cfg)
                vals.append(min(p["fA"], p["fI"], p["fP"], p["fP0"]))
            vals = np.array(vals); kk = int(vals.argmax())
            if best is None or vals[kk] < best[0]:
                best = (vals[kk], t, rs[kk])
        Rv, t, r = best
        p = pieces(d, q, S, t, r, cfg)
        rows.append(dict(delta=round(d, 4), x=x, R=round(Rv, 5), val=round(float(val[i, j]), 6), t=round(t, 4), r=round(r, 4),
                         m=round(p["m"], 4), zI=round(p["zI"], 4), zP=round(p["zP"], 4),
                         fA=round(p["fA"], 4), fI=round(p["fI"], 4), fP=round(p["fP"], 4)))
    return sg, geom, rows

if __name__ == "__main__":
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    geom = json.loads(sys.argv[2]) if len(sys.argv) > 2 else None
    sg, geom, rows = anatomy(over, geom)
    print("theta", sg, geom)
    for r in rows: print(r)
