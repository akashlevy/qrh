"""
Modified copy of the r4 model (03-six-sevenths/routes/r4/model.py) for the long-inverse agent.
Only the row count R_batch is changed; the LP is reused from the original module.

Extra cfg knobs (all default off => identical to the paper count):
  mix      : hypothetical 'fractional plain slot' mixed moment (coordinator idea):
             witness (r, m), r < 1:  f = 1 - delta r - delta zS - 2 q zP,
             zS = min(m, cap(r)), zP = min(S, cap(r) - zS), cap(r) = min((1-r)/2, (3-2r)/8)
  mixfull  : hypothetical FULL mixed moment |M_r S_m Q|^2 << U under r + 2m + 2z <= 1, 2r+8m+8z <= 3:
             f = 1 - delta (r+m) - 2 q zP  (only when m <= cap(r))
  lam,r1   : replace long count for r >= r1 by 1 - delta r + lam (r - r1)   (moment exponent 1+lam(r-r1));
             default lam = None (paper: e(r) = max(1,(1+5 kA r)/6))
  hb       : 'Heath-Brown split' count for long witnesses (rigorous, see notes): r>1 -> 1 - delta r/2
"""
import sys, numpy as np
sys.path.insert(0, "/home/user/qrh/03-six-sevenths/routes/r4")
import model as M4

def R_batch(delta, q, S, t, cfg):
    delta, q, t = np.broadcast_arrays(np.asarray(delta, float), np.asarray(q, float), np.asarray(t, float))
    kA, w1, w2, cP = cfg["kA"], cfg["w1"], cfg["w2"], cfg["cP"]
    mk, w2s = cfg["mk"], cfg["w2scale"]
    c1, c2 = (1 + w1) * mk, (3 + w2) * w2s
    lam = cfg.get("lam"); r1 = cfg.get("r1", 1.0)
    lo = np.maximum(t - 0.5, 0.0); hi = t
    one = np.ones_like(delta)
    pieces = []
    if lam is None:
        pieces.append((one, -delta)); pieces.append((one / 6, 5 * kA / 6 - delta))
    else:
        pieces.append((one, -delta)); pieces.append((one - lam * r1, lam - delta))
    pieces.append((one - 2 * q * c1, -delta + 2 * q * mk))
    pieces.append((one - 2 * q * c2, -delta + 4 * q * w2s))
    pieces.append((one - 2 * q * S, -delta))
    pieces.append((one - 2 * delta * t - 2 * q * cP * (1 - 2 * t), 2 * delta - 4 * q * cP))
    pieces.append((one - 2 * delta * t - 2 * q * S, 2 * delta))
    pieces.append((one - 2 * delta * t, 2 * delta))
    if cfg.get("mix") or cfg.get("mixfull"):
        # cap branches: mk(1+w1-r) = c1 - mk r ; w2s(3+w2-2r) = c2 - 2 w2s r
        for (ca, cb) in [(c1, -mk), (c2, -2 * w2s)]:
            pieces.append((one - delta * t - 2 * q * (ca - t), -2 * q * (cb + 1) * one))   # zS=m, zP=cap-m
            pieces.append((one - delta * ca, -delta - delta * cb))                     # zS=cap, zP=0
        pieces.append((one - delta * t - 2 * q * S, 0 * one))
        pieces.append((one - delta * t, 0 * one))
    if cfg.get("hb"):
        pieces.append((one, -delta / 2))
    cands = [lo, hi]
    for bp in [1 / kA * one, r1 * one, (c1 - c2) / (mk - 2 * w2s) * one if mk != 2 * w2s else one,
               (c1 - S) / mk * one, (c2 - S) / (2 * w2s) * one, c1 / mk * one, c2 / (2 * w2s) * one,
               (1 + w1) * one, t - (1 - S / cP) / 2, t - 0.5]:
        cands.append(bp)
    n = len(pieces)
    for i in range(n):
        for j in range(i + 1, n):
            a1, b1 = pieces[i]; a2, b2 = pieces[j]
            db = b1 - b2
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(np.abs(db) > 1e-14, (a2 - a1) / np.where(np.abs(db) > 1e-14, db, 1), lo)
            cands.append(r)
    C = np.stack(cands, axis=-1)
    C = np.clip(C, lo[..., None], hi[..., None])
    eps = 1e-12
    C = np.concatenate([C, np.clip((1 + w1) * np.ones_like(lo)[..., None] + eps, lo[..., None], hi[..., None]),
                        np.clip(np.ones_like(lo)[..., None] + eps, lo[..., None], hi[..., None])], axis=-1)
    r = C
    d_ = delta[..., None]; q_ = q[..., None]; t_ = t[..., None]
    if lam is None:
        fA = np.maximum(1.0, (1 + 5 * kA * r) / 6) - d_ * r
    else:
        fA = 1 + lam * np.clip(r - r1, 0, None) - d_ * r
    cap = np.clip(np.minimum(mk * (1 + w1 - r), w2s * (3 + w2 - 2 * r)), 0, None)
    zI = np.minimum(cap, S)
    fI = np.where(r < 1 + w1, 1 - d_ * r - 2 * q_ * zI, np.inf)
    m = np.clip(t_ - r, 0, None)
    zP = np.minimum(cP * np.clip(1 - 2 * m, 0, None), S)
    fP = np.where(m < 0.5, 1 - 2 * d_ * m - 2 * q_ * zP, np.inf)
    fP0 = 1 - 2 * d_ * m
    F = np.minimum(np.minimum(fA, fI), np.minimum(fP, fP0))
    if cfg.get("mix"):
        zS = np.minimum(m, cap); zPP = np.clip(np.minimum(S, cap - zS), 0, None)
        fM = np.where(r < 1 + w1, 1 - d_ * r - d_ * zS - 2 * q_ * zPP, np.inf)
        F = np.minimum(F, fM)
    if cfg.get("mixfull"):
        ok = (m <= cap + 1e-12) & (r < 1 + w1)
        zPP = np.clip(np.minimum(S, cap - m), 0, None)
        fM = np.where(ok, 1 - d_ * (r + m) - 2 * q_ * zPP, np.inf)
        F = np.minimum(F, fM)
    if cfg.get("hb"):
        F = np.minimum(F, np.where(r >= 1, 1 - d_ * r / 2, np.inf))
    return F.max(axis=-1)

M4.R_batch = R_batch
DEFAULT = dict(M4.DEFAULT)

def theta(over=None):
    return M4.theta(over)

if __name__ == "__main__":
    import json
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    sg, geom = theta(over)
    print(json.dumps(dict(cfg=over, theta=sg, geom=geom)))
