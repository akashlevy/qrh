"""
HEURISTIC 'cubic r4': the P78 Part II count machinery (copied from 03-six-sevenths/routes/r4/model.py,
imported as r4_model_copy) transplanted to the CUBIC architecture of NOTES.md sec 3:

  twist chi_A(m) = (m/A)_3,  B-row = completed Rankin-Selberg square theta(x)theta (T(A)^2),
  mu-absorption T(A)^2 gamma(psi_A) = mu(A) Hecke(A).

Changes w.r.t. r4 (all flagged HEURISTIC):
  copy exponent       1/6 -> c = 1/3   (low base, high-bin '-1/6', floor bin)
  Gram                gram_b * b/6 -> b/3 (cubic large sieve, sharp by DFDH); gram_q off by default
  energy              min over 4 orbit nodes (direct cubic LS / Poisson / RS-reflection with
                      conductor (QP)^cB / reflection+Poisson) of max(pieces), PLUS the normaliser
                      ('pole') term (2 beta' - 1)(1 + l')/3, beta' = abscissa for L(s,(alpha eta)^3)
  amplification       fA = max(1, 1/3 + (2/3) hmin kA r) - delta r ,  hmin = threshold of the cubic
                      P2-type mean square: min(3/2 [sextic-Kubota rep, cB=4],
                                               1 + (2 beta'-1)/3 [RS/A2 rep with normaliser]).
  marked inverse moment / plain fourth moment: kept exactly as in r4 (OPTIMISTIC: their cubic
                      analogues are not derived).
"""
import itertools, json, sys
import numpy as np
from scipy.optimize import linprog, minimize_scalar
import r4_model_copy as R4

DEFAULT = dict(R4.DEFAULT)
DEFAULT.update(c=1/3, cB=2, beta_p=0.875, gram_b=1.0, gram_q=0.0, caseA=None, hmin=None,
               pole=True, nodes="all")


def hmin_of(cfg):
    if cfg.get("hmin") is not None:
        return cfg["hmin"]
    return min(1.5, 1 + (2 * cfg["beta_p"] - 1) / 3) if cfg.get("pole") else 1.0


def R_batch_cubic(delta, q, S, t, cfg):
    """r4.R_batch with the amplification piece replaced by 1/3 + (2/3) hmin kA r."""
    # trick: r4 uses fA = max(1,(1+5 kA r)/6) - delta r.  We emulate (1/3 + a r) by kA' with
    # (1 + 5 kA' r)/6 = 1/3 + a r  is impossible (different intercept), so re-implement minimal copy.
    delta, q, t = np.broadcast_arrays(np.asarray(delta, float), np.asarray(q, float), np.asarray(t, float))
    kA, w1, w2, cP = cfg["kA"], cfg["w1"], cfg["w2"], cfg["cP"]
    mk, w2s = cfg["mk"], cfg["w2scale"]
    a_amp = (2 / 3) * hmin_of(cfg) * kA
    c1, c2 = (1 + w1) * mk, (3 + w2) * w2s
    lo = np.maximum(t - 0.5, 0.0); hi = t
    one = np.ones_like(delta)
    pieces = [(one, -delta), (one / 3, a_amp - delta),
              (one - 2 * q * c1, -delta + 2 * q * mk), (one - 2 * q * c2, -delta + 4 * q * w2s),
              (one - 2 * q * S, -delta),
              (one - 2 * delta * t - 2 * q * cP * (1 - 2 * t), 2 * delta - 4 * q * cP),
              (one - 2 * delta * t - 2 * q * S, 2 * delta), (one - 2 * delta * t, 2 * delta)]
    cands = [lo, hi, (2 / 3) / a_amp * one, (1 + w1) * one, t - 0.5,
             c1 / mk * one, c2 / (2 * w2s) * one, (c1 - S) / mk * one, (c2 - S) / (2 * w2s) * one]
    n = len(pieces)
    for i in range(n):
        for j in range(i + 1, n):
            a1, b1 = pieces[i]; a2, b2 = pieces[j]
            db = b1 - b2
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.where(np.abs(db) > 1e-14, (a2 - a1) / np.where(np.abs(db) > 1e-14, db, 1), lo)
            cands.append(r)
    C = np.clip(np.stack(cands, axis=-1), lo[..., None], hi[..., None])
    C = np.concatenate([C, np.clip((1 + w1) * np.ones_like(lo)[..., None] + 1e-12, lo[..., None], hi[..., None])], axis=-1)
    r = C
    d_ = delta[..., None]; q_ = q[..., None]; t_ = t[..., None]
    fA = np.maximum(1.0, 1 / 3 + a_amp * r) - d_ * r
    zI = np.clip(np.minimum(np.minimum(mk * (1 + w1 - r), w2s * (3 + w2 - 2 * r)), S), 0, None)
    fI = np.where(r < 1 + w1, 1 - d_ * r - 2 * q_ * zI, np.inf)
    m = np.clip(t_ - r, 0, None)
    zP = np.minimum(cP * np.clip(1 - 2 * m, 0, None), S)
    fP = np.where(m < 0.5, 1 - 2 * d_ * m - 2 * q_ * zP, np.inf)
    fP0 = 1 - 2 * d_ * m
    F = np.minimum(np.minimum(fA, fI), np.minimum(fP, fP0))
    return F.max(axis=-1)


def R_of(delta, x, S, cfg, nt=61):
    if cfg.get("perfect"):
        return 1 - np.asarray(delta, float)
    delta = np.asarray(delta, float); x = np.asarray(x, float)
    ts = np.linspace(1.0, cfg["tmax"], nt)
    vals = R_batch_cubic(delta[..., None], (x * delta)[..., None], S, ts, cfg)
    k = vals.argmin(axis=-1); best = vals.min(axis=-1)
    dt = (cfg["tmax"] - 1.0) / (nt - 1); tb = ts[k]
    for _ in range(25):
        dt *= 0.6
        for s in (-1, 1):
            tt = np.clip(tb + s * dt, 1.0, cfg["tmax"])
            v = R_batch_cubic(delta, x * delta, S, tt, cfg)
            better = v < best
            best = np.where(better, v, best); tb = np.where(better, tt, tb)
    return best


def lp_for_S(S, cfg, return_geom=False, Rcache=None):
    c = cfg["c"]; cB = cfg["cB"]; bp = cfg["beta_p"]
    s1 = S / (1 - S)
    def add(*fs, coef=None):
        cx = cy = c0 = 0.0
        for k, f in zip(coef, fs):
            cx += k * f[0]; cy += k * f[1]; c0 += k * f[2]
        return (cx, cy, c0)
    LX = (1.0, 0.0, 0.0); LY = (0.0, 1.0, 0.0); ONE = (0.0, 0.0, 1.0)
    ell = (-s1, 0.0, s1); h = (-1 / (1 - S), 0.0, 1 / (1 - S))
    M = add(LX, LY, coef=(1, 1)); B = add(LY, LX, coef=(1, -1))
    # high side rows (independent of node choice)
    high = []
    d0 = cfg["delta0"]
    if cfg.get("grh"):
        high.append(add(ONE, LY, ell, h, coef=(0.5, -0.5, -0.5, 1 - c)))
    else:
        deltas = np.linspace(d0, cfg["dmax"], cfg["dgrid"]); xs = np.linspace(0, 0.5, cfg["xgrid"])
        DD, XX = np.meshgrid(deltas, xs, indexing="ij")
        if Rcache is not None and S in Rcache:
            RR = Rcache[S]
        else:
            RR = R_of(DD, XX, S, cfg)
            if Rcache is not None: Rcache[S] = RR
        for dl, xx, R in zip(DD.ravel(), XX.ravel(), RR.ravel()):
            high.append(add(ONE, LY, ell, h, coef=((1 + dl) / 2, -(1 + dl) / 2, -0.5 + xx * dl, R + dl / 2 - c)))
        high.append(add(ONE, LY, ell, h, coef=((1 + d0) / 2, -(1 + d0) / 2, -0.5 + d0 / 2, 1 + d0 / 2 - c)))

    def nodes(Mp, lp):
        f = Mp; l = add(ONE, lp, coef=(1, 1)); lR = add(Mp, lp, ONE, coef=(cB, cB - 1, -1))
        Z = (0.0, 0.0, 0.0)
        nd = [[Z, add(l, f, coef=(1, -1)), add(l, f, coef=(2 / 3, -1 / 3))],
              [Z, add(l, f, coef=(1, -2 / 3))],
              [Z, add(lR, f, coef=(1, -1)), add(lR, f, coef=(2 / 3, -1 / 3))],
              [Z, add(lR, f, coef=(1, -2 / 3))]]
        if cfg.get("pole"):
            pole = add(ONE, lp, coef=((2 * bp - 1) / 3, (2 * bp - 1) / 3))
            nd = [n + [pole] for n in nd]
        return nd
    best = (9.0, None)
    rng = range(4) if cfg["nodes"] == "all" else cfg["nodes"]
    for choice in itertools.product(rng, repeat=2):
        A, b = [], []
        def ge_sigma(cx, cy, c0):
            A.append([cx, cy, -1.0]); b.append(-c0)
        for lam, ch in zip((0.0, 1.0), choice):
            d = add(ell, coef=(lam,))
            Mp = add(M, ell, coef=(1, -2 * lam)); lp = add(ell, coef=(1 - lam,))
            base = add(ONE, LX, ell, d, coef=(1 - c, c, -c, -1))
            gp = [(0.0, 0.0, 0.0)]
            if cfg["gram_b"] > 0: gp.append(add(B, coef=(cfg["gram_b"] * c,)))
            if cfg["gram_q"] > 0: gp.append(add(B, LY, d, coef=(2 * cfg["gram_q"], -cfg["gram_q"], cfg["gram_q"])))
            ep = nodes(Mp, lp)[ch]
            if cfg["caseA"] == "paper": ep = ep + [add(ONE, lp, Mp, coef=(1 / 4, 3 / 4, -1 / 2))]
            if cfg.get("third"): ep = ep + [add(Mp, lp, ONE, coef=(1 / 3, 4 / 3, -2 / 3))]
            for g in gp:
                for e in ep:
                    ge_sigma(*add(base, g, e, coef=(1, 0.5, 0.5)))
        for f in high:
            ge_sigma(*f)
        A.append([-1 - 2 * s1, -1, 0]); b.append(-0.05 - 2 * s1)
        A.append([-s1, -1, 0]); b.append(-0.05 - s1)
        A.append([-1 - s1, 0, 0]); b.append(-s1)
        res = linprog([0, 0, 1], A_ub=np.array(A), b_ub=np.array(b),
                      bounds=[(0, 1), (0, 1.5), (None, None)], method="highs")
        if res.status == 0 and res.fun < best[0]:
            lx, ly, sg = res.x
            best = (sg, dict(lx=lx, ly=ly, ell=s1 * (1 - lx), h=(1 - lx) / (1 - S), S=S, nodes=choice))
    return best if return_geom else best[0]


def theta(cfg_over=None, Sgrid=None):
    cfg = dict(DEFAULT); cfg.update(cfg_over or {})
    if Sgrid is None:
        Sgrid = np.linspace(0.02, 0.6, 30)
    cache = {}
    vals = [lp_for_S(S, cfg, Rcache=cache) for S in Sgrid]
    k = int(np.argmin(vals))
    a = Sgrid[max(k - 1, 0)]; c = Sgrid[min(k + 1, len(Sgrid) - 1)]
    res = minimize_scalar(lambda S: lp_for_S(S, cfg), bounds=(a, c), method="bounded", options=dict(xatol=1e-5))
    Sb = res.x if res.fun < vals[k] else Sgrid[k]
    return lp_for_S(Sb, cfg, return_geom=True)


if __name__ == "__main__":
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    sg, geom = theta(over)
    print(json.dumps(dict(cfg=over, theta=sg, geom=geom)))
