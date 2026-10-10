"""
slack-audit: extended copy of 03-six-sevenths/routes/r4/model.py (P78 Part II exponent model).

Unchanged semantics for the original knobs (mk, w2scale, caseA, wE, gram_b, gram_q, kA, w1, w2, cP,
perfect, delta0, dmax, tmax, xgrid, dgrid, grh).  Extra knobs (all default to the original behaviour):

  k3       : add a three-copy plain moment  sum_u |S_m|^6 |Q|^2 << U  for 3m + z/cP3 <= 1 (m <= 1/3)
             -> count fP3(m) = 1 - 3 delta m - 2 q min(cP3 (1-3m), S)      (needs a plain SIXTH moment)
  cP3      : slot capacity factor for k3 (default = cP)
  interm   : add the paper's intermediate-row constraint (d <= 1/2, t = 1, no slots, R = 1 - 2delta/3,
             z-contour 17/50, q = delta/2), eq. (adaptive-intermediate)
  tmin     : lower end of the detector cutoff range (paper: 1)
  drop_dl  : drop the d = ell low branch (diagnostic)
  side     : keep the side conditions M-2ell>=.05, ly-ell>=.05, lx>=ell (default True)
  dl_only  : (diagnostic) only d = 0 low branch -> same as drop_dl
  ramp     : marked capacity uses z <= mk*(m - r) only below r <= rcap (diagnostic)
Return of lp_for_S(..., return_duals=True) gives labelled LP multipliers.
"""
import numpy as np
from scipy.optimize import linprog, minimize_scalar

DEFAULT = dict(mk=0.5, w2scale=0.125, caseA="paper", wE=0.0, gram_b=1.0, gram_q=1.0,
               kA=1.0, w1=0.0, w2=0.0, cP=2/9, perfect=False, delta0=1/50,
               dmax=5/6, tmax=1.5, xgrid=11, dgrid=161, grh=False,
               k3=False, cP3=None, interm=False, tmin=1.0, drop_dl=False, side=True)

def R_batch(delta, q, S, t, cfg):
    delta, q, t = np.broadcast_arrays(np.asarray(delta, float), np.asarray(q, float), np.asarray(t, float))
    kA, w1, w2, cP = cfg["kA"], cfg["w1"], cfg["w2"], cfg["cP"]
    mk, w2s = cfg["mk"], cfg["w2scale"]
    k3 = cfg.get("k3", False); cP3 = cfg.get("cP3") if cfg.get("cP3") is not None else cP
    c1, c2 = (1 + w1) * mk, (3 + w2) * w2s
    lo = np.maximum(t - 0.5, 0.0); hi = t
    pieces = []
    one = np.ones_like(delta)
    pieces.append((one, -delta))
    pieces.append((one / 6, 5 * kA / 6 - delta))
    pieces.append((one - 2 * q * c1, -delta + 2 * q * mk))
    pieces.append((one - 2 * q * c2, -delta + 4 * q * w2s))
    pieces.append((one - 2 * q * S, -delta))
    pieces.append((one - 2 * delta * t - 2 * q * cP * (1 - 2 * t), 2 * delta - 4 * q * cP))
    pieces.append((one - 2 * delta * t - 2 * q * S, 2 * delta))
    pieces.append((one - 2 * delta * t, 2 * delta))
    if k3:
        # fP3 with m = t - r: 1 - 3 delta (t-r) - 2 q cP3 (1 - 3(t-r))
        pieces.append((one - 3 * delta * t - 2 * q * cP3 * (1 - 3 * t), 3 * delta - 6 * q * cP3))
        pieces.append((one - 3 * delta * t - 2 * q * S, 3 * delta))
        pieces.append((one - 3 * delta * t, 3 * delta))
    cands = [lo, hi]
    bps = [1 / kA * one, (c1 - c2) / (mk - 2 * w2s) * one if mk != 2 * w2s else one, (c1 - S) / mk * one,
           (c2 - S) / (2 * w2s) * one, c1 / mk * one, c2 / (2 * w2s) * one,
           (1 + w1) * one, t - (1 - S / cP) / 2, t - 0.5]
    if k3:
        bps += [t - 1 / 3, t - (1 - S / cP3) / 3]
    cands += bps
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
    extra = [np.clip((1 + w1) * np.ones_like(lo)[..., None] + eps, lo[..., None], hi[..., None])]
    if k3:
        extra.append(np.clip((t - 1 / 3)[..., None] - eps, lo[..., None], hi[..., None]))
    C = np.concatenate([C] + extra, axis=-1)
    r = C
    d_ = delta[..., None]; q_ = q[..., None]; t_ = t[..., None]
    fA = np.maximum(1.0, (1 + 5 * kA * r) / 6) - d_ * r
    zI = np.clip(np.minimum(np.minimum(mk * (1 + w1 - r), w2s * (3 + w2 - 2 * r)), S), 0, None)
    fI = np.where(r < 1 + w1, 1 - d_ * r - 2 * q_ * zI, np.inf)
    m = np.clip(t_ - r, 0, None)
    zP = np.minimum(cP * np.clip(1 - 2 * m, 0, None), S)
    fP = np.where(m < 0.5, 1 - 2 * d_ * m - 2 * q_ * zP, np.inf)
    fP0 = 1 - 2 * d_ * m
    F = np.minimum(np.minimum(fA, fI), np.minimum(fP, fP0))
    if k3:
        zP3 = np.minimum(cP3 * np.clip(1 - 3 * m, 0, None), S)
        fP3 = np.where(m <= 1 / 3 + 1e-15, 1 - 3 * d_ * m - 2 * q_ * zP3, np.inf)
        F = np.minimum(F, fP3)
    return F.max(axis=-1)

def R_of(delta, x, S, cfg, nt=61):
    if cfg.get("perfect"):
        return 1 - np.asarray(delta, float)
    delta = np.asarray(delta, float); x = np.asarray(x, float)
    tmin = cfg.get("tmin", 1.0)
    ts = np.linspace(tmin, cfg["tmax"], nt)
    vals = R_batch(delta[..., None], (x * delta)[..., None], S, ts, cfg)
    k = vals.argmin(axis=-1)
    best = vals.min(axis=-1)
    dt = (cfg["tmax"] - tmin) / (nt - 1)
    tb = ts[k]
    for _ in range(25):
        dt *= 0.6
        for s in (-1, 1):
            tt = np.clip(tb + s * dt, tmin, cfg["tmax"])
            v = R_batch(delta, x * delta, S, tt, cfg)
            better = v < best
            best = np.where(better, v, best); tb = np.where(better, tt, tb)
    return best

def lp_for_S(S, cfg, return_geom=False, Rcache=None, return_duals=False):
    s1 = S / (1 - S)
    A, b, lab = [], [], []
    def ge_sigma(cx, cy, c0, name):
        A.append([cx, cy, -1.0]); b.append(-c0); lab.append(name)
    ell = (-s1, 0.0, s1)
    h = (-1 / (1 - S), 0.0, 1 / (1 - S))
    def add(*fs, coef=None):
        cx = cy = c0 = 0.0
        for k, f in zip(coef, fs):
            cx += k * f[0]; cy += k * f[1]; c0 += k * f[2]
        return (cx, cy, c0)
    LX = (1.0, 0.0, 0.0); LY = (0.0, 1.0, 0.0); ONE = (0.0, 0.0, 1.0)
    M = add(LX, LY, coef=(1, 1)); B = add(LY, LX, coef=(1, -1))
    wE = cfg["wE"]
    lams = (0.0,) if cfg.get("drop_dl") else (0.0, 1.0)
    for lam in lams:
        d = add(ell, coef=(lam,))
        Mp = add(M, ell, coef=(1, -2 * lam)); lp = add(ell, coef=(1 - lam,))
        base = add(ONE, LX, ell, d, coef=(5 / 6, 1 / 6, -1 / 6, -1))
        gpieces = [((0.0, 0.0, 0.0), "g0")]
        if cfg["gram_b"] > 0: gpieces.append((add(B, coef=(cfg["gram_b"] / 6,)), "gram_b/6"))
        if cfg["gram_q"] > 0: gpieces.append((add(B, LY, d, coef=(2 * cfg["gram_q"], -cfg["gram_q"], cfg["gram_q"])), "gram_Pa2/Y"))
        epieces = [((0.0, 0.0, 0.0), "e0")]
        if cfg["caseA"] == "paper": epieces.append((add(ONE, lp, Mp, coef=(1 / 4, 3 / 4, -1 / 2)), "caseA_paper"))
        elif cfg["caseA"] == "refined": epieces.append((add(ONE, lp, Mp, coef=(1 / 3, 2 / 3, -2 / 3)), "caseA_ref"))
        if wE is not None and wE < 50:
            epieces.append((add(Mp, lp, ONE, coef=(1, 1, -1 - wE)), "wall"))
            if not cfg.get("drop_avg"):
                epieces.append((add(Mp, lp, ONE, coef=(1 / 3, 4 / 3, -2 / 3 - wE)), "cubicLS"))
        for g, gn in gpieces:
            for e, en in epieces:
                f = add(base, g, e, coef=(1, 0.5, 0.5))
                ge_sigma(*f, name=f"low d={'ell' if lam else 0} {gn} {en}")
    d0 = cfg["delta0"]
    if cfg.get("grh"):
        deltas = np.array([0.0]); xs = np.array([0.0])
        DD, XX = np.meshgrid(deltas, xs, indexing="ij"); RR = np.ones_like(DD)
    else:
        deltas = np.linspace(d0, cfg["dmax"], cfg["dgrid"]); xs = np.linspace(0, 0.5, cfg["xgrid"])
        DD, XX = np.meshgrid(deltas, xs, indexing="ij")
        if Rcache is not None and S in Rcache:
            RR = Rcache[S]
        else:
            RR = R_of(DD, XX, S, cfg)
            if Rcache is not None: Rcache[S] = RR
    for dl, xx, R in zip(DD.ravel(), XX.ravel(), RR.ravel()):
        f = add(ONE, LY, ell, h, coef=((1 + dl) / 2, -(1 + dl) / 2, -0.5 + xx * dl, R + dl / 2 - 1 / 6))
        ge_sigma(*f, name=f"high delta={dl:.4f} x={xx:.2f} R={R:.5f}")
    if not cfg.get("grh"):
        f = add(ONE, LY, ell, h, coef=((1 + d0) / 2, -(1 + d0) / 2, -0.5 + d0 / 2, 1 + d0 / 2 - 1 / 6))
        ge_sigma(*f, name="floor")
    if cfg.get("interm"):
        # E(1/2) with q = delta/2, z-contour 17/50:
        # sigma >= (1+dl)(1-ly)/2 - ell/2 + dl ell/2 + h(17/50-1/6) + (1/2)(R + dl/2 - 17/50)
        # interm = True/"paper": R = 1 - 2 delta/3 (paper: t = 1, no slots);
        # interm = "R0": no-slot count with optimized t (x = 0 count of Prop detector-counts), R0 = 1 - 2 d a/(3a - d)
        al = 5 * cfg["kA"] / 6
        for dl in np.linspace(d0, cfg["dmax"], 41):
            if cfg.get("interm") == "R0":
                R = 1 - 2 * dl * al / (3 * al - dl) if dl > d0 else 1.0
            else:
                R = 1 - 2 * dl / 3 if dl > d0 else 1.0
            f = add(ONE, LY, ell, h, coef=((1 + dl) / 2 + 0.5 * (R + dl / 2 - 17 / 50), -(1 + dl) / 2,
                                           -0.5 + dl / 2, 17 / 50 - 1 / 6))
            ge_sigma(*f, name=f"interm delta={dl:.3f}")
    if cfg.get("side", True):
        A.append([-1 - 2 * s1, -1, 0]); b.append(-0.05 - 2 * s1); lab.append("side M-2ell>=.05")
        A.append([-s1, -1, 0]); b.append(-0.05 - s1); lab.append("side ly-ell>=.05")
        A.append([-1 - s1, 0, 0]); b.append(-s1); lab.append("side lx>=ell")
    res = linprog([0, 0, 1], A_ub=np.array(A), b_ub=np.array(b),
                  bounds=[(0, 1), (0, 1.5), (None, None)], method="highs")
    if res.status != 0:
        return (9.0, None) if return_geom else 9.0
    if return_geom:
        lx, ly, sg = res.x
        geom = dict(lx=lx, ly=ly, ell=s1 * (1 - lx), h=(1 - lx) / (1 - S), S=S)
        if return_duals:
            mult = -res.ineqlin.marginals
            duals = [(lab[i], float(mult[i]), float(np.dot(A[i], res.x) - b[i])) for i in range(len(lab)) if mult[i] > 1e-9]
            return sg, geom, duals
        return sg, geom
    return res.fun

def theta(cfg_over=None, Sgrid=None, verbose=False, duals=False, Rcache=None):
    cfg = dict(DEFAULT); cfg.update(cfg_over or {})
    if Sgrid is None:
        Sgrid = np.linspace(0.05, 0.6, 45)
    if Rcache is None: Rcache = {}
    f = lambda S: lp_for_S(S, cfg, Rcache=Rcache)
    vals = [f(S) for S in Sgrid]
    k = int(np.argmin(vals))
    a = Sgrid[max(k - 1, 0)]; c = Sgrid[min(k + 1, len(Sgrid) - 1)]
    res = minimize_scalar(f, bounds=(a, c), method="bounded", options=dict(xatol=1e-6))
    Sb = res.x if res.fun < vals[k] else Sgrid[k]
    if duals:
        return lp_for_S(Sb, cfg, return_geom=True, Rcache=Rcache, return_duals=True)
    sg, geom = lp_for_S(Sb, cfg, return_geom=True, Rcache=Rcache)
    return sg, geom

def theta_sc(cfg_over=None, iters=3, kmin=0.75, sc_cP=True, sc_dmax=True, **kw):
    """self-consistent: plain capacity cP = 1/(6 kappa), kappa = max(kmin, 2 theta - 1); bins delta <= 2 theta - 1."""
    over = dict(cfg_over or {})
    sg, geom = theta(over)
    hist = [sg]
    for _ in range(iters):
        o2 = dict(over)
        kap = max(kmin, 2 * sg - 1)
        if sc_cP and "cP" not in over: o2["cP"] = 1 / (6 * kap)
        if sc_dmax: o2["dmax"] = 2 * sg - 1
        sg, geom = theta(o2)
        hist.append(sg)
    return sg, geom, hist

if __name__ == "__main__":
    import json, sys
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    if over.pop("_sc", False):
        kmin = over.pop("_kmin", 0.75)
        sg, geom, hist = theta_sc(over, kmin=kmin)
        print(json.dumps(dict(cfg=over, theta=sg, geom=geom, hist=hist, sc=True, kmin=kmin)))
    elif over.pop("_duals", False):
        sg, geom, du = theta(over, duals=True)
        print(json.dumps(dict(cfg=over, theta=sg, geom=geom)))
        for d in sorted(du, key=lambda z: -z[1]):
            print(f"  {d[1]:.5f}  {d[0]}")
    else:
        sg, geom = theta(over)
        print(json.dumps(dict(cfg=over, theta=sg, geom=geom)))
