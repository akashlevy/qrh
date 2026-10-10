"""
counts-new: generalized row-count laboratory for P78 Part II.

R(delta,x) = min_{t in [tmin,tmax]} max_{witness (r,m)} min_{bounds} f(r,m)
with the paper's pieces (fA, fI, fP, fP0; identical to slack-audit/model_sa.py R_batch) plus
optional NEW pieces (all exact: every piece is a max/min of affine functions of r, m = t - r):

  gm   : Guth-Maynard large values estimate transplanted to the sextic family, WITH sixth-power copies,
         |W| <= N^2 V^-2 + N^{18/5} V^-4 + H N^{12/5} V^-4 on rows w <= H, primitive count = (U/H)^{1/6}|W|,
         optimized over H >= U.  Closed form (U-exponents, N = U^s, V^2 = N^{1+delta}):
           cGM(s) = max( 1 + s(2/5 - 2d),  1/6 + s(9/10 - 7d/6),  1/6 + s(7/5 - 2d) ).
         Applied to M_r (s = r) and to plain powers S_m^k (s = k m, k = 1..kmax), coefficient-blind.
  hux  : Huxley large values  |W| <= N^2 V^-2 + H N^4 V^-6, with copies:
           cHux(s) = max( 1 + s(1 - 3d),  1/6 + s(1 - 4d/3) ).
  hm   : Halasz-Montgomery with GLH-sized off-diagonal Gram entries D^{1/2}: |W| <= N^{1-d} if d >= 1/2.
           cHM(s) = s(1-d) for d >= 1/2 (else no information).  (Copies cancel in HM: see notes.)
  p6   : plain sixth moment BEYOND the edge, sum_u |S_m|^6 << U for m <= m6:  1 - 3 d m.
  p2k  : plain 2k-th moment beyond the edge for m <= mk2: 1 - k d m (k = p2k).
  hk   : 2k-th moment of the Moebius polynomial over the full family, diagonal for H >= max(U, D^k), with copies:
           max(1 - k d r, 1/6 + k r (5/6 - d))   (k = hk; HYPOTHETICAL even within the edge).
Joint (M_r S_m product) pieces need the witness freedom r + m >= t (not = t); see joint_R().
"""
import sys, numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/slack-audit")
import model_sa as MS

NEWDEF = dict(gm=False, hux=False, hm=False, p6=False, m6=0.5, p2k=0, mk2=0.5, hk=0, kmax=4)


def _affines_new(delta, q, t, cfg, one):
    """list of affine components (a, b) in r (value a + b r) of the NEW pieces, for breakpoints."""
    comps = []
    d = delta
    def s_r(c0, c1):          # piece c0 + c1 * s with s = r
        comps.append((c0 * one, c1 * one))
    def s_km(c0, c1, k):      # piece c0 + c1 * s with s = k m = k (t - r)
        comps.append((c0 + c1 * k * t, -c1 * k * one))
    if cfg.get("gm"):
        for (c0, c1) in [(1.0, 0.4 - 2 * d), (1 / 6, 0.9 - 7 * d / 6), (1 / 6, 1.4 - 2 * d)]:
            s_r(c0, c1)
            for k in range(1, cfg["kmax"] + 1):
                s_km(c0, c1, k)
    if cfg.get("hux"):
        for (c0, c1) in [(1.0, 1 - 3 * d), (1 / 6, 1 - 4 * d / 3)]:
            s_r(c0, c1)
            for k in range(1, cfg["kmax"] + 1):
                s_km(c0, c1, k)
    if cfg.get("hm"):
        s_r(0.0, 1 - d)
        for k in range(1, cfg["kmax"] + 1):
            s_km(0.0, 1 - d, k)
    if cfg.get("p6"):
        s_km(1.0, -1.5 * d, 2)        # 1 - 3 d m
    if cfg.get("p2k"):
        k = cfg["p2k"]; s_km(1.0, -d * k / 2, 2)
    if cfg.get("hk"):
        k = cfg["hk"]; s_r(1 / 6, k * (5 / 6 - d)); s_r(1.0, -k * d)
    return comps


def _eval_new(r, d_, t_, cfg):
    """value of the NEW bounds at witness r (m = t - r); returns array (inf where no info)."""
    m = np.clip(t_ - r, 0, None)
    out = np.full(np.broadcast(r, d_).shape, np.inf)
    def cgm(s):
        return np.maximum(np.maximum(1 + s * (0.4 - 2 * d_), 1 / 6 + s * (0.9 - 7 * d_ / 6)), 1 / 6 + s * (1.4 - 2 * d_))
    def chux(s):
        return np.maximum(1 + s * (1 - 3 * d_), 1 / 6 + s * (1 - 4 * d_ / 3))
    def chm(s):
        return np.where(d_ >= 0.5, s * (1 - d_), np.inf)
    for fn, key in [(cgm, "gm"), (chux, "hux"), (chm, "hm")]:
        if cfg.get(key):
            out = np.minimum(out, fn(r))
            # plain: best power k (coefficient-blind estimate applied to S_m^k); m>0 only
            best = np.full_like(out, np.inf)
            for k in range(1, cfg["kmax"] + 1):
                best = np.minimum(best, fn(k * m))
            out = np.minimum(out, np.where(m > 1e-12, best, np.inf))
    if cfg.get("p6"):
        out = np.minimum(out, np.where(m <= cfg["m6"] + 1e-15, 1 - 3 * d_ * m, np.inf))
    if cfg.get("p2k"):
        k = cfg["p2k"]
        out = np.minimum(out, np.where(m <= cfg["mk2"] + 1e-15, 1 - k * d_ * m, np.inf))
    if cfg.get("hk"):
        k = cfg["hk"]
        out = np.minimum(out, np.maximum(1 - k * d_ * r, 1 / 6 + k * r * (5 / 6 - d_)))
    return out


def R_gen(delta, q, S, t, cfg):
    """exact max over r in [t-1/2, t] of min of all bounds (paper + new); chunked along axis 0 to bound memory."""
    delta, q, t = np.broadcast_arrays(np.asarray(delta, float), np.asarray(q, float), np.asarray(t, float))
    ncomp = 9 + len(_affines_new(delta.ravel()[:1], q.ravel()[:1], t.ravel()[:1], cfg, np.ones(1)))
    work = delta.size * (ncomp * (ncomp - 1) // 2 + 16)
    if delta.ndim >= 1 and delta.shape[0] > 1 and work > 3e7:
        nch = int(np.ceil(work / 3e7))
        idx = np.array_split(np.arange(delta.shape[0]), min(nch, delta.shape[0]))
        return np.concatenate([_R_gen(delta[i], q[i], S, t[i], cfg) for i in idx], axis=0)
    return _R_gen(delta, q, S, t, cfg)


def _R_gen(delta, q, S, t, cfg):
    kA, w1, w2, cP = cfg["kA"], cfg["w1"], cfg["w2"], cfg["cP"]
    mk, w2s = cfg["mk"], cfg["w2scale"]
    c1, c2 = (1 + w1) * mk, (3 + w2) * w2s
    lo = np.maximum(t - 0.5, 0.0); hi = t
    one = np.ones_like(delta)
    pieces = [(one, -delta), (one / 6, 5 * kA / 6 - delta),
              (one - 2 * q * c1, -delta + 2 * q * mk), (one - 2 * q * c2, -delta + 4 * q * w2s),
              (one - 2 * q * S, -delta),
              (one - 2 * delta * t - 2 * q * cP * (1 - 2 * t), 2 * delta - 4 * q * cP),
              (one - 2 * delta * t - 2 * q * S, 2 * delta), (one - 2 * delta * t, 2 * delta)]
    pieces += _affines_new(delta, q, t, cfg, one)
    cands = [lo, hi]
    bps = [1 / kA * one, (c1 - S) / mk * one, (c2 - S) / (2 * w2s) * one, c1 / mk * one,
           c2 / (2 * w2s) * one, (1 + w1) * one, t - (1 - S / cP) / 2, t - 0.5, one]
    if mk != 2 * w2s:
        bps.append((c1 - c2) / (mk - 2 * w2s) * one)
    if cfg.get("p6"):
        bps.append(t - cfg["m6"])
    if cfg.get("p2k"):
        bps.append(t - cfg["mk2"])
    if cfg.get("hk"):
        bps.append(one / cfg["hk"])
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
    extra = [np.clip(((1 + w1) * one)[..., None] + eps, lo[..., None], hi[..., None])]
    for e in (cfg.get("m6") if cfg.get("p6") else None, cfg.get("mk2") if cfg.get("p2k") else None):
        if e is not None:
            extra.append(np.clip((t - e)[..., None] - eps, lo[..., None], hi[..., None]))
    if cfg.get("hk"):
        extra.append(np.clip((one / cfg["hk"])[..., None] + eps, lo[..., None], hi[..., None]))
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
    F = np.minimum(F, _eval_new(r, d_, t_, cfg))
    return F.max(axis=-1)


def install():
    """route model_sa's count through R_gen (cfg keys of NEWDEF are passed through)."""
    MS.R_batch = R_gen
    for k, v in NEWDEF.items():
        MS.DEFAULT.setdefault(k, v)


# ---------------------------------------------------------------- joint (product) pieces, brute force
def joint_R(delta, x, S, cfg, lamJ=None, nt=101, nr=801, nm=201, sfree=True):
    """count with a JOINT moment for the product witness M_r S_m (total length s = r + m):
         joint(s) = 1 - delta s                       (s <= 1)
                  = 1 + lamJ (s - 1) - delta s        (s > 1)
       The adversary may take r + m = s anywhere in [t, t + 1/2] (detector Prop: r + m >= t).
       Brute force grid (for checks, x scalar)."""
    d = float(delta); q = x * d
    kA, w1, cP, mk, w2s = cfg["kA"], cfg["w1"], cfg["cP"], cfg["mk"], cfg["w2scale"]
    best = np.inf; tbest = None
    for t in np.linspace(cfg.get("tmin", 1.0), cfg["tmax"], nt):
        r = np.linspace(max(t - 0.5, 0), t, nr)[:, None]
        if sfree:
            m = np.linspace(0, 0.5, nm)[None, :]
            ok = (r + m >= t - 1e-12)
        else:
            m = np.clip(t - r, 0, 0.5); ok = np.ones_like(m, bool)
        fA = np.maximum(1.0, (1 + 5 * kA * r) / 6) - d * r
        zI = np.clip(np.minimum(np.minimum(mk * (1 + w1 - r), w2s * (3 - 2 * r)), S), 0, None)
        fI = np.where(r < 1 + w1, 1 - d * r - 2 * q * zI, np.inf)
        zP = np.minimum(cP * np.clip(1 - 2 * m, 0, None), S)
        fP = np.where(m < 0.5 + 1e-12, 1 - 2 * d * m - 2 * q * zP, np.inf)
        F = np.minimum(np.minimum(fA, fI), fP)
        F = np.minimum(F, 1 - 2 * d * m)
        if lamJ is not None:
            s = r + m
            F = np.minimum(F, np.where(s > 1, 1 + lamJ * (s - 1) - d * s, 1 - d * s))
        F = np.where(ok, F, -np.inf)
        v = F.max()
        if v < best:
            best, tbest = v, t
    return best, tbest


def R0(delta, alpha=5 / 6):
    return 1 - 2 * delta * alpha / (3 * alpha - delta)


if __name__ == "__main__":
    import json
    install()
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    sg, geom = MS.theta(over)
    print(json.dumps(dict(cfg=over, theta=sg, geom=geom)))
