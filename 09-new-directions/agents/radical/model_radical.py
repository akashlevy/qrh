"""
[hls copy of wall-spectral/model_ws.py] knob hls_a (exponent a in a partial hybrid large sieve D << KL + N L^a),
  combined with Hoff (H2 saving) or, if Hoff is None, with Cauchy-Schwarz over the marks. Use with bias_relax=False.
[wall-spectral copy of energy-wall/model_wall.py] knob Hoff (offset-only energy piece Hoff*l') and knob bias_relax (default True): if False the
cubic-sieve/Patterson-bias piece (M'+4l'-2)/3 is NOT relaxed by wE (it is a genuine main term, DFH).
[energy-wall copy] adds knob wall_factor f: wall piece f*(M'+l'-1) - wE (f=1: paper/W(wE)).
R4: parametrized model of the Part II bookkeeping of the 7/8 paper (OpenAI, 30 Sep 2026).

theta(cfg) = min over geometry (l_x, l_y, ell) of sigma0 subject to

LOW  (Prop probe-low, Prop probe-gram, Lemma probe-row-norm/reflected-energy), for d in {0, ell}:
     sigma0 >= 5/6 + l_x/6 - ell/6 - d + g(d)/2 + e(d)/2
     g(d) = max(0, gram_b * b/6, gram_q*(2b - l_y + d))           b = l_y - l_x
     e(d) = max(0, caseA(M',l'), M'+l'-1-w_E, (M'+4l'-2)/3 - w_E)  M' = M-2d, l' = ell-d
            caseA = (1+3l'-2M')/4  (paper count O/2 of frozen powerful row parts)
                  = (1+2l'-2M')/3  (refined count O/6 + 2A/3; elementary, see notes)
                  = none           (hypothetical)
HIGH (Lemma high-bin, eq. common-high-exponent at d = h, Prop detector-counts):
     sigma0 >= (1+delta)(1-l_y)/2 - ell/2 + x delta ell + h (R(delta,x) + delta/2 - 1/6)
     for delta in [delta0, dmax], x in [0,1/2];  h = 1 - l_x + ell
FLOOR: sigma0 >= (1+delta0)(1-l_y)/2 - ell/2 + delta0 ell/2 + h(1 + delta0/2 - 1/6)

Row count R(delta, x) = min_{t in [1, tmax]} max_{r in [t-1/2, t]} min{ fA(r), fI(r), fP(t-r), fP0(t-r) }
   fA(r)  = max(1, (1 + 5 kA r)/6) - delta r              (sixth-power amplification; kA=1 paper)
   fI(r)  = 1 - delta r - 2 q zI(r),  r < 1 + w1          (Lemma marked)
            zI = clip(min((1+w1-r)/2, (3+w2-2r)/8, S), 0)
   fP(m)  = 1 - 2 delta m - 2 q min(cP (1-2m), S), m<1/2   (Lemma plain, two copies; cP = 2/9 paper)
   fP0(m) = 1 - 2 delta m
   q = x delta,  S = ell/h (supply at base U = Z^h)
Perfect counts option: R = 1 - delta.

For fixed supply S, ell = S(1-l_x)/(1-S) and h = (1-l_x)/(1-S); then every constraint is linear
in (l_x, l_y, sigma0): solve an LP for each S on a grid and minimize over S.
"""
import numpy as np
from scipy.optimize import linprog, minimize_scalar

DEFAULT = dict(mk=0.5, w2scale=0.125, caseA="paper", wE=0.0, gram_b=1.0, gram_q=1.0,
               kA=1.0, wall_factor=1.0, w1=0.0, w2=0.0, cP=2/9, perfect=False, delta0=1/50,
               dmax=5/6, tmax=1.5, xgrid=11, dgrid=161, grh=False)

# ------------------------------------------------------------------ row count
def R_batch(delta, q, S, t, cfg):
    """exact max over r in [t-1/2, t] of min of the count bounds; arrays broadcast (delta,q,t)."""
    delta, q, t = np.broadcast_arrays(np.asarray(delta, float), np.asarray(q, float), np.asarray(t, float))
    kA, w1, w2, cP = cfg["kA"], cfg["w1"], cfg["w2"], cfg["cP"]
    mk, w2s = cfg["mk"], cfg["w2scale"]
    c1, c2 = (1 + w1) * mk, (3 + w2) * w2s
    lo = np.maximum(t - 0.5, 0.0); hi = t
    # affine pieces in r: value = a + b r
    pieces = []
    one = np.ones_like(delta)
    pieces.append((one, -delta))                                    # fA piece 1
    pieces.append((one / 6, 5 * kA / 6 - delta))                    # fA piece 2
    pieces.append((one - 2 * q * c1, -delta + 2 * q * mk))            # fI with mk(1+w1-r)
    pieces.append((one - 2 * q * c2, -delta + 4 * q * w2s))          # fI with w2s(3+w2-2r)
    pieces.append((one - 2 * q * S, -delta))                        # fI with supply
    # fP with m = t - r: 1 - 2 delta (t-r) - 2 q cP (1 - 2(t-r))
    pieces.append((one - 2 * delta * t - 2 * q * cP * (1 - 2 * t), 2 * delta - 4 * q * cP))
    pieces.append((one - 2 * delta * t - 2 * q * S, 2 * delta))
    pieces.append((one - 2 * delta * t, 2 * delta))                 # fP0
    cands = [lo, hi]
    # breakpoints
    for bp in [1 / kA * one, (c1 - c2) / (mk - 2 * w2s) * one if mk != 2 * w2s else one, (c1 - S) / mk * one, (c2 - S) / (2 * w2s) * one, c1 / mk * one, c2 / (2 * w2s) * one,
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
    C = np.stack(cands, axis=-1)                                   # (..., K)
    C = np.clip(C, lo[..., None], hi[..., None])
    eps = 1e-12
    # also right-limits at the end of the marked range
    C = np.concatenate([C, np.clip((1 + w1) * np.ones_like(lo)[..., None] + eps, lo[..., None], hi[..., None])], axis=-1)
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
    return F.max(axis=-1)

def R_of(delta, x, S, cfg, nt=61):
    """min over t of R_t; delta, x arrays (same shape)."""
    if cfg.get("perfect"):
        return 1 - np.asarray(delta, float)
    delta = np.asarray(delta, float); x = np.asarray(x, float)
    ts = np.linspace(1.0, cfg["tmax"], nt)
    vals = R_batch(delta[..., None], (x * delta)[..., None], S, ts, cfg)    # (..., nt)
    k = vals.argmin(axis=-1)
    best = vals.min(axis=-1)
    # local refinement around best grid t
    dt = (cfg["tmax"] - 1.0) / (nt - 1)
    tb = ts[k]
    for _ in range(25):
        dt *= 0.6
        for s in (-1, 1):
            tt = np.clip(tb + s * dt, 1.0, cfg["tmax"])
            v = R_batch(delta, x * delta, S, tt, cfg)
            better = v < best
            best = np.where(better, v, best); tb = np.where(better, tt, tb)
    return best

# ------------------------------------------------------------------ LP for fixed supply
def lp_for_S(S, cfg, return_geom=False, Rcache=None):
    s1 = S / (1 - S)
    # ell = s1 (1 - lx) ; h = (1 - lx)/(1 - S)
    # variables v = (lx, ly, sigma); constraint rows a.v <= rhs  (we write sigma >= lin -> lin - sigma <= 0)
    A, b = [], []
    def ge_sigma(cx, cy, c0):           # sigma >= cx lx + cy ly + c0
        A.append([cx, cy, -1.0]); b.append(-c0)
    # helper: linear forms as (cx, cy, c0)
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
    # ---- low side at d = 0 and d = ell
    for lam in (0.0, 1.0):
        d = add(ell, coef=(lam,))
        Mp = add(M, ell, coef=(1, -2 * lam)); lp = add(ell, coef=(1 - lam,))
        base = add(ONE, LX, ell, d, coef=(5 / 6, 1 / 6, -1 / 6, -1))
        gpieces = [(0.0, 0.0, 0.0)]
        if cfg["gram_b"] > 0: gpieces.append(add(B, coef=(cfg["gram_b"] / 6,)))
        if cfg["gram_q"] > 0: gpieces.append(add(B, LY, d, coef=(2 * cfg["gram_q"], -cfg["gram_q"], cfg["gram_q"])))
        epieces = [(0.0, 0.0, 0.0)]
        if cfg["caseA"] == "paper": epieces.append(add(ONE, lp, Mp, coef=(1 / 4, 3 / 4, -1 / 2)))
        elif cfg["caseA"] == "refined": epieces.append(add(ONE, lp, Mp, coef=(1 / 3, 2 / 3, -2 / 3)))
        if wE is not None and wE < 50:
            fw = cfg.get("wall_factor", 1.0)
            if cfg.get('hls_f') is not None:
                # [hls] block-wise bound D << N + K + L*beta with square-root block norms beta ~ sqrt(K N)
                #   (unproved 'twisted quadratic sieve'): diagonal-in-P excess max(M'-1, (M'+l'-1)/2) = hls_f*(M'+l'-1)
                #   with hls_f = 1/2; off-diagonal by (H2) (Hoff).
                f_ = cfg['hls_f']
                epieces.append(add(Mp, lp, ONE, coef=(f_, f_, -f_)))
                epieces.append(add(Mp, lp, ONE, coef=(1.0, 0.0, -1.0)))
                epieces.append(add(Mp, lp, ONE, coef=(0.0, cfg.get('Hoff', 0.0), 0.0)))
            elif cfg.get('hls_a') is not None:
                # [hls] diagonal-in-P energy from a (partial) hybrid large sieve D << KL + N L^a:
                #   sum_P |c_P|^2 G_PP << (K + N L^{a-1}) N  -> excess M' + a l' - 1;
                # off-diagonal: (H2) with saving L^{-theta2} (Hoff = 1 - theta2) if given, else Cauchy-Schwarz
                #   over P (off <= L * diag) -> excess M' + (1+a) l' - 1.
                a_ = cfg['hls_a']
                if cfg.get('Hoff') is None:
                    epieces.append(add(Mp, lp, ONE, coef=(1.0, 1.0 + a_, -1.0)))
                else:
                    epieces.append(add(Mp, lp, ONE, coef=(1.0, a_, -1.0)))
                    epieces.append(add(Mp, lp, ONE, coef=(0.0, cfg['Hoff'], 0.0)))
            elif cfg.get('Hoff') is None:
                epieces.append(add(Mp, lp, ONE, coef=(fw, fw, -fw - wE)))
            else:   # (H1)+(H2_theta2): E << K N L'^{Hoff}, Hoff = 1 - theta2, valid for all N <= K L
                epieces.append(add(Mp, lp, ONE, coef=(0.0, cfg['Hoff'], 0.0)))
            epieces.append(add(Mp, lp, ONE, coef=(1 / 3, 4 / 3, -2 / 3 - (wE if cfg.get('bias_relax', True) else 0.0))))
        for g in gpieces:
            for e in epieces:
                f = add(base, g, e, coef=(1, 0.5, 0.5))
                ge_sigma(*f)
    # ---- high side
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
        # (1+dl)(1-ly)/2 - ell/2 + xx dl ell + h (R + dl/2 - 1/6)
        f = add(ONE, LY, ell, h, coef=((1 + dl) / 2, -(1 + dl) / 2, -0.5 + xx * dl, R + dl / 2 - 1 / 6))
        ge_sigma(*f)
    # ---- floor bin (R=1, q <= delta0/2)
    if not cfg.get("grh"):
        f = add(ONE, LY, ell, h, coef=((1 + d0) / 2, -(1 + d0) / 2, -0.5 + d0 / 2, 1 + d0 / 2 - 1 / 6))
        ge_sigma(*f)
    # ---- side conditions: M - 2 ell >= 0.05, l_y - ell >= 0.05, lx >= ell (d <= l_x), lx,ly >= 0
    A.append([-1, -1, 0]); b.append(-0.05); A[-1] = list(np.array(A[-1]) + np.array([-2 * s1, 0, 0]) * -1)
    # rewrite: M - 2ell >= 0.05  <=>  -(lx + ly) + 2 s1 (1 - lx) <= -0.05  <=> (-1 - 2 s1) lx - ly <= -0.05 - 2 s1
    A[-1] = [-1 - 2 * s1, -1, 0]; b[-1] = -0.05 - 2 * s1
    A.append([-s1, -1, 0]); b.append(-0.05 - s1)        # ly - ell >= 0.05 : -ly + s1(1-lx) <= -0.05
    A.append([-1 - s1, 0, 0]); b.append(-s1)            # lx >= ell : s1(1-lx) - lx <= 0
    res = linprog([0, 0, 1], A_ub=np.array(A), b_ub=np.array(b),
                  bounds=[(0, 1), (0, 1.5), (None, None)], method="highs")
    if res.status != 0:
        return (9.0, None) if return_geom else 9.0
    if return_geom:
        lx, ly, sg = res.x
        return sg, dict(lx=lx, ly=ly, ell=s1 * (1 - lx), h=(1 - lx) / (1 - S), S=S)
    return res.fun

def theta(cfg_over=None, Sgrid=None, verbose=False):
    cfg = dict(DEFAULT); cfg.update(cfg_over or {})
    if Sgrid is None:
        Sgrid = np.linspace(0.05, 0.6, 45)
    vals = [lp_for_S(S, cfg) for S in Sgrid]
    k = int(np.argmin(vals))
    a = Sgrid[max(k - 1, 0)]; c = Sgrid[min(k + 1, len(Sgrid) - 1)]
    res = minimize_scalar(lambda S: lp_for_S(S, cfg), bounds=(a, c), method="bounded", options=dict(xatol=1e-6))
    Sb = res.x if res.fun < vals[k] else Sgrid[k]
    sg, geom = lp_for_S(Sb, cfg, return_geom=True)
    return sg, geom

if __name__ == "__main__":
    import json, sys
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    sg, geom = theta(over)
    print(json.dumps(dict(cfg=over, theta=sg, geom=geom)))
