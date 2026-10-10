"""combo_model.py -- consolidated exponent model of the WEIGHTED high side (agent `combo`, 10 Oct 2026).

One model with every knob found in session 10, defaults = PROVED inputs.

Geometry (base Z): slot length ell, M = l_x + l_y, b = l_y - l_x, h = 1 + ell - l_x, U = Z^d, d = h, L = ell/h.
Low side (P78 + 08, Gram factor g, lossless wall relaxation w):
    B_low = 5/6 + M/12 - ell/6 - b(1-g)/12 + (M + ell - 1 - w)_+/6 ,   default M = 1 + w - ell (relaxed wall)
    side conditions (08 certificate): Gram domination l_y - ell - (2 - g/6) b > 0, l_y < 1/2, 2 l_y > h, ...
High side, bin (delta, x), a = (1+delta)/2, q = delta x:
    E + B = a + h(z_c - 1/6) - ell/2 + q ell - l_y/2 + d (R + sigma - z_c),     z_c = 17/50
    sigma = min( sigma_pw = delta n0/2                           [P78 pointwise numerator == E_old],
                 sigma_4  = (1-R)/4 - q zN/2                     [08 weighted numerator, Hoelder + 4th moment],
                 sigma_6  = (1+eta-R)/6                          [HYP sixth moment, excess eta],
                 sigma_LS = (max(R, l_y/d) v ((1+R)/2 - q zN - sLS) - R)/2   [HYP restricted large sieve on C] )
    n0 = 1 - l_y/d, zN = min(L, cN (1 - 2 n0 + 6 eP)).
Count R(delta,x) = max(1-delta, min_{t in [1,3/2]} max_{r+m=t, m<=1/2} min(W_wit(r), W_plain(m))):
    inverse-marked (r<1): 1 - delta r - 2 q zM,  zM = min(supply_M, mk (1 + sM_eff - r), [(3-2r)/8 if cap82])
                          sM_eff = max(sM, w/2) (wall feeds lem:marked as r + 2z <= 1 + w/2)
    long inverse (r>=1):  1 - delta + (alpha - delta)(r - 1)
    plain (m):            1 - 2 delta m - 2 q zP,  zP = min(L, max(cP (1 - 2m + 6 eP), [1/2 - 2m if plainLS]))
kappa = min(3/4, 2 theta - 1) ('dynamic', Collet variable kappa in [7/10,3/4]); cP = cN = 1/(6 kappa) unless set.

The inner max over the adversary's split is computed EXACTLY by bisection (monotonicity of both witness
lines in r is checked: needs mk <= 1, cP <= 1); the outer min over t is a grid + golden refinement (any t is
a valid choice for us, so this can only overestimate R).
"""
import numpy as np
from dataclasses import dataclass, replace, asdict, field
from typing import Optional

ZC = 17 / 50
TARGET = 41 / 47
INF = 1e9


@dataclass(frozen=True)
class Lev:
    # counts: inverse-marked witness (lem:marked)
    mk: float = 0.5            # capacity z <= mk (1 + sM_eff - r)          [P78: 1/2]
    sM: float = 0.0            # cross-mark saving: r + 2z <= 1 + sM and 2r + 8z <= 3 + 4 sM   [P78: 0]
    cap82: bool = True         # cubic branch 2r + 8z <= 3 kept               [P78: kept]
    mark_supply: str = 'unlimited'   # 'unlimited' = slot-power lemma (PROVED, agent amplitude); 'physical' = L
    # counts + numerator: lem:plain
    cP: Optional[float] = None # plain capacity rate (None -> 1/(6 kappa))
    cN: Optional[float] = None # numerator capacity rate (None -> cP if joint else 1/(6 kappa))
    joint: bool = True         # improved lem:plain improves plain witness and numerator together (Theorem G)
    eP: float = 0.0            # power saving at lem:plain's edge: capacity c (1 - 2n + 6 eP)
    plainLS: bool = True       # sextic large sieve for short plain witnesses, 2m + z <= 1/2 (PROVED, plain-cap P2)
    # long inverse witness (lem:inverse-amplification)
    alpha: float = 5 / 6
    eL: float = 0.0            # power saving beyond the diagonal: Sum|M_r|^2 << U^{max(1, 1 + alpha(r-1) - eL)}
    # numerator
    eta: Optional[float] = None   # sixth moment excess (None = not available)
    sLS: float = 0.0           # restricted large sieve saving over the 4th-moment bound on bad rows
    Eold: bool = True          # pointwise numerator branch (P78)
    # low side
    w: float = 0.0             # lossless wall relaxation
    wall_marked: bool = True   # wall feeds lem:marked (r + 2z <= 1 + w/2)
    g: float = 1.0             # Gram factor: Gram term P_a^{g/6}
    gram_dom: bool = True      # enforce Gram domination l_y - ell - (2 - g/6) b >= 0
    # kappa
    kappa: object = 'dynamic'  # 'dynamic' -> min(3/4, 2 theta - 1); or a number (08: 3/4)
    DH: bool = False


PROVED = Lev()
PAPER08 = Lev(kappa=0.75, mark_supply='physical', plainLS=False)


def caps(lev, kappa):
    cP = lev.cP if lev.cP is not None else 1 / (6 * kappa)
    if lev.cN is not None:
        cN = lev.cN
    else:
        cN = cP if (lev.joint and lev.cP is not None) else 1 / (6 * kappa)
    return cP, cN


# ----------------------------------------------------------------------------------------------- counts
def _maxmin(f_dec, g_inc, lo, hi, it=40):
    """max_{r in [lo,hi]} min(f_dec(r), g_inc(r)), f nonincreasing, g nondecreasing (vectorized)."""
    flo, glo = f_dec(lo), g_inc(lo)
    fhi, ghi = f_dec(hi), g_inc(hi)
    a = lo.copy(); b = hi.copy()
    for _ in range(it):
        mid = 0.5 * (a + b)
        hm = f_dec(mid) - g_inc(mid)
        pos = hm >= 0
        a = np.where(pos, mid, a); b = np.where(pos, b, mid)
    cross = np.maximum(np.minimum(f_dec(a), g_inc(a)), np.minimum(f_dec(b), g_inc(b)))
    out = np.where(flo <= glo, flo, np.where(fhi >= ghi, ghi, cross))
    return out


def count_R(dl, x, L, lev, kappa, nt=41, golden=28, detail=False):
    """Row-count exponent R for arrays dl, x (same shape (n,)); L = slot supply (U-units)."""
    dl = np.asarray(dl, float)[:, None]; x = np.asarray(x, float)[:, None]
    if lev.DH:
        return (1 - dl)[:, 0]
    assert lev.mk <= 1 + 1e-12, 'monotonicity needs mk <= 1'
    q = dl * x
    cP, _ = caps(lev, kappa)
    assert cP <= 1 + 1e-12
    sMe = max(lev.sM, lev.w / 2 if lev.wall_marked else 0.0)
    LM = INF if lev.mark_supply == 'unlimited' else L
    mk, alpha, eP = lev.mk, lev.alpha, lev.eP

    def Winv(r):
        z = np.minimum(LM, mk * (1 + sMe - r))
        if lev.cap82:   # cubic branch; a cross-mark saving sM lowers P78's excess e = max(0, r+2z-1, (2r+8z-3)/4) by sM
            z = np.minimum(z, (3 + 4 * lev.sM - 2 * r) / 8)
        return 1 - dl * r - 2 * q * np.maximum(z, 0)

    def Wpl(m):
        cap = cP * (1 - 2 * m + 6 * eP)
        if lev.plainLS:
            cap = np.maximum(cap, 0.5 - 2 * m)
        return 1 - 2 * dl * m - 2 * q * np.clip(cap, 0, L)

    eL = lev.eL
    rk = 1 + (eL / alpha if alpha > 0 else 0.0)        # kink: e(r) = max(1, 1 + alpha (r-1) - eL)

    def Wlong(r):
        return np.maximum(1.0, 1 + alpha * (r - 1) - eL) - dl * r

    def F(t):
        # short branch r in [t-1/2, 1)  (m = t - r in (t-1, 1/2])
        lo = (t - 0.5) * np.ones_like(dl * t); hi = np.ones_like(lo) * (1 - 1e-12)
        Gs = _maxmin(Winv, lambda r: Wpl(t - r), lo, hi)
        # long branch r in [1, t]; Wlong decreasing on [1, rk], slope alpha - delta on [rk, t]
        lo2 = np.ones_like(lo); tt = t * np.ones_like(lo)
        k2 = np.minimum(rk, tt)
        G1 = _maxmin(Wlong, lambda r: Wpl(t - r), lo2, k2)          # decreasing part (empty if eL = 0)
        inc = (alpha - dl) >= 0
        G2_inc = np.minimum(Wlong(tt), Wpl(t - tt))                  # both nondecreasing -> r = t
        G2_dec = _maxmin(Wlong, lambda r: Wpl(t - r), k2, tt)
        G2 = np.where(inc, G2_inc, G2_dec)
        Gl = np.maximum(G1, G2)
        return np.maximum(Gs, Gl)

    ts = np.linspace(1, 1.5, nt)[None, :]
    Ft = F(ts)
    i = np.argmin(Ft, axis=1)
    best = Ft[np.arange(len(i)), i]
    tbest = ts[0, i]
    if golden:
        h = 0.5 / (nt - 1)
        a = np.clip(tbest - h, 1, 1.5)[:, None]; b = np.clip(tbest + h, 1, 1.5)[:, None]
        gr = (np.sqrt(5) - 1) / 2
        c = b - gr * (b - a); d = a + gr * (b - a)
        Fc = F(c); Fd = F(d)
        for _ in range(golden):
            left = Fc <= Fd
            b = np.where(left, d, b); a = np.where(left, a, c)
            d_new = np.where(left, c, a + gr * (b - a)); c_new = np.where(left, b - gr * (b - a), d)
            c, d = c_new, d_new
            Fc = F(c); Fd = F(d)
        cand = np.minimum(Fc, Fd)[:, 0]
        tc = np.where(Fc <= Fd, c, d)[:, 0]
        upd = cand < best
        best = np.where(upd, cand, best); tbest = np.where(upd, tc, tbest)
    R = np.maximum(best, (1 - dl)[:, 0])
    if detail:
        return R, tbest
    return R


# ----------------------------------------------------------------------------------------------- geometry
def geom(ell, b, lev, M=None):
    if M is None:
        M = 1 + lev.w - ell
    lx = (M - b) / 2; ly = (M + b) / 2; h = 1 - lx + ell
    return M, lx, ly, h


def B_low(ell, b, lev, M=None):
    M, lx, ly, h = geom(ell, b, lev, M)
    return 5 / 6 + M / 12 - ell / 6 - b * (1 - lev.g) / 12 + max(0.0, M + ell - 1 - lev.w) / 6


def side_ok(ell, b, lev, M=None, kappa=0.75, report=False):
    """08-certificate side conditions that do not involve the bin optimisation."""
    M, lx, ly, h = geom(ell, b, lev, M)
    c = {
        'ell<=1/5': 0.2 - ell,
        'lx>ell': lx - ell, 'ly>ell': ly - ell, '1-3ell>0': 1 - 3 * ell,
        'Gram domination': (ly - ell - (2 - lev.g / 6) * b) if lev.gram_dom else 1.0,
        'b>=0': b,
        'row extension 5ell>h': 5 * ell - h,
        'reflected numerator short 2ly>h': 2 * ly - h,
        'ly<1/2': 0.5 - ly,
        'small rows': -(h * (ZC - 1 / 6) - ly / 2 + 63 / 5000),
        'kappa range': (kappa - 0.7) if lev.kappa == 'dynamic' else 1.0,
    }
    ok = all(v >= -1e-12 for v in c.values())
    if report:
        return ok, {k: v for k, v in c.items() if v < 1e-12}
    return ok


# ----------------------------------------------------------------------------------------------- bins
def make_bins(dstep=0.01, xs=(0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5), dmin=0.02, dmax=0.75):
    ds = np.arange(dmin, dmax + 1e-9, dstep)
    D, X = np.meshgrid(ds, np.asarray(xs, float), indexing='ij')
    return D.ravel(), X.ravel()


BINS_SEARCH = make_bins(0.01)
BINS_FINE = make_bins(0.0025, xs=np.linspace(0, 0.5, 21))

_RC = {}
_LFREE = {}
L_FREE_LO = 0.19


def R_cached(bins, L, lev, kappa, Lq=1e-3, kq=1e-5):
    """R on bins; kappa rounded UP (conservative).  If R(L_FREE_LO) == R(inf) on all bins, R does not
    depend on L for L >= L_FREE_LO (R is nonincreasing in L) and a single array is cached; otherwise L is
    rounded DOWN to a multiple of Lq (conservative)."""
    kr = np.ceil(kappa / kq - 1e-9) * kq
    keyl = replace(lev, eta=None, sLS=0.0, Eold=True, g=1.0, gram_dom=True,
                   w=(lev.w if lev.wall_marked else 0.0))
    base = (id(bins[0]), len(bins[0]), round(kr, 7), keyl)
    if base not in _LFREE:
        Rinf = count_R(bins[0], bins[1], 10.0, lev, kr)
        Rlo = count_R(bins[0], bins[1], L_FREE_LO, lev, kr)
        _LFREE[base] = Rinf if np.max(np.abs(Rinf - Rlo)) < 1e-12 else None
    if _LFREE[base] is not None and L >= L_FREE_LO:
        return _LFREE[base], kr
    Lr = np.floor(L / Lq) * Lq if L < 10 else L
    key = base + (round(Lr, 6),)
    if key not in _RC:
        _RC[key] = count_R(bins[0], bins[1], Lr, lev, kr)
    return _RC[key], kr


FLOOR = (np.array([0.02, 0.02]), np.array([0.0, 0.5]))   # floor bin: delta <= 1/50, R = 1, q <= 1/100


def _EB(dl, x, R, ell, ly, h, d, L, lev, kap):
    _, cN = caps(lev, kap)
    q = dl * x; a = (1 + dl) / 2; n0 = max(0.0, 1 - ly / d)
    zN = min(L, max(0.0, cN * (1 - 2 * n0 + 6 * lev.eP)))
    base = a + h * (ZC - 1 / 6) - ell / 2 + q * ell - ly / 2
    s4 = (1 - R) / 4 - q * zN / 2
    sig = s4.copy()
    if lev.Eold:
        sig = np.minimum(sig, dl * n0 / 2)
    if lev.eta is not None:
        sig = np.minimum(sig, (1 + lev.eta - R) / 6)
    if lev.sLS > 0:
        rho = np.maximum(np.maximum(R, ly / d), (1 + R) / 2 - q * zN - lev.sLS)
        sig = np.minimum(sig, (rho - R) / 2)
    EB = base + d * (R + sig - ZC)
    return EB, sig, s4, n0, zN


def high_terms(ell, b, lev, kappa, bins=None, M=None, R=None, exactL=False, d=None, floor=True):
    """E + B on all bins (row dyad U = Z^d, default d = h), followed by the two floor-bin entries."""
    if bins is None:
        bins = BINS_SEARCH
    M, lx, ly, h = geom(ell, b, lev, M)
    if d is None:
        d = h
    L = ell / d
    dl, x = bins
    if R is None:
        if exactL or d != h:
            R = count_R(dl, x, L, lev, kappa); kap = kappa
        else:
            R, kap = R_cached(bins, L, lev, kappa)
    else:
        kap = kappa
    EB, sig, s4, n0, zN = _EB(dl, x, R, ell, ly, h, d, L, lev, kap)
    if floor:
        Rf = np.ones(2)
        EBf, sigf, s4f, _, _ = _EB(FLOOR[0], FLOOR[1], Rf, ell, ly, h, d, L, lev, kap)
        EB = np.concatenate([EB, EBf]); R = np.concatenate([R, Rf])
        sig = np.concatenate([sig, sigf]); s4 = np.concatenate([s4, s4f])
    return EB, dict(R=R, sig=sig, s4=s4, n0=n0, zN=zN, L=L, h=h, ly=ly, lx=lx, M=M)


def B_high(ell, b, lev, kappa, bins=None, M=None, detail=False, exactL=False, d=None):
    EB, info = high_terms(ell, b, lev, kappa, bins, M, exactL=exactL, d=d)
    i = int(np.argmax(EB))
    if detail:
        dl, x = bins if bins is not None else BINS_SEARCH
        dl = np.concatenate([dl, FLOOR[0]]); x = np.concatenate([x, FLOOR[1]])
        return EB[i], dict(delta=float(dl[i]), x=float(x[i]), R=float(info['R'][i]), sigma=float(info['sig'][i]),
                           sigma4=float(info['s4'][i]), n0=info['n0'], zN=info['zN'], L=info['L'],
                           floor=bool(i >= len(EB) - 2))
    return EB[i]


def dcheck(ell, b, lev, kappa, bins=None, fracs=(1.0, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5)):
    """max over bins of E+B at row dyads d = f*h (d = h should bind; Collet: convex in d)."""
    M, lx, ly, h = geom(ell, b, lev)
    return [(f, float(np.max(high_terms(ell, b, lev, kappa, bins, d=f * h)[0]))) for f in fracs]


def worstE(ell, b, B, lev, kappa=0.75, bins=None, M=None, detail=False, exactL=False):
    """agent-style: worst high-side exponent at fixed (ell, b, B)."""
    r = B_high(ell, b, lev, kappa, bins, M, detail, exactL)
    if detail:
        return r[0] - B, r[1]
    return r - B


# ----------------------------------------------------------------------------------------------- theta
def kappa_of(theta, lev):
    if lev.kappa == 'dynamic':
        return min(0.75, 2 * theta - 1)
    return float(lev.kappa)


def theta_at_b(b, lev, kappa, bins=None, ell_lo=0.13, ell_hi=0.23, it=30):
    """min over ell of max(B_low, B_high) at fixed b (B_high increasing, B_low decreasing in ell)."""
    def f(ell):
        return B_high(ell, b, lev, kappa, bins) - B_low(ell, b, lev)
    lo, hi = ell_lo, ell_hi
    if f(hi) < 0:
        return B_low(hi, b, lev), hi
    if f(lo) > 0:
        return B_high(lo, b, lev, kappa, bins), lo
    for _ in range(it):
        mid = (lo + hi) / 2
        if f(mid) <= 0:
            lo = mid
        else:
            hi = mid
    ell = lo
    return max(B_low(ell, b, lev), B_high(ell, b, lev, kappa, bins)), ell


def theta_fixed_kappa(lev, kappa, bins=None, b_lo=0.0, b_hi=0.17, nb=12, golden=18, side=True):
    bs = np.linspace(b_lo, b_hi, nb)
    vals = []
    for b in bs:
        th, ell = theta_at_b(b, lev, kappa, bins)
        if side and not side_ok(ell, b, lev, kappa=kappa):
            th = 9.0
        vals.append((th, b, ell))
    j = int(np.argmin([v[0] for v in vals]))
    a = bs[max(j - 1, 0)]; c = bs[min(j + 1, nb - 1)]
    gr = (np.sqrt(5) - 1) / 2

    def fb(b):
        th, ell = theta_at_b(b, lev, kappa, bins)
        if side and not side_ok(ell, b, lev, kappa=kappa):
            th = 9.0
        return th, ell
    x1 = c - gr * (c - a); x2 = a + gr * (c - a)
    f1 = fb(x1); f2 = fb(x2)
    for _ in range(golden):
        if f1[0] <= f2[0]:
            c, x2, f2 = x2, x1, f1; x1 = c - gr * (c - a); f1 = fb(x1)
        else:
            a, x1, f1 = x1, x2, f2; x2 = a + gr * (c - a); f2 = fb(x2)
    cands = vals + [(f1[0], x1, f1[1]), (f2[0], x2, f2[1])]
    return min(cands)


def theta_best(lev, bins=None, verbose=False, **kw):
    """best theta; kappa dynamic handled by fixed-point iteration kappa = min(3/4, 2 theta - 1)."""
    if lev.kappa == 'dynamic':
        kappa = kw.pop('kappa0', 0.7498)
        for _ in range(8):
            th, b, ell = theta_fixed_kappa(lev, kappa, bins, **kw)
            k2 = kappa_of(th, lev)
            if verbose:
                print(f'  kappa={kappa:.6f} -> theta={th:.7f} (b={b:.4f}, ell={ell:.5f})')
            if abs(k2 - kappa) < 2e-6:
                kappa = k2
                break
            kappa = k2
        th, b, ell = theta_fixed_kappa(lev, kappa, bins, **kw)
        return dict(theta=th, b=b, ell=ell, kappa=kappa)
    kappa = float(lev.kappa)
    th, b, ell = theta_fixed_kappa(lev, kappa, bins, **kw)
    return dict(theta=th, b=b, ell=ell, kappa=kappa)


def verify(res, lev, bins=None):
    """re-evaluate a theta-optimum on the fine bins with exact L; returns margin B_low - B_high (>=0 good)."""
    bins = bins or BINS_FINE
    ell, b, kappa = res['ell'], res['b'], res['kappa']
    bh, det = B_high(ell, b, lev, kappa, bins, detail=True, exactL=True)
    return B_low(ell, b, lev) - bh, det


def closes(lev, bins=None, target=TARGET):
    """does theta < target?  kappa dynamic evaluated at the target (kappa = 2*target - 1)."""
    if lev.kappa == 'dynamic':
        kappa = min(0.75, 2 * target - 1)
    else:
        kappa = float(lev.kappa)
    th, b, ell = theta_fixed_kappa(lev, kappa, bins)
    return th, b, ell


BINS_MID = make_bins(0.0025, xs=(0, 0.2, 0.3, 0.4, 0.45, 0.5))


def best_over_b(ell, B, lev, kappa, bins=None, bs=np.linspace(0, 0.15, 31), M_of=None, refine=True):
    """agent-style criterion: min over b of the worst high exponent at fixed (ell, B)."""
    vals = [(worstE(ell, b, B, lev, kappa, bins), b) for b in bs]
    e, b0 = min(vals)
    if refine:
        st = bs[1] - bs[0]
        for b in np.linspace(max(0, b0 - st), b0 + st, 9):
            vals.append((worstE(ell, b, B, lev, kappa, bins), b))
        e, b0 = min(vals)
    return e, b0


def bisect_lever(fun, lo, hi, it=14):
    """fun(v) = worst exponent / margin, monotone in v, with opposite signs at lo and hi.
    Returns the threshold on the CLOSING side (fun < 0) of the final bracket."""
    flo = fun(lo); fhi = fun(hi)
    if (flo < 0) == (fhi < 0):
        return None
    for _ in range(it):
        mid = (lo + hi) / 2
        if (fun(mid) < 0) == (fhi < 0):
            hi = mid
        else:
            lo = mid
    return hi if fhi < 0 else lo


# ----------------------------------------------------------------------------------------------- target test
def ell_T(b, lev, T):
    """slot length at which B_low = T on the relaxed wall M = 1 + w - ell."""
    return 4 * (11 / 12 + lev.w / 12 - b * (1 - lev.g) / 12 - T)


def margin_T(lev, T, bins=None, kappa=None, bs=np.linspace(0, 0.2, 21), golden=16, detail=False):
    """min over b of [B_high(ell_T(b), b) - T]  (< 0  <=>  theta < T is proved in the model).
    B_high is increasing in ell (checked), so the smallest admissible ell = ell_T(b) is optimal."""
    if kappa is None:
        kappa = kappa_of(T, lev)

    def f(b):
        ell = ell_T(b, lev, T)
        if not side_ok(ell, b, lev, kappa=kappa):
            return 9.0
        return B_high(ell, b, lev, kappa, bins) - T
    vals = [(f(b), b) for b in bs]
    v0, b0 = min(vals)
    j = list(bs).index(b0)
    a = bs[max(j - 1, 0)]; c = bs[min(j + 1, len(bs) - 1)]
    gr = (np.sqrt(5) - 1) / 2
    x1 = c - gr * (c - a); x2 = a + gr * (c - a); f1 = f(x1); f2 = f(x2)
    for _ in range(golden):
        if f1 <= f2:
            c, x2, f2 = x2, x1, f1; x1 = c - gr * (c - a); f1 = f(x1)
        else:
            a, x1, f1 = x1, x2, f2; x2 = a + gr * (c - a); f2 = f(x2)
    best = min([(v0, b0), (f1, x1), (f2, x2)])
    if detail:
        b = best[1]; ell = ell_T(b, lev, T)
        bh, det = B_high(ell, b, lev, kappa, bins, detail=True)
        det.update(b=b, ell=ell, kappa=kappa, margin=best[0])
        return best[0], det
    return best[0]


def theta_star(lev, bins=None, lo=0.84, hi=0.879, tol=2e-8, verbose=False):
    """best theta: kappa fixed-point (kappa = min(3/4, 2 theta - 1)) around a bisection in T at fixed kappa."""
    def solve(kappa):
        a, c = lo, hi
        if margin_T(lev, c, bins, kappa) > 0:
            return None
        while c - a > tol:
            m = (a + c) / 2
            if margin_T(lev, m, bins, kappa) < 0:
                c = m
            else:
                a = m
        return c
    if lev.kappa != 'dynamic':
        return solve(float(lev.kappa))
    kappa = 0.75
    th = None
    for _ in range(6):
        th = solve(kappa)
        if th is None:
            return None
        k2 = min(0.75, 2 * th - 1)
        if verbose:
            print(f'   kappa {kappa:.6f} -> theta {th:.8f}')
        if abs(k2 - kappa) < 1e-6:
            break
        kappa = k2
    return th


def closes41(lev, bins=None):
    return margin_T(lev, TARGET, bins)
