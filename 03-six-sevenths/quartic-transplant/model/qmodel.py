"""
Parametrized exponent model of Part II of the 7/8 paper (sextic, Q(sqrt-3)) and of its quartic transplant (Q(i)).

theta(cfg) = min over geometry (l_x, l_y, ell) of sigma0 subject to (every line = one lemma or its analogue):

LOW  (Prop probe-low; Lemma probe-row-norm/reflected-energy; Prop probe-gram), for d in [0, ell]:
     sigma0 >= 1 - h/n - d + g(d)/2 + e(M - 2d, ell - d)/2,     h = 1 - l_x + ell,  M = l_x + l_y,  b = l_y - l_x
       1/n  : pole of zeta_F(nz) from the Poisson rows u a^n  (C(s) = s + l_x/2 - 1 + h/n),
       -d   : rescaled slots (q_p^{-3/2} compensation weight, Z^d tuples, X'^{1/2}),
       g(d) = max(0, gb*b/n, 2b - l_y + d)  (Gram: diagonal, exceptional frequencies k = a^n e, generic),
       e    : reflected-energy loss (cubic theta / quartic Kubota-Diaconu geometry + large sieves).
HIGH (Lemma high-bin, eq. common-high-exponent at d = h, Prop detector-counts):
     sigma0 >= (1+delta)(1-l_y)/2 - ell/2 + x delta ell + h (R(delta,x) + delta/2 - 1/n),
     delta in [delta0, dmax], x in [0, 1/2]
FLOOR: R = 1 at delta0, q <= delta0/2.
R(delta,x) = min_{t in [1,3/2]} max_{r in [t-1/2,t]} min{ counts }, m = t - r, q = x delta, supply S = ell/h:
     fA  = eamp(r) [min edir(r)] - delta r                       (n-th power amplification [+ direct sieve])
     fI  = 1 - delta r - 2q zI(r),  r < rI                        (marked inverse moment)
     fP  = 1 - 2 delta m - 2q min(cP (1-2m), S), m < 1/2          (plain fourth moment, two copies)
     fP0 = 1 - 2 delta m                                          (plain, no slots)
     fDI = min_z edir(r+z) - delta r - 2qz;  fDP = min_z edir(2m+z) - 2 delta m - 2qz   (direct de Faveri)
"""
import numpy as np, itertools, json, sys
from scipy.optimize import minimize
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from energy_closed import e_cubic_paper, e_quartic

# ----------------------------------------------------------------------------- configurations
DF4 = [(1, 0), (0, 1), (0.75, 0.5), (0.5, 0.75)]       # de Faveri n=4: A + B + A^{3/4}B^{1/2} + A^{1/2}B^{3/4}
DF6 = [(1, 0), (0, 1), (5 / 6, 1 / 3), (1 / 3, 5 / 6)]  # de Faveri n=6
PERF = [(1, 0), (0, 1)]

def cfg_sextic_paper(tmin=1.0, tmax=1.5, dirs=None):
    return dict(name="sextic paper (calibration)", n=6, energy=lambda M, l: e_cubic_paper(M, l), gb=1.0,
                amp=[(1.0, 0.0), (1 / 6, 5 / 6)], dirs_inv=dirs, dirs_plain=dirs,
                zI=[[(0.5, -0.5), (3 / 8, -0.25)]], rI=1.0, cP=2 / 9, kappa_dyn=False,
                delta0=1 / 50, dmax=5 / 6, slots_dir=dirs is not None, Z0small=17 / 50, tmin=tmin, tmax=tmax)

def cfg_quartic(scenario="A", kappa=None, dmax=None, tmin=0.3, tmax=3.0):
    """Quartic transplant (Q(i)) with the DERIVED Kubota geometry.
    Lemma-marked capacity = union over the two terminal evaluations of the quartic canonical family
    (rows M = 2(r+z)-1, columns N = r+z, marks z0 = z, in base U):
       direct (unreflected) sieve on columns N,   reflected: dual length T_d = 3M + 3z0 - N, rows sieve x quad. marks."""
    if scenario == "A":            # unconditional: de Faveri n=4 everywhere
        zI = [[(2 / 3, -1.0)],                                           # direct dF: M <= N/2  <=> r + z <= 2/3
              [(1.0, -1.0), (3 / 7, -4 / 7), (11 / 24, -15 / 24)]]       # reflected dF: M<N, 3M+3z<2N, 11M+9z<7N
        rI = 11 / 15; energy = "dF"; dinv = DF4; dpl = DF4
    elif scenario == "C":          # Lindelof on average (cube-shift aspect) for the Kubota partial sums
        zI = [[(1.0, -1.0)]]; rI = 1.0; energy = "LoA"; dinv = PERF; dpl = DF4
    elif scenario == "D":          # hypothetical perfect quartic large sieve A+B for all vectors (false by dFDH)
        zI = [[(1.0, -1.0)]]; rI = 1.0; energy = "perfect"; dinv = PERF; dpl = PERF
    elif scenario == "BGL":        # BGL (AB)^{2/3} instead of de Faveri (older proved input)
        zI = [[(2 / 3 * 0 + 0.6, -1.0)],                                 # direct BGL: 2(M+N)/3 <= N <=> M <= N/2 -> r+z<=2/3
              [(1.0, -1.0), (3 / 7, -4 / 7), (8 / 17, -11 / 17)]]
        zI[0] = [(2 / 3, -1.0)]
        rI = 8 / 11; energy = "BGL"; dinv = [(1, 0), (0, 1), (2 / 3, 2 / 3)]; dpl = dinv
    else:
        raise ValueError(scenario)
    k1 = 1.0 / rI
    amp = [(1.0, 0.0), (0.25, 0.75 * k1)]                    # e(r) = max(1, (1 + 3 k1 r)/4): u -> u a^4 amplification
    return dict(name=f"quartic derived ({scenario})", n=4,
                energy=(lambda M, l, s=energy: e_quartic(M, l, s, True)), gb=1.0,
                amp=amp, dirs_inv=dinv, dirs_plain=dpl, zI=zI, rI=rI, cP=None, kappa_dyn=True, kappa=kappa,
                delta0=1 / 50, dmax=dmax, slots_dir=True, Z0small=0.26, tmin=tmin, tmax=tmax)

def cfg_theta_route(kappa=None, dmax=None, tmin=0.3, tmax=3.0):
    """R5's ASSUMED quartic 'theta route' (conductor N(h)^2, quadratic dual twist), with every other lemma transplanted
    6 -> 4 (n=4 poles/Gram/amplification, plain capacity 1/(4 kappa)).  NOT the derived geometry; reference only."""
    return dict(name="hypothetical quartic theta route (N^2, quadratic dual)", n=4,
                energy=lambda M, l: e_cubic_paper(M, l), gb=1.0,
                amp=[(1.0, 0.0), (0.25, 0.75)], dirs_inv=DF4, dirs_plain=DF4,
                zI=[[(0.5, -0.5), (3 / 8, -0.25)]], rI=1.0, cP=None, kappa_dyn=True, kappa=kappa,
                delta0=1 / 50, dmax=dmax, slots_dir=True, Z0small=0.26, tmin=tmin, tmax=tmax)

# ----------------------------------------------------------------------------- exact row count
def _breaks(pcs):
    out = set()
    for (a0, a1), (b0, b1) in itertools.combinations(pcs, 2):
        if a1 != b1:
            r = (b0 - a0) / (a1 - b1)
            if r > 0 and abs(a0 + a1 * r - max(c0 + c1 * r for c0, c1 in pcs)) < 1e-12:
                out.add(r)
    return sorted(out)

def _emax(pcs, r):
    return np.max(np.stack([c0 + c1 * r for c0, c1 in pcs]), axis=0)

def _lines(cfg, d, q, t, S, cP):
    """all affine pieces (a + b r) of all count functions (candidate generation for the exact max-min)."""
    L = []
    one = np.ones_like(d)
    for c0, c1 in cfg["amp"]:
        L.append((c0 * one, c1 - d))
    L.append((one, -d))
    for alt in cfg["zI"]:
        for c0, c1 in alt:
            L.append((1 - 2 * q * c0, -d - 2 * q * c1))
    L.append((1 - 2 * q * S, -d))
    L.append((1 - 2 * d * t, 2 * d))                                           # plain, m = t - r
    L.append((1 - 2 * d * t - 2 * q * cP * (1 - 2 * t), 2 * d - 4 * q * cP))
    L.append((1 - 2 * d * t - 2 * q * S, 2 * d))
    if cfg["dirs_inv"]:
        P = cfg["dirs_inv"]
        for c0, c1 in P:
            L.append((c0 * one, c1 - d))                                       # fDI z=0 / direct inverse
            L.append((c0 + c1 * S - 2 * q * S, c1 - d))                        # fDI z=S
        for bp in _breaks(P):
            L.append((_emax(P, bp) - 2 * q * bp, -d + 2 * q))                  # fDI r+z=bp
    if cfg["dirs_plain"]:
        P = cfg["dirs_plain"]
        for c0, c1 in P:
            L.append((c0 + c1 * 2 * t - 2 * d * t, -2 * c1 + 2 * d))           # fDP z=0
            L.append((c0 + c1 * (2 * t + S) - 2 * d * t - 2 * q * S, -2 * c1 + 2 * d))
        for bp in _breaks(P):
            L.append((_emax(P, bp) - 2 * d * t - 2 * q * (bp - 2 * t), 2 * d - 4 * q))
    return L

def _zI(cfg, r, S):
    alts = [np.min(np.stack([c0 + c1 * r for c0, c1 in alt] + [S + 0 * r]), axis=0) for alt in cfg["zI"]]
    return np.clip(np.max(np.stack(alts), axis=0), 0, None)

def _F(cfg, r, d, q, t, S, cP):
    """min of all count exponents at inverse length r, plain length m = t - r (arrays broadcast)."""
    m = t - r
    eamp = _emax(cfg["amp"], r)
    if cfg["dirs_inv"]:
        eamp = np.minimum(eamp, _emax(cfg["dirs_inv"], r))
    fA = eamp - d * r
    fI = np.where(r < cfg["rI"] - 1e-12, 1 - d * r - 2 * q * _zI(cfg, r, S), np.inf)
    zP = np.minimum(cP * np.clip(1 - 2 * m, 0, None), S)
    fP = np.where(m < 0.5, 1 - 2 * d * m - 2 * q * zP, np.inf)
    fP0 = 1 - 2 * d * m
    F = np.minimum(np.minimum(fA, fI), np.minimum(fP, fP0))
    if cfg["slots_dir"]:
        if cfg["dirs_inv"]:
            P = cfg["dirs_inv"]
            c = [_emax(P, r + S) - d * r - 2 * q * S]
            for bp in _breaks(P):
                z = bp - r
                c.append(np.where((z >= 0) & (z <= S), _emax(P, bp + 0 * r) - d * r - 2 * q * z, np.inf))
            F = np.minimum(F, np.min(np.stack(c), axis=0))
        if cfg["dirs_plain"]:
            P = cfg["dirs_plain"]
            mm = np.clip(m, 0, None)
            c = [_emax(P, 2 * mm) - 2 * d * mm, _emax(P, 2 * mm + S) - 2 * d * mm - 2 * q * S]
            for bp in _breaks(P):
                z2 = bp - 2 * mm
                c.append(np.where((z2 >= 0) & (z2 <= S), _emax(P, bp + 0 * r) - 2 * d * mm - 2 * q * z2, np.inf))
            F = np.minimum(F, np.min(np.stack(c), axis=0))
    return F

def R_t(cfg, d, q, t, S, cP):
    t = np.asarray(t, float) + 0 * d
    lo = np.maximum(t - 0.5, 0.0); hi = t
    L = _lines(cfg, d, q, t, S, cP)
    cands = [lo, hi, cfg["rI"] + 0 * d, (cfg["rI"] - 1e-9) + 0 * d]
    for alt in cfg["zI"]:
        for c0, c1 in alt:
            if c1 != 0:
                cands.append(-c0 / c1 + 0 * d); cands.append((S - c0) / c1 + 0 * d)
    for P in (cfg["dirs_inv"], cfg["dirs_plain"]):
        if P:
            for bp in _breaks(P):
                for v in (bp, bp - S, t - bp / 2, t - (bp - S) / 2):
                    cands.append(v + 0 * d)
    cands.append(t - (1 - S / max(cP, 1e-9)) / 2 + 0 * d)
    for (a1, b1), (a2, b2) in itertools.combinations(L, 2):
        db = b1 - b2
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(np.abs(db) > 1e-13, (a2 - a1) / np.where(np.abs(db) > 1e-13, db, 1.0), lo)
        cands.append(r)
    C = np.clip(np.stack(cands, axis=-1), lo[..., None], hi[..., None])
    vals = _F(cfg, C, d[..., None], q[..., None], t[..., None], S, cP)
    return vals.max(axis=-1)

def R_of(cfg, d, x, S, cP, nt=41, refine=18):
    d = np.asarray(d, float); x = np.asarray(x, float); q = x * d
    tmin, tmax = cfg.get("tmin", 1.0), cfg.get("tmax", 1.5)
    ts = np.linspace(tmin, tmax, nt)
    vals = np.stack([R_t(cfg, d, q, t, S, cP) for t in ts], axis=-1)
    k = vals.argmin(axis=-1); best = vals.min(axis=-1); tb = ts[k]
    dt = (tmax - tmin) / (nt - 1)
    for _ in range(refine):
        dt *= 0.6
        for s in (-1, 1):
            tt = np.clip(tb + s * dt, tmin, tmax)
            v = R_t(cfg, d, q, tt, S, cP)
            better = v < best
            best = np.where(better, v, best); tb = np.where(better, tt, tb)
    return np.clip(best, 0.0, None)

# ----------------------------------------------------------------------------- theta
class Model:
    def __init__(self, cfg, kappa=None, dgrid=81, xgrid=11, Sgrid=None, nt=41, refine=10):
        self.cfg = cfg; self.n = cfg["n"]
        if cfg["cP"] is not None:
            self.cP = cfg["cP"]
        else:
            self.cP = 1.0 / (4.0 * kappa)            # quartic plain capacity 1/(4 kappa)  [sextic: 1/(6 kappa)]
        dmax = cfg["dmax"] if cfg["dmax"] is not None else (kappa if kappa is not None else 1.0)
        self.deltas = np.linspace(cfg["delta0"], dmax, dgrid)
        self.xs = np.linspace(0, 0.5, xgrid)
        self.Sgrid = np.linspace(0.0, 0.6, 61) if Sgrid is None else Sgrid
        DD, XX = np.meshgrid(self.deltas, self.xs, indexing="ij")
        self.DD, self.XX = DD, XX
        self.Rtab = np.stack([R_of(cfg, DD, XX, S, self.cP, nt=nt, refine=refine) for S in self.Sgrid])

    def R_at(self, S, exact=False):
        if exact:
            return R_of(self.cfg, self.DD, self.XX, S, self.cP, nt=81, refine=30)
        S = min(max(S, self.Sgrid[0]), self.Sgrid[-1])
        k = np.searchsorted(self.Sgrid, S)
        if k == 0: return self.Rtab[0]
        if k >= len(self.Sgrid): return self.Rtab[-1]
        s0, s1 = self.Sgrid[k - 1], self.Sgrid[k]
        w = (S - s0) / (s1 - s0)
        return self.Rtab[k - 1] * (1 - w) + self.Rtab[k] * w

    def parts(self, p, nd=9, exact=False):
        lx, ly, ell = p; n = self.n; cfg = self.cfg
        h = 1 - lx + ell; b = ly - lx; M = lx + ly
        low = -9
        for d in np.linspace(0, ell, nd):
            g = max(0.0, cfg["gb"] * b / n, 2 * b - ly + d) if b > 0 else 0.0
            e = cfg["energy"](M - 2 * d, ell - d)
            low = max(low, 1 - h / n - d + g / 2 + e / 2)
        S = ell / h
        R = self.R_at(S, exact)
        hi = ((1 + self.DD) * (1 - ly) / 2 - ell / 2 + self.XX * self.DD * ell + h * (R + self.DD / 2 - 1 / n)).max()
        d0 = cfg["delta0"]
        fl = (1 + d0) * (1 - ly) / 2 - ell / 2 + d0 * ell / 2 + h * (1 + d0 / 2 - 1 / n)
        small = h * (cfg["Z0small"] - 1 / n) - ly / 2 + 0.0126
        return low, hi, fl, small

    def theta_at(self, p):
        lx, ly, ell = p
        mg = self.cfg.get("margin", 0.05)
        if min(lx, ly, ell) < 0 or ell > lx + 1e-12 or ly - ell < mg or lx + ly - 2 * ell < mg or lx > 1.2 or ly > 1.2:
            return 9.0
        if ell > self.cfg.get("ell_frac", 1.0) * lx + 1e-12:
            return 9.0
        low, hi, fl, small = self.parts(p)
        return max(low, hi, fl) + 10 * max(0.0, small)

    def optimize(self, starts=None, seed=0):
        rng = np.random.default_rng(seed)
        if starts is None:
            starts = [(0.355, 0.478, 0.167), (0.3, 0.4, 0.0), (0.25, 0.35, 0.05), (0.2, 0.3, 0.1), (0.3, 0.3, 0.1),
                      (0.15, 0.3, 0.02), (0.35, 0.25, 0.05), (0.4, 0.45, 0.15), (0.1, 0.2, 0.0)]
            starts += [tuple(rng.uniform([0.05, 0.05, 0.0], [0.6, 0.6, 0.25])) for _ in range(12)]
        best = None
        for s in starts:
            r = minimize(self.theta_at, s, method="Nelder-Mead",
                         options=dict(xatol=1e-7, fatol=1e-9, maxiter=4000, maxfev=8000))
            if best is None or r.fun < best.fun:
                best = r
        # polish
        r = minimize(self.theta_at, best.x, method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-11, maxiter=8000))
        if r.fun < best.fun: best = r
        p = best.x
        low, hi, fl, small = self.parts(p, nd=41, exact=True)
        return dict(theta=float(max(low, hi, fl)), theta_interp=float(best.fun), lx=float(p[0]), ly=float(p[1]), ell=float(p[2]), h=float(1 - p[0] + p[2]),
                    M=float(p[0] + p[1]), low=float(low), high=float(hi), floor=float(fl), small=float(small))
