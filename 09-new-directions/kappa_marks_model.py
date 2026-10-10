"""
Exponent model for Part II of the 7/8 paper with GENERALISED MARK RESCALING.

Paper: marked term  eta-bar(p) I_{eta;p}(X, Y, Z p)        minus  p^{-3/2} I_eta(X/p, Y/p, Z)
Here:  marked term  eta-bar(p) I_{eta;p}(X, Y, Z p^kappa)  minus  p^{-3/2} I_eta(X/p, Y/p, Z p^(kappa-1))
(kappa = 1 is the paper).  The local multiplier becomes Q^{(kappa-1)(x+z-1)} * Q^{z-1} B_p, B_p as in the paper.

Derived consequences (see NOTES.md):
  signal      C_k(s) = l_x/3 + s - 5/6 + ell*(1 + (kappa-1)s - 5 kappa/6)       (slope D = 1+(kappa-1)ell)
  row range   h_k    = 1 + kappa*ell - l_x
  low side    D*sigma >= 5/6 + l_x/6 - ell + 5 kappa ell/6 - d + g(d)/2 + e(d)/2
              energy  e(d) = max(0, (1 + 3l' + (kappa-1)ell - 2M')/4 , M' + l' - 1 + (1-kappa) ell)
  high side   (1+delta)(1-l_y)/2 + h_k (R + delta/2 - 1/6) - ell/2 + x delta ell
              + (kappa-1) ell (1+delta)/2  <  D*sigma
  floor bin   same with delta0, R=1, x delta = delta0/2
Row counts R(delta,x,S) reused from ../03-six-sevenths/routes/r4/model.py with supply S = ell/h_k.
For fixed (S, sigma) every constraint is linear in (l_x, l_y): bisection on sigma.
"""
import sys, json
import numpy as np
from scipy.optimize import linprog, minimize_scalar
sys.path.insert(0, "/home/user/qrh/03-six-sevenths/routes/r4")
import model as M4

def feasible(S, sigma, kap, cfg, RR, DD, XX):
    # ell = S h, h = 1 + kap ell - lx  => ell = S(1-lx)/(1-kap S) ; affine in lx
    a1 = S / (1 - kap * S)                       # ell = a1 (1 - lx)
    ELL = (-a1, 0.0, a1)                         # (cx, cy, c0) affine forms in (lx, ly)
    LX = (1.0, 0.0, 0.0); LY = (0.0, 1.0, 0.0); ONE = (0.0, 0.0, 1.0)
    def add(*pairs):
        cx = cy = c0 = 0.0
        for k, f in pairs:
            cx += k * f[0]; cy += k * f[1]; c0 += k * f[2]
        return (cx, cy, c0)
    H = add((1, ONE), (kap, ELL), (-1, LX))
    Dslope = add((1, ONE), (kap - 1, ELL))       # D = 1 + (kap-1) ell
    A, b = [], []
    def le(f):          # f(lx,ly) <= 0
        A.append([f[0], f[1]]); b.append(-f[2])
    def ge_sigma(f):    # f <= sigma * D
        le(add((1, f), (-sigma, Dslope)))
    Mf = add((1, LX), (1, LY)); B = add((1, LY), (-1, LX))
    wE = cfg.get("wE", 0.0)
    for lam in (0.0, 1.0):
        d = add((lam, ELL))
        Mp = add((1, Mf), (-2 * lam, ELL)); lp = add((1 - lam, ELL))
        base = add((5/6, ONE), (1/6, LX), (-1 + 5 * kap / 6, ELL), (-1, d))
        gp = [(0, 0, 0), add((cfg.get("gram_b", 1.0) / 6, B))]
        if cfg.get("gram_q", 1.0) > 0:
            gp.append(add((2, B), (-1, LY), (1, d)))
        ep = [(0, 0, 0)]
        ep.append(add((1/4, ONE), (3/4, lp), ((kap - 1) / 4, ELL), (-1/2, Mp)))       # case A
        ep.append(add((1, Mp), (1, lp), (1 - kap, ELL), (-1 - wE, ONE)))              # wall
        for g in gp:
            for e in ep:
                ge_sigma(add((1, base), (0.5, g), (0.5, e)))
    for dl, xx, R in zip(DD.ravel(), XX.ravel(), RR.ravel()):
        f = add(((1 + dl) / 2, ONE), (-(1 + dl) / 2, LY), (R + dl / 2 - 1/6, H),
                (-0.5 + xx * dl + (kap - 1) * (1 + dl) / 2, ELL))
        ge_sigma(f)
    d0 = cfg["delta0"]
    f = add(((1 + d0) / 2, ONE), (-(1 + d0) / 2, LY), (1 + d0 / 2 - 1/6, H),
            (-0.5 + d0 / 2 + (kap - 1) * (1 + d0) / 2, ELL))
    ge_sigma(f)
    # side conditions (paper): M - 2 ell >= 0.05, l_y - ell >= 0.05, l_x >= ell
    le(add((-1, Mf), (2, ELL), (0.05, ONE)))
    le(add((-1, LY), (1, ELL), (0.05, ONE)))
    le(add((-1, LX), (1, ELL)))
    res = linprog([0, 0], A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, 1), (0, 1.5)], method="highs")
    return res.status == 0, (res.x if res.status == 0 else None)

def theta_for_S(S, kap, cfg, cache):
    if S not in cache:
        deltas = np.linspace(cfg["delta0"], cfg["dmax"], cfg["dgrid"]); xs = np.linspace(0, 0.5, cfg["xgrid"])
        DD, XX = np.meshgrid(deltas, xs, indexing="ij")
        RR = M4.R_of(DD, XX, S, cfg) if not cfg.get("perfect") else 1 - DD
        if cfg.get("grh"):
            DD = np.array([[0.0]]); XX = np.array([[0.0]]); RR = np.array([[1.0]])
        cache[S] = (RR, DD, XX)
    RR, DD, XX = cache[S]
    lo, hi = 0.80, 0.95
    ok, _ = feasible(S, hi, kap, cfg, RR, DD, XX)
    if not ok:
        return 9.0, None
    for _ in range(40):
        mid = (lo + hi) / 2
        ok, x = feasible(S, mid, kap, cfg, RR, DD, XX)
        if ok: hi = mid
        else: lo = mid
    ok, x = feasible(S, hi, kap, cfg, RR, DD, XX)
    return hi, x

def theta(kap, cfg_over=None, Sgrid=None):
    cfg = dict(M4.DEFAULT); cfg.update(cfg_over or {})
    cache = {}
    if Sgrid is None:
        Sgrid = np.linspace(0.02, 0.45, 44)
    vals = [theta_for_S(S, kap, cfg, cache)[0] for S in Sgrid]
    k = int(np.argmin(vals))
    a = Sgrid[max(k - 1, 0)]; c = Sgrid[min(k + 1, len(Sgrid) - 1)]
    res = minimize_scalar(lambda S: theta_for_S(S, kap, cfg, cache)[0], bounds=(a, c), method="bounded",
                          options=dict(xatol=1e-6))
    S = res.x if res.fun < vals[k] else Sgrid[k]
    sg, x = theta_for_S(S, kap, cfg, cache)
    geom = None
    if x is not None:
        lx, ly = x
        a1 = S / (1 - kap * S); ell = a1 * (1 - lx)
        geom = dict(lx=lx, ly=ly, ell=ell, h=1 + kap * ell - lx, S=S)
    return sg, geom

if __name__ == "__main__":
    kap = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    over = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    sg, geom = theta(kap, over)
    print(json.dumps(dict(kappa=kap, cfg=over, theta=sg, geom=geom)))
