"""
Idea 1 (Gram P_a^{1/6} split) in the R4 Part II LP model.

Low side Gram exponent g(d) is replaced by
    g = max(0, b/6 - 2 s(M'), 2b - l_y + d),     M' = M - 2d,
where s(M') = min((1-M')/2, 1 - 7M'/6, ...) is the power saving (in Z) of a direct
(non-Cauchy-Schwarz) bound for the theta-pole component <A^sig, B> over its CS bound.
  'LS'      : saving from Voronoi for theta_6 x quadratic char + Heath-Brown cubic LS
              s = -1/2 max(M'-1, 2(M'-1), 7M'/3 - 2)
  'random'  : s = (1-M')/2  (ideal large sieve; heuristic only)
  'none'    : s = 0 (paper)
Implemented by monkeypatching the low-side g pieces in model_r4.lp_for_S.
"""
import sys, json
import numpy as np
import model_r4 as M

def make_lp(mode):
    orig = M.lp_for_S
    def lp_for_S(S, cfg, return_geom=False, Rcache=None):
        # rebuild low side with modified gram pieces by temporarily replacing gram via cfg flags
        return orig(S, cfg, return_geom=return_geom, Rcache=Rcache)
    return lp_for_S

# --- direct re-implementation of the low side with theta-split pieces (copy of model_r4 logic) ---
from scipy.optimize import linprog

def lp_for_S(S, cfg, mode, return_geom=False, Rcache=None):
    s1 = S / (1 - S)
    A, b = [], []
    def ge_sigma(cx, cy, c0):
        A.append([cx, cy, -1.0]); b.append(-c0)
    ell = (-s1, 0.0, s1)
    h = (-1 / (1 - S), 0.0, 1 / (1 - S))
    def add(*fs, coef=None):
        cx = cy = c0 = 0.0
        for k, f in zip(coef, fs):
            cx += k * f[0]; cy += k * f[1]; c0 += k * f[2]
        return (cx, cy, c0)
    LX = (1.0, 0.0, 0.0); LY = (0.0, 1.0, 0.0); ONE = (0.0, 0.0, 1.0)
    Mx = add(LX, LY, coef=(1, 1)); B = add(LY, LX, coef=(1, -1))
    wE = cfg["wE"]
    for lam in (0.0, 1.0):
        d = add(ell, coef=(lam,))
        Mp = add(Mx, ell, coef=(1, -2 * lam)); lp = add(ell, coef=(1 - lam,))
        base = add(ONE, LX, ell, d, coef=(5 / 6, 1 / 6, -1 / 6, -1))
        gpieces = [(0.0, 0.0, 0.0)]
        # P^2/Y' piece (nonexceptional frequencies), unchanged
        if not cfg.get("dropq"): gpieces.append(add(B, LY, d, coef=(2, -1, 1)))
        # theta piece: b/6 - 2 s(M')
        if mode == "none":
            gpieces.append(add(B, coef=(1 / 6,)))
        elif mode == "LS":
            # -2s = max(M'-1, 2(M'-1), 7M'/3 - 2)
            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 1, -1)))
            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 2, -2)))
            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 7 / 3, -2)))
        elif mode == "LSmarked":
            # reflected length with marked slots: L = Q^2 Z^{l'-1};  -2s = max(M'+l'-1, 7M'/3+2l'-2, -M'/6)
            gpieces.append(add(B, Mp, lp, ONE, coef=(1 / 6, 1, 1, -1)))
            gpieces.append(add(B, Mp, lp, ONE, coef=(1 / 6, 7 / 3, 2, -2)))
            gpieces.append(add(B, Mp, coef=(1 / 6, -1 / 6)))
        elif mode == "randommarked":
            gpieces.append(add(B, Mp, lp, ONE, coef=(1 / 6, 1, 1, -1)))
            gpieces.append(add(B, Mp, coef=(1 / 6, -1 / 6)))
        elif mode == "random":
            gpieces.append(add(B, Mp, ONE, coef=(1 / 6, 1, -1)))
        epieces = [(0.0, 0.0, 0.0)]
        if cfg["caseA"] == "paper": epieces.append(add(ONE, lp, Mp, coef=(1 / 4, 3 / 4, -1 / 2)))
        elif cfg["caseA"] == "refined": epieces.append(add(ONE, lp, Mp, coef=(1 / 3, 2 / 3, -2 / 3)))
        epieces.append(add(Mp, lp, ONE, coef=(1, 1, -1 - wE)))
        epieces.append(add(Mp, lp, ONE, coef=(1 / 3, 4 / 3, -2 / 3 - wE)))
        for g in gpieces:
            for e in epieces:
                ge_sigma(*add(base, g, e, coef=(1, 0.5, 0.5)))
    d0 = cfg["delta0"]
    deltas = np.linspace(d0, cfg["dmax"], cfg["dgrid"]); xs = np.linspace(0, 0.5, cfg["xgrid"])
    DD, XX = np.meshgrid(deltas, xs, indexing="ij")
    key = (S,)
    if Rcache is not None and key in Rcache:
        RR = Rcache[key]
    else:
        RR = M.R_of(DD, XX, S, cfg)
        if Rcache is not None: Rcache[key] = RR
    for dl, xx, R in zip(DD.ravel(), XX.ravel(), RR.ravel()):
        ge_sigma(*add(ONE, LY, ell, h, coef=((1 + dl) / 2, -(1 + dl) / 2, -0.5 + xx * dl, R + dl / 2 - 1 / 6)))
    ge_sigma(*add(ONE, LY, ell, h, coef=((1 + d0) / 2, -(1 + d0) / 2, -0.5 + d0 / 2, 1 + d0 / 2 - 1 / 6)))
    A.append([-1 - 2 * s1, -1, 0]); b.append(-0.05 - 2 * s1)
    A.append([-s1, -1, 0]); b.append(-0.05 - s1)
    A.append([-1 - s1, 0, 0]); b.append(-s1)
    res = linprog([0, 0, 1], A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, 1), (0, 1.5), (None, None)], method="highs")
    if res.status != 0:
        return (9.0, None) if return_geom else 9.0
    if return_geom:
        lx, ly, sg = res.x
        return sg, dict(lx=lx, ly=ly, ell=s1 * (1 - lx), h=(1 - lx) / (1 - S), S=S, M=lx + ly, b=ly - lx)
    return res.fun

def theta(mode, over=None):
    from scipy.optimize import minimize_scalar
    cfg = dict(M.DEFAULT); cfg.update(over or {})
    cache = {}
    Sgrid = np.linspace(0.05, 0.6, 45)
    vals = [lp_for_S(S, cfg, mode, Rcache=cache) for S in Sgrid]
    k = int(np.argmin(vals))
    a = Sgrid[max(k - 1, 0)]; c = Sgrid[min(k + 1, len(Sgrid) - 1)]
    res = minimize_scalar(lambda S: lp_for_S(S, cfg, mode), bounds=(a, c), method="bounded", options=dict(xatol=1e-7))
    Sb = res.x if res.fun < vals[k] else Sgrid[k]
    return lp_for_S(Sb, cfg, mode, return_geom=True)

if __name__ == "__main__":
    for mode in ["none", "LS", "LSmarked", "random", "randommarked"]:
        for over in [{}, {"caseA": "refined"}]:
            sg, g = theta(mode, over)
            print(mode, over, round(sg, 6), {k: round(v, 4) for k, v in g.items()})
