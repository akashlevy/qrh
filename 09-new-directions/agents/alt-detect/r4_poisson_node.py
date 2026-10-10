"""Sextic r4 model (paper counts) + the trivial 'Poisson node' for the theta-row energy:
   ||B||^2 via Poisson in m -> Mobius mean square -> BGL sextic large sieve:
   excess = max(0, (1+l') - 2M'/3)   [diag once (ZP')^3 <= Q'^2].
Purpose: detect whether r4/template exponents are being used outside their validity range
(huge Q = XY >> Z geometry)."""
import itertools, sys, json
import numpy as np
from scipy.optimize import linprog, minimize_scalar
import r4_model_copy as R4

def lp_for_S(S, cfg, box, Rcache):
    s1 = S / (1 - S)
    def add(*fs, coef=None):
        cx = cy = c0 = 0.0
        for k, f in zip(coef, fs):
            cx += k * f[0]; cy += k * f[1]; c0 += k * f[2]
        return (cx, cy, c0)
    LX = (1.0, 0.0, 0.0); LY = (0.0, 1.0, 0.0); ONE = (0.0, 0.0, 1.0)
    ell = (-s1, 0.0, s1); h = (-1 / (1 - S), 0.0, 1 / (1 - S))
    M = add(LX, LY, coef=(1, 1)); B = add(LY, LX, coef=(1, -1))
    d0 = cfg["delta0"]
    deltas = np.linspace(d0, cfg["dmax"], cfg["dgrid"]); xs = np.linspace(0, 0.5, cfg["xgrid"])
    DD, XX = np.meshgrid(deltas, xs, indexing="ij")
    if S not in Rcache: Rcache[S] = R4.R_of(DD, XX, S, cfg)
    RR = Rcache[S]
    high = [add(ONE, LY, ell, h, coef=((1 + dl) / 2, -(1 + dl) / 2, -0.5 + xx * dl, R + dl / 2 - 1 / 6))
            for dl, xx, R in zip(DD.ravel(), XX.ravel(), RR.ravel())]
    high.append(add(ONE, LY, ell, h, coef=((1 + d0) / 2, -(1 + d0) / 2, -0.5 + d0 / 2, 1 + d0 / 2 - 1 / 6)))
    best = (9.0, None)
    for choice in itertools.product((0, 1), repeat=2):
        A, b = [], []
        def ge(cx, cy, c0): A.append([cx, cy, -1.0]); b.append(-c0)
        for lam, ch in zip((0.0, 1.0), choice):
            d = add(ell, coef=(lam,))
            Mp = add(M, ell, coef=(1, -2 * lam)); lp = add(ell, coef=(1 - lam,))
            base = add(ONE, LX, ell, d, coef=(5 / 6, 1 / 6, -1 / 6, -1))
            gp = [(0.0, 0.0, 0.0), add(B, coef=(1 / 6,)), add(B, LY, d, coef=(2, -1, 1))]
            if ch == 0:   # paper: caseA + wall + third
                ep = [(0.0, 0.0, 0.0), add(ONE, lp, Mp, coef=(1 / 4, 3 / 4, -1 / 2)),
                      add(Mp, lp, ONE, coef=(1, 1, -1)), add(Mp, lp, ONE, coef=(1 / 3, 4 / 3, -2 / 3))]
            else:         # Poisson node
                ep = [(0.0, 0.0, 0.0), add(ONE, lp, Mp, coef=(1, 1, -2 / 3))]
            for g in gp:
                for e in ep:
                    ge(*add(base, g, e, coef=(1, 0.5, 0.5)))
        for f in high: ge(*f)
        A.append([-1 - 2 * s1, -1, 0]); b.append(-0.05 - 2 * s1)
        A.append([-s1, -1, 0]); b.append(-0.05 - s1)
        A.append([-1 - s1, 0, 0]); b.append(-s1)
        res = linprog([0, 0, 1], A_ub=np.array(A), b_ub=np.array(b), bounds=box, method="highs")
        if res.status == 0 and res.fun < best[0]:
            lx, ly, sg = res.x
            best = (sg, dict(lx=lx, ly=ly, ell=s1 * (1 - lx), h=(1 - lx) / (1 - S), S=S, nodes=choice))
    return best

if __name__ == "__main__":
    cfg = dict(R4.DEFAULT)
    for box, Smax in [([(0, 1), (0, 1.5), (None, None)], 0.6), ([(0, 2), (0, 3), (None, None)], 0.8),
                      ([(0, 4), (0, 6), (None, None)], 0.9)]:
        cache = {}
        Sg = np.linspace(0.05, Smax, 24)
        vals = [lp_for_S(S, cfg, box, cache) for S in Sg]
        k = int(np.argmin([v[0] for v in vals]))
        print("box", [bb[1] for bb in box[:2]], "Smax", Smax, "theta=%.6f" % vals[k][0],
              {kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in vals[k][1].items()}, flush=True)
