"""Item 2 (cont.): first-order correlation of rho and pi with Gauss-sum sequences h = gamma_j(R)^e alpha(m)^k,
j in {1,2,3}, e in {-1,+1}, |k|<=6, classes mod f in {3,12}; squarefree rows and prime rows; Q = 1500, 3000, 6000.
Reports cells with z > 3 vs a +-1-twist null and checks replication across Q (z > 2 at all three)."""
import numpy as np, itertools
from a2lib import *
rng = np.random.default_rng(3)
fs = [3, 12]; ks = range(-6, 7); allz = {}
for Q in (1500, 3000, 6000):
    R = Rows(Q); W = R.w; K = len(R.scales)
    for rname, rows in (("sf", np.where(R.sf)[0]), ("prime", np.where(R.sf & (R.nfac == 1))[0])):
        Om = R.Om[rows]; n = len(rows)
        labs = {f: np.unique((R.ma[rows] % f) * f + R.mb[rows] % f, return_inverse=True)[1] for f in fs}
        for xname, zf in (("rho", np.ones(n)), ("pi", R.zetaZ[rows])):
            X = R.V[:, rows] * zf[None, :]; D = cross_est(X.T, W)
            def rat(h, lab):
                nc = lab.max() + 1; Y = X * (Om * h)[None, :]
                T = np.stack([np.bincount(lab, y.real, nc) + 1j * np.bincount(lab, y.imag, nc) for y in Y], 1)
                return cross_est(T, W).sum() / (Om ** 2 * np.abs(h) ** 2 * D).sum()
            null = {f: np.array([rat(rng.choice([-1., 1.], n), labs[f]) for _ in range(40)]) for f in fs}
            for j in (1, 2, 3):
                for e in (-1, 1):
                    for k in ks:
                        h = R.G[j][rows] ** e * R.alpha_m[rows] ** k
                        for f in fs:
                            v = rat(h, labs[f]); allz[(rname, xname, j, e, k, f, Q)] = (v, (v - null[f].mean()) / null[f].std())
    print("Q=%d done" % Q, flush=True)
for rname in ("sf", "prime"):
    for xname in ("rho", "pi"):
        cells = [k[:6] for k in allz if k[0] == rname and k[1] == xname and k[6] == 6000]
        big = [(c, allz[c + (6000,)]) for c in cells if allz[c + (6000,)][1] > 3]
        rep = [c for c in cells if all(allz[c + (q,)][1] > 2 for q in (1500, 3000, 6000))]
        print("%-5s x=%-3s: %d cells; at Q=6000 z>3 in %d; replicated (z>2 at all Q): %s" % (rname, xname, len(cells), len(big),
              ["j=%d e=%+d k=%+d f=%d: %s" % (c[2], c[3], c[4], c[5], [round(allz[c + (q,)][0], 2) for q in (1500, 3000, 6000)]) for c in rep]))
        top = sorted(cells, key=lambda c: -allz[c + (6000,)][1])[:3]
        print("        top-3 z at Q=6000: " + "; ".join("j=%d e=%+d k=%+d f=%d Rat=%.2f z=%.1f" % (c[2], c[3], c[4], c[5], *allz[c + (6000,)]) for c in top))
