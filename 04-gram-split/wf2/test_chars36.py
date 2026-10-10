"""Per-character coherence for characters psi of G = (O/36)^x (648 characters, built as joint eigenvectors of the
multiplication operators).  Rat(psi) = E|sum_m Om x h psi|^2 / sum Om^2 |x h|^2 (distinct-scale estimate).
For each case, list the top characters at each Q and test whether the same characters are large at Q = 1500, 3000, 6000.
Control: x = rho, h = 1 (theta constant term) should give reproducible top characters."""
import numpy as np
from a2lib import *
f = 36
els = [(a, b) for a in range(f) for b in range(f) if norm(a, b) % 2 and norm(a, b) % 3]
idx = {e: i for i, e in enumerate(els)}; nG = len(els)
def mulmod(x, y): a, b = mul(x[0], x[1], y[0], y[1]); return (a % f, b % f)
rng = np.random.default_rng(0)
A = np.zeros((nG, nG), complex)
for g in [els[i] for i in rng.choice(nG, 12, replace=False)] + [(1, 1), (2, 1), (1, 3), (5, 0), (7, 0), (1, 6), (0, 1)]:
    if g not in idx: continue
    Pg = np.zeros((nG, nG))
    for i, e in enumerate(els): Pg[idx[mulmod(e, g)], i] = 1
    A += (rng.normal() + 1j * rng.normal()) * Pg
ev, Vc = np.linalg.eig(A)
CH = Vc / Vc[idx[(1, 0)]][None, :]           # characters: CH[class, psi], value 1 at identity
CH = CH / np.abs(CH)
assert np.allclose(np.abs(CH.conj().T @ CH) / nG, np.eye(nG), atol=1e-6), "characters not orthogonal"
print("built %d characters of (O/36)^x" % nG)
cases = {"rho*1": lambda R, r: (np.ones(len(r)), np.ones(len(r))),
         "pi*1": lambda R, r: (R.zetaZ[r], np.ones(len(r))),
         "pi*mu*alpha^5": lambda R, r: (R.zetaZ[r], (-1.0) ** R.nfac[r] * R.alpha_m[r] ** 5),
         "pi*alpha^4": lambda R, r: (R.zetaZ[r], R.alpha_m[r] ** 4),
         "pi*alpha^-2": lambda R, r: (R.zetaZ[r], R.alpha_m[r] ** -2)}
res = {}
for Q in (1500, 3000, 6000):
    R = Rows(Q); W = R.w; K = len(R.scales); r = np.where(R.sf)[0]; Om = R.Om[r]
    lab = np.array([idx[(a % f, b % f)] for a, b in zip(R.ma[r], R.mb[r])])
    for name, fn in cases.items():
        zf, h = fn(R, r)
        X = R.V[:, r] * zf[None, :]; D = (Om ** 2 * np.abs(h) ** 2 * cross_est(X.T, W)).sum()
        Y = X * (Om * h)[None, :]
        S = np.stack([np.bincount(lab, y.real, nG) + 1j * np.bincount(lab, y.imag, nG) for y in Y], 1)
        C = CH.T @ S                                       # (psi, K): sum_m x h psi(m)
        res[(name, Q)] = cross_est(C, W) / D
for name in cases:
    r1, r3, r6 = res[(name, 1500)], res[(name, 3000)], res[(name, 6000)]
    print("\n%s: mean Rat over psi: %.2f / %.2f / %.2f ; max: %.1f / %.1f / %.1f" % (name, r1.mean(), r3.mean(), r6.mean(), r1.max(), r3.max(), r6.max()))
    top6 = np.argsort(-r6)[:5]
    print("  top-5 characters at Q=6000 and their Rat at Q=1500, 3000:", "; ".join("%.1f [%.1f, %.1f]" % (r6[i], r1[i], r3[i]) for i in top6))
    top3 = np.argsort(-r3)[:5]
    print("  top-5 characters at Q=3000 and their Rat at Q=1500, 6000:", "; ".join("%.1f [%.1f, %.1f]" % (r3[i], r1[i], r6[i]) for i in top3))
    cc = np.corrcoef(np.stack([r1, r3, r6]))
    print("  correlation of the 648-vectors Rat(psi) between Q (1500-3000, 1500-6000, 3000-6000): %.2f %.2f %.2f" % (cc[0, 1], cc[0, 2], cc[1, 2]))
    both = np.where((r3 > 8) & (r6 > 8))[0]
    print("  characters with Rat > 8 at both Q=3000 and 6000: %d  (expected by chance if Rat ~ Exp(1): %.2f)" % (len(both), nG * np.exp(-8) ** 2))
