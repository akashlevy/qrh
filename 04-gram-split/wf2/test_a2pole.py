"""A2-residue fingerprint: if sum_m rho(m) gamma_j(m)^e alpha(m)^k psi(m) N(m)^{-z} (on squarefree m) had a pole beyond the
random-level exponent, the class-resolved untwisted sums would show a coherent character psi mod 36 reproducible across Q.
x = rho * gamma_j(R)^e, j in {1,2,3}, e = +-1, k in [-6,6]; cross-Q correlation of the 648-profile and #psi with Rat>8 at both Q.
(gamma_3 is a ray function: j=3 must reproduce rho's own k=0 coherence -> control.)"""
import numpy as np
from a2lib import *
exec(open("test_chars36.py").read().split("cases = {")[0])
out = {}
for Q in (3000, 6000):
    R = Rows(Q); W = R.w; r = np.where(R.sf)[0]; Om = R.Om[r]
    lab = np.array([idx[(a % f, b % f)] for a, b in zip(R.ma[r], R.mb[r])])
    d = cross_est(R.V[:, r].T, W)
    for j in (1, 2, 3):
        for e in (-1, 1):
            X = R.V[:, r] * (R.G[j][r] ** e)[None, :]
            for k in range(-6, 7):
                Y = X * (Om * R.alpha_m[r] ** k)[None, :]
                S = np.stack([np.bincount(lab, y.real, nG) + 1j * np.bincount(lab, y.imag, nG) for y in Y], 1)
                out[(j, e, k, Q)] = cross_est(CH.T @ S, W) / (Om ** 2 * d).sum()
for j in (1, 2, 3):
    for e in (-1, 1):
        line = []
        for k in range(-6, 7):
            a, b = out[(j, e, k, 3000)], out[(j, e, k, 6000)]
            line.append("%+d:%.2f/%d" % (k, np.corrcoef(a, b)[0, 1], int(((a > 8) & (b > 8)).sum())))
        print("x = rho*gamma_%d^%+d : " % (j, e) + " ".join(line))
