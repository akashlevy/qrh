"""Sweep of the reproducibility fingerprint: for x in {rho, pi}, h in {alpha^k, mu alpha^k}, |k| <= 6, rows in {sf, prime},
the 648-vector Rat(psi) over characters psi of (O/36)^x is computed at Q = 3000 and 6000; a genuine coherent direction
shows up as cross-Q correlation of the profiles and as characters large at both Q."""
import numpy as np
from a2lib import *
exec(open("test_chars36.py").read().split("cases = {")[0])         # reuse the character construction
out = {}
for Q in (3000, 6000):
    R = Rows(Q); W = R.w
    for rname, r in (("sf", np.where(R.sf)[0]), ("prime", np.where(R.sf & (R.nfac == 1))[0])):
        Om = R.Om[r]; lab = np.array([idx[(a % f, b % f)] for a, b in zip(R.ma[r], R.mb[r])]); mu = (-1.0) ** R.nfac[r]
        for xname, zf in (("rho", np.ones(len(r))), ("pi", R.zetaZ[r])):
            X = R.V[:, r] * zf[None, :]; d = cross_est(X.T, W)
            for hname, hb in (("alpha^k", np.ones(len(r))), ("mu*alpha^k", mu)):
                if rname == "prime" and hname == "mu*alpha^k": continue
                for k in range(-6, 7):
                    h = hb * R.alpha_m[r] ** k
                    Y = X * (Om * h)[None, :]
                    S = np.stack([np.bincount(lab, y.real, nG) + 1j * np.bincount(lab, y.imag, nG) for y in Y], 1)
                    out[(rname, xname, hname, k, Q)] = cross_est(CH.T @ S, W) / (Om ** 2 * d).sum()
print("%-6s %-4s %-11s  k : corr(Rat_3000, Rat_6000)  #psi with Rat>8 at both   max_6000" % ("rows", "x", "h"))
for rname in ("sf", "prime"):
    for xname in ("rho", "pi"):
        for hname in ("alpha^k", "mu*alpha^k"):
            if (rname, xname, hname, 0, 3000) not in out: continue
            line = []
            for k in range(-6, 7):
                a, b = out[(rname, xname, hname, k, 3000)], out[(rname, xname, hname, k, 6000)]
                line.append("%+d:%.2f/%d" % (k, np.corrcoef(a, b)[0, 1], int(((a > 8) & (b > 8)).sum())))
            print("%-6s %-4s %-11s  " % (rname, xname, hname) + " ".join(line))
