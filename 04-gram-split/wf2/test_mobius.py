"""Follow-up of a weak replicated cell in test_corr: rho vs conj(gamma_2) alpha^4 (f=12).  Since conj(zeta) conj(gamma_2) alpha^4
 = conj(psi_1) gamma_2^3 alpha^4 = conj(psi_1) mu alpha^5 (squarefree), this is the same as pi vs mu(m) alpha(m)^5 x ray.
k-profiles of Rat for x=pi, h = mu(R) alpha(m)^k [and x = rho, h = conj(gamma_2) alpha^k], f in {3,12,36}, three Q, with null z."""
import numpy as np
from a2lib import *
rng = np.random.default_rng(7)
fs = [3, 12, 36]; ks = list(range(-6, 7)); out = {}
for Q in (1500, 3000, 6000):
    R = Rows(Q); W = R.w; rows = np.where(R.sf)[0]; Om = R.Om[rows]; n = len(rows)
    mu = (-1.0) ** R.nfac[rows]
    labs = {f: np.unique((R.ma[rows] % f) * f + R.mb[rows] % f, return_inverse=True)[1] for f in fs}
    for xname, zf, hb in (("pi*mu", R.zetaZ[rows], mu), ("pi*mu*psi1", R.zetaZ[rows], mu * R.psi1[rows]),
                          ("rho*conj(g2)", np.ones(n), np.conj(R.G[2][rows])), ("rho*g2", np.ones(n), R.G[2][rows])):
        X = R.V[:, rows] * zf[None, :]; D = cross_est(X.T, W)
        def rat(h, lab):
            nc = lab.max() + 1; Y = X * (Om * h)[None, :]
            T = np.stack([np.bincount(lab, y.real, nc) + 1j * np.bincount(lab, y.imag, nc) for y in Y], 1)
            return cross_est(T, W).sum() / (Om ** 2 * np.abs(h) ** 2 * D).sum()
        null = {f: np.array([rat(rng.choice([-1., 1.], n), labs[f]) for _ in range(40)]) for f in fs}
        for k in ks:
            for f in fs:
                v = rat(hb * R.alpha_m[rows] ** k, labs[f]); out[(xname, k, f, Q)] = (v, (v - null[f].mean()) / null[f].std())
for xname in ("pi*mu", "pi*mu*psi1", "rho*conj(g2)", "rho*g2"):
    print("x*h = %s * alpha^k :   Rat (z) at Q=1500 | 3000 | 6000" % xname)
    for f in fs:
        print("  f=%d" % f)
        for k in ks:
            print("    k=%+d  " % k + " | ".join("%6.2f (%5.1f)" % out[(xname, k, f, Q)] for Q in (1500, 3000, 6000)))
