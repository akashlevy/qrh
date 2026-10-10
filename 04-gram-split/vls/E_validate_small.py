"""Validate E_build at small size (KMAX = 13, s < 2^14, Q = 600) against a direct expA-style computation with zomega.
Also checks (zeta6/p)_6 = zeta6^{(Np-1)/6}."""
import numpy as np, sys
GS = "/Users/akashlevy/Documents/riemann/gram_split"
sys.path.insert(0, GS)
from zomega import *
X = 2 ** 14; Q = 600
P = primes_upto(X)
gt = np.array([Pr.gauss() / np.sqrt(Pr.N) for Pr in P])
# Euler check for units
bad = sum(int(Pr.symbol(np.array([1]), np.array([1]))[0]) != ((Pr.N - 1) // 6) % 6 for Pr in P)
print("primes with (zeta6/p)_6 != zeta6^{(Np-1)/6}:", bad, "of", len(P))
S = squarefree_products(P, X - 1)
ma, mb, mN = elements_in_annulus(Q / 2, Q)
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1.0 / (u * (1 - u)) + 4.0); return out
Om = bump(mN / Q, .5, 1); keep = Om > 1e-12; ma, mb, mN = ma[keep], mb[keep], mN[keep]
J = np.stack([Pr.symbol(-ma, -mb) for Pr in P])
K = 6; A = np.zeros((K, len(ma)), complex)
for (sa, sb, sN, idx) in S:
    if sN < 256: continue
    i = int(np.floor(np.log2(sN / 256)))
    if i >= K: continue
    rr = sN / (256 * 2 ** i); w = bump(np.array([rr]), 1, 2)[0] * rr ** -0.5 / (256 * 2 ** i)
    if w == 0: continue
    g = 1 + 0j
    for t, ii in enumerate(idx):
        g *= gt[ii]
        for jj in idx[t + 1:]:
            e1 = P[ii].symbol(np.array([P[jj].a]), np.array([P[jj].b]))[0]; e2 = P[jj].symbol(np.array([P[ii].a]), np.array([P[ii].b]))[0]
            g *= Z6[(e1 + e2) % 6]
    rows = J[list(idx)]; badr = (rows < 0).any(0); e = rows.sum(0) % 6
    v = np.conj(Z6[e]); v[badr] = 0
    A[i] += w * g * v
A2 = np.load(GS + "/vls/data/A_Q%d_E13.npy" % Q); m2 = np.load(GS + "/vls/data/m_Q%d_E13.npy" % Q)
key = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(m2[0], m2[1]))}
ix = np.array([key[(int(a), int(b))] for a, b in zip(ma, mb)])
print("rows", len(ma), len(m2[0]), " max |A_direct - A_E| =", np.max(np.abs(A - A2[:, ix])), " max|A| =", np.max(np.abs(A)))
