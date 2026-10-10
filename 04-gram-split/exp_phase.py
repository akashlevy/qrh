"""Do the twisted theta sums F_j(c) = sum_m Om rho(m) (m/c)_6^j share a common phase across primes c?
Unbiased cross-scale estimator of E[F(c) conj F(c')] for c != c'. Compare mean off-diagonal with diagonal."""
import numpy as np, sys
from zomega import *
Q = 3000; i0 = 4
A = np.load("A_Q%d.npy" % Q); ma, mb, mN = np.load("m_Q%d.npy" % Q)
K = A.shape[0]; Ys = np.array([256 * 2 ** i for i in range(K)], float)
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (u * (1 - u)) + 4); return out
Om = bump(mN / Q, .5, 1)
V = A[i0:] * Ys[i0:, None] ** (1 / 3); nV = len(V)
P = [Pr for Pr in primes_upto(int(sys.argv[1]) if len(sys.argv) > 1 else 3000) if Pr.N > 30]
for j in (1, 2, 3, 4, 5):
    T = np.zeros((len(P), nV), complex)
    for k, Pr in enumerate(P):
        e = Pr.symbol(ma, mb); ok = e >= 0
        psi = np.where(ok, Z6[(j * e) % 6], 0)
        T[k] = (V * Om * psi).sum(axis=1)
    # unbiased cross: sum over scale pairs i != i' of T[c,i] conj T[c',i']
    S = T.sum(axis=1); S2 = (T * T.conj()).sum(axis=1)       # per c
    G = np.outer(S, S.conj()) - (T @ T.conj().T)              # sum_{i != i'} T[c,i] conj T[c',i']
    npairs = nV * (nV - 1)
    G /= npairs
    diag = np.real(np.diag(G)); off = G[~np.eye(len(P), dtype=bool)]
    print("j=%d: mean diag |F|^2 = %.4f ; mean offdiag Re F(c)conjF(c') = %.4f (se %.4f) ; ratio %.3f" %
          (j, diag.mean(), off.real.mean(), off.real.std() / np.sqrt(len(off)), off.real.mean() / diag.mean()))
