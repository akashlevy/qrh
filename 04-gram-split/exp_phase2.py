"""Common phase test of F_j(c) * conj(gamma_j-ish(c)) across primes c, using cross-scale unbiased products."""
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
P = [Pr for Pr in primes_upto(int(sys.argv[1]) if len(sys.argv) > 1 else 2000) if Pr.N > 30]
def gam(Pr, j):
    if Pr.kind == "split":
        p = Pr.q; x = np.arange(1, p); jj = Pr.symbol(x, np.zeros_like(x))
        return np.sum(Z6[(j * jj) % 6] * np.exp(-2j * np.pi * ((Pr.b * x) % p) / p)) / np.sqrt(Pr.N)
    q = Pr.q; c, d = np.meshgrid(np.arange(q), np.arange(q), indexing="ij"); c = c.ravel(); d = d.ravel()
    nz = (c != 0) | (d != 0); c = c[nz]; d = d[nz]; jj = Pr.symbol(c, d); sgn = -1 if Pr.a < 0 else 1
    return np.sum(Z6[(j * jj) % 6] * np.exp(-2j * np.pi * sgn * d / q)) / np.sqrt(Pr.N)
for j in (1, 2, 3, 4, 5):
    T = np.zeros((len(P), nV), complex)
    for k, Pr in enumerate(P):
        e = Pr.symbol(ma, mb); ok = e >= 0
        psi = np.where(ok, Z6[(j * e) % 6], 0)
        T[k] = (V * Om * psi).sum(axis=1)
    for jg in range(6):
        gj = np.array([gam(Pr, jg) if jg % 6 else 1.0 for Pr in P])
        for conjflag in (False, True):
            ph = np.conj(gj) if conjflag else gj
            Tp = T * ph[:, None]
            S = Tp.sum(axis=1)
            G = np.outer(S, S.conj()) - (Tp @ Tp.conj().T)
            G /= nV * (nV - 1)
            diag = np.real(np.diag(G)).mean(); off = G[~np.eye(len(P), dtype=bool)]
            r = off.real.mean() / diag; se = off.real.std() / np.sqrt(len(off)) / diag
            if abs(r) > 4 * se and abs(r) > 0.05:
                print("j=%d  phase %s gamma_%d : offdiag/diag = %.3f (se %.3f)  <-- common phase" % (j, "conj" if conjflag else "", jg, r, se))
    print("j=%d done" % j)
