"""Sextic: cross-scale estimate of <|F_j(mu)|^2> / sum Om^2 |What rho|^2 for all j = 0..5, in broad bins of N(mu)^2/Q.
Voronoi prediction (this file's orientation, rho from g(s,-m)): coherent (main term) only for j = 0 and j = 1;
j = 2..5 have primitive dual characters chi_mu^{1-j} and no n = 0 term.  (Cross-scale estimator; noisy.)"""
import numpy as np, sys
from zomega import *
Q = int(sys.argv[1]) if len(sys.argv) > 1 else 3000; i0 = int(sys.argv[2]) if len(sys.argv) > 2 else 4
A = np.load("A_Q%d.npy" % Q); ma, mb, mN = np.load("m_Q%d.npy" % Q)
K = A.shape[0]; Ys = np.array([256 * 2 ** i for i in range(K)], float)
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (u * (1 - u)) + 4); return out
Om = bump(mN / Q, .5, 1)
V = A[i0:] * Ys[i0:, None] ** (1 / 3)
pairs = [(i, j) for i in range(len(V)) for j in range(i + 1, len(V))]
bn = np.mean([np.sum(Om ** 2 * V[i] * np.conj(V[j])).real for i, j in pairs])
edges = [0, 0.25, 4, 64, 1e9]
out = {}
for Pr in primes_upto(4000):
    e = Pr.symbol(ma, mb); ok = e >= 0; x = Pr.N ** 2 / Q
    b = int(np.searchsorted(edges, x) - 1)
    for j in range(6):
        psi = np.where(ok, Z6[(j * e) % 6], 0)
        T = [np.sum(Om * V[i] * psi) for i in range(len(V))]
        out.setdefault((j, b), []).append(np.mean([(T[i] * np.conj(T[k])).real for i, k in pairs]) / bn)
print("Q=%d, scales i>=%d; entries: mean ratio (s.e.) [#mu]" % (Q, i0))
for b in range(len(edges) - 1):
    print("N(mu)^2/Q in [%g,%g): " % (edges[b], edges[b + 1]) + "  ".join(
        "j=%d %6.2f(%4.2f)" % (j, np.mean(out[(j, b)]), np.std(out[(j, b)]) / np.sqrt(len(out[(j, b)]))) for j in range(6))
        + "  [%d]" % len(out[(0, b)]))
