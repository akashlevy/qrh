"""Voronoi test: coherence ratio of the quadratic twist (j=3) and other j as a function of N(mu)^2/Q."""
import numpy as np, sys
from zomega import *
Q = int(sys.argv[1]); i0 = int(sys.argv[2]) if len(sys.argv) > 2 else 4
A = np.load("A_Q%d.npy" % Q); ma, mb, mN = np.load("m_Q%d.npy" % Q)
K = A.shape[0]; Ys = np.array([256 * 2 ** i for i in range(K)], float)
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (u * (1 - u)) + 4); return out
Om = bump(mN / Q, .5, 1)
V = A[i0:] * Ys[i0:, None] ** (1 / 3)
pairs = [(i, j) for i in range(len(V)) for j in range(i + 1, len(V))]
bn = np.mean([np.sum(Om ** 2 * V[i] * np.conj(V[j])).real for i, j in pairs])
P = [Pr for Pr in primes_upto(4000)]
bins = {}
for Pr in P:
    e = Pr.symbol(ma, mb); ok = e >= 0
    x = Pr.N ** 2 / Q
    b = int(np.floor(np.log2(x))) if x > 0 else -99
    for j in (1, 2, 3):
        psi = np.where(ok, Z6[(j * e) % 6], 0)
        T = [np.sum(Om * V[i] * psi) for i in range(len(V))]
        r = np.mean([(T[i] * np.conj(T[k])).real for i, k in pairs]) / bn
        bins.setdefault((j, b), []).append(r)
print("Q=%d  sum Om^2|rho|^2 (cross-scale est) = %.4f" % (Q, bn))
print("ratio = <|F(mu)|^2>/sum Om^2|rho|^2 ; random twist -> 1; Voronoi predicts << 1 when N(mu)^2/Q << 1")
for j in (1, 2, 3):
    print("j=%d:" % j)
    for b in sorted(set(k[1] for k in bins if k[0] == j)):
        r = np.array(bins[(j, b)])
        print("   N(mu)^2/Q ~ 2^%3d : n=%3d  mean ratio %7.3f  (se %.3f)" % (b, len(r), r.mean(), r.std() / np.sqrt(len(r))))
