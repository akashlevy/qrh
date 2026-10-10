"""
Cubic analogue, dual-length test: for fixed small mu, R_j(mu) = |F_j(mu)|^2 / sum|aW|^2 as Q grows.
Voronoi: for j not in {0, k} (k = multiplier exponent, = 4 here) the n = 0 dual term vanishes and the
dual sum has length ~ C N(mu)^2 / Q, so R_j(mu) -> 0 (rapidly) once Q >> C N(mu)^2.
Random-twist behaviour would keep R_j ~ 1 for every Q.
Gauss sums of primes are cached in cubic_gauss_cache.npz.
"""
import numpy as np, sys, os, time
from zomega import *
from cubic_theta import cubic_gauss_echeck, squarefree_with_gauss, primary_odd_upto

Qmax = int(float(sys.argv[1])) if len(sys.argv) > 1 else 400000
Qs = [Qmax // 16, Qmax // 8, Qmax // 4, Qmax // 2, Qmax]
t0 = time.time()
primes = primes_upto(Qmax)
cache = "cubic_gauss_cache_%d.npz" % Qmax
if os.path.exists(cache):
    gp = np.load(cache)["gp"]
else:
    gp = np.array([cubic_gauss_echeck(P) for P in primes]); np.savez(cache, gp=gp)
print("primes", len(primes), "gauss done t=%.0f" % (time.time() - t0), flush=True)
A, B, N, G = squarefree_with_gauss(Qmax, primes, gp)
print("sf n:", len(N), "t=%.0f" % (time.time() - t0), flush=True)
bl = [b for b in primary_odd_upto(int(round(Qmax ** (1 / 3))) + 1) if b[2] ** 3 <= Qmax]

def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (u * (1 - u)) + 4); return out

mus = [P for P in primes if P.N <= 400]
res = {}
for Q in Qs:
    xs = []
    for (ba, bb, bn) in bl:
        keep = (N * bn ** 3 > Q / 2) & (N * bn ** 3 <= Q)
        if not keep.any(): continue
        ca, cb = mul(ba, bb, ba, bb); ca, cb = mul(ca, cb, ba, bb)
        xa, xb = mul(A[keep], B[keep], ca, cb)
        xs.append((xa, xb, N[keep] * bn ** 3, np.sqrt(bn) * G[keep]))
    xa = np.concatenate([t[0] for t in xs]); xb = np.concatenate([t[1] for t in xs])
    xN = np.concatenate([t[2] for t in xs]); ax = np.concatenate([t[3] for t in xs])
    w = bump(xN / Q, 0.5, 1.0) * ax; den = np.sum(np.abs(w) ** 2)
    for P in mus:
        e = P.symbol(xa, xb); ok = e >= 0
        for j in range(6):
            F = np.sum(np.where(ok, w * Z6[(j * e) % 6], 0))
            res[(Q, P.N, P.a, P.b, j)] = abs(F) ** 2 / den
    print("Q=%d done (#x=%d) t=%.0f" % (Q, len(xN), time.time() - t0), flush=True)

norms = sorted(set(P.N for P in mus))
for j in range(6):
    print("\nj=%d   R_j averaged over the primes of each norm;  columns Q = %s" % (j, Qs))
    for nn in norms:
        ps = [P for P in mus if P.N == nn]
        row = [np.mean([res[(Q, nn, P.a, P.b, j)] for P in ps]) for Q in Qs]
        print("  N(mu)=%4d (%d primes)  " % (nn, len(ps)) + "  ".join("%8.3f" % r for r in row)
              + "    N^2/Q at Qmax = %.4f" % (nn ** 2 / Qmax))
