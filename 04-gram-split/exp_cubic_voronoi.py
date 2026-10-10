"""
Exact cubic analogue of lemma V (Patterson's cubic theta; coefficients a(x) from cubic_theta.py).

(1) Multiplicative twists  F_j(mu) = sum_x a(x) W(N x/Q) chi_mu(x)^j,  chi_mu = (./mu)_6,  mu a good prime.
    R_j(mu) = |F_j(mu)|^2 / sum_x |a(x) W|^2   (random unimodular twist => R ~ 1).
    Voronoi prediction: the a-sum combines chi_mu(a)^{-j} with the multiplier chi_mu(a)^k of the cusp a/mu.
      * j == k (mod 6): Ramanujan sum, n = 0 term survives -> main term (R >> 1 for small mu);
      * j != k, j != 0:  primitive character chi_mu^{k-j} -> NO n = 0 term, dual length ~ C N(mu)^2/Q,
        so R ~ 0 when N(mu)^2 << Q and R ~ 1 on average when N(mu)^2 >> Q.
(2) Additive twists A(a/mu) = sum_x a(x) W e(x a/mu), small mu: A(a/mu) ~ kappa(a) * (constant-term main term),
    so A(a/mu)/A(1/mu) should equal chi_mu(a)^k for one k: a direct measurement of the cusp multiplier.
"""
import numpy as np, sys, time
from zomega import *
from cubic_theta import cubic_theta_coeffs

Q = float(sys.argv[1]) if len(sys.argv) > 1 else 4e5
Lmax = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
t0 = time.time()

def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (u * (1 - u)) + 4); return out

xa, xb, xN, ax, primes, gp = cubic_theta_coeffs(Q / 2, Q)
w = bump(xN / Q, 0.5, 1.0) * ax
den = np.sum(np.abs(w) ** 2)
print("Q=%g  #x=%d  sum|aW|^2=%.1f  t=%.1f" % (Q, len(xN), den, time.time() - t0), flush=True)

# ---- (2) additive twists at small split primes: read off the multiplier -------------------------
print("\n(2) additive twists A(a/mu), a = 1..p-1 (integers), compared with chi_mu(a)^k")
for P in [p for p in primes if p.kind == "split" and p.N <= 80][:10]:
    p = P.q
    pa, pb = P.a - P.b, -P.b                       # conj(mu)
    X, Yc = mul(xa, xb, pa, pb)                    # x * conj(mu) = X + Yc w ; Im = Yc sqrt3/2
    Yr = Yc % p
    A = np.array([np.sum(w * np.exp(2j * np.pi * ((a * Yr) % p) / p)) for a in range(p)])
    a_s = np.arange(1, p); ea = P.symbol(a_s, np.zeros_like(a_s))
    rat = A[1:] / A[1]
    fits = []
    for k in range(6):
        pred = Z6[(k * ea) % 6] / Z6[(k * ea[0]) % 6]
        fits.append(np.mean(np.abs(rat - pred)))
    kbest = int(np.argmin(fits))
    print("  N(mu)=%3d  |A(0)|=%8.1f  mean|A(a)|,a!=0 =%8.1f  sd|A(a)|=%6.1f   mean|ratio - chi^k|, k=0..5: %s  -> k=%d"
          % (P.N, abs(A[0]), np.mean(np.abs(A[1:])), np.std(np.abs(A[1:])), " ".join("%.3f" % f for f in fits), kbest))

# ---- (1) multiplicative twists ---------------------------------------------------------------
print("\n(1) multiplicative twists: mean R_j by bins of N(mu)^2/Q")
mus = [P for P in primes_upto(Lmax)]
rng = np.random.default_rng(1)
if len(mus) > 700:
    small = [P for P in mus if P.N <= 1500]
    big = [P for P in mus if P.N > 1500]
    mus = small + list(rng.choice(big, 700 - len(small), replace=False))
bins = {}
for P in mus:
    e = P.symbol(xa, xb); ok = e >= 0
    x2 = P.N ** 2 / Q
    b = int(np.floor(np.log2(x2)))
    for j in range(6):
        F = np.sum(np.where(ok, w * Z6[(j * e) % 6], 0))
        bins.setdefault((j, b), []).append(abs(F) ** 2 / den)
print("t=%.1f" % (time.time() - t0))
bs = sorted(set(k[1] for k in bins))
print("log2(N(mu)^2/Q) " + " ".join("%5d" % b for b in bs))
print("#mu             " + " ".join("%5d" % len(bins[(0, b)]) for b in bs))
for j in range(6):
    print("j=%d  mean R     " % j + " ".join("%5.2f" % np.mean(bins[(j, b)]) for b in bs))
for j in range(6):
    print("j=%d  median R   " % j + " ".join("%5.2f" % np.median(bins[(j, b)]) for b in bs))
