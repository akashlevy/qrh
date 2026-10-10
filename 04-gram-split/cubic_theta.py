"""
Exact coefficients of (a translate-combination of) Patterson's cubic theta function, for the
cubic analogue of lemma V.

Source formula (paper78, eq. reflection-infinity-coefficient, from Dunn-Radziwill (5.7)):
the coefficient of  theta-bar  at the index lambda^{-3} x, x = n b^3 == 1 (mod 3), n squarefree, b primary, is
        3^{5/2} |b| gt3(n),     gt3(n) = q_n^{-1/2} sum_{v mod n} (v/n)_3 echeck(v/n),
with echeck(z) = exp(2 pi i (z + zbar)) and (v/n)_3 = chi_n(v)^2.  Restricting the frequency to the class
lambda^3 mu == 1 (mod 3) is an additive (congruence) condition, so the sequence
        a(x) = |b| gt3(n)   (x = n b^3 primary, n sf),   a(x) = 0 otherwise
is the coefficient sequence of a finite combination of translates of theta-bar: it obeys the same
Voronoi formula, with the same multiplier, as theta-bar itself.  We also drop x divisible by 2
(another congruence mask), so that every prime involved is a Prime of zomega.
"""
import numpy as np
from math import isqrt
from zomega import *

def cubic_gauss_echeck(P):
    """gt3(P) = N^{-1/2} sum_x (x/P)_3 echeck(x/P)."""
    if P.kind == "split":
        p = P.q
        x = np.arange(1, p, dtype=np.int64)
        j = P.symbol(x, np.zeros_like(x))
        ph = np.exp(2j * np.pi * ((x * (2 * P.a - P.b)) % p) / p)     # echeck(x/pi), pi = a + b w
        return np.sum(Z6[(2 * j) % 6] * ph) / np.sqrt(P.N)
    q = P.q
    c, d = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    c = c.ravel(); d = d.ravel(); nz = (c != 0) | (d != 0); c = c[nz]; d = d[nz]
    j = P.symbol(c, d)
    sgn = -1 if P.a < 0 else 1                                            # pi = sgn*q
    ph = np.exp(2j * np.pi * sgn * ((2 * c - d) % q) / q)
    return np.sum(Z6[(2 * j) % 6] * ph) / np.sqrt(P.N)

def squarefree_with_gauss(Qmax, primes, gp):
    """all squarefree primary n (prime to 6) with N(n) <= Qmax, with gt3(n) by twisted multiplicativity
       gt3(n p) = gt3(n) gt3(p) (n/p)_3 (p/n)_3 = gt3(n) gt3(p) (n/p)_3^2   (cubic reciprocity, primary)."""
    rt = isqrt(int(Qmax))
    small = [i for i, P in enumerate(primes) if P.N <= rt]
    A = np.array([1], np.int64); B = np.array([0], np.int64); N = np.array([1], np.int64); G = np.array([1 + 0j])
    for i in small:
        P = primes[i]
        ok = N * P.N <= Qmax
        if not ok.any(): continue
        a0, b0, n0, g0 = A[ok], B[ok], N[ok], G[ok]
        e = P.symbol(a0, b0)                    # chi_P(n), n coprime to P
        na, nb = mul(a0, b0, P.a, P.b)
        A = np.concatenate([A, na]); B = np.concatenate([B, nb]); N = np.concatenate([N, n0 * P.N])
        G = np.concatenate([G, g0 * gp[i] * Z6[(4 * e) % 6]])
    o = np.argsort(N, kind="stable"); A, B, N, G = A[o], B[o], N[o], G[o]
    outA, outB, outN, outG = [A], [B], [N], [G]
    for i, P in enumerate(primes):
        if P.N <= rt: continue
        k = np.searchsorted(N, Qmax // P.N, side="right")
        if k == 0: continue
        a0, b0, n0, g0 = A[:k], B[:k], N[:k], G[:k]
        e = P.symbol(a0, b0)
        na, nb = mul(a0, b0, P.a, P.b)
        outA.append(na); outB.append(nb); outN.append(n0 * P.N); outG.append(g0 * gp[i] * Z6[(4 * e) % 6])
    return np.concatenate(outA), np.concatenate(outB), np.concatenate(outN), np.concatenate(outG)

def primary_odd_upto(Nmax):
    R = isqrt(int(4 * Nmax / 3)) + 2
    out = []
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            n = norm(a, b)
            if 1 <= n <= Nmax and a % 3 == 1 and b % 3 == 0 and n % 2 == 1:
                out.append((a, b, n))
    return out

def cubic_theta_coeffs(Qlo, Qhi):
    """x = n b^3 in the annulus Qlo < N(x) <= Qhi, x primary and odd; returns xa, xb, N, a(x)."""
    Qhi = int(Qhi); Qlo = int(Qlo)
    primes = primes_upto(Qhi)
    gp =np.array([cubic_gauss_echeck(P) for P in primes])
    A, B, N, G = squarefree_with_gauss(Qhi, primes, gp)
    xs = []
    for (ba, bb, bn) in primary_odd_upto(int(round(Qhi ** (1 / 3))) + 1):
        if bn ** 3 > Qhi: continue
        ca, cb = mul(ba, bb, ba, bb); ca, cb = mul(ca, cb, ba, bb)          # b^3
        keep = (N * bn ** 3 > Qlo) & (N * bn ** 3 <= Qhi)
        xa, xb = mul(A[keep], B[keep], ca, cb)
        xs.append((xa, xb, N[keep] * bn ** 3, np.sqrt(bn) * G[keep]))
    xa = np.concatenate([t[0] for t in xs]); xb = np.concatenate([t[1] for t in xs])
    xN = np.concatenate([t[2] for t in xs]); ax = np.concatenate([t[3] for t in xs])
    return xa, xb, xN, ax, primes, gp
