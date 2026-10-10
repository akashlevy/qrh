"""Check the row scalar of the completed cubic reflection (paper78 eq. reflection-row-phase).
For a good prime p active with j_p = 1 the scalar carries omega_{p,1} chi_p(sigma_p)^{-2},
omega_{p,1} = tau^-_{p,1} tau^+_{p,3} chi_p(eps_p)^{-3}.  tau^-_{p,1} = chi_p(-1) gamma_1(p), tau^+_{p,3} = gamma_3(p).
Claim: gamma_1(p) gamma_3(p) = conj(chi_p(4)) chi_p(-1) gamma_2(p)^2   (from eq. signal),
so the row scalar is (ray/Hecke character of p) x gamma_2(p)^2 = (Hecke) x conj(gamma_2(p)),
which is NOT a function of p modulo any fixed modulus (cubic Gauss sums at primes equidistribute).
"""
import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *

def gam(P, j):
    p = P.q
    x = np.arange(1, p, dtype=np.int64)
    s = P.symbol(x, np.zeros_like(x))
    ph = np.exp(-2j*np.pi*((P.b*x) % p)/p)
    return np.sum(Z6[(j*s) % 6]*ph)/np.sqrt(p)

def sym(P, a, b):
    return int(P.symbol(np.array([a]), np.array([b]))[0])

P = [Q for Q in primes_upto(20000) if Q.kind == "split"]
err = 0.0
rows = []
for Q in P:
    g1, g2, g3 = gam(Q,1), gam(Q,2), gam(Q,3)
    c4 = Z6[sym(Q,4,0)]; cm1 = Z6[sym(Q,-1,0)]
    lhs = g1*g3
    rhs = np.conj(c4)*cm1*g2**2
    err = max(err, abs(lhs-rhs))
    # class of the prime mod 72 (a,b mod 72)
    rows.append(((Q.a % 72, Q.b % 72), lhs, g2**3))
print("primes tested:", len(P), " max |g1 g3 - conj(chi(4)) chi(-1) g2^2| =", err)
# Is g1*g3 * alpha^k a function of the class mod 72 for some k?  alpha = -g2^3.
from collections import defaultdict
for k in range(6):
    cls = defaultdict(list)
    for c, v, a3 in rows:
        cls[c].append(v*(-a3)**k)
    spread = [np.ptp(np.angle(np.array(L)*np.conj(L[0]))) for L in cls.values() if len(L) >= 3]
    print("k=%d: classes with >=3 primes: %d, median within-class phase spread of g1g3*alpha^k: %.3f rad" % (k, len(spread), np.median(spread)))
# sanity: a genuine ray quantity (conj(chi(4))chi(-1)) has zero spread
cls = defaultdict(list)
for Q in P:
    cls[(Q.a % 72, Q.b % 72)].append(np.conj(Z6[sym(Q,4,0)])*Z6[sym(Q,-1,0)])
spread = [np.ptp(np.angle(np.array(L)*np.conj(L[0]))) for L in cls.values() if len(L) >= 3]
print("control conj(chi_p(4))chi_p(-1): median within-class spread %.3g" % np.median(spread))
