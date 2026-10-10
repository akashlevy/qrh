"""Composite squarefree R = pi1*pi2 (split primes, distinct norms): check
 (a) prod_{p|R} chi_p(R/p)^{-2} gamma_1(p) gamma_3(p) = gamma_1(R) gamma_3(R)   (cross terms of the row scalar cancel)
 (b) gamma_1(R) gamma_3(R) = conj(chi_R(4)) chi_R(-1) gamma_2(R)^2."""
import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
P = [Q for Q in primes_upto(400) if Q.kind == "split"]
def symexp(Q, a, b):
    return Q.symbol(np.atleast_1d(np.asarray(a)), np.atleast_1d(np.asarray(b)))
def gam_mod(primes, Ra, Rb, j):
    N = 1
    for Q in primes: N *= Q.q
    x = np.arange(N, dtype=np.int64)
    tot = np.zeros(N, np.int64); ok = np.ones(N, bool)
    for Q in primes:
        s = symexp(Q, x, np.zeros_like(x)); ok &= s >= 0; tot += np.where(s >= 0, s, 0)
    ph = np.exp(-2j*np.pi*((Rb*x) % N)/N)
    return np.sum(np.where(ok, Z6[(j*tot) % 6], 0)*ph)/np.sqrt(N)
err_a = err_b = 0.0; n = 0
for i in range(len(P)):
    for k in range(i+1, len(P)):
        p1, p2 = P[i], P[k]
        if p1.q == p2.q: continue
        Ra, Rb = mul(p1.a, p1.b, p2.a, p2.b)
        g = {j: gam_mod([p1, p2], Ra, Rb, j) for j in (1, 2, 3)}
        gp = {(Q, j): gam_mod([Q], Q.a, Q.b, j) for Q in (p1, p2) for j in (1, 3)}
        c12 = Z6[int(symexp(p1, p2.a, p2.b)[0])]; c21 = Z6[int(symexp(p2, p1.a, p1.b)[0])]
        lhs = c12**-2 * c21**-2 * gp[(p1,1)]*gp[(p1,3)]*gp[(p2,1)]*gp[(p2,3)]
        err_a = max(err_a, abs(lhs - g[1]*g[3]))
        chi4 = Z6[(int(symexp(p1,4,0)[0]) + int(symexp(p2,4,0)[0])) % 6]
        chim1 = Z6[(int(symexp(p1,-1,0)[0]) + int(symexp(p2,-1,0)[0])) % 6]
        err_b = max(err_b, abs(g[1]*g[3] - np.conj(chi4)*chim1*g[2]**2))
        n += 1
        if n >= 60: break
    if n >= 60: break
print("pairs:", n, " max err (a):", err_a, " max err (b):", err_b)
