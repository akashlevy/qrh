"""Attack D, check (N3): generalised Lemma Z for an arbitrary uniform twist exponent j.

Paper78 eq. reflection-row-phase / reflection-cross-prime-phase give, for an active prime p with exponent j
(j != 0, 4) in the completed reflection of the cubic-theta row (theta-bar, multiplier shift +2):
     chi_p(sigma_p)^{-2} omega_{p,j} = xi_{p,j} chi_p(c/p)^{2j+2},
     omega_{p,j} = tau^-_{p,j} tau^+_{p,j+2} chi_p(eps_p)^{-j-2},   tau^-_{p,m} = chi_p(-1)^m gamma_m(p),
with xi_{p,j} depending on p, j and the fixed data only.  For a squarefree R = p1 p2 (c = c_F R, take c_F = 1)
the claim is that the moving cross terms cancel against the twisted multiplicativity of gamma_j, gamma_{j+2}:
     prod_{p | R} [ chi_p(R/p)^{2j+2} gamma_j(p) gamma_{j+2}(p) ]  =  gamma_j(R) gamma_{j+2}(R) * (ray class phase)
so that zeta_R ~ abar(R)^2 gamma_j(R) gamma_{j+2}(R) * (ray), of Gauss class 2j+2 mod 3.
We test: ratio = LHS / (gamma_j(R) gamma_{j+2}(R)) is a function of the classes of (p1, p2) mod 12.
"""
import sys, numpy as np
from collections import defaultdict
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
from zomega import *          # noqa
from a2lib import gamma_prime, sym1, gauss_bruteforce

P = [Q for Q in primes_upto(400) if Q.kind == "split"]
pairs = [(P[i], P[k]) for i in range(len(P)) for k in range(i + 1, len(P)) if P[i].N != P[k].N and P[i].N * P[k].N < 20000]
print("pairs R = p1 p2:", len(pairs))
for j in [1, 2, 3, 5]:
    jj = (j + 2) % 6
    d = defaultdict(list); mod_err = 0.0
    for (A, B) in pairs:
        Ra, Rb = mul(A.a, A.b, B.a, B.b)
        lhs = 1.0 + 0j
        for (X, Y) in [(A, B), (B, A)]:
            lhs *= Z6[(sym1(X, Y.a, Y.b) * (2 * j + 2)) % 6] * gamma_prime(X, j) * (gamma_prime(X, jj) if jj else 1.0)
        rhs = gauss_bruteforce(Ra, Rb, [A, B], j) * (gauss_bruteforce(Ra, Rb, [A, B], jj) if jj else 1.0)
        mod_err = max(mod_err, abs(abs(rhs) - 1))
        d[(A.a % 12, A.b % 12, B.a % 12, B.b % 12)].append(lhs / rhs)
    sp = max(np.ptp(np.angle(np.array(L) * np.conj(L[0]))) for L in d.values() if len(L) >= 2)
    ncls = sum(1 for L in d.values() if len(L) >= 2)
    print("j=%d: zeta ~ gamma_%d gamma_%d, class 2j+2 = %d mod 3; |gamma_j(R)gamma_j+2(R)| - 1 <= %.1e; "
          "within-class spread of ratio over %d classes: %.2e" % (j, j, jj, (2 * j + 2) % 3, mod_err, ncls, sp))
