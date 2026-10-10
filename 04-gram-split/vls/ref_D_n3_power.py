# Referee check: does the N3 test (D_lemmaZ_general.py) discriminate the exponent 2j+2 from alternatives?
import sys, numpy as np
from collections import defaultdict
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
from zomega import *
from a2lib import gamma_prime, sym1, gauss_bruteforce
P = [Q for Q in primes_upto(300) if Q.kind == "split"]
pairs = [(P[i], P[k]) for i in range(len(P)) for k in range(i + 1, len(P)) if P[i].N != P[k].N and P[i].N * P[k].N < 12000]
print("pairs", len(pairs))
cache = {}
for j in [1, 2]:
    jj = (j + 2) % 6
    rhsl = []
    for (A, B) in pairs:
        Ra, Rb = mul(A.a, A.b, B.a, B.b)
        rhsl.append(gauss_bruteforce(Ra, Rb, [A, B], j) * (gauss_bruteforce(Ra, Rb, [A, B], jj) if jj else 1.0))
    for e in [2*j+2, -(2*j+2), 2*j+1, 2*j+3, 0]:
        d = defaultdict(list); allr = []
        for (A, B), rhs in zip(pairs, rhsl):
            lhs = 1.0 + 0j
            for (X, Y) in [(A, B), (B, A)]:
                lhs *= Z6[(sym1(X, Y.a, Y.b) * e) % 6] * gamma_prime(X, j) * (gamma_prime(X, jj) if jj else 1.0)
            r = lhs / rhs; d[(A.a % 12, A.b % 12, B.a % 12, B.b % 12)].append(r); allr.append(r)
        sp = max(np.ptp(np.angle(np.array(L) * np.conj(L[0]))) for L in d.values() if len(L) >= 2)
        allr = np.array(allr)
        print("j=%d exponent %+d: within-class spread %.2e ; distinct global values %d" % (j, e, sp, len(set(np.round(np.angle(allr), 6)))))
