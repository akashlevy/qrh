import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split"); sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
from zomega import *
from a2lib import gamma_prime, sym1
Ps = [P for P in primes_upto(6000) if P.kind == "split"]
e1 = e2 = e3 = 0.0; cls = {}
for P in Ps:
    g = {j: gamma_prime(P, j) for j in range(1, 6)}
    al = cval(P.a, P.b) / abs(cval(P.a, P.b))
    e1 = max(e1, abs(g[2] ** 3 - (-1) * al))                     # gamma2^3 = mu*alpha, mu(p) = -1
    chi4 = Z6[sym1(P, 4, 0) % 6]
    G = np.conj(chi4) * g[3]
    e2 = max(e2, abs(g[1] * g[2] - (-1) * al * G))               # gamma1 gamma2 = mu alpha G
    e3 = max(e3, abs(np.conj(g[2]) ** 4 - np.conj(al) ** 2 * np.conj(G) * g[1]))   # form 3: zeta ~ abar^2 Gbar gamma1
    key = (P.a % 72, P.b % 72)
    cls.setdefault(key, []).append(g[3])
spread = max(max(abs(v - vs[0]) for v in vs) for vs in cls.values())
print("n split primes", len(Ps))
print("gamma2^3 = mu alpha  err", e1)
print("gamma1 gamma2 = mu alpha G (G = conj chi(4) gamma3) err", e2)
print("conj(gamma2)^4 = abar^2 conj(G) gamma1 err", e3)
print("gamma3 constant on classes (a,b) mod 72: max spread", spread, "#classes", len(cls))
