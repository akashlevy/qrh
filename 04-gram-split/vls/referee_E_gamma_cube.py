import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2"); sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from a2lib import *
from cubic_theta import cubic_gauss_echeck
P = primes_upto(3000)
out = []
for Pr in P[:60]:
    g2 = gamma_prime(Pr, 2); al = cval(Pr.a, Pr.b); al /= abs(al)
    g3c = cubic_gauss_echeck(Pr)
    out.append((Pr.N, Pr.kind, np.round(g2 ** 3 / al, 6), np.round(g3c / g2, 4)))
r = np.array([o[2] for o in out])
print("gamma_2(p)^3 / alpha(p) over 60 primes: real parts set", sorted(set(np.round(r.real, 6))), "max imag", np.abs(r.imag).max())
print("examples:", out[:6])
