"""(c') prime rows R = p (split primary, N(p) <= 30000): compare check_row_phase.py's
   w(p) = conj(alpha(p))^2 tau^-_{p,1} tau^+_{p,3}   with   zetaZ(p) = psi_1(p) conj(gamma_2(p))^4 :
   the paper's zeta at a prime row differs from w only by chi_p(sigma_p)^{-2} chi_p(eps_p)^{-3}, residue symbols of
   fixed elements, so w/zetaZ must be a ray class function (constant on classes mod 72 suffices here)."""
import numpy as np
from collections import defaultdict
from a2lib import *
P = [Pr for Pr in primes_upto(30000) if Pr.kind == "split"]
cls = defaultdict(list); cls_w = defaultdict(list)
for Pr in P:
    al = cval(Pr.a, Pr.b); al /= abs(al)
    g2 = gamma_prime(Pr, 2); g1m = gamma_prime(Pr, 1, -1); g3 = gamma_prime(Pr, 3)
    psi1 = Z6[(-sym1(Pr, *LAM) - sym1(Pr, 4, 0) + sym1(Pr, -1, 0)) % 6]
    zZ = psi1 * np.conj(g2) ** 4
    w = np.conj(al) ** 2 * g1m * g3
    key = (Pr.a % 72, Pr.b % 72)
    cls[key].append(w / zZ); cls_w[key].append(w)
nc = sum(len(L) >= 2 for L in cls.values())
bad = sum(np.max(np.abs(np.array(L) - L[0])) > 1e-9 for L in cls.values() if len(L) >= 2)
spr = np.median([np.max(np.abs(np.array(L) - L[0])) for L in cls_w.values() if len(L) >= 2])
vals = sorted(set(np.round(np.angle(np.concatenate([np.array(L) for L in cls.values()])) / (np.pi / 6)).astype(int) % 12))
print("split primes:", len(P), " classes mod 72 with >=2 distinct primes:", nc)
print("classes where w/zetaZ is NOT constant: %d   (values of w/zetaZ, in units of pi/6: %s)" % (bad, vals))
print("for contrast, median within-class spread max|w - w0| of w itself: %.2f" % spr)
