"""Patterson bias on the principal row k = 1: B(nu) = sum_{n sf monic deg nu} chi6bar(disc n) chi3bar(Res(P,n))
(= c_nu^{-1} sum beta(n) chi_P(n)^{-2}); compare |B|^2/q^nu (square-root size = 1) and |B|/q^{5nu/6}."""
import sys, math, itertools
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ff_lib import Fq, resultant, deriv
q = int(sys.argv[1]); numax = int(sys.argv[2])
F = Fq(q)
for P, lam in [((1,), 0), ((0, 1), 1)]:
    for nu in range(1, numax + 1):
        B = 0; cnt = 0
        for c in itertools.product(range(q), repeat=nu):
            n = tuple(c) + (1,)
            r = resultant(n, deriv(n, q), q)
            if r == 0: continue
            rp = resultant(P, n, q) if lam else 1
            if rp == 0: continue
            B += F.chi(r, 5) * F.chi(rp, 4); cnt += 1
        print(f"q={q} lam={lam} nu={nu}: count={cnt:7d} |B|^2/q^nu={abs(B)**2/q**nu:9.3f}  |B|/q^(5nu/6)={abs(B)/q**(5*nu/6):.4f}", flush=True)
