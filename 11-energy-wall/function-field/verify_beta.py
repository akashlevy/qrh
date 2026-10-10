"""(G5') beta(n) = gamma_2(n) conj(alpha(n)) equals c_{deg n} * conj(chi6(disc n)) with |c_d| = 1 depending only on d."""
import sys, math, cmath
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ff_lib import *
for q, maxd in [(7, 4), (13, 3), (19, 2)]:
    F = Fq(q); tau2 = F.tau(2)
    eps0 = -((-(-tau2 / math.sqrt(q))) ** 3)
    vals = {}
    for d in range(1, maxd + 1):
        for n in monic_polys(d, q):
            rnd = resultant(n, deriv(n, q), q)
            if rnd == 0: continue
            mu = (-1) ** len(factor(n, q))
            g2 = mu * (-tau2 / math.sqrt(q)) ** d * F.chi(rnd, 2)
            beta = g2 * eps0.conjugate() ** d
            disc = (-1) ** (d * (d - 1) // 2) * rnd % q
            r = beta / F.chi(disc, 5)
            vals.setdefault(d, set()).add((round(r.real, 9), round(r.imag, 9)))
    print(q, {d: len(v) for d, v in vals.items()}, "  c_d:", {d: sorted(v)[0] for d, v in vals.items()})
