"""Brute-force check of ffenergy.c on tiny cases: S_{k,P} via direct resultants (Python)."""
import sys, math
import numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ff_lib import *
from ffdriver import run, choose_marks

for (q, kappa, nu, lam, nm) in [(7, 2, 3, 1, 3), (7, 1, 3, 2, 2), (13, 2, 2, 1, 2), (7, 3, 3, 1, 2), (7, 3, 1, 1, 2), (7, 4, 2, 2, 2), (7, 2, 0, 1, 1)]:
    F = Fq(q)
    marks = choose_marks(q, lam, nm)
    R = run(q, kappa, nu, lam, marks, threads=2, tag="test")
    ks = [tuple((int(ki) // q**i) % q for i in range(kappa)) + (1,) for ki in R["kidx"]]   # C row order
    ns = [n for n in monic_polys(nu, q) if resultant(n, deriv(n, q), q) != 0]
    maxerr = 0.0
    assert len(ks) == sum(1 for k in monic_polys(kappa, q) if resultant(k, deriv(k, q), q) != 0)
    for i, k in enumerate(ks):
        for j, P in enumerate(marks):
            s = 0
            for n in ns:
                rk = resultant(k, n, q); rp = resultant(P, n, q)
                if rk == 0 or rp == 0:
                    continue
                s += F.chi(resultant(n, deriv(n, q), q), 5) * F.chi(rk, 3) * F.chi(rp, -2 % 6)
            maxerr = max(maxerr, abs(s - R["S"][i, j]))
    print(q, kappa, nu, lam, "rows", len(ks), "max |S_C - S_py| =", maxerr)
