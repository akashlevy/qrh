"""
Brute-force checks over F_q[T] (q = 7, 13), squarefree monic n of degree <= 3:
 (G1) g_j(n) := sum_{a mod n} (a/n)_6^j e(a/n) = mu(n) (-tau_j)^{deg n} (n'/n)_6^j   (j = 1..5)
 (G2) mu(n) = (-1)^{deg n} chi2(disc n),  disc n = (-1)^{d(d-1)/2} Res(n, n')
 (G3) sextic reciprocity (a/b)_6 = (-1)^{((q-1)/6) deg a deg b} (b/a)_6, a,b monic coprime
 (G4) gamma_2(P)^3 = -alpha(P), alpha(n) := eps0^{deg n}, eps0 = -(tau_2/sqrt q)^3 ... (constant)
 (G5) beta(n) := gamma_2(n) conj(alpha(n)) = C^{deg n} chi6(disc n) with |C| = 1 constant
"""
import sys, cmath, math, itertools
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ff_lib import *


def gauss_brute(n, j, F):
    q = F.q
    d = deg(n)
    tot = 0
    for a in itertools.product(range(q), repeat=d):
        a = trim(a)
        if not a:
            continue
        r = resultant(n, a, q)
        if r == 0:
            continue
        c1 = laurent(a, n, q, 1)[0]
        tot += F.chi(r, j) * cmath.exp(2j * math.pi * c1 / q)
    return tot


def run(q, maxd):
    F = Fq(q)
    tau = {j: F.tau(j) for j in range(1, 6)}
    bad = {"G1": 0, "G2": 0, "G3": 0, "G5": 0}
    cnt = 0
    eps_vals, C_vals = set(), set()
    for d in range(1, maxd + 1):
        for n in monic_polys(d, q):
            if not is_squarefree(n, q):
                continue
            fac = factor(n, q)
            mu = (-1) ** len(fac)
            dn = deriv(n, q)
            rnd = resultant(n, dn, q)
            disc = (-1) ** (d * (d - 1) // 2) * rnd % q
            # G2
            chi2disc = F.chi(disc, 3)
            if abs(mu - (-1) ** d * chi2disc) > 1e-9:
                bad["G2"] += 1
            for j in range(1, 6):
                gb = gauss_brute(n, j, F)
                gc = mu * (-tau[j]) ** d * F.chi(rnd, j)
                if abs(gb - gc) > 1e-6:
                    bad["G1"] += 1
            # gamma_2(P)^3 for primes
            g2 = mu * (-tau[2] / math.sqrt(q)) ** d * F.chi(rnd, 2)
            if len(fac) == 1:
                eps = -(g2 ** 3)
                eps_vals.add((d, round((eps ** (1.0 / 1)).real, 6), round(eps.imag, 6)))
            cnt += 1
    # eps0: gamma_2(P)^3 = -eps0^{deg P}
    eps0 = -((-tau[2] / math.sqrt(q)) ** 3) * (-1)   # check below numerically
    # compute eps0 from degree-1 primes directly
    e1 = [-(((-1) * (-tau[2] / math.sqrt(q)) * F.chi(1, 2)) ** 3)]
    # G3 reciprocity on random small pairs
    for da in range(1, 3):
        for db in range(1, 3):
            for a in monic_polys(da, q):
                for b in itertools.islice(monic_polys(db, q), 0, 40):
                    if resultant(a, b, q) == 0:
                        continue
                    lhs = F.chi(resultant(b, a, q))      # (a/b)_6
                    rhs = F.chi(resultant(a, b, q))      # (b/a)_6
                    sgn = (-1) ** (((q - 1) // 6) * da * db)
                    if abs(lhs - sgn * rhs) > 1e-9:
                        bad["G3"] += 1
    # G5: beta(n) = gamma2(n) conj(alpha(n)), alpha = eps0^deg with gamma2(P)^3 = -alpha(P)
    # degree-1 prime P = T: n' = 1, Res = 1 -> gamma2(T) = -(-tau2/sqrt q)
    g2T = -(-tau[2] / math.sqrt(q))
    eps0 = -(g2T ** 3)
    for d in range(1, maxd + 1):
        for n in monic_polys(d, q):
            if not is_squarefree(n, q):
                continue
            fac = factor(n, q)
            mu = (-1) ** len(fac)
            rnd = resultant(n, deriv(n, q), q)
            disc = (-1) ** (d * (d - 1) // 2) * rnd % q
            g2 = mu * (-tau[2] / math.sqrt(q)) ** d * F.chi(rnd, 2)
            beta = g2 * (eps0.conjugate()) ** d
            ratio = beta / F.chi(disc, 1)
            C_vals.add((d, round(ratio.real, 6), round(ratio.imag, 6)))
            if len(fac) == 1 and abs(g2 ** 3 + eps0 ** d) > 1e-9:
                bad["G5"] += 1
    print(f"q={q}: {cnt} squarefree n (deg<= {maxd}); failures {bad}")
    print("   eps0 =", eps0, "|eps0| =", abs(eps0), " (gamma_2(P)^3 = -eps0^deg P checked as G5)")
    byd = {}
    for d, re, im in sorted(C_vals):
        byd.setdefault(d, set()).add((re, im))
    for d in byd:
        vals = byd[d]
        print(f"   deg {d}: beta(n)/chi6(disc n) takes {len(vals)} value(s): {sorted(vals)[:3]}")
    C = complex(*sorted(byd[1])[0])
    print("   C =", C, " check C^d:", [(d, complex(*sorted(byd[d])[0]), C ** d) for d in byd])


if __name__ == "__main__":
    run(7, 3)
    run(13, 2)
