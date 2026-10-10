"""
Cheap Z[omega] sanity checks (cB-theory).  NUMERICAL, exact up to float rounding.
 (A) p-parts of SEXTIC Gauss sums: |g(p^l, p^k)| / p^{k/2} for l = 0..5, k = 0..6 at split p = 7, 13.
     Expect: k=0 -> 1; 0<k<=l with 6|k -> p^{k/2}(1-1/p); k = l+1 -> p^{l/2}; else 0.
 (B) Orientation of B_h.  P78 notation: chi_n = (./n)_6, gamma_j(n;a) = N^{-1/2} sum_v chi_n(v)^j e(av/n).
     Task's B_h coefficient  mu(n) g~(h,n) = mu(n) gamma_2(n;h)  (cubic Gauss sum, numerator h).
     Claim: R4(n,h) = mu(n) gamma_2(n;h) / conj(gamma_1(n;h^4))  does not depend on h  (sextic index h^4),
            R2(n,h) = mu(n) gamma_2(n;h) / conj(gamma_1(n;h^2))  does depend on h      (index h^2 = other orientation).
     Also R4(n)/(R4(p1)R4(p2)) must be a ray quantity (quadratic reciprocity sign).
"""
import numpy as np
from math import isqrt
from itertools import combinations

Z6 = np.exp(2j * np.pi / 6)


def norm(a, b): return a * a - a * b + b * b
def mul(a1, b1, a2, b2): return a1 * a2 - b1 * b2, a1 * b2 + a2 * b1 - b1 * b2


def primary(a, b):
    for (ua, ub) in [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]:
        x, y = mul(a, b, ua, ub)
        if x % 3 == 2 and y % 3 == 0:   # P78-style primary: = -1 mod 3 or 1 mod 3; either is a fixed choice
            return x, y
    raise ValueError


def split_primes(pmax):
    out = []
    for p in range(7, pmax + 1, 6):
        if any(p % d == 0 for d in range(2, isqrt(p) + 1)): continue
        for b in range(1, p):
            d = 4 * p - 3 * b * b
            if d < 0: break
            r = isqrt(d)
            if r * r == d and norm((b + r) // 2, b) == p:
                a = (b + r) // 2; break
        a, b = primary(a, b)
        for (x, y) in {(a, b), primary(a - b, -b)}:
            r = (-x * pow(y, -1, p)) % p   # omega mod pi:  x + y w = 0 mod pi
            out.append(dict(a=x, b=y, p=p, r=r))
    return out


def hensel_omega(r, p, k):
    """lift root of w^2 + w + 1 from mod p to mod p^k."""
    m = p
    for _ in range(k):
        m2 = m * p
        f = (r * r + r + 1) % m2
        r = (r - f * pow(2 * r + 1, -1, m2)) % m2
        m = m2
    return r % p ** k


def sext_exp(x, P):
    """exponent j with x^{(p-1)/6} = (1+omega)^j mod pi; -1 if pi | x."""
    p, r = P['p'], P['r']
    x = x % p
    if x == 0: return -1
    v = pow(int(x), (p - 1) // 6, p)
    z = 1
    for j in range(6):
        if z == v: return j
        z = z * (1 + r) % p
    raise ValueError


# ---------------- (A) p-parts of sextic Gauss sums ----------------
print("(A) |g(p^l,p^k)| / p^{k/2},   sextic, split p (rows l = 0..5, columns k = 0..6)")
for P in [next(P for P in split_primes(13) if P['p'] == pp) for pp in (7, 13)]:
    p = P['p']
    tab = np.zeros((6, 7))
    exps = np.array([sext_exp(x, P) for x in range(p)])
    for k in range(7):
        if k == 0:
            tab[:, 0] = 1; continue
        M = p ** k
        x = np.arange(M)
        j = exps[x % p]
        ok = j >= 0
        ch = np.where(ok, Z6 ** ((k * j) % 6), 0)
        # primitive additive character of O/pi^k = Z/p^k: x -> exp(2 pi i x t / p^k), t a unit
        for l in range(6):
            num = p ** l
            g = np.sum(ch * np.exp(2j * np.pi * ((x * num) % M) / M))
            tab[l, k] = abs(g) / p ** (k / 2)
    print(f"  p = {p}:")
    for l in range(6):
        print("   l=%d: " % l + " ".join(f"{v:9.3f}" for v in tab[l]))
    print(f"   (p^0.5={p**0.5:.3f}, p={p}, p^1.5={p**1.5:.3f}, p^2={p**2}, p^2.5={p**2.5:.1f}, p^3(1-1/p)={p**3*(1-1/p):.1f})")

# ---------------- (B) orientation of B_h ----------------
print("\n(B) orientation:  R4 = mu*gamma_2(n;h)/conj gamma_1(n;h^4),  R2 = ... h^2  (several h per n)")
primes = split_primes(61)
small = primes[:10]


def gamma_comp(Ps, j, a_mod):
    """gamma_j(n; a) for n = prod of distinct split primes Ps, by direct sum over O/n = Z/N (CRT)."""
    N = 1
    for P in Ps: N *= P['p']
    # element n of Z[omega] and its residue: O/n ~ Z/N with omega -> rho (CRT of the r's)
    from sympy.ntheory.modular import crt
    rho = int(crt([P['p'] for P in Ps], [P['r'] for P in Ps])[0]) if len(Ps) > 1 else Ps[0]['r']
    na, nb = 1, 0
    for P in Ps: na, nb = mul(na, nb, P['a'], P['b'])
    x = np.arange(N)
    jj = np.zeros(N, dtype=int); ok = np.ones(N, bool)
    for P in Ps:
        e = np.array([sext_exp(t, P) for t in range(P['p'])])[x % P['p']]
        ok &= e >= 0; jj += np.where(e >= 0, e, 0)
    ch = np.where(ok, Z6 ** ((j * jj) % 6), 0)
    # additive character e(a x / n): Tr-form; for x in Z, a in Z: Tr(a x conj(n)/N(n)) = a x (2 na - nb)/N
    tr = (2 * na - nb) % N
    return np.sum(ch * np.exp(2j * np.pi * ((x * (a_mod % N) * tr) % N) / N)) / np.sqrt(N)


def elt_mod(a, b, Ps):
    """residue in Z/N of a + b omega."""
    from sympy.ntheory.modular import crt
    N = 1
    for P in Ps: N *= P['p']
    rho = int(crt([P['p'] for P in Ps], [P['r'] for P in Ps])[0]) if len(Ps) > 1 else Ps[0]['r']
    return (a + b * rho) % N


hs = [(2, 1), (3, 1), (5, 2), (4, -1), (7, 3)]   # assorted h in Z[omega]
R4single = {}
for Ps in [[small[0]], [small[1]], [small[2]], [small[0], small[2]], [small[1], small[3]], [small[0], small[3], small[4]]]:
    mu = (-1) ** len(Ps)
    out4, out2 = [], []
    for (ha, hb) in hs:
        hm = elt_mod(ha, hb, Ps)
        if any(hm % P['p'] == 0 for P in Ps): continue
        N = np.prod([P['p'] for P in Ps])
        g2 = gamma_comp(Ps, 2, hm)
        g1h4 = gamma_comp(Ps, 1, pow(hm, 4, int(N)))
        g1h2 = gamma_comp(Ps, 1, pow(hm, 2, int(N)))
        out4.append(mu * g2 / np.conj(g1h4)); out2.append(mu * g2 / np.conj(g1h2))
    lab = "*".join(str(P['p']) for P in Ps)
    print(f"  n = {lab:10s} R4: " + " ".join(f"{np.angle(z)/(2*np.pi)*6:+.3f}" for z in out4)
          + "   | R2: " + " ".join(f"{np.angle(z)/(2*np.pi)*6:+.3f}" for z in out2) + "   (args in units of 2pi/6)")
    print(f"               |R4| = {np.round(np.abs(out4), 6)}")
