"""
function-field: exact arithmetic in A = F_q[T], q prime, q = 1 (mod 6).

Polynomials are tuples of ints mod q, low degree first (p[i] = coeff of T^i).
Characters: fix a generator g of F_q^*, zeta6 := g^((q-1)/6) in F_q, and identify
mu_6(F_q) -> mu_6(C) by zeta6^j -> exp(2 pi i j/6).  chi6(c) := c^((q-1)/6) read in C,
so chi6(g) = e^{2 pi i/6}.  Residue symbols for monic irreducible P of degree d:
(a/P)_6 = a^((q^d-1)/6) mod P = chi6(Res(P,a)), extended multiplicatively in the
denominator, so (a/m)_6 = chi6(Res(m,a)) for monic m (Res(m,a) = prod_{m(t)=0} a(t)).
Additive character: e(f) = psi0(coefficient of T^{-1} in f), psi0(x) = exp(2 pi i x/q).
"""
import itertools, cmath, math
from functools import lru_cache


def primroot(q):
    for g in range(2, q):
        if all(pow(g, (q - 1) // r, q) != 1 for r in prime_factors(q - 1)):
            return g
    raise ValueError


def prime_factors(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


class Fq:
    def __init__(self, q):
        assert q % 6 == 1
        self.q = q
        self.g = primroot(q)
        self.dlog = [None] * q
        x = 1
        for j in range(q - 1):
            self.dlog[x] = j
            x = x * self.g % q
        self.z6 = [cmath.exp(2j * math.pi * t / 6) for t in range(6)]

    # sextic exponent of c in F_q^* (chi6(c) = z6[e6(c)]), None for c = 0
    def e6(self, c):
        c %= self.q
        return None if c == 0 else self.dlog[c] % 6

    def chi(self, c, j=1):
        e = self.e6(c)
        return 0 if e is None else self.z6[(j * e) % 6]

    def tau(self, j):
        """Gauss sum over F_q of chi6^j with psi0."""
        q = self.q
        return sum(self.chi(x, j) * cmath.exp(2j * math.pi * x / q) for x in range(1, q))


# ---------------- polynomial arithmetic over F_q ----------------
def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def deg(p):
    return len(p) - 1


def padd(a, b, q):
    n = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % q for i in range(n)])


def pmul(a, b, q):
    if not a or not b:
        return ()
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % q
    return trim(r)


def pdivmod(a, b, q):
    a = list(a)
    inv = pow(b[-1], q - 2, q)
    db = len(b) - 1
    quo = [0] * max(len(a) - db, 1)
    for i in range(len(a) - 1, db - 1, -1):
        c = a[i] * inv % q
        if c:
            quo[i - db] = c
            for j in range(db + 1):
                a[i - db + j] = (a[i - db + j] - c * b[j]) % q
    return trim(quo), trim(a[:db] if db > 0 else [])


def pmod(a, b, q):
    return pdivmod(a, b, q)[1]


def deriv(p, q):
    return trim([(i * p[i]) % q for i in range(1, len(p))])


def resultant(a, b, q):
    """Res(a,b) = lc(a)^deg b prod_{a(t)=0} b(t)."""
    a, b = trim(a), trim(b)
    if not a or not b:
        return 0
    res = 1
    while True:
        da, db = deg(a), deg(b)
        if db == 0:
            return res * pow(b[0], da, q) % q
        if da == 0:
            return res * pow(a[0], db, q) % q
        r = pmod(a, b, q)
        if not r:
            return 0
        # Res(a,b) = (-1)^{da db} lc(b)^{da - dr} Res(b, r)
        dr = deg(r)
        if (da * db) % 2:
            res = -res
        res = res * pow(b[-1], da - dr, q) % q
        a, b = b, r


def monic_polys(d, q):
    for c in itertools.product(range(q), repeat=d):
        yield tuple(c) + (1,)


def is_squarefree(p, q):
    return resultant(p, deriv(p, q), q) != 0 if deg(p) > 0 else True


@lru_cache(maxsize=None)
def irreducibles(d, q):
    """monic irreducibles of degree d (brute force by sieving products)."""
    if d == 0:
        return []
    allp = set(monic_polys(d, q))
    red = set()
    for e in range(1, d // 2 + 1):
        for a in monic_polys(e, q):
            for b in monic_polys(d - e, q):
                red.add(pmul(a, b, q))
    return sorted(p for p in allp if p not in red)


def factor(p, q, maxd=None):
    """factor monic p into monic irreducibles by trial division (small degrees)."""
    out = []
    rem = trim(p)
    d = 1
    while deg(rem) > 0:
        if 2 * d > deg(rem):
            out.append(rem)
            break
        found = False
        for P in irreducibles(d, q):
            qq, r = pdivmod(rem, P, q)
            while not r:
                out.append(P)
                rem = qq
                qq, r = pdivmod(rem, P, q)
                found = True
        d += 1
    return out


def laurent(a, k, q, nterms):
    """coefficients c_1..c_nterms of T^{-1},...,T^{-nterms} in a/k (deg a < deg k, k monic)."""
    # a/k = sum c_j T^{-j}; multiply: a = k * sum c_j T^{-j}  -> long division of a*T^nterms by k
    kk = len(k) - 1
    num = [0] * nterms + list(a) + [0] * (kk + 1)
    num = trim(num)
    quo, _ = pdivmod(num, k, q) if num else ((), ())
    quo = list(quo) + [0] * (nterms + 1)
    # quotient of a T^nterms / k = sum_j c_j T^{nterms - j}
    return [quo[nterms - j] % q for j in range(1, nterms + 1)]


def e_res(f_num, f_den, F, nterms=1):
    """e(f) for f = f_num/f_den with deg f_num < deg f_den."""
    c1 = laurent(pmod(f_num, f_den, F.q), f_den, F.q, 1)[0]
    return cmath.exp(2j * math.pi * c1 / F.q)


# ---------------- random irreducibles (Ben-Or test) ----------------
def pgcd(a, b, q):
    a, b = trim(a), trim(b)
    while b:
        a, b = b, pmod(a, b, q)
    if not a:
        return a
    inv = pow(a[-1], q - 2, q)
    return tuple(c * inv % q for c in a)


def ppowmod(base, e, m, q):
    result = (1,)
    base = pmod(base, m, q)
    while e:
        if e & 1:
            result = pmod(pmul(result, base, q), m, q)
        base = pmod(pmul(base, base, q), m, q)
        e >>= 1
    return result


def is_irreducible(P, q):
    d = deg(P)
    if d <= 0:
        return False
    if d == 1:
        return True
    x = (0, 1)
    xp = x
    for i in range(1, d // 2 + 1):
        xp = ppowmod(xp, q, P, q)           # T^{q^i} mod P
        g = pgcd(P, padd(xp, (0, q - 1), q), q)
        if deg(g) > 0:
            return False
    return True


def random_irreducibles(d, q, count, seed=1):
    import random as _r
    rnd = _r.Random(seed * 7919 + q * 31 + d)
    out, seen = [], set()
    tries = 0
    while len(out) < count and tries < 200000:
        tries += 1
        P = tuple(rnd.randrange(q) for _ in range(d)) + (1,)
        if P in seen:
            continue
        seen.add(P)
        if is_irreducible(P, q):
            out.append(P)
    return out
