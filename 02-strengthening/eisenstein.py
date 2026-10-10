"""
Arithmetic in Z[w], w = exp(2 pi i/3):  primes, sextic residue symbols, cubic Gauss sums.
Elements are pairs (a, b) <-> a + b w.  Norm N = a^2 - a b + b^2.
Conventions follow the 7/8 paper: primary = congruent to 1 mod 3;
e(z) = exp(4 pi i Im(z)/sqrt3);  chi_n(u) = (u/n)_6;  gamma_j(n) = N(n)^{-1/2} sum_x chi_n(x)^j e(x/n).
"""
import cmath, math
import numpy as np

W = complex(-0.5, math.sqrt(3) / 2)
Z6 = [cmath.exp(2j * math.pi * k / 6) for k in range(6)]

def cx(e): return e[0] + e[1] * W
def norm(e): a, b = e; return a * a - a * b + b * b
def mul(x, y):
    a, b = x; c, d = y           # (a+bw)(c+dw) = ac + (ad+bc)w + bd w^2,  w^2 = -1-w
    return (a * c - b * d, a * d + b * c - b * d)
UNITS = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]   # 1, w, w^2, -1, -w, -w^2
def primary(e):
    for u in UNITS:
        a, b = mul(u, e)
        if a % 3 == 1 and b % 3 == 0:
            return (a, b)
    return None

def is_prime_int(n):
    return n > 1 and all(n % i for i in range(2, int(n ** 0.5) + 1))

def eisenstein_primes(maxnorm):
    """primary primes pi (not above 3) with N(pi) <= maxnorm, with residue-field data."""
    out = []
    for q in range(2, maxnorm + 1):
        if not is_prime_int(q) or q == 3:
            continue
        if q % 3 == 1:
            # find a+bw of norm q
            found = None
            for a in range(0, int(2 * math.sqrt(q)) + 2):
                for b in range(-int(2 * math.sqrt(q)) - 2, int(2 * math.sqrt(q)) + 2):
                    if a * a - a * b + b * b == q:
                        found = (a, b); break
                if found: break
            p1 = primary(found)
            p2 = primary((found[0] - found[1], -found[1]))   # conjugate: a + b w-bar = (a-b) - b w
            for p in {p1, p2}:
                out.append(("split", p, q))
        elif q % 3 == 2 and q * q <= maxnorm:
            out.append(("inert", primary((q, 0)), q))
    return out

class ResidueField:
    """Z[w]/(pi) for a prime pi; maps elements to F_q (split) or F_{q^2} (inert)."""
    def __init__(self, kind, p, q):
        self.kind, self.p, self.q = kind, p, q
        if kind == "split":
            # find r with r^2 + r + 1 = 0 mod q and pi | (w - r): a + b r = 0 mod q
            a, b = p
            for r in range(q):
                if (r * r + r + 1) % q == 0 and (a + b * r) % q == 0:
                    self.r = r; break
            self.size = q
        else:
            self.size = q * q
    def red(self, e):
        a, b = e
        if self.kind == "split":
            return (a + b * self.r) % self.q
        return (a % self.q, b % self.q)
    # field ops for inert case: elements (x, y) <-> x + y w over F_q
    def fmul(self, x, y):
        q = self.q
        a, b = x; c, d = y
        return ((a * c - b * d) % q, (a * d + b * c - b * d) % q)
    def fpow(self, x, n):
        res = (1, 0); base = x
        while n:
            if n & 1: res = self.fmul(res, base)
            base = self.fmul(base, base); n >>= 1
        return res
    def sextic_exp(self, e):
        """k in 0..5 with (e/pi)_6 = exp(2 pi i k/6), or None if pi | e."""
        x = self.red(e)
        if self.kind == "split":
            if x == 0: return None
            v = pow(x, (self.q - 1) // 6, self.q)
            # identify v with w6^k where w6 = image of -w^2 = 1 + w (a primitive 6th root): -w^2 -> -(r^2)
            w6 = (-(self.r * self.r)) % self.q
            t = 1
            for k in range(6):
                if t == v: return k
                t = t * w6 % self.q
            raise ValueError("not a 6th root")
        else:
            if x == (0, 0): return None
            v = self.fpow(x, (self.size - 1) // 6)
            w6 = (1, 1)                      # 1 + w = -w^2, primitive 6th root
            t = (1, 0)
            for k in range(6):
                if t == v: return k
                t = self.fmul(t, w6)
            raise ValueError("not a 6th root")

def e_add(z):
    return cmath.exp(4j * math.pi * z.imag / math.sqrt(3))

def gauss_prime(F, j):
    """gamma_j(pi) by direct summation over residues."""
    p = F.p; pc = cx(p); q = F.q
    s = 0
    if F.kind == "split":
        # residues represented by integers t = 0..q-1 (since Z -> F_q surjective)
        for t in range(1, q):
            k = F.sextic_exp((t, 0))
            s += Z6[(j * k) % 6] * e_add(t / pc)
    else:
        for x in range(q):
            for y in range(q):
                if x == 0 and y == 0: continue
                k = F.sextic_exp((x, y))
                s += Z6[(j * k) % 6] * e_add(cx((x, y)) / pc)
    return s / math.sqrt(F.size)

if __name__ == "__main__":
    pr = eisenstein_primes(400)
    print("primes:", len(pr), pr[:5])
    bad = 0
    for kind, p, q in pr:
        F = ResidueField(kind, p, q)
        g1, g2, g3 = gauss_prime(F, 1), gauss_prime(F, 2), gauss_prime(F, 3)
        # paper identity (eq. gj): gamma_2^3 = mu * alpha = -alpha  for a prime
        alpha = cx(p) / abs(cx(p))
        bad += abs(g2 ** 3 + alpha) > 1e-8
        bad += abs(abs(g1) - 1) > 1e-8
    print("checks gamma_2(p)^3 = -alpha(p), |gamma_1|=1  failures:", bad)
