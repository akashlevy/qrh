"""
Arithmetic in O = Z[omega] for numerical experiments on the sextic additive polynomial
A_m(Y) = Y^{-1} sum_s W(N s/Y) (N s/Y)^{-1/2} gtilde(s,-m)   (s primary, squarefree, (s,6)=1)

Conventions follow the 7/8 paper:
  e(z) = exp(4 pi i Im(z)/sqrt3);  primary = congruent to 1 mod 3;
  chi_s(x) = (x/s)_6, the sextic residue symbol, = zeta6^j with zeta6 = 1+omega = e^{i pi/3};
  g(s,h) = sum_{d mod s} chi_s(d) e(h d/s).
Elements are pairs (a,b) <-> a + b*omega, omega^2 = -1 - omega.
"""
import numpy as np
from math import isqrt

SQ3 = np.sqrt(3.0)
Z6 = np.exp(1j * np.pi / 3 * np.arange(6))          # zeta6^j

def norm(a, b):
    return a * a - a * b + b * b

def cval(a, b):
    return (a - b / 2.0) + 1j * (b * SQ3 / 2.0)

def mul(a1, b1, a2, b2):
    # (a1 + b1 w)(a2 + b2 w) = a1a2 + (a1b2 + a2b1) w + b1b2 w^2,  w^2 = -1 - w
    return a1 * a2 - b1 * b2, a1 * b2 + a2 * b1 - b1 * b2

def primary(a, b):
    """unique associate congruent to 1 mod 3 (a ≡ 1, b ≡ 0 mod 3); assumes coprime to 3."""
    units = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]   # 1,-1,w,-w,1+w=-w^2,-(1+w)=w^2
    for (ua, ub) in units:
        x, y = mul(a, b, ua, ub)
        if x % 3 == 1 and y % 3 == 0:
            return x, y
    raise ValueError("not coprime to 3: %d %d" % (a, b))

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def find_split(p):
    """primary pi = a + b w with norm p (p ≡ 1 mod 3)."""
    # solve a^2 - ab + b^2 = p by search over b
    for b in range(1, isqrt(4 * p // 3) + 2):
        disc = 4 * p - 3 * b * b
        if disc < 0: break
        r = isqrt(disc)
        if r * r == disc:
            for a in ((b + r) // 2, (b - r) // 2):
                if norm(a, b) == p:
                    return primary(a, b)
    raise ValueError(p)

class Prime:
    """a primary prime pi of O outside S = {lambda, 2}"""
    def __init__(self, a, b, N, kind, q):
        self.a, self.b, self.N, self.kind, self.q = a, b, N, kind, q
        self.e6 = (N - 1) // 6
        if kind == "split":
            p = q
            # omega ≡ r mod pi:  a + b r ≡ 0 mod p
            self.r = (-a * pow(b, -1, p)) % p
            assert (self.r * self.r + self.r + 1) % p == 0
            z = (1 + self.r) % p            # zeta6 = 1 + omega
            self.roots = [pow(z, j, p) for j in range(6)]
            assert len(set(self.roots)) == 6
        else:
            # O/q = F_q[w]/(w^2+w+1); zeta6 = 1 + w
            self.roots = []
            x = (1, 0)
            for j in range(6):
                self.roots.append(x)
                x = self._fmul(x[0], x[1], 1, 1)
            assert len(set(self.roots)) == 6

    def _fmul(self, a1, b1, a2, b2):
        q = self.q
        x, y = mul(a1, b1, a2, b2)
        return x % q, y % q

    def symbol(self, xa, xb):
        """exponent j in Z/6 with (x/pi)_6 = zeta6^j; -1 where pi | x. Vectorized over int arrays."""
        xa = np.asarray(xa, np.int64); xb = np.asarray(xb, np.int64)
        if self.kind == "split":
            p = self.q
            v = (xa + xb * self.r) % p
            out = np.full(v.shape, -1, np.int64)
            nz = v != 0
            y = _powmod(v[nz], self.e6, p)
            res = np.full(y.shape, -1, np.int64)
            for j, rt in enumerate(self.roots):
                res[y == rt] = j
            assert (res >= 0).all()
            out[nz] = res
            return out
        else:
            q = self.q
            va = xa % q; vb = xb % q
            out = np.full(va.shape, -1, np.int64)
            nz = (va != 0) | (vb != 0)
            ya, yb = _powF(va[nz], vb[nz], self.e6, q)
            res = np.full(ya.shape, -1, np.int64)
            for j, (ra, rb) in enumerate(self.roots):
                res[(ya == ra) & (yb == rb)] = j
            assert (res >= 0).all()
            out[nz] = res
            return out

    def gauss(self):
        """g(pi,1) = sum_x chi(x) e(x/pi)."""
        if self.kind == "split":
            p = self.q
            x = np.arange(1, p, dtype=np.int64)
            j = self.symbol(x, np.zeros_like(x))
            # e(x/pi) = exp(-2 pi i b x / p)
            ph = np.exp(-2j * np.pi * ((self.b * x) % p) / p)
            return np.sum(Z6[j] * ph)
        else:
            q = self.q
            c, d = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
            c = c.ravel(); d = d.ravel()
            nz = (c != 0) | (d != 0)
            c = c[nz]; d = d[nz]
            j = self.symbol(c, d)
            # pi = -q:  Im(x/pi) = -d sqrt3/(2q)  -> e = exp(-2 pi i d/q)
            sgn = -1 if self.a < 0 else 1
            ph = np.exp(-2j * np.pi * sgn * d / q)
            return np.sum(Z6[j] * ph)

def _powmod(v, e, p):
    r = np.ones_like(v); b = v % p
    while e:
        if e & 1: r = (r * b) % p
        b = (b * b) % p; e >>= 1
    return r

def _powF(va, vb, e, q):
    ra = np.ones_like(va); rb = np.zeros_like(vb)
    ba, bb = va % q, vb % q
    while e:
        if e & 1:
            ra, rb = (ra * ba - rb * bb) % q, (ra * bb + rb * ba - rb * bb) % q
        ba, bb = (ba * ba - bb * bb) % q, (2 * ba * bb - bb * bb) % q
        e >>= 1
    return ra, rb

def primes_upto(Nmax):
    """primary primes of O with norm <= Nmax, excluding lambda and 2."""
    out = []
    for p in sieve(Nmax):
        p = int(p)
        if p == 3 or p == 2: continue
        if p % 3 == 1:
            a, b = find_split(p)
            out.append(Prime(a, b, p, "split", p))
            a2, b2 = primary(a - b, -b)   # conjugate of a+bw is a + b w^2 = (a-b) - b w
            out.append(Prime(a2, b2, p, "split", p))
        elif p * p <= Nmax:
            a, b = primary(p, 0)
            out.append(Prime(a, b, p * p, "inert", p))
    out.sort(key=lambda P: P.N)
    return out

def squarefree_products(primes, Nmax):
    """all squarefree primary s = product of distinct primes, N(s) <= Nmax. Returns list of (a,b,N,idx tuple)."""
    res = [(1, 0, 1, ())]
    def rec(start, a, b, N, idx):
        for i in range(start, len(primes)):
            P = primes[i]
            if N * P.N > Nmax: break
            na, nb = mul(a, b, P.a, P.b)
            res.append((na, nb, N * P.N, idx + (i,)))
            rec(i + 1, na, nb, N * P.N, idx + (i,))
    rec(0, 1, 0, 1, ())
    return res

def elements_in_annulus(lo, hi):
    """all m = a + b w with lo < N(m) <= hi and gcd(N(m),6)=1."""
    R = isqrt(int(4 * hi / 3)) + 2
    a, b = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
    a = a.ravel(); b = b.ravel()
    N = norm(a, b)
    keep = (N > lo) & (N <= hi) & (N % 2 != 0) & (N % 3 != 0)
    return a[keep].astype(np.int64), b[keep].astype(np.int64), N[keep]
