"""cB-numerics shared library: loading the C output, cubic symbols mod h, coefficient families.

Conventions as in cubic_data.c (= alt-detect/conductor_test.py):
  primary n = A + B w, A = 1 mod 3, B = 0 mod 3;  G~(n) = g(1,n)/sqrt(N n), e(Tr(x/n)).
  chi_h(n) := (n/h)_3 = (h/n)_3 (cubic reciprocity, both primary; checked in verify.py).
"""
import numpy as np, math, os

W3 = np.exp(2j * np.pi / 3)
SQ3 = math.sqrt(3.0)
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

PRIME_DT = np.dtype([("a", "<i8"), ("b", "<i8"), ("N", "<i8"), ("kind", "<i4"), ("p", "<i4"),
                     ("r", "<i8"), ("gre", "<f8"), ("gim", "<f8")])
SF_DT = np.dtype([("A", "<i4"), ("B", "<i4"), ("N", "<i4"), ("mu", "<i4"), ("gre", "<f8"), ("gim", "<f8")])


def load(prefix):
    P = np.fromfile(prefix + "_primes.bin", dtype=PRIME_DT)
    S = np.fromfile(prefix + "_sf.bin", dtype=SF_DT)
    o = np.argsort(S["N"], kind="stable")
    return P, S[o]


def cval(a, b):
    return (np.asarray(a, float) - np.asarray(b, float) / 2) + 1j * np.asarray(b, float) * SQ3 / 2


def mulw(a1, b1, a2, b2):
    return a1 * a2 - b1 * b2, a1 * b2 + a2 * b1 - b1 * b2


def norm(a, b):
    return a * a - a * b + b * b


def ffmul(u1, v1, u2, v2, q):
    return (u1 * u2 - v1 * v2) % q, (u1 * v2 + u2 * v1 - v1 * v2) % q


def ffpow(u, v, e, q):
    ru, rv = 1, 0
    while e:
        if e & 1:
            ru, rv = ffmul(ru, rv, u, v, q)
        u, v = ffmul(u, v, u, v, q)
        e >>= 1
    return ru, rv


class Mod:
    """cubic character x -> (x / pi)_3 for a primary prime pi (row of the primes table).
    __call__(A, B) returns exponent j in {0,1,2} or 3 if pi | x (vectorised)."""

    def __init__(self, row):
        self.kind = int(row["kind"]); self.p = int(row["p"]); self.r = int(row["r"])
        self.a = int(row["a"]); self.b = int(row["b"]); self.N = int(row["N"])
        p = self.p
        if self.kind == 0:
            x = np.arange(p, dtype=np.int64)
            e = (p - 1) // 3
            v = np.ones(p, dtype=np.int64); base = x.copy(); ee = e
            while ee:
                if ee & 1:
                    v = v * base % p
                base = base * base % p; ee >>= 1
            lut = np.full(p, 3, np.int8)
            lut[v == 1] = 0; lut[v == self.r % p] = 1; lut[v == (self.r * self.r) % p] = 2
            lut[0] = 3
            self.lut = lut
        else:
            q = p
            u, w = np.meshgrid(np.arange(q, dtype=np.int64), np.arange(q, dtype=np.int64), indexing="ij")
            u = u.ravel(); w = w.ravel()
            e = (q * q - 1) // 3
            ru, rv = np.ones_like(u), np.zeros_like(w); bu, bv = u.copy(), w.copy()
            while e:
                if e & 1:
                    ru, rv = (ru * bu - rv * bv) % q, (ru * bv + rv * bu - rv * bv) % q
                bu, bv = (bu * bu - bv * bv) % q, (2 * bu * bv - bv * bv) % q
                e >>= 1
            lut = np.full(q * q, 3, np.int8)
            lut[(ru == 1) & (rv == 0)] = 0; lut[(ru == 0) & (rv == 1)] = 1
            lut[(ru == q - 1) & (rv == q - 1)] = 2
            lut[0] = 3
            self.lut = lut

    def __call__(self, A, B):
        p = self.p
        if self.kind == 0:
            return self.lut[(A.astype(np.int64) % p + (B.astype(np.int64) % p) * self.r) % p]
        return self.lut[(A.astype(np.int64) % p) * p + (B.astype(np.int64) % p)]


def chi_exp(mods, A, B):
    """exponent of chi_h(n) = (n/h)_3 for h = prod of the given primes; 3+ means not coprime."""
    j = np.zeros(len(A), np.int64); bad = np.zeros(len(A), bool)
    for m in mods:
        t = m(A, B).astype(np.int64)
        bad |= t == 3; j += t
    j %= 3
    j[bad] = -1
    return j


def nu_char(A, B, u, v):
    """ray-class character mod 9 on primary n = 1 + 3z: w^{u (A-1)/3 + v B/3}."""
    z1 = ((A.astype(np.int64) - 1) // 3) % 3; z2 = (B.astype(np.int64) // 3) % 3
    return W3 ** ((u * z1 + v * z2) % 3)
