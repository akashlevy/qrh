"""
numerics: fast shared arithmetic for the energy-wall step-1 numerics (agent `numerics`).

Built on the edge-theta conventions (edge_lib.py), validated against them in conv_check.py:
  chi_n(x) = (x/n)_6, gamma_j(n) = N(n)^{-1/2} sum_x chi_n(x)^j e(x/n), e(z) = exp(4 pi i Im z/sqrt3),
  beta(n) = gamma_2(n) conj(alpha(n)) on squarefree primary n (P78 theta coefficient),
  psi_k(n) = chi_k(n)^3 = prod_{p|k} (n/p)_2 (quadratic rows), marks chi_P(n)^{-2} = conj (n/P)_3.
Exact data come from cubic_data.c (cB-numerics, copied): G~(n) = sum_x (x/n)_3 e(Tr(x/n))/sqrt(N n)
for every squarefree primary n, by twisted multiplicativity.  Since Tr(z) = Tr((lambda z)/sqrt(-3)) with
lambda = 1 + 2w = sqrt(-3), the two additive characters differ by x -> lambda x, so
      gamma_2(n) = (lambda/n)_3 * G~(n)        (checked in conv_check.py at primes and composites).
"""
import math, os, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W = complex(-0.5, math.sqrt(3) / 2)
PRIME_DT = np.dtype([("a", "<i8"), ("b", "<i8"), ("N", "<i8"), ("kind", "<i4"), ("p", "<i4"), ("r", "<i8"),
                     ("gre", "<f8"), ("gim", "<f8")])
SF_DT = np.dtype([("A", "<i4"), ("B", "<i4"), ("N", "<i4"), ("mu", "<i4"), ("gre", "<f8"), ("gim", "<f8")])


def load_cd(nmax):
    pre = os.path.join(HERE, "data", f"cd{nmax}")
    if not os.path.exists(pre + "_sf.bin"):
        subprocess.run([os.path.join(HERE, "cubic_data"), str(nmax), pre, "2"], check=True)
    P = np.fromfile(pre + "_primes.bin", dtype=PRIME_DT)
    S = np.fromfile(pre + "_sf.bin", dtype=SF_DT)
    return P, S


def powmod_vec(x, e, q):
    x = np.asarray(x, np.int64) % q
    res = np.ones_like(x); base = x.copy()
    while e:
        if e & 1: res = res * base % q
        base = base * base % q; e >>= 1
    return res


class Prime:
    """a primary prime pi of Z[w] (not above 2, 3): kind 0 split (O/pi = F_q via w -> r), 1 inert (q = p)."""
    __slots__ = ("a", "b", "N", "kind", "p", "r", "g2t")
    def __init__(self, rec):
        self.a, self.b, self.N = int(rec["a"]), int(rec["b"]), int(rec["N"])
        self.kind, self.p, self.r = int(rec["kind"]), int(rec["p"]), int(rec["r"])
        self.g2t = complex(rec["gre"], rec["gim"])      # cubic_data G~(pi)

    def red(self, A, B):
        """reduction of A + B w mod pi: split -> int in [0, q); inert -> (A mod q, B mod q)."""
        if self.kind == 0:
            return (np.asarray(A, np.int64) + np.asarray(B, np.int64) * self.r) % self.p
        return np.asarray(A, np.int64) % self.p, np.asarray(B, np.int64) % self.p

    def sextic_exp(self, A, B):
        """k in 0..5 with (x/pi)_6 = zeta6^k, zeta6 = 1 + w;  -1 where pi | x.  Vectorised."""
        A = np.atleast_1d(np.asarray(A, np.int64)); B = np.atleast_1d(np.asarray(B, np.int64))
        q = self.p
        out = np.full(A.shape, -1, np.int8)
        if self.kind == 0:
            x = self.red(A, B)
            v = powmod_vec(x, (q - 1) // 6, q)
            w6 = (-(self.r * self.r)) % q
            t = 1
            for k in range(6):
                out[(v == t) & (x != 0)] = k; t = t * w6 % q
            return out
        # inert: F_{q^2} = F_q[w]
        u, v = A % q, B % q
        e = (q * q - 1) // 6
        ru, rv = np.ones_like(u), np.zeros_like(v); bu, bv = u.copy(), v.copy()
        while e:
            if e & 1:
                ru, rv = (ru * bu - rv * bv) % q, (ru * bv + rv * bu - rv * bv) % q
            bu, bv = (bu * bu - bv * bv) % q, (2 * bu * bv - bv * bv) % q
            e >>= 1
        tu, tv = 1, 0
        for k in range(6):
            out[(ru == tu) & (rv == tv) & ((u != 0) | (v != 0))] = k
            tu, tv = (tu - tv) % q, (tu + tv - tv) % q        # times (1 + w)
        return out

    def quad(self, A, B, qr=None):
        """quadratic symbol (x/pi)_2 = (x/pi)_6^3 in {-1, 0, 1} (int8), vectorised; qr: optional QR table."""
        q = self.p
        if qr is None: qr = qr_table(q)
        if self.kind == 0:
            x = self.red(A, B)
        else:   # (x/q)_2 = x^{(q^2-1)/2} = N(x)^{(q-1)/2} mod q : Legendre of the norm
            A = np.asarray(A, np.int64); B = np.asarray(B, np.int64)
            x = (A * A - A * B + B * B) % q
        return qr[x]


def qr_table(q):
    """int8 table t[x] = Legendre (x/q), x in [0, q)."""
    t = -np.ones(q, np.int8)
    i = np.arange(1, (q - 1) // 2 + 1, dtype=np.int64)
    t[(i * i) % q] = 1
    t[0] = 0
    return t


def cubic_class_table(pr):
    """for a split prime pi of norm q: table c[x] in {0,1,2,-1} with (x/pi)_3 = w^c (x in F_q)."""
    q = pr.p
    x = np.arange(q, dtype=np.int64)
    v = powmod_vec(x, (q - 1) // 3, q)
    c = np.full(q, -1, np.int8)
    c[v == 1] = 0; c[v == pr.r % q] = 1; c[v == (pr.r * pr.r) % q] = 2
    c[0] = -1
    return c


# ---------------------------------------------------------------- columns, rows, factorisation
import ctypes, time
_lib = ctypes.CDLL(os.path.join(HERE, "libpsi.so"))
_i64p = np.ctypeslib.ndpointer(dtype=np.int64, flags="C_CONTIGUOUS")
_u8p = np.ctypeslib.ndpointer(dtype=np.uint8, flags="C_CONTIGUOUS")
_f64p = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
_lib.build_qr.argtypes = [ctypes.c_int64, _i64p, _i64p, _u8p]
_lib.build_psi.argtypes = [ctypes.c_int64, _i64p, _i64p, _i64p, ctypes.c_int64, _i64p, _i64p, _i64p, _i64p,
                           _u8p, _i64p, _f64p, ctypes.c_int]


class Arith:
    """prime table (cubic_data records) + factorisation of squarefree primary elements."""
    def __init__(self, nmax):
        P, S = load_cd(nmax)
        self.nmax = nmax
        self.P = [Prime(r) for r in P]
        self.S = S
        self.byrat = {}
        for i, pr in enumerate(self.P):
            self.byrat.setdefault(pr.p, []).append(i)
        # cubic exponent of (lambda/pi)_3, lambda = 1 + 2w
        self.lam3 = np.array([int(pr.sextic_exp(1, 2)[0]) % 3 for pr in self.P], np.int64)
        M = int(S["N"].max()) + 1
        spf = np.zeros(M, np.int64)
        for p in range(2, int(M ** 0.5) + 1):
            if spf[p] == 0:
                blk = spf[p * p::p]; blk[blk == 0] = p
        self.spf = spf

    def factor(self, A, B, N):
        """list of prime indices dividing each element (squarefree primary, prime to 6 assumed)."""
        out = []
        spf = self.spf
        for a, b, n in zip(A.tolist(), B.tolist(), N.tolist()):
            fac = []; m = n
            while m > 1:
                p = int(spf[m]) or m
                e = 0
                while m % p == 0: m //= p; e += 1
                ids = self.byrat[p]
                pr0 = self.P[ids[0]]
                if pr0.kind == 1:
                    fac.append(ids[0])
                elif e == 2:
                    fac.extend(ids)
                else:
                    fac.append(ids[0] if (a + b * pr0.r) % p == 0 else ids[1])
            out.append(fac)
        return out

    def elements(self, lo, hi):
        """squarefree primary odd elements with lo < N <= hi (cubic_data order -> sorted by N)."""
        S = self.S
        sel = (S["N"] > lo) & (S["N"] <= hi) & ~((S["A"] % 2 == 0) & (S["B"] % 2 == 0))
        idx = np.nonzero(sel)[0]
        o = np.argsort(S["N"][idx], kind="stable"); idx = idx[o]
        A = S["A"][idx].astype(np.int64); B = S["B"][idx].astype(np.int64); Nn = S["N"][idx].astype(np.int64)
        mu = S["mu"][idx].astype(np.int64)
        Gt = S["gre"][idx] + 1j * S["gim"][idx]
        fac = self.factor(A, B, Nn)
        return A, B, Nn, mu, Gt, fac


class Cols:
    """columns n sf primary odd, N0 < N(n) <= 2 N0: beta = gamma_2 conj(alpha), w = sin^2(pi(N/N0 - 1))."""
    def __init__(self, ar, N0):
        A, B, Nn, mu, Gt, fac = ar.elements(N0, 2 * N0)
        self.A, self.B, self.N, self.mu, self.fac = A, B, Nn, mu, fac
        self.n = len(A)
        lam = np.array([sum(int(ar.lam3[i]) for i in f) % 3 for f in fac], np.int64)
        self.gamma2 = np.exp(2j * np.pi * lam / 3) * Gt
        self.alpha = (A + B * W) / np.sqrt(Nn)
        self.beta = self.gamma2 * np.conj(self.alpha)
        self.w = np.sin(np.pi * (Nn / N0 - 1)) ** 2
        self.N0 = N0


class Rows:
    """rows k sf primary odd with N(k) <= KMAX (k = 1 included), sorted by norm, with factor lists."""
    def __init__(self, ar, KMAX):
        A, B, Nn, mu, Gt, fac = ar.elements(0, KMAX)
        self.A, self.B, self.N, self.mu, self.fac = A, B, Nn, mu, fac
        self.n = len(A)
        lam = np.array([sum(int(ar.lam3[i]) for i in f) % 3 for f in fac], np.int64)
        self.gamma2 = np.exp(2j * np.pi * lam / 3) * Gt
        self.alpha = (A + B * W) / np.sqrt(Nn)
        self.arg = np.angle(A + B * W)
        # psi factor lists: (q, r) for single split primes, (q, -1) for inert or conjugate pairs
        fq, fr = [], []; ptr = [0]
        for f in fac:
            seen = {}
            for i in f:
                pr = ar.P[i]
                seen.setdefault(pr.p, []).append(pr)
            for p, lst in seen.items():
                if len(lst) == 1 and lst[0].kind == 0:
                    fq.append(p); fr.append(lst[0].r)
                else:
                    fq.append(p); fr.append(-1)
            ptr.append(len(fq))
        self.fq = np.array(fq, np.int64); self.fr = np.array(fr, np.int64); self.ptr = np.array(ptr, np.int64)
        qs = np.unique(self.fq)
        self.qs = qs
        self.fqi = np.searchsorted(qs, self.fq).astype(np.int64)
        off = np.zeros(len(qs) + 1, np.int64)
        off[1:] = np.cumsum((qs >> 3) + 1)
        self.off = off
        self.bits = np.zeros(int(off[-1]), np.uint8)
        _lib.build_qr(len(qs), qs, off[:-1].copy(), self.bits)

    def psi_block(self, cols, i0, i1, nth=2):
        out = np.empty((i1 - i0, cols.n), np.float64)
        ptr = (self.ptr[i0:i1 + 1] - self.ptr[i0]).astype(np.int64)
        s0, s1 = self.ptr[i0], self.ptr[i1]
        _lib.build_psi(cols.n, cols.A, cols.B, cols.N, i1 - i0, ptr,
                       np.ascontiguousarray(self.fq[s0:s1]), np.ascontiguousarray(self.fr[s0:s1]),
                       np.ascontiguousarray(self.fqi[s0:s1]), self.bits, self.off[:-1].copy(), out, nth)
        return out


def rowsums(rows, cols, Ccoef, Dabs, nrows=None, block=384, nth=2, verbose=True):
    """S[k] = sum_n Ccoef[n] psi_k(n)  (complex, nrows x m);  D[k] = sum_n Dabs[n] 1_{(n,k)=1}  (nrows x md)."""
    t0 = time.time()
    nr = rows.n if nrows is None else nrows
    Cs = np.ascontiguousarray(np.hstack([Ccoef.real, Ccoef.imag]))
    Da = np.ascontiguousarray(Dabs)
    m = Ccoef.shape[1]
    S = np.empty((nr, m), complex); D = np.empty((nr, Da.shape[1]))
    for b0 in range(0, nr, block):
        b1 = min(nr, b0 + block)
        Psi = rows.psi_block(cols, b0, b1, nth)
        R = Psi @ Cs
        S[b0:b1] = R[:, :m] + 1j * R[:, m:]
        np.abs(Psi, out=Psi)
        D[b0:b1] = Psi @ Da
        if verbose and (b0 // block) % 20 == 0:
            print(f"   rows {b1}/{nr}  {time.time()-t0:.1f}s", flush=True)
    return S, D


def mark_chi_m2(pr, cols):
    """chi_P(n)^{-2} = conj (n/P)_3 on the columns (0 where P | n)."""
    k = pr.sextic_exp(cols.A, cols.B).astype(np.int64)
    return np.where(k < 0, 0, np.exp(-2j * np.pi * 2 * np.where(k < 0, 0, k) / 6))


def gamma_prime(pr, j):
    """gamma_j(pi) in the edge-theta convention, direct O(q) sum (split primes only)."""
    q = pr.p
    t = np.arange(1, q, dtype=np.int64)
    k = pr.sextic_exp(t, np.zeros_like(t)).astype(np.int64)
    # e(x/pi) = exp(4 pi i Im(x/pi)/sqrt3), x/pi = x conj(pi)/q
    pc = complex(pr.a, 0) + pr.b * W
    z = t * np.conj(pc) / q
    return np.sum(np.exp(2j * np.pi * j * k / 6) * np.exp(4j * np.pi * z.imag / math.sqrt(3))) / math.sqrt(q)
