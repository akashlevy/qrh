"""
edge-theta: shared arithmetic for the truth-check numerics (Z[w], w = exp(2 pi i/3)).

Vectorised sextic residue symbols via per-prime lookup tables, theta coefficients
beta(n) = gamma_2(n) conj(alpha(n)) on squarefree primary n (twisted multiplicativity),
sextic Gauss sums gamma_1, the ray-class factor G(n) of the identity
gamma_1 gamma_2 = mu alpha G (P78 eq. signal; wall-spectral identity_check.py), and the
Moebius partner a(n) = mu(n) conj G(n) = conj(gamma_1(n) beta(n)).

Re-uses energy-wall/eisenstein.py and energy-wall/dual_energy.py (read only).
"""
import sys, math, cmath
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/energy-wall")
from eisenstein import eisenstein_primes, ResidueField, cx, mul, norm, Z6, primary   # noqa
import dual_energy as de                                                              # noqa

W = complex(-0.5, math.sqrt(3) / 2)


def setup(maxnorm):
    primes = eisenstein_primes(maxnorm + 10)
    fields = {p[1]: ResidueField(*p) for p in primes}
    pnorm = {p[1]: (p[2] if p[0] == "split" else p[2] ** 2) for p in primes}
    plist = sorted(fields, key=lambda p: pnorm[p])
    return fields, pnorm, plist


class Symb:
    """sextic exponent lookup for (x / pi)_6, x = A + B w given as int arrays; -1 if pi | x."""
    def __init__(self, F):
        self.F = F
        q = F.q
        if F.kind == "split":
            x = np.arange(q, dtype=np.int64)
            v = de.powmod_vec(x, (q - 1) // 6, q)
            w6 = (-(F.r * F.r)) % q
            lut = np.full(q, -1, dtype=np.int8); t = 1
            for k in range(6):
                lut[v == t] = k; t = t * w6 % q
            lut[0] = -1
            self.lut = lut
        else:
            # F_{q^2} = F_q[w]; table over (x mod q, y mod q)
            xs, ys = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
            xs = xs.ravel().astype(np.int64); ys = ys.ravel().astype(np.int64)
            e = (q * q - 1) // 6
            rx, ry = np.ones_like(xs), np.zeros_like(ys)
            bx, by = xs.copy(), ys.copy()
            while e:
                if e & 1:
                    rx, ry = (rx * bx - ry * by) % q, (rx * by + ry * bx - ry * by) % q
                bx, by = (bx * bx - by * by) % q, (2 * bx * by - by * by) % q
                e >>= 1
            lut = np.full(q * q, -1, dtype=np.int8)
            tx, ty = 1, 0
            for k in range(6):
                lut[(rx == tx) & (ry == ty)] = k
                tx, ty = (tx * 1 - ty * 1) % q, (tx * 1 + ty * 1 - ty * 1) % q   # times (1 + w)
            lut[0] = -1
            self.lut = lut

    def __call__(self, A, B):
        F = self.F; q = F.q
        if F.kind == "split":
            return self.lut[(A + B * F.r) % q]
        return self.lut[(A % q) * q + (B % q)]


def gauss_vec(F, j):
    """gamma_j(pi) = N(pi)^{-1/2} sum_x chi_pi(x)^j e(x/pi), e(z) = exp(4 pi i Im z / sqrt3)."""
    if j == 2:
        return de.gauss2_vec(F)
    pc = cx(F.p)
    if F.kind == "split":
        t = np.arange(1, F.q, dtype=np.int64)
        k = de.sextic_exp_vec(F, t, np.zeros_like(t)).astype(np.int64)
        z = t / pc
        return np.sum(np.exp(2j * np.pi * j * k / 6) * np.exp(4j * np.pi * z.imag / math.sqrt(3))) / math.sqrt(F.q)
    q = F.q
    xs, ys = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    xs = xs.ravel(); ys = ys.ravel(); keep = (xs != 0) | (ys != 0)
    xs, ys = xs[keep], ys[keep]
    k = de.sextic_exp_vec(F, xs, ys).astype(np.int64)
    z = (xs + ys * W) / pc
    return np.sum(np.exp(2j * np.pi * j * k / 6) * np.exp(4j * np.pi * z.imag / math.sqrt(3))) / q


def Gamma(c):
    a, b = c
    return (1 + 1j ** (-b % 4) + 1j ** (a % 4) + 1j ** ((b - a) % 4)) / 2


class Columns:
    """squarefree primary n with lo < N(n) <= hi; theta beta(n), Moebius partner a(n), etc."""
    def __init__(self, lo, hi, fields, pnorm, plist, symb):
        A, B, Nn, fac = de.primary_sf_elements(lo, hi, fields, plist, pnorm)
        odd = np.array([all(pnorm[p] != 4 for p in f) for f in fac], bool)   # (n, 6) = 1
        keep = np.nonzero(odd)[0]
        A, B, Nn = A[keep], B[keep], Nn[keep]; fac = [fac[i] for i in keep]
        self.A, self.B, self.N, self.fac = A, B, Nn, fac
        self.n = len(A)
        self.mu = np.array([(-1) ** len(f) for f in fac], float)
        need = sorted({p for f in fac for p in f}, key=lambda p: pnorm[p])
        self.primes = need
        g2 = {p: gauss_vec(fields[p], 2) for p in need}
        beta = np.empty(self.n, complex)
        chi4 = np.empty(self.n, complex)
        for i, f in enumerate(fac):
            g = 1 + 0j; acc = (1, 0); c4 = 1 + 0j
            for p in f:
                k = fields[p].sextic_exp(acc) if acc != (1, 0) else 0
                g *= g2[p] * Z6[(4 * k) % 6]; acc = mul(acc, p)
                k4 = fields[p].sextic_exp((4, 0))
                c4 *= Z6[k4]
            nc = complex(A[i], 0) + B[i] * W
            beta[i] = g * np.conj(nc) / abs(nc)
            chi4[i] = c4
        self.beta = beta
        self.alpha = (A + B * W) / np.sqrt(Nn)
        self.G = np.conj(chi4) * np.array([Gamma((int(a), int(b))) for a, b in zip(A, B)])
        # Moebius partner: a(n) = mu(n) conj G(n)  (= conj(gamma_1 beta) by gamma_1 gamma_2 = mu alpha G)
        self.amob = self.mu * np.conj(self.G)
        self.w = np.sin(np.pi * (Nn / lo - 1)) ** 2 if lo > 0 else np.ones(self.n)


def primary_all(H):
    """all primary h (a = 1 mod 3, b = 0 mod 3) with 0 < N(h) <= H, sorted by norm."""
    R = int(2 * math.sqrt(H)) + 3
    a = np.arange(-R, R + 1)
    A, B = np.meshgrid(a, a, indexing="ij"); A = A.ravel(); B = B.ravel()
    sel = (A % 3 == 1) & (B % 3 == 0)
    A, B = A[sel], B[sel]
    Nn = A * A - A * B + B * B
    sel = (Nn > 0) & (Nn <= H)
    A, B, Nn = A[sel].astype(np.int64), B[sel].astype(np.int64), Nn[sel].astype(np.int64)
    o = np.argsort(Nn, kind="stable")
    return A[o], B[o], Nn[o]


def check_identity(cols, fields, nmax=60):
    """verify gamma_1 gamma_2 = mu alpha G at the first primes, i.e. G(p) = -gamma_1 gamma_2/alpha."""
    bad = 0
    for p in cols.primes[:nmax]:
        F = fields[p]
        g1, g2 = gauss_vec(F, 1), gauss_vec(F, 2)
        al = cx(p) / abs(cx(p))
        k4 = F.sextic_exp((4, 0))
        G = np.conj(Z6[k4]) * Gamma(p)
        bad += abs(g1 * g2 + al * G) > 1e-6
    return bad
