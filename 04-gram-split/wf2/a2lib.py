"""Shared helpers for the A2 / pi(m) = rho(m) zeta_m numerics (wf2/).

Conventions are those of gram_split/zomega.py and paper78 (Q(omega), primary = 1 mod 3,
chi_c(x) = (x/c)_6 with the denominator c, e(z) = exp(4 pi i Im z / sqrt 3),
gamma_j(c) = q_c^{-1/2} sum_{v mod c} chi_c(v)^j e(v/c)).

Data: D_numerics/data/A_Q{Q}_K12.npy, A[i, m] = A_m(Y_i), Y_i = 256 * 2^i, built from gtilde(s,-m),
so V_i(m) = Y_i^{1/3} A_m(Y_i) / What = rho(m) + w_i(m), rho = Res_{w=2/3} sum_s gtilde(s,-m) N(s)^{-w}
(squarefree s only, as in the data).  Noise w_i is (heuristically) independent across scales.
"""
import numpy as np, sys, os
GS = "/Users/akashlevy/Documents/riemann/gram_split"
sys.path.insert(0, GS)
from zomega import *          # noqa

LAM = (1, 2)                  # lambda = 1 + 2 omega
UNITS = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]


def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    uu = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (uu * (1 - uu)) + 4); return out


_r = np.linspace(1, 2, 200001)
WHAT = np.trapezoid(bump(_r, 1, 2) * _r ** (-5 / 6), _r)


# ---------------------------------------------------------------- Gauss sums at primes
def gamma_prime(P, j, sign=+1):
    """q^{-1/2} sum_{v mod p} chi_p(v)^j e(sign * v/p) for a primary prime P (split or inert)."""
    if P.kind == "split":
        p = P.q; x = np.arange(1, p, dtype=np.int64); jj = P.symbol(x, np.zeros_like(x))
        ph = np.exp(sign * -2j * np.pi * ((P.b * x) % p) / p)
        return np.sum(Z6[(j * jj) % 6] * ph) / np.sqrt(P.N)
    q = P.q
    c, d = np.meshgrid(np.arange(q), np.arange(q), indexing="ij"); c = c.ravel(); d = d.ravel()
    nz = (c != 0) | (d != 0); c = c[nz]; d = d[nz]; jj = P.symbol(c, d)
    # P = -q (primary).  x/P = -(c + d w)/q, Im = -d sqrt3/(2q)  ->  e(x/P) = exp(-2 pi i d / q)
    s = -1 if P.a < 0 else 1
    ph = np.exp(sign * 2j * np.pi * (-s) * d / q)
    return np.sum(Z6[(j * jj) % 6] * ph) / np.sqrt(P.N)


def sym1(P, a, b):
    """exponent of (x/P)_6 for a single element x = a + b w (-1 if P | x)."""
    return int(P.symbol(np.array([a], np.int64), np.array([b], np.int64))[0])


# ---------------------------------------------------------------- brute-force Gauss sum mod R (validation)
def gauss_bruteforce(Ra, Rb, primes_of_R, j):
    """gamma_j(R) for squarefree primary R = Ra + Rb w, by direct summation over O/R.
    chi_R(v) = prod_{p | R} chi_p(v)."""
    N = norm(Ra, Rb)
    g = np.gcd(Ra, Rb); g = abs(int(g)) if g else 1
    x, y = np.meshgrid(np.arange(N // g, dtype=np.int64), np.arange(g, dtype=np.int64), indexing="ij")
    x = x.ravel(); y = y.ravel()
    tot = np.zeros(len(x), np.int64); ok = np.ones(len(x), bool)
    for P in primes_of_R:
        s = P.symbol(x, y); ok &= s >= 0; tot += np.where(s >= 0, s, 0)
    # v / R = v * conj(R) / N;  conj(R) = (Ra - Rb) - Rb w
    ca, cb = mul(x, y, Ra - Rb, -Rb)
    ph = np.exp(2j * np.pi * (cb % N) / N)
    return np.sum(np.where(ok, Z6[(j * tot) % 6], 0) * ph) / np.sqrt(N)


# ---------------------------------------------------------------- rows
class Rows:
    """Rows m of the Q-data with factorisation into primary primes, unit, primary part R,
    and the exact row scalar zeta_m (squarefree rows)."""

    def __init__(self, Q, data="K12", scales=None):
        self.Q = Q
        if data == "K12":
            A = np.load(GS + "/D_numerics/data/A_Q%d_K12.npy" % Q)
            ma, mb, mN = np.load(GS + "/D_numerics/data/m_Q%d_K12.npy" % Q)
        else:
            A = np.load(GS + "/A_Q%d.npy" % Q); ma, mb, mN = np.load(GS + "/m_Q%d.npy" % Q)
        K = A.shape[0]
        self.Ys = 256.0 * 2 ** np.arange(K)
        self.Vall = A * self.Ys[:, None] ** (1 / 3) / WHAT
        self.ma, self.mb, self.mN = ma.astype(np.int64), mb.astype(np.int64), mN.astype(np.int64)
        self.Om = bump(mN / Q, .5, 1)
        self.n = len(ma)
        self.P = primes_upto(int(mN.max()))
        self.byq = {}
        for k, Pr in enumerate(self.P):
            self.byq.setdefault(Pr.q, []).append(k)
        self._factor()
        self._gauss()
        self._zeta()
        if scales is None:
            scales = list(range(4, K))
        self.set_scales(scales)

    # -- factorisation
    def _factor(self):
        n = self.n; Nmax = int(self.mN.max())
        spf = np.zeros(Nmax + 1, np.int64)
        for p in sieve(Nmax):
            p = int(p)
            spf[p::p][spf[p::p] == 0] = p
        self.fac = []; self.sf = np.zeros(n, bool); self.u = np.zeros(n, int)
        self.Ra = np.zeros(n, np.int64); self.Rb = np.zeros(n, np.int64)
        for t in range(n):
            a, b, N = int(self.ma[t]), int(self.mb[t]), int(self.mN[t])
            rem = N; fl = []; sq = True; ca, cb = a, b
            while rem > 1:
                p = int(spf[rem]); e = 0
                while rem % p == 0: rem //= p; e += 1
                for k in self.byq[p]:
                    Pr = self.P[k]
                    while True:
                        if Pr.kind == "split":
                            if (ca + cb * Pr.r) % p: break
                        else:
                            if ca % p or cb % p: break
                        # divide by Pr: x / Pr = x * conj(Pr) / N(Pr)
                        xa, xb = mul(ca, cb, Pr.a - Pr.b, -Pr.b)
                        assert xa % Pr.N == 0 and xb % Pr.N == 0
                        ca, cb = xa // Pr.N, xb // Pr.N
                        if k in fl: sq = False
                        fl.append(k)
            # (ca, cb) is now a unit
            assert norm(ca, cb) == 1, (a, b, ca, cb)
            self.u[t] = UNITS.index((ca, cb))
            Ra, Rb = 1, 0
            for k in fl: Ra, Rb = mul(Ra, Rb, self.P[k].a, self.P[k].b)
            assert Ra % 3 == 1 and Rb % 3 == 0
            ua, ub = UNITS[self.u[t]]
            assert mul(ua, ub, Ra, Rb) == (a, b)
            self.Ra[t], self.Rb[t] = Ra, Rb
            self.fac.append(sorted(fl)); self.sf[t] = sq
        self.nfac = np.array([len(f) for f in self.fac])
        # incidence matrix (rows x primes), squarefree rows only meaningful
        self.Inc = np.zeros((n, len(self.P)), np.int8)
        for t, fl in enumerate(self.fac):
            for k in fl: self.Inc[t, k] += 1

    # -- Gauss sums and symbols at the primes that divide rows
    def _gauss(self):
        used = sorted(set(k for f in self.fac for k in f))
        self.used = used
        self.g = {j: np.zeros(len(self.P), complex) for j in (1, 2, 3)}
        self.gm1 = np.zeros(len(self.P), complex)          # tau^-_{p,1}
        for k in used:
            Pr = self.P[k]
            for j in (1, 2, 3):
                self.g[j][k] = gamma_prime(Pr, j)
            self.gm1[k] = gamma_prime(Pr, 1, sign=-1)
        self.e_lam = np.zeros(len(self.P), int); self.e_4 = np.zeros(len(self.P), int)
        self.e_m1 = np.zeros(len(self.P), int); self.e_2 = np.zeros(len(self.P), int)
        for k in used:
            Pr = self.P[k]
            self.e_lam[k] = sym1(Pr, *LAM); self.e_4[k] = sym1(Pr, 4, 0)
            self.e_m1[k] = sym1(Pr, -1, 0); self.e_2[k] = sym1(Pr, 2, 0)
        self.cross = {}
        for f in self.fac:
            for x in range(len(f)):
                for y in range(x + 1, len(f)):
                    i, k = f[x], f[y]
                    if (i, k) not in self.cross:
                        self.cross[(i, k)] = (sym1(self.P[i], self.P[k].a, self.P[k].b),
                                              sym1(self.P[k], self.P[i].a, self.P[i].b))

    def gamma_R(self, t, j):
        """gamma_j(R_t) by twisted multiplicativity: gamma_j(ab) = (chi_a(b) chi_b(a))^j gamma_j(a) gamma_j(b)."""
        f = self.fac[t]; g = 1.0 + 0j
        for x in range(len(f)):
            g *= self.g[j][f[x]]
            for y in range(x + 1, len(f)):
                e1, e2 = self.cross[(f[x], f[y])]
                g *= Z6[(j * (e1 + e2)) % 6]
        return g

    def chiR(self, t, e_arr):
        return Z6[int(sum(e_arr[k] for k in self.fac[t])) % 6]

    # -- the row scalar
    def _zeta(self):
        n = self.n
        self.alpha = np.array([cval(a, b) / abs(cval(a, b)) for a, b in zip(self.Ra, self.Rb)])   # alpha(R)
        self.alpha_m = np.array([cval(a, b) / abs(cval(a, b)) for a, b in zip(self.ma, self.mb)])
        self.G = {j: np.array([self.gamma_R(t, j) if self.sf[t] else 0 for t in range(n)]) for j in (1, 2, 3)}
        # Lemma Z closed form (c_F = 1, kappa_0 dropped):  psi_1(R) conj(gamma_2(R))^4,
        # psi_1(R) = chi_R(lambda)^{-1} conj(chi_R(4)) chi_R(-1)
        self.psi1 = np.array([np.conj(self.chiR(t, self.e_lam)) * np.conj(self.chiR(t, self.e_4)) * self.chiR(t, self.e_m1)
                              for t in range(n)])
        self.zetaZ = np.where(self.sf, self.psi1 * np.conj(self.G[2]) ** 4, 0)
        # paper's product formula (eq. reflection-row-phase), c_F = 1, all row primes active with j_p = 1,
        # constants -(i/81) phihat(h0) conj(kappa_F) dropped:
        #   zeta = conj(alpha(c))^2 prod_p chi_p(sigma_p)^{-2} tau^-_{p,1} tau^+_{p,3} chi_p(eps_p)^{-3},
        #   sigma_p = lambda^2 R/p,  eps_p = -((lambda^3 R/p) sigma_p)^{-1}
        zP = np.zeros(n, complex)
        for t in range(n):
            if not self.sf[t]: continue
            f = self.fac[t]; z = np.conj(self.alpha[t]) ** 2
            for k in f:
                Pr = self.P[k]
                ca, cb = 1, 0
                for k2 in f:
                    if k2 != k: ca, cb = mul(ca, cb, self.P[k2].a, self.P[k2].b)
                l2a, l2b = mul(LAM[0], LAM[1], LAM[0], LAM[1])
                sa, sb = mul(l2a, l2b, ca, cb)                       # sigma_p = lambda^2 R/p
                l3a, l3b = mul(l2a, l2b, LAM[0], LAM[1])
                xa, xb = mul(*mul(l3a, l3b, ca, cb), sa, sb)         # (lambda^3 R/p) sigma_p
                e_sig = sym1(Pr, sa, sb); e_x = sym1(Pr, xa, xb)
                # chi(eps) = chi(-1) * conj(chi(x))  ->  chi(eps)^{-3} = conj(chi(-1))^3 chi(x)^3
                e_eps_m3 = (-3 * self.e_m1[k] + 3 * e_x) % 6
                z *= Z6[(-2 * e_sig) % 6] * self.gm1[k] * self.g[3][k] * Z6[e_eps_m3]
            zP[t] = z
        self.zetaP = zP

    # -- scales / estimators
    def set_scales(self, scales):
        self.scales = list(scales)
        self.V = self.Vall[self.scales]
        C = ((self.V * self.Om) @ self.V.conj().T).real
        K = len(self.scales)
        off = np.array([np.mean([C[i, j] for j in range(K) if j != i]) for i in range(K)])
        self.sig2 = np.mean(C[~np.eye(K, dtype=bool)])          # ||rho||^2_Omega estimate
        self.noise = np.diag(C) - self.sig2
        self.w = 1 / self.noise; self.w /= self.w.sum()


def cross_est(T, w):
    """Unbiased estimate of |F|^2 from per-scale estimates T[..., i] = F + noise_i (independent noises):
    sum_{i != j} w_i w_j T_i conj(T_j) / sum_{i != j} w_i w_j.   T: array (..., K)."""
    S = (T * w).sum(-1)
    num = np.abs(S) ** 2 - (np.abs(T) ** 2 * w ** 2).sum(-1)
    return num / (w.sum() ** 2 - (w ** 2).sum())


def jack_cross(T, w):
    """estimate and leave-one-scale-out jackknife s.e."""
    K = T.shape[-1]
    est = cross_est(T, w)
    reps = []
    for i in range(K):
        keep = [j for j in range(K) if j != i]
        reps.append(cross_est(T[..., keep], w[keep]))
    reps = np.array(reps)
    se = np.sqrt((K - 1) / K * ((reps - reps.mean(0)) ** 2).sum(0))
    return est, se, reps
