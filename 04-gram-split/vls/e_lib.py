"""Attack E analysis helpers.  Rows for the E data (scales Y = 2^k), sparse incidence, quadratic-twist matrices,
cubic theta columns.  Builds on wf2/a2lib.py (exact zeta_m via Lemma Z, cross-scale estimators)."""
import numpy as np, sys, time
import scipy.sparse as sps
GS = "/Users/akashlevy/Documents/riemann/gram_split"
sys.path.insert(0, GS); sys.path.insert(0, GS + "/wf2")
import a2lib
from a2lib import *            # noqa  (zomega, bump, WHAT, cross_est, jack_cross, gamma_prime, sym1, UNITS, LAM)

DATA = GS + "/vls/data/"
KMIN = 8
UNITS6 = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]     # zeta6^c


class RowsE(a2lib.Rows):
    """Same as a2lib.Rows but for vls/data/A_Q{Q}_E{KMAX}.npy (rows: 6 associate blocks of the primary reps),
    sparse incidence, optional product-formula check on a subset only."""

    def __init__(self, Q, KMAX=21, ymin_over_Q=4.0, zeta_check=300, src="E"):
        self.Q = Q
        if src == "E":
            A = np.load(DATA + "A_Q%d_E%d.npy" % (Q, KMAX)); ma, mb, mN = np.load(DATA + "m_Q%d_E%d.npy" % (Q, KMAX))
            self.Ys = 2.0 ** (KMIN + np.arange(A.shape[0]))
        else:   # K12 data
            A = np.load(GS + "/D_numerics/data/A_Q%d_K12.npy" % Q); ma, mb, mN = np.load(GS + "/D_numerics/data/m_Q%d_K12.npy" % Q)
            self.Ys = 256.0 * 2 ** np.arange(A.shape[0])
        self.Vall = A * self.Ys[:, None] ** (1 / 3) / WHAT
        self.ma, self.mb, self.mN = ma.astype(np.int64), mb.astype(np.int64), mN.astype(np.int64)
        self.Om = bump(mN / Q, .5, 1)
        self.n = len(ma)
        self.P = primes_upto(int(mN.max()))
        self.byq = {}
        for k, Pr in enumerate(self.P):
            self.byq.setdefault(Pr.q, []).append(k)
        t0 = time.time()
        self._factor()
        self._gauss()
        self._zeta_fast(zeta_check)
        self.tinit = time.time() - t0
        scales = [i for i in range(len(self.Ys)) if self.Ys[i] >= ymin_over_Q * Q]
        self.set_scales(scales)
        # primary representative id (ideal) and unit class of each row
        self.rep_key = {}
        self.ideal = np.zeros(self.n, np.int64)
        keys = list(zip(self.Ra.tolist(), self.Rb.tolist()))
        for t, k in enumerate(keys):
            self.ideal[t] = self.rep_key.setdefault(k, len(self.rep_key))
        self.nideal = len(self.rep_key)
        # unit class c with m = zeta6^c R
        uidx = {u: c for c, u in enumerate(UNITS6)}
        self.ucls = np.array([uidx[UNITS[u]] for u in self.u])

    def _factor(self):
        n = self.n; Nmax = int(self.mN.max())
        spf = np.zeros(Nmax + 1, np.int64)
        for p in sieve(Nmax):
            p = int(p)
            spf[p::p][spf[p::p] == 0] = p
        self.fac = []; self.sf = np.zeros(n, bool); self.u = np.zeros(n, int)
        self.Ra = np.zeros(n, np.int64); self.Rb = np.zeros(n, np.int64)
        cache = {}
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
                        xa, xb = mul(ca, cb, Pr.a - Pr.b, -Pr.b)
                        ca, cb = xa // Pr.N, xb // Pr.N
                        if k in fl: sq = False
                        fl.append(k)
            self.u[t] = UNITS.index((ca, cb))
            Ra, Rb = 1, 0
            for k in fl: Ra, Rb = mul(Ra, Rb, self.P[k].a, self.P[k].b)
            self.Ra[t], self.Rb[t] = Ra, Rb
            self.fac.append(sorted(fl)); self.sf[t] = sq
        self.nfac = np.array([len(f) for f in self.fac])
        rows = np.repeat(np.arange(n), self.nfac); cols = np.concatenate([np.array(f, np.int64) for f in self.fac])
        self.IncS = sps.csr_matrix((np.ones(len(rows), np.int64), (rows, cols)), shape=(n, len(self.P)))

    def _zeta_fast(self, ncheck):
        n = self.n
        self.alpha = np.array([cval(a, b) / abs(cval(a, b)) for a, b in zip(self.Ra, self.Rb)])
        self.alpha_m = np.array([cval(a, b) / abs(cval(a, b)) for a, b in zip(self.ma, self.mb)])
        self.G = {j: np.array([self.gamma_R(t, j) if self.sf[t] else 0 for t in range(n)]) for j in (1, 2, 3)}
        self.psi1 = np.array([np.conj(self.chiR(t, self.e_lam)) * np.conj(self.chiR(t, self.e_4)) * self.chiR(t, self.e_m1)
                              for t in range(n)])
        self.zetaZ = np.where(self.sf, self.psi1 * np.conj(self.G[2]) ** 4, 0)
        # product-formula check on a random subset of squarefree rows
        rng = np.random.default_rng(1)
        sub = rng.choice(np.where(self.sf)[0], size=min(ncheck, self.sf.sum()), replace=False)
        err = 0.0
        for t in sub:
            f = self.fac[t]; z = np.conj(self.alpha[t]) ** 2
            for k in f:
                Pr = self.P[k]; ca, cb = 1, 0
                for k2 in f:
                    if k2 != k: ca, cb = mul(ca, cb, self.P[k2].a, self.P[k2].b)
                l2a, l2b = mul(LAM[0], LAM[1], LAM[0], LAM[1])
                sa, sb = mul(l2a, l2b, ca, cb)
                l3a, l3b = mul(l2a, l2b, LAM[0], LAM[1])
                xa, xb = mul(*mul(l3a, l3b, ca, cb), sa, sb)
                e_sig = sym1(Pr, sa, sb); e_x = sym1(Pr, xa, xb)
                e_eps_m3 = (-3 * self.e_m1[k] + 3 * e_x) % 6
                z *= Z6[(-2 * e_sig) % 6] * self.gm1[k] * self.g[3][k] * Z6[e_eps_m3]
            err = max(err, abs(z / self.zetaZ[t] - 1))
        self.zeta_check_err = err
        self.zeta_check_n = len(sub)


def quad_twist(R, rows, nu_a, nu_b, block=3000):
    """TW[row, nu] = (nu / m_row)_2 = prod_{p | m} chi_p(nu)^3  (0 if (nu, m) != 1), int8, for squarefree rows."""
    used = np.unique(R.IncS[rows].indices)
    nn = len(nu_a)
    Ep = np.zeros((len(R.P), nn), np.int8); Zp = np.zeros((len(R.P), nn), np.int8)
    for k in used:
        e = R.P[k].symbol(nu_a, nu_b)
        Zp[k] = e < 0; Ep[k] = np.where(e >= 0, e, 0) % 2
    Inc = R.IncS[rows]
    TW = np.zeros((len(rows), nn), np.int8)
    for c0 in range(0, nn, block):
        par = (Inc @ Ep[:, c0:c0 + block].astype(np.int64)) % 2
        zer = (Inc @ Zp[:, c0:c0 + block].astype(np.int64)) > 0
        TW[:, c0:c0 + block] = np.where(zer, 0, 1 - 2 * par)
    return TW


def nu_columns(numax, numin=7):
    Pn = primes_upto(numax)
    nus = [(a, b, N) for (a, b, N, idx) in squarefree_products(Pn, numax) if N >= numin]
    nu_a = np.array([x[0] for x in nus], np.int64); nu_b = np.array([x[1] for x in nus], np.int64)
    nu_N = np.array([x[2] for x in nus], np.int64)
    o = np.argsort(nu_N, kind="stable")
    return nu_a[o], nu_b[o], nu_N[o]


def pi_layout(R):
    """squarefree ideals x 6 unit classes.  Returns IX (nI,6) row indices, OmI, Vp (K,nI,6) per-scale pi, pi2 (nI,6) |pi|^2 est."""
    sfi = np.unique(R.ideal[R.sf])
    pos = {v: k for k, v in enumerate(sfi)}
    nI = len(sfi)
    IX = -np.ones((nI, 6), np.int64)
    for t in np.where(R.sf)[0]:
        IX[pos[R.ideal[t]], R.ucls[t]] = t
    assert (IX >= 0).all()
    OmI = R.Om[IX[:, 0]]
    Vp = R.V[:, IX] * R.zetaZ[IX][None]
    pi2 = cross_est(np.moveaxis(Vp, 0, -1), R.w)
    return IX, OmI, Vp, pi2


def divides(da, db, dN, xa, xb):
    """boolean array: d | x for elements x (vectorised); d = da + db w of norm dN."""
    ya, yb = mul(xa, xb, da - db, -db)          # x * conj(d)
    return (ya % dN == 0) & (yb % dN == 0)
