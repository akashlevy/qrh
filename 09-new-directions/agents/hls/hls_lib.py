"""
hls: numerics for the hybrid large sieve (HLS) over Z[w], w = exp(2 pi i/3).

Family (rows):  xi_{k,P}(n) = psi_k(n) * eta_P(n),
    psi_k(n) = chi_k(n)^3 = prod_{pi | k} (n/pi)_6^3      (quadratic, k squarefree primary odd, N(k) <= K)
    eta_P(n) = chi_P(n)^{-2} = (n/P)_6^{-2}               (cubic, P split primary prime, L0 <= N(P) < 2 L0)
Columns: n squarefree primary odd, N0 < N(n) <= 2 N0.

(HLS) asks  sum_{k,P} |sum_n a_n xi_{k,P}(n)|^2 << (#rows + #cols) ||a||^2, i.e.
    lambda_max(S S^*) << #rows + #cols,  S = (xi_{k,P}(n))  (#rows x #cols).
Structure used for fast mat-vecs:  S a = Psi (diag(a) H^T)   (Psi: #k x #cols real, H: #P x #cols).
S^* S = (Psi^T Psi) o (H^* H)  (Schur/Hadamard product of the quadratic and cubic Gram matrices).

Reuses edge-theta/edge_lib.py and energy-wall/{eisenstein,dual_energy}.py (read only).
"""
import sys, math, time
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/energy-wall")
from edge_lib import setup, Symb, Gamma          # noqa: E402
import dual_energy as de                         # noqa: E402
from eisenstein import cx, Z6                    # noqa: E402

TWOPI = 2 * math.pi


def gauss_all(F):
    """gamma_j(pi) for j = 0..5 (index j), e(z) = exp(4 pi i Im z / sqrt3), normalised by N(pi)^{-1/2}."""
    pc = cx(F.p)
    if F.kind == "split":
        t = np.arange(1, F.q, dtype=np.int64)
        k = de.sextic_exp_vec(F, t, np.zeros_like(t)).astype(np.int64)
        z = t / pc
        add = np.exp(4j * np.pi * z.imag / math.sqrt(3))
        out = np.array([np.sum(np.exp(TWOPI * 1j * j * k / 6) * add) for j in range(6)]) / math.sqrt(F.q)
        return out
    q = F.q
    xs, ys = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    xs = xs.ravel(); ys = ys.ravel(); keep = (xs != 0) | (ys != 0)
    xs, ys = xs[keep], ys[keep]
    k = Symb(F)(xs, ys).astype(np.int64)
    z = (xs + ys * complex(-0.5, math.sqrt(3) / 2)) / pc
    add = np.exp(4j * np.pi * z.imag / math.sqrt(3))
    return np.array([np.sum(np.exp(TWOPI * 1j * j * k / 6) * add) for j in range(6)]) / q


def elements(lo, hi, fields, pnorm, plist, mod4=False):
    A, B, Nn, fac = de.primary_sf_elements(lo, hi, fields, plist, pnorm)
    odd = np.array([all(pnorm[p] != 4 for p in f) for f in fac], bool)
    if mod4:
        odd &= (A % 4 == 1) & (B % 4 == 0)
    keep = np.nonzero(odd)[0]
    A, B, Nn = A[keep], B[keep], Nn[keep]; fac = [fac[i] for i in keep]
    o = np.argsort(Nn, kind="stable")
    return A[o], B[o], Nn[o], [fac[i] for i in o]


class Family:
    def __init__(self, N0, K, L0, nP=None, Kmin=0, mod4=False, gauss=True, inert_marks=False,
                 verbose=True, dtype=np.float64):
        t0 = time.time()
        self.N0, self.K, self.L0 = N0, K, L0
        fields, pnorm, plist = setup(max(2 * N0, K, 2 * L0) + 10)
        self.fields, self.pnorm, self.plist = fields, pnorm, plist
        # columns
        A, B, Nn, fac = elements(N0, 2 * N0, fields, pnorm, plist, mod4=mod4)
        self.A, self.B, self.Nn, self.fac = A, B, Nn, fac
        self.nc = len(A)
        self.mu = np.array([(-1) ** len(f) for f in fac], float)
        self.alpha = (A + B * complex(-0.5, math.sqrt(3) / 2)) / np.sqrt(Nn)
        # rows k
        kA, kB, kN, kfac = elements(Kmin, K, fields, pnorm, plist)
        if Kmin == 0 and (len(kN) == 0 or kN[0] != 1):
            raise RuntimeError("k = 1 missing")
        self.kN, self.kfac = kN, kfac
        self.nk = len(kN)
        # marks P
        marks = [p for p in plist if L0 <= pnorm[p] < 2 * L0 and (inert_marks or fields[p].kind == "split")]
        if nP: marks = marks[:nP]
        self.marks = marks
        self.nP = len(marks)
        # sextic exponents of columns at the needed primes
        need = sorted({p for f in kfac for p in f} | set(marks), key=lambda p: pnorm[p])
        self.E = {p: Symb(fields[p])(A, B) for p in need}
        Psi = np.empty((self.nk, self.nc), dtype=dtype)
        for i, f in enumerate(kfac):
            v = np.ones(self.nc, dtype=dtype)
            for p in f:
                e = self.E[p]
                v *= np.where(e < 0, 0, np.where(e % 2 == 1, -1, 1)).astype(dtype)
            Psi[i] = v
        self.Psi = Psi
        H = np.empty((self.nP, self.nc), complex)
        for j, p in enumerate(marks):
            e = self.E[p].astype(np.int64)
            H[j] = np.where(e < 0, 0, np.exp(-TWOPI * 1j * 2 * e / 6))
        self.H = H
        self.nrows = self.nk * self.nP
        if verbose:
            print(f"[Family] N0={N0} cols={self.nc}  K={K} k's={self.nk}  L0={L0} P's={self.nP}  "
                  f"rows={self.nrows}  ({time.time()-t0:.1f}s)", flush=True)
        if gauss:
            self._gauss(verbose, t0)

    # ---------------------------------------------------------------- Gauss-sum (theta-type) vectors
    def _gauss(self, verbose, t0):
        fields = self.fields
        primes = sorted({p for f in self.fac for p in f}, key=lambda p: self.pnorm[p])
        g = {p: gauss_all(fields[p]) for p in primes}
        G = np.ones((self.nc, 6), complex)
        for i, f in enumerate(self.fac):
            v = np.ones(6, complex)
            for p in f:
                v = v * g[p]
            # cross symbols  prod_{i != l} chi_{p_i}(p_l)^j   (pure CRT, no reciprocity used)
            ex = 0
            for a in f:
                Fa = fields[a]
                for b in f:
                    if a == b: continue
                    k = Fa.sextic_exp(b)
                    ex += k
            v = v * np.array([Z6[(j * ex) % 6] for j in range(6)])
            G[i] = v
        self.gam = G             # gam[:, j] = gamma_j(n)
        self.beta = G[:, 2] * np.conj(self.alpha)
        if verbose:
            print(f"[Family] Gauss sums at {len(primes)} primes ({time.time()-t0:.1f}s)", flush=True)

    # ---------------------------------------------------------------- operators
    def S(self, a):
        """(S a)_{k,P} = sum_n a_n psi_k(n) eta_P(n)  ->  array (nk, nP)."""
        M = a[:, None] * self.H.T                      # nc x nP
        dt = self.Psi.dtype
        return (self.Psi @ M.real.astype(dt)).astype(np.float64) + 1j * (self.Psi @ M.imag.astype(dt))

    def Sh(self, Y):
        """(S^* Y)_n = sum_{k,P} conj(psi_k(n) eta_P(n)) Y_{k,P}."""
        dt = self.Psi.dtype
        Z = (self.Psi.T @ Y.real.astype(dt)).astype(np.float64) + 1j * (self.Psi.T @ Y.imag.astype(dt))
        return np.sum(np.conj(self.H.T) * Z, axis=1)

    def R(self, a, split_k1=False):
        """Rayleigh quotient sum_{k,P}|S a|^2 / ||a||^2 (optionally without the k = 1 rows)."""
        Y = self.S(a)
        tot = np.sum(np.abs(Y) ** 2) / np.sum(np.abs(a) ** 2)
        if not split_k1:
            return float(tot)
        k1 = np.sum(np.abs(Y[0]) ** 2) / np.sum(np.abs(a) ** 2)
        return float(tot), float(k1)

    def top(self, k=4, tol=1e-6, maxiter=None, side="auto", seed=0):
        """largest eigenvalues of S^* S (= of S S^*), via block Lanczos (scipy eigsh)."""
        from scipy.sparse.linalg import LinearOperator, eigsh
        if side == "auto":
            side = "col" if self.nc <= self.nrows else "row"
        rng = np.random.default_rng(seed)
        if side == "col":
            op = LinearOperator((self.nc, self.nc), matvec=lambda x: self.Sh(self.S(np.ravel(x))),
                                dtype=complex)
            v0 = rng.standard_normal(self.nc) + 1j * rng.standard_normal(self.nc)
        else:
            sh = (self.nk, self.nP)
            op = LinearOperator((self.nrows, self.nrows),
                                matvec=lambda y: self.S(self.Sh(np.ravel(y).reshape(sh))).ravel(),
                                dtype=complex)
            v0 = rng.standard_normal(self.nrows) + 1j * rng.standard_normal(self.nrows)
        vals, vecs = eigsh(op, k=k, which="LA", tol=tol, maxiter=maxiter, v0=v0)
        o = np.argsort(vals)[::-1]
        vals, vecs = vals[o], vecs[:, o]
        if side == "row":
            # convert to column-side eigenvectors a = S^* y / ||S^* y||
            cv = np.stack([self.Sh(vecs[:, i].reshape(self.nk, self.nP)) for i in range(vecs.shape[1])], 1)
            cv /= np.linalg.norm(cv, axis=0)
            vecs = cv
        return vals, vecs


def rand_control(fam, kind, seed=1):
    """replace the arithmetic matrices by random ones with the same zero pattern.
       kind = 'phase6': every entry of S an independent random 6th root of unity (needs explicit S; small only)
              'prod'  : Psi -> random signs, H -> random cube roots (product structure kept)
              'psi'   : only Psi random (quadratic part random, cubic part arithmetic)
              'H'     : only H random"""
    rng = np.random.default_rng(seed)
    import copy
    g = copy.copy(fam)
    if kind in ("prod", "psi"):
        sgn = rng.choice(np.array([-1, 1], dtype=fam.Psi.dtype), size=fam.Psi.shape)
        g.Psi = np.where(fam.Psi == 0, 0, sgn).astype(fam.Psi.dtype)
    if kind in ("prod", "H"):
        ph = np.exp(TWOPI * 1j * rng.integers(0, 3, size=fam.H.shape) / 3)
        g.H = np.where(fam.H == 0, 0, ph)
    return g
