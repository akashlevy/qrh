"""
R1 numerics: where the column term reappears in paper2's dual (B-level) for H < D, and whether the
signed (coprime) form cancels it.

B-level column sums (paper2 eq. intro-gauss-family, xi = 1, f = 1):
    B_h(X) = sum_{n sf primary, (n,6)=1, X < N n <= 2X} a(n) chi_n(h) W(N n / X),
    a(n) = conj(alpha(n)) gamma_2(n)  (|a(n)| = 1),   W(x) = sin^2(pi (x-1)).
For every squarefree primary v put
    Q_v = sum_{0 < N h <= calH} | sum_{n : v | n} a(n) chi_n(h) W(N n/X) |^2   (>= 0).
Then
    Q_1                  = the positive dual form (paper2's E(calH, X) * X),
    sum_v |mu(v)| Q_v    = what paper2's positivity step (Lemma smooth-mean-square applied with
                           coefficient mu(f)) effectively bounds,
    Q_cop = sum_v mu(v) Q_v = sum_h sum_{(n1,n2)=1} a(n1) conj(a(n2)) chi_{n1}(h) conj(chi_{n2}(h)) W W
                         (the TRUE coprime structure produced by Poisson in u; no diagonal at all).
Heuristic (GRH): Q_1 ~ diag ~ c calH X ;  Q_cop = O(X^2).   Target needed for M_kappa: X^2.
calH = X^{2 - kappa} corresponds to original rows H = X^kappa.

Also reports the A-level mean square  S(H) = sum_{0<N u<=X^kappa} |A_u(X)|^2  (mu coefficients, same
columns and weight) divided by its diagonal, as a check that M_kappa itself looks true numerically.
"""
import sys, math, time, itertools
import numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/strengthening")
from eisenstein import eisenstein_primes, ResidueField, gauss_prime, primary, norm, cx, mul, Z6
from wall_experiment import elements_upto

def inert_pow_table(q, ra, rb, e):
    """vectorised x^e in F_q[w] (w^2 = -1 - w) for x = a + b w."""
    a = ra % q; b = rb % q
    ra_, rb_ = np.ones_like(a), np.zeros_like(b)
    ba, bb = a.copy(), b.copy()
    while e:
        if e & 1:
            ra_, rb_ = (ra_ * ba - rb_ * bb) % q, (ra_ * bb + rb_ * ba - rb_ * bb) % q
        ba, bb = (ba * ba - bb * bb) % q, (2 * ba * bb - bb * bb) % q
        e >>= 1
    return ra_, rb_

def sextic_exponents(F, ra, rb):
    """exponent k in 0..5 of (h/pi)_6 for every row h = a + b w; -1 where pi | h."""
    q = F.q
    if F.kind == "split":
        x = (ra + rb * F.r) % q
        e = (q - 1) // 6
        res = np.ones_like(x); base = x.copy()
        while e:
            if e & 1: res = res * base % q
            base = base * base % q; e >>= 1
        w6 = (-(F.r * F.r)) % q
        lut = np.full(q, -1, dtype=np.int8); t = 1
        for k in range(6):
            lut[t] = k; t = t * w6 % q
        out = np.full(len(x), -1, dtype=np.int8)
        nz = x != 0
        out[nz] = lut[res[nz]]
        return out
    # inert: field F_{q^2}
    za = (ra % q == 0) & (rb % q == 0)
    pa, pb = inert_pow_table(q, ra, rb, (q * q - 1) // 6)
    out = np.full(len(ra), -1, dtype=np.int8)
    t = (1, 0); w6 = (1, 1)
    for k in range(6):
        m = (pa == t[0]) & (pb == t[1]) & ~za
        out[m] = k
        t = ((t[0] * w6[0] - t[1] * w6[1]) % q, (t[0] * w6[1] + t[1] * w6[0] - t[1] * w6[1]) % q)
    return out

def setup_columns(X):
    primes = [p for p in eisenstein_primes(2 * X + 10) if p[2] != 2]
    fields = {pr[1]: ResidueField(*pr) for pr in primes}
    pnorm = {pr[1]: (pr[2] if pr[0] == "split" else pr[2] ** 2) for pr in primes}
    plist = sorted(fields, key=lambda p: pnorm[p])
    cols = []
    for e in elements_upto(2 * X):
        if norm(e) <= X: continue
        pe = primary(e)
        if pe is None or pe != e: continue
        rem = e; fac = []; ok = True; Nn = norm(e)
        for p in plist:
            if pnorm[p] > Nn: break
            Fp = fields[p]; red = Fp.red(rem)
            if (red == 0 if Fp.kind == "split" else red == (0, 0)):
                pc = (p[0] - p[1], -p[1]); prod = mul(rem, pc)
                rem = (prod[0] // pnorm[p], prod[1] // pnorm[p]); fac.append(p); Nn //= pnorm[p]
                red = Fp.red(rem)
                if (red == 0 if Fp.kind == "split" else red == (0, 0)): ok = False; break
                if Nn == 1: break
        if ok and Nn == 1: cols.append((e, tuple(fac)))
    need = sorted({p for _, f in cols for p in f}, key=lambda p: pnorm[p])
    g2 = {p: gauss_prime(fields[p], 2) for p in need}
    coef = []
    for e, fac in cols:
        g = 1 + 0j; acc = (1, 0)
        for p in fac:
            k = fields[p].sextic_exp(acc) if acc != (1, 0) else 0
            g *= g2[p] * Z6[(4 * k) % 6]                 # gamma_2(ab) = gamma_2(a) gamma_2(b) chi_b(a)^4
            acc = mul(acc, p)
        al = cx(e) / abs(cx(e))
        coef.append(g * al.conjugate())
    coef = np.array(coef)
    Ncol = np.array([norm(e) for e, _ in cols], dtype=float)
    W = np.sin(math.pi * (Ncol / X - 1)) ** 2
    mu = np.array([(-1) ** len(f) for _, f in cols], dtype=float)
    return cols, fields, pnorm, need, coef, W, mu

def run(X, kappas, chunk=8000, seed_check=True):
    t0 = time.time()
    cols, fields, pnorm, need, coef, W, mu = setup_columns(X)
    nc = len(cols)
    # divisor structure: v (squarefree, as a sorted tuple of primes) -> list of column indices
    vmap = {}
    for j, (_, fac) in enumerate(cols):
        for r in range(len(fac) + 1):
            for sub in itertools.combinations(fac, r):
                vmap.setdefault(sub, []).append(j)
    vkeys = list(vmap)
    vidx = [np.array(vmap[v]) for v in vkeys]
    vmu = np.array([(-1) ** len(v) for v in vkeys], dtype=float)
    print(f"X={X}: {nc} columns, {len(need)} primes in support, {len(vkeys)} divisors v; setup {time.time()-t0:.1f}s",
          flush=True)
    calHs = {k: X ** (2 - k) for k in kappas}
    Hmax = max(calHs.values())
    rows = elements_upto(int(Hmax) + 1)
    ra = np.array([r[0] for r in rows], dtype=np.int64); rb = np.array([r[1] for r in rows], dtype=np.int64)
    rN = ra * ra - ra * rb + rb * rb
    o = np.argsort(rN, kind="stable"); ra, rb, rN = ra[o], rb[o], rN[o]
    nrows = len(rN)
    cuts = {k: int(np.searchsorted(rN, calHs[k], side="right")) for k in kappas}
    cw = coef * W
    # accumulators per kappa
    Qv = {k: np.zeros(len(vkeys)) for k in kappas}
    diag = {k: 0.0 for k in kappas}
    for s0 in range(0, nrows, chunk):
        s1 = min(nrows, s0 + chunk)
        ca, cb = ra[s0:s1], rb[s0:s1]
        tab = {p: sextic_exponents(fields[p], ca, cb) for p in need}
        C = np.zeros((s1 - s0, nc), dtype=np.complex128)
        Mk = np.zeros((s1 - s0, nc), dtype=bool)
        for j, (_, fac) in enumerate(cols):
            ex = np.zeros(s1 - s0, dtype=np.int16); m = np.ones(s1 - s0, dtype=bool)
            for p in fac:
                t = tab[p]; m &= (t >= 0); ex += t
            C[:, j] = np.where(m, np.exp(1j * math.pi / 3 * ex), 0)
            Mk[:, j] = m
        C *= cw[None, :]
        dcum = np.cumsum(Mk.astype(np.float64) @ (np.abs(cw) ** 2))
        nks = {k: min(cuts[k] - s0, s1 - s0) for k in kappas}
        for k in kappas:
            if nks[k] > 0: diag[k] += float(dcum[nks[k] - 1])
        for i, idx in enumerate(vidx):
            S = C[:, idx].sum(axis=1)
            cs = np.cumsum(S.real ** 2 + S.imag ** 2)
            for k in kappas:
                if nks[k] > 0: Qv[k][i] += float(cs[nks[k] - 1])
        print(f"  rows {s1}/{nrows}  {time.time()-t0:.0f}s", flush=True)
    out = {}
    for k in kappas:
        calH = calHs[k]
        q1 = Qv[k][vkeys.index(())]
        qabs = float(np.sum(Qv[k]))
        qcop = float(np.sum(vmu * Qv[k]))
        out[k] = dict(calH=calH, rows=cuts[k], Q1=q1, diag=diag[k], Qabs=qabs, Qcop=qcop)
        print(f"  kappa={k:.4f} calH=X^{2-k:.3f}={calH:9.0f} rows={cuts[k]:7d} | "
              f"Q1/diag={q1/diag[k]:.3f}  Q1/X^2={q1/X**2:8.2f}  sum|Q_v|/X^2={qabs/X**2:8.2f}  "
              f"Q_cop/X^2={qcop/X**2:+.3f}  Q_cop/Q1={qcop/q1:+.4f}", flush=True)
    # A-level check with the same columns: S(H) = sum_{0<N u<=X^kappa} |A_u(X)|^2, A_u = sum mu(n) chi_n(u) W
    Hmax_A = max(X ** k for k in kappas)
    rowsA = elements_upto(int(Hmax_A) + 1)
    ua = np.array([r[0] for r in rowsA], dtype=np.int64); ub = np.array([r[1] for r in rowsA], dtype=np.int64)
    uN = ua * ua - ua * ub + ub * ub; o = np.argsort(uN, kind="stable"); ua, ub, uN = ua[o], ub[o], uN[o]
    tab = {p: sextic_exponents(fields[p], ua, ub) for p in need}
    A = np.zeros(len(uN), dtype=complex); dA = np.zeros(len(uN))
    for j, (_, fac) in enumerate(cols):
        ex = np.zeros(len(uN), dtype=np.int16); m = np.ones(len(uN), dtype=bool)
        for p in fac:
            t = tab[p]; m &= (t >= 0); ex += t
        A += np.where(m, mu[j] * W[j] * np.exp(1j * math.pi / 3 * ex), 0)
        dA += np.where(m, W[j] ** 2, 0)
    cumA = np.cumsum(np.abs(A) ** 2); cumD = np.cumsum(dA)
    for k in kappas:
        i = int(np.searchsorted(uN, X ** k, side="right"))
        out[k].update(SA=float(cumA[i - 1]), SAdiag=float(cumD[i - 1]), rowsA=i)
        print(f"  A-level kappa={k:.4f}: rows={i:6d}  S(H)/diag = {cumA[i-1]/cumD[i-1]:.3f}   "
              f"S(H)/X^2 = {cumA[i-1]/X**2:.3f}   S(H)/(H X) = {cumA[i-1]/(X**k*X):.3f}", flush=True)
    return out

if __name__ == "__main__":
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    kappas = [float(a) for a in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.5, 0.7, 6 / 7, 1.0]
    run(X, kappas)
