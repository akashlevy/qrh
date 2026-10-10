"""Attack E, stage 2: the sextic additive polynomial at large Q, scales Y = 2^k (KMIN <= k <= KMAX), s up to X = 2^(KMAX+1).

  A_k(m) = Y^{-1} sum_{s sf primary, (s,6)=1, N s in (Y,2Y)} W(Ns/Y) (Ns/Y)^{-1/2} gtilde(s,-m),   Y = 2^k,
  gtilde(s,-m) = conj(chi_s(-m)) gtilde(s),  gtilde(s1 s2) = chi_{s1}(s2) chi_{s2}(s1) gtilde(s1) gtilde(s2)

(identical to D_numerics/expA.py, whose output A_Q*_K12 is reproduced in E_validate.py).

Associates.  For m = zeta6^c m0 (m0 primary), chi_s(-m) = zeta6^{c t(s)} chi_s(-m0), t(s) = sum_{p|s} (N p - 1)/6 mod 6,
since (zeta6/p)_6 = zeta6^{(Np-1)/6} (Euler).  So we accumulate the six components
  A^{(t)}_k(m0) = sum_{s: t(s)=t} coef_s conj(chi_s(-m0)),   A_k(zeta6^c m0) = sum_t conj(zeta6^{c t}) A^{(t)}_k(m0),
which needs the symbol table only on primary representatives (1/6 of the rows).

usage: python E_build.py KMAX Q1 Q2 ...       (KMIN = 8)
output: vls/data/A_Q{Q}_E{KMAX}.npy  (K x nrows, K12 layout: row k-KMIN <-> Y = 2^k)  and  m_Q{Q}_E{KMAX}.npy
"""
import numpy as np, sys, os, time
from math import isqrt
from multiprocessing import Pool, shared_memory
GS = "/Users/akashlevy/Documents/riemann/gram_split"
sys.path.insert(0, GS); sys.path.insert(0, GS + "/vls")
from zomega import Prime, Z6, elements_in_annulus, mul
from E_primes import powmod_vec

DATA = GS + "/vls/data/"
KMIN = 8
CZ6 = np.conj(Z6)
UNITS6 = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]       # zeta6^c, c = 0..5


def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1.0 / (u * (1 - u)) + 4.0); return out


def load_primes(X):
    d = np.load(DATA + "primes_X%d.npz" % X)
    return {k: d[k] for k in d.files}


def match_roots(y, roots):
    """y: (n,) or (nb, n) residues; roots: (n, 6) or (nb, 1, 6) -> index j with y == roots[..., j]; -1 if none."""
    out = np.full(y.shape, -1, np.int8)
    for j in range(6):
        out[y == roots[..., j]] = j
    return out


def symbols_of_elements_mod_primes(xa, xb, P, idx):
    """chi_{P[idx_k]}(x_k), elementwise (split primes only).  -1 where divisible."""
    p = P["q"][idx]; r = P["r"][idx]
    v = (xa % p + (xb % p) * r) % p
    y = powmod_vec(v, (p - 1) // 6, p)
    out = match_roots(y, P["roots"][idx])
    out[v == 0] = -1
    assert (out[v != 0] >= 0).all()
    return out


def zprime(P, i):
    return Prime(int(P["a"][i]), int(P["b"][i]), int(P["N"][i]), "inert" if P["kind"][i] else "split", int(P["q"][i]))


# ------------------------------------------------------------------ pair-symbol tables for small primes
def pair_tables(P, nsmall):
    """T1T[x, f] = chi_f(x) (prime x mod small prime f),  T2T[x, f] = chi_x(f) (small f mod prime x).  int8, col nsmall = 0 (pad)."""
    nP = len(P["N"])
    T1T = np.zeros((nP, nsmall + 1), np.int8); T2T = np.zeros((nP, nsmall + 1), np.int8)
    split = np.where(P["kind"] == 0)[0]; inert = np.where(P["kind"] == 1)[0]
    for f in range(nsmall):
        Pf = zprime(P, f)
        e = Pf.symbol(P["a"], P["b"]); e[e < 0] = 0          # only x = f itself; never used (squarefree)
        T1T[:, f] = e
        e2 = np.zeros(nP, np.int8)
        e2[split] = symbols_of_elements_mod_primes(np.full(len(split), P["a"][f]), np.full(len(split), P["b"][f]), P, split)
        for x in inert:
            e2[x] = zprime(P, x).symbol(np.array([P["a"][f]]), np.array([P["b"][f]]))[0]
        e2[e2 < 0] = 0
        T2T[:, f] = e2
    return T1T, T2T


# ------------------------------------------------------------------ squarefree s with gtilde(s), t(s), factor lists
def build_s(P, X):
    t0 = time.time()
    N = P["N"]; gt = P["gt"]; nP = len(N)
    nsmall = int(np.searchsorted(N, isqrt(X), side="right"))
    T1T, T2T = pair_tables(P, nsmall)
    print("pair tables (%d small primes) [%.1fs]" % (nsmall, time.time() - t0), flush=True)
    tp = ((N - 1) // 6) % 6
    # phase 1: products of small primes
    sN = np.array([1], np.int64); sg = np.array([1 + 0j]); st = np.array([0], np.int8)
    sf = np.full((1, 6), -1, np.int32); sw = np.array([0], np.int8)
    for f in range(nsmall):
        ok = sN * N[f] < X
        if not ok.any(): continue
        fac = sf[ok]; fpad = np.where(fac < 0, nsmall, fac)
        e = T1T[f][fpad].astype(np.int64).sum(1) + T2T[f][fpad].astype(np.int64).sum(1)
        nf = fac.copy(); nf[np.arange(len(nf)), sw[ok]] = f
        sN = np.concatenate([sN, sN[ok] * N[f]]); sg = np.concatenate([sg, sg[ok] * gt[f] * Z6[e % 6]])
        st = np.concatenate([st, ((st[ok] + tp[f]) % 6).astype(np.int8)])
        sf = np.concatenate([sf, nf]); sw = np.concatenate([sw, sw[ok] + 1])
    o = np.argsort(sN, kind="stable"); sN, sg, st, sf, sw = sN[o], sg[o], st[o], sf[o], sw[o]
    print("small-smooth s: %d [%.1fs]" % (len(sN), time.time() - t0), flush=True)
    fpad_all = np.where(sf < 0, nsmall, sf)
    LN, LG, LT, LF, LW = [sN], [sg], [st], [sf], [sw]
    for x in range(nsmall, nP):
        k = int(np.searchsorted(sN, (X - 1) // N[x], side="right"))
        if k == 0: break
        fp = fpad_all[:k]
        e = T1T[x][fp].astype(np.int64).sum(1) + T2T[x][fp].astype(np.int64).sum(1)
        nf = sf[:k].copy(); nf[np.arange(k), sw[:k]] = x
        LN.append(sN[:k] * N[x]); LG.append(sg[:k] * gt[x] * Z6[e % 6]); LT.append(((st[:k] + tp[x]) % 6).astype(np.int8))
        LF.append(nf); LW.append(sw[:k] + 1)
    sN = np.concatenate(LN); sg = np.concatenate(LG); st = np.concatenate(LT); sf = np.concatenate(LF); sw = np.concatenate(LW)
    print("all squarefree s < X: %d [%.1fs]" % (len(sN), time.time() - t0), flush=True)
    return sN, sg, st, sf, sw


# ------------------------------------------------------------------ symbol table on rows (shared memory)
def _jblock(args):
    shm_name, shape, i0, i1, a, b, q, kind, r, roots, xa, xb = args
    shm = shared_memory.SharedMemory(name=shm_name)
    J = np.ndarray(shape, np.int8, buffer=shm.buf)
    sp = np.where(kind == 0)[0]
    if len(sp):
        p = q[sp][:, None]; v = (xa[None, :] % p + (xb[None, :] % p) * r[sp][:, None]) % p
        e = ((q[sp] - 1) // 6)[:, None]
        res = np.ones_like(v); bb = v.copy()
        for bit in range(int(e.max()).bit_length()):
            m = ((e >> bit) & 1).astype(bool)
            res = np.where(m, (res * bb) % p, res); bb = (bb * bb) % p
        out = match_roots(res, roots[sp][:, None, :])
        out[v == 0] = -1
        J[i0 + sp] = out
    for k in np.where(kind == 1)[0]:
        Pr = Prime(int(a[k]), int(b[k]), int(q[k]) ** 2, "inert", int(q[k]))
        J[i0 + k] = Pr.symbol(xa, xb).astype(np.int8)
    shm.close()
    return i1 - i0


def _acc(args):
    shm_name, shape, tasks, nK = args
    shm = shared_memory.SharedMemory(name=shm_name)
    J = np.ndarray(shape, np.int8, buffer=shm.buf)
    nrows = shape[1]
    A = np.zeros((nK, 6, nrows), complex)
    for (k, idx, tc, coef) in tasks:
        Jg = J[idx]                                   # (C, w, nrows)
        bad = (Jg < 0).any(axis=1)
        e = Jg.sum(axis=1, dtype=np.int16) % 6
        vals = CZ6[e]; vals[bad] = 0
        Cm = np.zeros((6, len(coef)), complex); Cm[tc, np.arange(len(coef))] = coef
        A[k] += Cm @ vals
    shm.close()
    return A


if __name__ == "__main__":
    KMAX = int(sys.argv[1]); Qs = [int(x) for x in sys.argv[2:]]
    X = int(os.environ.get("PRIMES_X", 2 ** 22))   # prime / Gauss-sum table; s restricted to < 2^(KMAX+1)
    Xs = 2 ** (KMAX + 1)
    nK = KMAX - KMIN + 1
    t0 = time.time()
    P = load_primes(X)
    keepP = P["N"] < Xs
    P = {k: (v[keepP] if (hasattr(v, "shape") and v.shape and v.shape[0] == len(keepP)) else v) for k, v in P.items()}
    nP = len(P["N"])
    print("primes with N < %d: %d [%.1fs]" % (Xs, nP, time.time() - t0), flush=True)
    sfile = DATA + "slist_X%d.npz" % Xs
    if os.path.exists(sfile):
        d = np.load(sfile); sN, sg, st, sf, sw = d["sN"], d["sg"], d["st"], d["sf"], d["sw"]
    else:
        sN, sg, st, sf, sw = build_s(P, Xs)
        np.savez(sfile, sN=sN, sg=sg, st=st, sf=sf, sw=sw)
    # scales
    ks = np.floor(np.log2(sN.astype(float))).astype(int)
    ok = (ks >= KMIN) & (ks <= KMAX)
    sN, sg, st, sf, sw, ks = sN[ok], sg[ok], st[ok], sf[ok], sw[ok], ks[ok]
    rr = sN / 2.0 ** ks
    w = bump(rr, 1, 2) * rr ** -0.5 / 2.0 ** ks
    ok = w > 0
    sN, sg, st, sf, sw, ks, coef = sN[ok], sg[ok], st[ok], sf[ok], sw[ok], ks[ok], (w * sg)[ok]
    print("s in scales 2^%d..2^%d: %d  (max omega %d) [%.1fs]" % (KMIN, KMAX, len(sN), sw.max(), time.time() - t0), flush=True)
    o = np.lexsort((sw, ks)); sN, st, sf, sw, ks, coef = sN[o], st[o], sf[o], sw[o], ks[o], coef[o]
    for Q in Qs:
        ma, mb, mN = elements_in_annulus(Q / 2, Q)
        Om = bump(mN / Q, 0.5, 1.0); keep = Om > 1e-12
        ma, mb, mN = ma[keep].astype(np.int64), mb[keep].astype(np.int64), mN[keep].astype(np.int64)
        rep = (ma % 3 == 1) & (mb % 3 == 0)
        ra, rb, rN = ma[rep], mb[rep], mN[rep]; nr = len(ra)
        xa, xb = -ra, -rb
        shape = (nP, nr)
        shm = shared_memory.SharedMemory(create=True, size=nP * nr)
        try:
            J = np.ndarray(shape, np.int8, buffer=shm.buf)
            B = max(64, int(2e7 // max(nr, 1)))
            args = [(shm.name, shape, i0, min(i0 + B, nP), P["a"][i0:i0 + B], P["b"][i0:i0 + B], P["q"][i0:i0 + B],
                     P["kind"][i0:i0 + B], P["r"][i0:i0 + B], P["roots"][i0:i0 + B], xa, xb) for i0 in range(0, nP, B)]
            with Pool(11) as pool:
                pool.map(_jblock, args)
            print("Q=%d: %d primary reps (%d rows); symbol table %s [%.1fs]" % (Q, nr, len(ma), shape, time.time() - t0), flush=True)
            # chunk tasks
            C = max(256, int(1.5e7 // nr))
            tasks = []
            bounds = np.flatnonzero(np.diff(ks * 8 + sw)) + 1
            starts = np.concatenate([[0], bounds]); ends = np.concatenate([bounds, [len(ks)]])
            for s0, s1 in zip(starts, ends):
                wv = int(sw[s0]); k = int(ks[s0]) - KMIN
                for c0 in range(s0, s1, C):
                    c1 = min(c0 + C, s1)
                    tasks.append((k, sf[c0:c1, :wv].astype(np.int64), st[c0:c1].astype(np.int64), coef[c0:c1]))
            ntask = 44
            groups = [tasks[i::ntask] for i in range(ntask)]
            with Pool(11) as pool:
                parts = pool.map(_acc, [(shm.name, shape, g, nK) for g in groups])
            At = sum(parts)                       # (nK, 6, nr)  t-components on primary reps
        finally:
            shm.close(); shm.unlink()
        print("Q=%d accumulated (%d tasks) [%.1fs]" % (Q, len(tasks), time.time() - t0), flush=True)
        # expand to all associates
        Afull = np.zeros((nK, 6 * nr), complex); MA = np.zeros(6 * nr, np.int64); MB = np.zeros(6 * nr, np.int64)
        tt = np.arange(6)
        for c, (ua, ub) in enumerate(UNITS6):
            ph = CZ6[(c * tt) % 6]                 # conj(zeta6^{c t})
            Afull[:, c * nr:(c + 1) * nr] = np.einsum("t,ktr->kr", ph, At)
            MA[c * nr:(c + 1) * nr], MB[c * nr:(c + 1) * nr] = mul(np.int64(ua), np.int64(ub), ra, rb)
        MN = np.tile(rN, 6)
        np.save(DATA + "A_Q%d_E%d.npy" % (Q, KMAX), Afull)
        np.save(DATA + "m_Q%d_E%d.npy" % (Q, KMAX), np.stack([MA, MB, MN]))
        np.save(DATA + "At_Q%d_E%d.npy" % (Q, KMAX), At)
        print("Q=%d saved [%.1fs]" % (Q, time.time() - t0), flush=True)
