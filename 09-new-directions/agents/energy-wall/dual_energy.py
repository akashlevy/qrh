"""
Dual-side (reflected) marked energy beyond the wall, numerically, in Z[w].

E(K) = sum_{k sf primary, N(k) <= K} | sum_{n sf primary, N0 < N(n) <= 2 N0} c(n) psi_k(n) w(N(n)/N0) |^2
psi_k(n) = chi_k(n)^3 = prod_{pi | k} (n/pi)_6^3     (quadratic, as in P78 Lemma hybrid-energy)
c(n) candidates:
  theta      : gamma_2(n) conj(alpha(n))           (normalised cubic Gauss sum = cubic-theta type coefficient)
  theta*X    : theta(n) * prod_slots sum_{p in slot} N(p)^{-1/2} chi_p(n)^{-2}      (P78 dual marks)
  theta*chiP : theta(n) chi_P(n)^{-2}, single prime P (the per-P, "level P" energy)
  rand*|X|   : random unit phases * |X(n)|          (control: generic coefficients)
  adv        : |theta X|(n) psi_{k0}(n)              (aligned with one row: attains the column term)
Diagonal prediction (W holds): E(K)/#rows(K) = sum_n |c(n)|^2 w^2 (times the coprimality density).
Large-sieve worst case: E(K)/#rows ~ N0/K * diagonal.
Reflection proves diagonality only for K <= N0/P (P = norm of the active mark product).
"""
import math, sys, time, json
import numpy as np
from eisenstein import eisenstein_primes, ResidueField, norm, cx, mul, Z6

def powmod_vec(x, e, q):
    x = x.astype(np.int64) % q
    res = np.ones_like(x); base = x.copy()
    while e:
        if e & 1: res = res * base % q
        base = base * base % q; e >>= 1
    return res

def sextic_exp_vec(F, A, B):
    """exponent k of (a+bw / pi)_6 for arrays A,B; -1 where pi | element."""
    q = F.q
    if F.kind == "split":
        x = (A + B * F.r) % q
        v = powmod_vec(x, (q - 1) // 6, q)
        w6 = (-(F.r * F.r)) % q
        lut = np.full(q, -1, dtype=np.int8); t = 1
        for k in range(6):
            lut[t] = k; t = t * w6 % q
        out = lut[v]; out[x == 0] = -1
        return out
    out = np.empty(len(A), dtype=np.int8)
    for i, (a, b) in enumerate(zip(A, B)):
        k = F.sextic_exp((int(a), int(b)))
        out[i] = -1 if k is None else k
    return out

def gauss2_vec(F):
    """gamma_2(pi) = N(pi)^{-1/2} sum_x chi_pi(x)^2 e(x/pi), e(z)=exp(4 pi i Im z/sqrt3)."""
    pc = cx(F.p)
    if F.kind == "split":
        t = np.arange(1, F.q, dtype=np.int64)
        k = sextic_exp_vec(F, t, np.zeros_like(t)).astype(np.int64)
        z = t / pc
        return np.sum(np.exp(2j * np.pi * 2 * k / 6) * np.exp(4j * np.pi * z.imag / math.sqrt(3))) / math.sqrt(F.q)
    q = F.q
    xs, ys = np.meshgrid(np.arange(q), np.arange(q), indexing="ij")
    xs = xs.ravel(); ys = ys.ravel(); keep = (xs != 0) | (ys != 0)
    xs, ys = xs[keep], ys[keep]
    k = sextic_exp_vec(F, xs, ys).astype(np.int64)
    z = (xs + ys * complex(-0.5, math.sqrt(3) / 2)) / pc
    return np.sum(np.exp(2j * np.pi * 2 * k / 6) * np.exp(4j * np.pi * z.imag / math.sqrt(3))) / q

def primary_sf_elements(lo, hi, fields, plist, pnorm):
    """primary (a=1 mod 3, b=0 mod 3) elements with lo < N <= hi, with squarefree factorisation."""
    R = int(2 * math.sqrt(hi)) + 3
    a = np.arange(-R, R + 1); A, B = np.meshgrid(a, a, indexing="ij"); A = A.ravel(); B = B.ravel()
    Nn = A * A - A * B + B * B
    sel = (Nn > lo) & (Nn <= hi) & (A % 3 == 1) & (B % 3 == 0)
    A, B, Nn = A[sel].astype(np.int64), B[sel].astype(np.int64), Nn[sel].astype(np.int64)
    fac = [[] for _ in range(len(A))]
    sf = np.ones(len(A), bool)
    rem = Nn.copy()
    for p in plist:
        q = pnorm[p]
        if q > hi: break
        F = fields[p]
        if F.kind == "split":
            div = (A + B * F.r) % F.q == 0
        else:
            div = (A % F.q == 0) & (B % F.q == 0)
        idx = np.nonzero(div)[0]
        if len(idx) == 0: continue
        # quotient n / p = n * conj(p) / N(p); conj(a+bw) = (a-b) - b w
        pa, pb = p; ca, cb = pa - pb, -pb
        Aq = A[idx] * ca - B[idx] * cb; Bq = A[idx] * cb + B[idx] * ca - B[idx] * cb
        Aq //= q; Bq //= q
        if F.kind == "split":
            div2 = (Aq + Bq * F.r) % F.q == 0
        else:
            div2 = (Aq % F.q == 0) & (Bq % F.q == 0)
        sf[idx[div2]] = False
        for i in idx: fac[i].append(p)
        rem[idx] //= q
    ok = sf & (rem == 1)
    keep = np.nonzero(ok)[0]
    return A[keep], B[keep], Nn[keep], [fac[i] for i in keep]

def main(N0=20000, KMAX=20000, slots=((31,),), singleP=(7, 31, 97), seed=1, out=None):
    t0 = time.time()
    hi = max(2 * N0, KMAX)
    primes = [p for p in eisenstein_primes(hi + 10)]
    fields = {pr[1]: ResidueField(*pr) for pr in primes}
    pnorm = {pr[1]: (pr[2] if pr[0] == "split" else pr[2] ** 2) for pr in primes}
    plist = sorted(fields, key=lambda p: pnorm[p])
    cA, cB, cN, cfac = primary_sf_elements(N0, 2 * N0, fields, plist, pnorm)
    rA, rB, rN, rfac = primary_sf_elements(0, KMAX, fields, plist, pnorm)
    o = np.argsort(rN, kind="stable"); rA, rB, rN = rA[o], rB[o], rN[o]; rfac = [rfac[i] for i in o]
    print(f"cols {len(cA)} rows {len(rA)} primes {len(plist)}  {time.time()-t0:.1f}s", flush=True)
    # gamma_2 at primes dividing columns, then twisted multiplicativity gamma_2(ab)=gamma_2(a)gamma_2(b)chi_b(a)^4
    need = sorted({p for f in cfac for p in f}, key=lambda p: pnorm[p])
    g2 = {p: gauss2_vec(fields[p]) for p in need}
    # check gamma_2(p)^3 = -alpha(p)
    bad = sum(abs(g2[p] ** 3 + cx(p) / abs(cx(p))) > 1e-6 for p in need[:200])
    print(f"gamma_2 at {len(need)} primes; identity failures {bad}  {time.time()-t0:.1f}s", flush=True)
    theta = np.empty(len(cA), complex)
    for i, f in enumerate(cfac):
        g = 1 + 0j; acc = (1, 0)
        for p in f:
            k = fields[p].sextic_exp(acc) if acc != (1, 0) else 0
            g *= g2[p] * Z6[(4 * k) % 6]; acc = mul(acc, p)
        theta[i] = g * np.conj(complex(cA[i], 0) + cB[i] * complex(-0.5, math.sqrt(3) / 2)) / math.sqrt(cN[i])
    w = np.sin(np.pi * (cN / N0 - 1)) ** 2
    # sextic exponent tables of columns at all primes up to KMAX (for rows) and slot primes
    rowprimes = sorted({p for f in rfac for p in f}, key=lambda p: pnorm[p])
    par = {}; zer = {}
    for p in rowprimes:
        e = sextic_exp_vec(fields[p], cA, cB)
        zer[p] = e < 0; par[p] = (e % 2 == 1) & (e >= 0)        # (n/p)_6^3 = (-1)^e
    def chi_m2(p):   # chi_p(n)^{-2} on columns
        e = sextic_exp_vec(fields[p], cA, cB).astype(np.int64)
        return np.where(e >= 0, np.exp(-2j * np.pi * 2 * e / 6), 0)
    # mark polynomials
    split_by_norm = {}
    for p in plist:
        if fields[p].kind == "split": split_by_norm.setdefault(pnorm[p], []).append(p)
    def slot_poly(p0):
        ps = [p for nn in range(p0, 2 * p0) for p in split_by_norm.get(nn, [])]
        X = np.zeros(len(cA), complex)
        for p in ps: X += chi_m2(p) / math.sqrt(pnorm[p])
        return X, ps
    coefs = {}; info = {}
    coefs["theta"] = theta
    for sl in slots:
        X = np.ones(len(cA), complex); tot = 1; allps = []
        for p0 in sl:
            Xi, ps = slot_poly(p0); X *= Xi; tot *= p0; allps += ps
        name = "theta*X" + "x".join(map(str, sl))
        coefs[name] = theta * X; info[name] = dict(P=tot, nprimes=len(allps))
        rng = np.random.default_rng(seed)
        coefs["rand*|X|" + "x".join(map(str, sl))] = np.exp(2j * np.pi * rng.random(len(cA))) * np.abs(X)
    for P in singleP:
        p = split_by_norm[min(nn for nn in split_by_norm if nn >= P)][0]
        coefs[f"theta*chiP{pnorm[p]}"] = theta * chi_m2(p); info[f"theta*chiP{pnorm[p]}"] = dict(P=pnorm[p])
    # adversarial: aligned with a row k0 of norm about KMAX/4
    i0 = np.searchsorted(rN, KMAX // 4)
    k0 = rfac[i0]
    v0 = np.ones(len(cA))
    for p in k0: v0 = v0 * np.where(zer[p], 0, np.where(par[p], -1, 1))
    Xm, _ = slot_poly(slots[0][0])
    coefs["adv"] = np.abs(theta * Xm) * v0
    names = list(coefs)
    C = np.stack([coefs[n] * w for n in names], axis=1)          # ncols x m
    diag = np.sum(np.abs(C) ** 2, axis=0)
    # rows
    S = np.zeros((len(rA), len(names)), complex)
    for i, f in enumerate(rfac):
        if not f:
            v = np.ones(len(cA))
        else:
            pr = np.zeros(len(cA), bool); zz = np.zeros(len(cA), bool)
            for p in f: pr ^= par[p]; zz |= zer[p]
            v = np.where(zz, 0.0, np.where(pr, -1.0, 1.0))
        S[i] = v @ C
    print(f"rows done {time.time()-t0:.1f}s", flush=True)
    E = np.cumsum(np.abs(S) ** 2, axis=0)
    res = {}
    for j, n in enumerate(names):
        line = []
        for f in (1 / 256, 1 / 128, 1 / 64, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1, 2):
            K = f * N0
            if K > KMAX: continue
            i = np.searchsorted(rN, K, side="right")
            if i < 20: continue
            line.append((f, float(E[i - 1, j] / i / diag[j])))
        res[n] = dict(per_row_over_diag=line, info=info.get(n, {}),
                      max_row_ratio=float(np.max(np.abs(S[:, j]) ** 2) / diag[j]))
        print(f"{n:16s} " + "  ".join(f"K/N0={a:g}:{b:.3f}" for a, b in line) +
              f"   max|S|^2/diag={res[n]['max_row_ratio']:.1f}", flush=True)
    if out:
        json.dump(dict(N0=N0, KMAX=KMAX, cols=len(cA), rows=len(rA), res=res), open(out, "w"), indent=1)
    return res

if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    main(**args)
