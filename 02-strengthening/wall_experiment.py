"""
Empirical test of the marked-energy wall (Conjecture W) in Z[w].

T(k) = X^{-1/2} sum_{n' sf primary, X < N(n') <= 2X}  a(n') chi_{n'}(k P^4) w(N(n')/X),
a(n) = conj(alpha(n)) gamma_2(n)          (the paper's coefficient; 'control' drops alpha)

Unmarked (P = 1): the paper's reflection bound is  sum_{N(k)<=K} |T(k)|^2 << K + K^2/X  (wall at K = X).
Marked by a fixed prime P of norm p (twist chi_{n'}(P)^4, conductor P): reflection gives K + K^2 p / X
(wall at K = X/p).  W(w) predicts no growth of the per-row energy past the wall.
We report  E(K)/#rows(K)  as K/wall grows.
"""
import math, cmath, sys, time
import numpy as np
from eisenstein import (eisenstein_primes, ResidueField, gauss_prime, primary, norm, cx, mul, Z6)

def elements_upto(K):
    R = int(2 * math.sqrt(K)) + 3
    out = []
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            n = a * a - a * b + b * b
            if 0 < n <= K:
                out.append((a, b))
    return out

def sextic_table(F, rows_a, rows_b):
    """vectorized exponent of (k/pi)_6 for all rows; -1 where pi | k."""
    q = F.q
    if F.kind == "split":
        x = (rows_a + rows_b * F.r) % q
        e = (q - 1) // 6
        res = np.ones_like(x); base = x.copy()
        while e:
            if e & 1: res = res * base % q
            base = base * base % q; e >>= 1
        w6 = (-(F.r * F.r)) % q
        lut = np.full(q, -1, dtype=np.int8)
        t = 1
        for k in range(6):
            lut[t] = k; t = t * w6 % q
        out = np.full(len(x), -1, dtype=np.int8)
        nz = x != 0
        out[nz] = lut[res[nz]]
        return out
    out = np.empty(len(rows_a), dtype=np.int8)
    for i, (a, b) in enumerate(zip(rows_a, rows_b)):
        k = F.sextic_exp((int(a), int(b)))
        out[i] = -1 if k is None else k
    return out

def main(X=3000, KMAX=30000, Ps=(1, 7, 13, 31, 61), control=False):
    t0 = time.time()
    primes = [p for p in eisenstein_primes(2 * X + 10) if p[2] != 2]
    fields = {pr[1]: ResidueField(*pr) for pr in primes}
    pnorm = {pr[1]: (pr[2] if pr[0] == "split" else pr[2] ** 2) for pr in primes}
    g2 = {p: gauss_prime(F, 2) for p, F in fields.items()}
    # squarefree primary columns with X < N <= 2X, factored over the prime list
    plist = sorted(fields, key=lambda p: pnorm[p])
    cols = []
    for e in elements_upto(2 * X):
        if norm(e) <= X: continue
        pe = primary(e)
        if pe is None or pe != e: continue
        # factor via norms
        rem = e; fac = []; ok = True; Nn = norm(e)
        for p in plist:
            if pnorm[p] > Nn: break
            F = fields[p]
            if F.red(rem) == 0 if F.kind == "split" else F.red(rem) == (0, 0):
                # divide rem by p exactly: rem * conj(p) / N(p)
                pc = (p[0] - p[1], -p[1])
                prod = mul(rem, pc)
                if prod[0] % pnorm[p] or prod[1] % pnorm[p]: continue
                rem = (prod[0] // pnorm[p], prod[1] // pnorm[p])
                fac.append(p)
                if (F.red(rem) == 0 if F.kind == "split" else F.red(rem) == (0, 0)):
                    ok = False; break       # not squarefree
                Nn //= pnorm[p]
                if Nn == 1: break
        if not ok or Nn != 1: continue
        cols.append((e, fac))
    # gamma_2 via twisted multiplicativity gamma_2(ab) = gamma_2(a) gamma_2(b) chi_b(a)^4  (validated in eisenstein.py)
    def chi_exp(p, x):
        k = fields[p].sextic_exp(x); return k
    coef = []
    for e, fac in cols:
        g = 1 + 0j; acc = (1, 0)
        for p in fac:
            k = chi_exp(p, acc) if acc != (1, 0) else 0     # chi_p(acc)^4 * ...
            g *= g2[p] * Z6[(4 * k) % 6]
            acc = mul(acc, p)
        al = cx(e) / abs(cx(e))
        coef.append(g if control else g * al.conjugate())
    coef = np.array(coef)
    Ncol = np.array([norm(e) for e, _ in cols], dtype=float)
    wts = np.sin(math.pi * (Ncol / X - 1)) ** 2
    print(f"columns: {len(cols)}  primes: {len(primes)}  setup {time.time()-t0:.0f}s", flush=True)
    rows = elements_upto(KMAX)
    ra = np.array([r[0] for r in rows], dtype=np.int64); rb = np.array([r[1] for r in rows], dtype=np.int64)
    rN = ra * ra - ra * rb + rb * rb
    order = np.argsort(rN); ra, rb, rN = ra[order], rb[order], rN[order]
    # exponent tables for all primes dividing some column
    need = sorted({p for _, fac in cols for p in fac}, key=lambda p: pnorm[p])
    tab = {p: sextic_table(fields[p], ra, rb) for p in need}
    print(f"rows: {len(rows)}  tables {time.time()-t0:.0f}s", flush=True)
    results = {}
    for P in Ps:
        if P == 1:
            Pel, pP = None, 1
        else:
            cand = [p for p in plist if pnorm[p] == P]
            Pel, pP = cand[0], P
        T = np.zeros(len(rows), dtype=complex)
        for j, (e, fac) in enumerate(cols):
            if Pel is not None and Pel in fac:
                continue                                    # keep (n', P) = 1
            ex = np.zeros(len(rows), dtype=np.int16); mask = np.ones(len(rows), dtype=bool)
            for p in fac:
                t = tab[p]; mask &= (t >= 0); ex += t
            phase = coef[j] * wts[j]
            if Pel is not None:
                # chi_{n'}(P)^4 : constant per column
                kk = 0
                for p in fac:
                    kk += fields[p].sextic_exp(Pel)
                phase *= Z6[(4 * kk) % 6]
            T += np.where(mask, phase * np.exp(1j * math.pi / 3 * ex), 0)
        T /= math.sqrt(X)
        E = np.cumsum(np.abs(T) ** 2)
        wall = X / pP
        rep = []
        for K in [wall * f for f in (0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)]:
            if K > KMAX or K < 1: continue
            i = np.searchsorted(rN, K, side="right")
            if i < 50: continue
            rep.append((K / wall, E[i - 1] / i))
        base = np.mean(np.abs(T[: max(50, np.searchsorted(rN, wall))]) ** 2)
        results[P] = rep
        print(f"P-norm {pP:4d} (wall K={wall:.0f}):  per-row energy at K/wall = " +
              ", ".join(f"{f:g}:{v:.3f}" for f, v in rep), flush=True)
    return results

if __name__ == "__main__":
    control = "--control" in sys.argv
    print("coefficients:", "gamma_2 (no alpha twist, CONTROL)" if control else "conj(alpha) gamma_2 (paper)")
    main(control=control)
