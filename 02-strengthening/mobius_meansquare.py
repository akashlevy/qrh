"""
Direct numerical test of the hypothesis that gives theta = 6/7:

   sum_{0 < N(u) <= H} |A_u(D)|^2  <<  D^{1+eps} H        for H = D^kappa with kappa < 1,
   A_u(D) = sum_{n sf primary, D < N(n) <= 2D} mu(n) chi_n(u) w(N(n)/D),

i.e. the companion paper's Prop thm:ms (= the 7/8 paper's raw inverse moment) pushed
BELOW the large-sieve diagonal H = D.  Fed into the sixth-power extraction it gives
theta = (6 + 5 kappa)/12, which equals 6/7 at kappa = 6/7.

We report the per-row energy  E(H) / (#rows * D)  (diagonal value ~ mean of w^2 * #cols / D)
as H ranges from D^0.4 to D^1.2.  Rows include ALL nonzero elements, in particular
the sixth powers u = a^6 for which A_u ~ A_1 (the extraction rows).
"""
import math, sys
import numpy as np
import wall_experiment as we
from eisenstein import Z6

def run(D=3000, kappas=(0.4, 0.5, 0.6, 0.7, 6/7, 0.9, 1.0, 1.1, 1.2)):
    HMAX = int(D ** max(kappas)) + 1
    primes = [p for p in we.eisenstein_primes(2 * D + 10) if p[2] != 2]
    fields = {pr[1]: we.ResidueField(*pr) for pr in primes}
    pnorm = {pr[1]: (pr[2] if pr[0] == "split" else pr[2] ** 2) for pr in primes}
    plist = sorted(fields, key=lambda p: pnorm[p])
    cols = []
    for e in we.elements_upto(2 * D):
        if we.norm(e) <= D: continue
        pe = we.primary(e)
        if pe is None or pe != e: continue
        rem = e; fac = []; ok = True; Nn = we.norm(e)
        for p in plist:
            if pnorm[p] > Nn: break
            F = fields[p]; red = F.red(rem)
            if (red == 0 if F.kind == "split" else red == (0, 0)):
                pc = (p[0] - p[1], -p[1]); prod = we.mul(rem, pc)
                rem = (prod[0] // pnorm[p], prod[1] // pnorm[p]); fac.append(p); Nn //= pnorm[p]
                red = F.red(rem)
                if (red == 0 if F.kind == "split" else red == (0, 0)): ok = False; break
                if Nn == 1: break
        if ok and Nn == 1: cols.append((e, fac))
    mu = np.array([(-1) ** len(f) for _, f in cols], float)
    Ncol = np.array([we.norm(e) for e, _ in cols], float)
    wts = np.sin(math.pi * (Ncol / D - 1)) ** 2
    diag = float(np.sum(wts ** 2)) / D
    rows = we.elements_upto(HMAX)
    ra = np.array([r[0] for r in rows], np.int64); rb = np.array([r[1] for r in rows], np.int64)
    rN = ra * ra - ra * rb + rb * rb; o = np.argsort(rN); ra, rb, rN = ra[o], rb[o], rN[o]
    need = sorted({p for _, f in cols for p in f}, key=lambda p: pnorm[p])
    tab = {p: we.sextic_table(fields[p], ra, rb) for p in need}
    A = np.zeros(len(rows), complex)
    for j, (e, fac) in enumerate(cols):
        ex = np.zeros(len(rows), np.int16); mask = np.ones(len(rows), bool)
        for p in fac:
            t = tab[p]; mask &= (t >= 0); ex += t
        A += np.where(mask, mu[j] * wts[j] * np.exp(1j * math.pi / 3 * ex), 0)
    A /= math.sqrt(D)
    E = np.cumsum(np.abs(A) ** 2)
    print(f"D={D}: {len(cols)} squarefree columns; diagonal per-row energy ~ {diag:.4f}; rows up to N(u) <= {HMAX}")
    for k in kappas:
        H = D ** k
        i = np.searchsorted(rN, H, side="right")
        print(f"  kappa={k:.4f}  H={H:9.1f}  rows={i:7d}  per-row energy / diagonal = {E[i-1]/i/diag:.3f}")
    # the extraction rows themselves: |A_1|^2 and |A_{u}|^2 for u = 1 (normalized)
    print(f"  |A_1(D)|^2 / diag = {abs(A[np.argmin(rN + (ra != 1) * 10**9)])**2 / diag:.3f}")

if __name__ == "__main__":
    run(D=int(sys.argv[1]) if len(sys.argv) > 1 else 3000)
