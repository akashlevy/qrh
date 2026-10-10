"""Item 4: unreflected theta pairing with the TRUE completed row (exp_pairing2.py 'nu' variant), K12 data,
weighted cross-scale estimator.
  B_m^{(q)} = sum_{c sf primary, N c in [Z,2Z]} W(Nc/Z) c_A chi_c(q m) / sqrt(N c),
  c_A = gamma_2(c) conj(alpha(c)) conj(G(c)) nu(c),   G(c) = conj(chi_c(4)) gamma_3(c),  nu(c) = conj(chi_c(-1)),
the extra fixed element q in chi_c(q m) = chi_c(q) chi_c(m) runs over 10 ray characters nu.  (paper78 l.3440ff;
wf/V_gaps.md 1.6).  This B already contains zeta_m implicitly: its reflection is zeta_m sum_mu beta_mu psi_m(mu).
  T = sum_m Om rho(m) B_m.   R_rand = E|T|^2 / sum Om^2 |rho|^2 |B|^2 ;  R_CS = E|T|^2 / (sum Om^2|rho|^2)(sum Om^2|B|^2).
For comparison: the 'zeta-free' pairing <rho, conj(zeta) B> over squarefree rows (the sum Lemma (V) would bound).
"""
import numpy as np, sys, time
from a2lib import *
t0 = time.time()
Q = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
Zs = [int(z) for z in sys.argv[2:]] or [3000, 10000, 20000]
R = Rows(Q); W = R.w; K = len(R.scales)
Om = R.Om; n = R.n
P = primes_upto(2 * max(Zs))
print("Q=%d rows %d; primes up to %d: %d  [%.0fs]" % (Q, n, 2 * max(Zs), len(P), time.time() - t0), flush=True)
g2 = np.array([gamma_prime(Pr, 2) for Pr in P]); g3 = np.array([gamma_prime(Pr, 3) for Pr in P])
e4 = np.array([sym1(Pr, 4, 0) for Pr in P]); em1 = np.array([sym1(Pr, -1, 0) for Pr in P])
print("gauss sums done [%.0fs]" % (time.time() - t0), flush=True)
twists = [(1, 0), (7, 3), (4, 3), (-2, 3), (13, 0), (10, 3), (-5, -3), (1, 3), (7, 0), (-11, -3)]
Jq = {}
for (qa, qb) in twists:
    xa, xb = mul(R.ma, R.mb, qa, qb)
    Jq[(qa, qb)] = np.stack([Pr.symbol(xa, xb) for Pr in P]).astype(np.int8)
print("symbols done [%.0fs]" % (time.time() - t0), flush=True)
S = squarefree_products(P, 2 * max(Zs))
crossc = {}
def cr(ii, jj):
    if (ii, jj) not in crossc:
        crossc[(ii, jj)] = (sym1(P[ii], P[jj].a, P[jj].b), sym1(P[jj], P[ii].a, P[ii].b))
    return crossc[(ii, jj)]
sfrows = R.sf; zb = np.where(sfrows, np.conj(R.zetaZ), 0)
dr = cross_est(R.V.T, W)                       # per-row |rho|^2 estimate
for Z in Zs:
    B = {k: np.zeros(n, complex) for k in twists}
    nc = 0
    for (sa, sb, sN, idx) in S:
        if not (Z <= sN <= 2 * Z): continue
        w = bump(np.array([sN / Z]), 1, 2)[0] / np.sqrt(sN)
        if w == 0: continue
        nc += 1
        g = 1.0 + 0j; g3c = 1.0 + 0j; ex4 = 0; exm1 = 0
        for t, ii in enumerate(idx):
            g *= g2[ii]; g3c *= g3[ii]; ex4 += e4[ii]; exm1 += em1[ii]
            for jj in idx[t + 1:]:
                e1, e2 = cr(ii, jj)
                g *= Z6[(2 * (e1 + e2)) % 6]; g3c *= Z6[(3 * (e1 + e2)) % 6]
        nu = Z6[ex4 % 6] * np.conj(g3c) * np.conj(Z6[exm1 % 6])      # conj(G(c) chi_c(-1))
        al = cval(sa, sb); al = np.conj(al / abs(al))
        cA = w * g * al * nu
        for k in twists:
            rows_ = Jq[k][list(idx)]
            bad = (rows_ < 0).any(axis=0); e = rows_.astype(np.int64).sum(axis=0) % 6
            val = Z6[e]; val[bad] = 0
            B[k] += cA * val
    print("\nZ=%d: %d moduli c; L = Q^2/Z = %.0f" % (Z, nc, Q ** 2 / Z))
    res = {"true <rho,B>": [], "zeta-free <rho, conj(zeta) B> (sf rows)": []}
    for k in twists:
        for lab, Bv, rowmask in (("true <rho,B>", B[k], np.ones(n, bool)),
                                 ("zeta-free <rho, conj(zeta) B> (sf rows)", B[k] * zb, sfrows)):
            T = ((R.V * (Om * Bv)[None, :]) [:, rowmask]).sum(1)
            est = cross_est(T, W)
            rnd = (Om ** 2 * dr * np.abs(Bv) ** 2)[rowmask].sum()
            cs2 = (Om ** 2 * dr)[rowmask].sum() * (Om ** 2 * np.abs(Bv) ** 2)[rowmask].sum()
            reps = [cross_est(T[[j for j in range(K) if j != i]], W[[j for j in range(K) if j != i]]) / rnd for i in range(K)]
            se = np.sqrt((K - 1) / K * ((np.array(reps) - np.mean(reps)) ** 2).sum())
            res[lab].append((est / rnd, se, est / cs2))
    for lab, L in res.items():
        L = np.array(L)
        print("  %-42s R_rand per twist: %s" % (lab, " ".join("%5.2f" % x for x in L[:, 0])))
        print("  %-42s mean R_rand = %.3f (s.e. over twists %.3f; typical jackknife s.e. %.2f);  mean R_CS = %.2e ;  1/#rows = %.2e"
              % ("", L[:, 0].mean(), L[:, 0].std() / np.sqrt(len(L)), np.median(L[:, 1]), L[:, 2].mean(), 1 / n))
    print("[%.0fs]" % (time.time() - t0), flush=True)
