"""Attack E calibration (exact, noise-free): quadratic twists of the row scalar zeta alone, at small conductor, as Q grows.

  F_cl(nu) = sum_{R in class cl, R primary sf, N R in (Q/2,Q]} Om(N R/Q) x(R) (nu/R)_2 ,  classes cl = R mod 4 (12),
  Rat(nu)  = sum_cl |F_cl(nu)|^2 / sum_R Om^2 |x|^2 |(nu/R)_2|^2       (random: 1; structurally small: << 1)
for x = zeta = psi1 conj(gamma2)^4 (Lemma Z), conj(gamma2) alpha-bar, gamma2, and "1" (degree 1).
A Voronoi/dual sum of length C N(nu)^k / Q makes Rat -> 0 at fixed nu as Q grows (k = 1 for "1"; k = 2 expected for the
cubic-Gauss-sum sequences, which are GL(2)-metaplectic).  Compare with rho and pi from E_ls.py.
usage: python E_zeta_twist.py NUMAX Q1 Q2 ...
"""
import numpy as np, sys, time, json
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *
import scipy.sparse as sps

t0 = time.time()
NUMAX = int(sys.argv[1]); Qs = [int(x) for x in sys.argv[2:]]
nu_a, nu_b, nu_N = nu_columns(NUMAX)
res = {}
for Q in Qs:
    P = primes_upto(Q)
    S = [t for t in squarefree_products(P, Q) if t[2] > Q / 2]
    Ra = np.array([t[0] for t in S], np.int64); Rb = np.array([t[1] for t in S], np.int64); RN = np.array([t[2] for t in S], np.int64)
    fac = [t[3] for t in S]
    used = sorted(set(k for f in fac for k in f))
    g2 = np.zeros(len(P), complex); e_lam = np.zeros(len(P), int); e_4 = np.zeros(len(P), int); e_m1 = np.zeros(len(P), int)
    for k in used:
        g2[k] = gamma_prime(P[k], 2); e_lam[k] = sym1(P[k], *LAM); e_4[k] = sym1(P[k], 4, 0); e_m1[k] = sym1(P[k], -1, 0)
    cross = {}; G2 = np.zeros(len(S), complex); psi1 = np.zeros(len(S), complex)
    for t, f in enumerate(fac):
        g = 1 + 0j
        for x in range(len(f)):
            g *= g2[f[x]]
            for y in range(x + 1, len(f)):
                key = (f[x], f[y])
                if key not in cross:
                    cross[key] = sym1(P[f[x]], P[f[y]].a, P[f[y]].b) + sym1(P[f[y]], P[f[x]].a, P[f[x]].b)
                g *= Z6[(2 * cross[key]) % 6]
        G2[t] = g
        psi1[t] = np.conj(Z6[sum(e_lam[k] for k in f) % 6]) * np.conj(Z6[sum(e_4[k] for k in f) % 6]) * Z6[sum(e_m1[k] for k in f) % 6]
    zeta = psi1 * np.conj(G2) ** 4
    Om = bump(RN / Q, .5, 1)
    alR = np.array([cval(a, b) for a, b in zip(Ra, Rb)]); alR /= np.abs(alR)
    rows = np.repeat(np.arange(len(S)), [len(f) for f in fac]); cols = np.concatenate([np.array(f) for f in fac])
    Inc = sps.csr_matrix((np.ones(len(rows), np.int64), (rows, cols)), shape=(len(S), len(P)))
    Ep = np.zeros((len(P), len(nu_N)), np.int8); Zp = np.zeros((len(P), len(nu_N)), np.int8)
    for k in used:
        e = P[k].symbol(nu_a, nu_b); Zp[k] = e < 0; Ep[k] = np.where(e >= 0, e, 0) % 2
    par = (Inc @ Ep.astype(np.int64)) % 2; zer = (Inc @ Zp.astype(np.int64)) > 0
    TW = np.where(zer, 0, 1 - 2 * par).astype(np.float64)
    lab = np.unique((Ra % 4) * 4 + (Rb % 4), return_inverse=True)[1]; ncl = lab.max() + 1
    out = {}
    for name, x in (("zeta", zeta), ("conj(g2)*conj(alpha)", np.conj(G2) * np.conj(alR)), ("g2", G2), ("1", np.ones(len(S)))):
        num = np.zeros(len(nu_N))
        for c in range(ncl):
            m = lab == c
            num += np.abs(TW[m].T @ (Om[m] * x[m])) ** 2
        den = (Om ** 2 * np.abs(x) ** 2) @ TW ** 2
        out[name] = num / den
    res[Q] = out
    print("Q=%d: %d ideals [%.0fs]" % (Q, len(S), time.time() - t0), flush=True)
edges = [7, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
edges = [e for e in edges if e <= NUMAX] + [NUMAX + 1]
for name in res[Qs[0]]:
    print("\n%s: mean Rat per bin of N(nu)  (columns: Q = %s)" % (name, Qs))
    for b in range(len(edges) - 1):
        sel = (nu_N >= edges[b]) & (nu_N < edges[b + 1])
        if not sel.any(): continue
        print("  [%5d,%5d) #nu=%4d  " % (edges[b], edges[b + 1], sel.sum()) + "  ".join("%7.3f" % res[Q][name][sel].mean() for Q in Qs) +
              "    N^2/Q at bin centre: " + " ".join("%.2g" % (edges[b] * edges[b + 1] / Q) for Q in Qs))
json.dump({str(Q): {k: v.tolist() for k, v in o.items()} for Q, o in res.items()} | {"nu_N": nu_N.tolist()}, open(DATA + "zeta_twist.json", "w"))
print("[%.0fs]" % (time.time() - t0))
