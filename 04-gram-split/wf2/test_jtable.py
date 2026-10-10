"""Positive control + automorphy fingerprint: sextic twists chi_nu(m)^j = (m/nu)_6^j, nu prime, j = 0..5,
for x = rho (known: j = 1 coherent, the multiplier of theta_6 at a/nu) and x = pi = rho zeta (unknown),
restricted to squarefree rows, unit classes separately (f = 3).
Also: degree-1 calibrations for the quadratic twist (nu/m)_2 on ALL rows vs squarefree rows (exact),
to show what 'structurally small' looks like and how much the squarefree restriction hides it.
"""
import numpy as np, sys, time
from a2lib import *
t0 = time.time()
Q = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
NUMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
R = Rows(Q); W = R.w; K = len(R.scales)
rows = np.where(R.sf)[0]; n = len(rows); Om = R.Om[rows]
lab = np.unique((R.ma[rows] % 3) * 3 + R.mb[rows] % 3, return_inverse=True)[1]; nc = lab.max() + 1
X = {"rho": R.V[:, rows], "pi": R.V[:, rows] * R.zetaZ[rows][None, :]}
drow = {k: cross_est(v.T, W) for k, v in X.items()}
Pn = [Pr for Pr in primes_upto(NUMAX) if Pr.N >= 7]
edges = [7, 32, 128, 512, 2048, NUMAX + 1]
res = {}
for Pr in Pn:
    e = Pr.symbol(R.ma[rows], R.mb[rows]); ok = e >= 0
    b = int(np.searchsorted(edges, Pr.N, side="right") - 1)
    for j in range(6):
        psi = np.where(ok, Z6[(j * np.where(ok, e, 0)) % 6], 0)
        for name, Xs in X.items():
            T = np.zeros((nc, K), complex)
            Y = Xs * (Om * psi)[None, :]
            for i in range(K):
                T[:, i] = np.bincount(lab, Y[i].real, nc) + 1j * np.bincount(lab, Y[i].imag, nc)
            num = cross_est(T, W).sum(); den = (Om ** 2 * np.abs(psi) ** 2 * drow[name]).sum()
            res.setdefault((name, j, b), []).append(num / den)
print("Q=%d, sf rows %d, unit classes separately.  sextic twist (m/nu)_6^j, nu prime.  mean ratio (s.e. over nu) [#nu]" % (Q, n))
for name in X:
    print("x = %s" % name)
    for b in range(len(edges) - 1):
        if (name, 0, b) not in res: continue
        print("  N(nu) in [%4d,%4d): " % (edges[b], edges[b + 1]) + "  ".join(
            "j=%d %6.2f(%4.2f)" % (j, np.mean(res[(name, j, b)]), np.std(res[(name, j, b)]) / np.sqrt(len(res[(name, j, b)])))
            for j in range(6)) + "  [%d]" % len(res[(name, 0, b)]))
print("[%.1fs]" % (time.time() - t0), flush=True)

# ---- degree-1 calibration of the QUADRATIC twist on all rows vs squarefree rows (exact sequences)
print("\nExact degree-1 calibrations, quadratic twist (nu/m)_2 = prod_p chi_p(nu)^{3 v_p(m)}, nu squarefree primary, unit classes separately")
allrows = np.arange(R.n); labA = np.unique((R.ma % 3) * 3 + R.mb % 3, return_inverse=True)[1]
nus = [(a, b, N) for (a, b, N, idx) in squarefree_products(primes_upto(4 * Q), 4 * Q) if N >= 7]
nu_a = np.array([x[0] for x in nus], np.int64); nu_b = np.array([x[1] for x in nus], np.int64); nu_N = np.array([x[2] for x in nus])
Ek = {}
TWall = np.ones((R.n, len(nus)))
for k in R.used:
    e = R.P[k].symbol(nu_a, nu_b)
    v = np.where(e < 0, 0.0, np.where((3 * e) % 6 == 0, 1.0, -1.0))
    mult = R.Inc[:, k]
    for t in np.where(mult > 0)[0]:
        TWall[t] *= v ** int(mult[t])
seqs = {"1": np.ones(R.n), "conj(alpha(m))^2": np.conj(R.alpha_m) ** 2,
        "gamma2(R) [sf]": np.where(R.sf, R.G[2], 0)}
edges2 = [7, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 4 * Q + 1]
for rset_name, rset in (("ALL rows", np.ones(R.n, bool)), ("SQUAREFREE rows", R.sf)):
    print("  %s:" % rset_name)
    for sname, x in seqs.items():
        if "sf" in sname and rset_name == "ALL rows": continue
        y = R.Om * x * rset
        S = np.zeros((len(nus), labA.max() + 1), complex)
        for c in range(labA.max() + 1):
            sel = labA == c
            S[:, c] = TWall[sel].T @ y[sel]
        rr = (np.abs(S) ** 2).sum(1) / ((np.abs(y) ** 2) @ (TWall ** 2))
        line = []
        for b in range(len(edges2) - 1):
            sel = (nu_N >= edges2[b]) & (nu_N < edges2[b + 1])
            if sel.sum(): line.append("%d:%.3f" % (edges2[b], rr[sel].mean()))
        print("    %-18s " % sname + "  ".join(line))
print("[%.1fs]" % (time.time() - t0))
