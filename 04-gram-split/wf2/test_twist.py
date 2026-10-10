"""Item 3 (structural test): quadratic twists of pi = rho * zeta versus rho alone.

  F'_c(nu) = sum_{m in class c, m squarefree} Om(m) pi(m) (nu/m)_2 ,   nu = lambda^4 mu  squarefree primary, (nu,6)=1,
  (lambda^4 mu / m)_2 = prod_{p | m} chi_p(lambda^4 mu)^3 = (nu/R_m)_2  (lambda^4 is a square).
  Rat(nu; f) = sum_{classes c mod f} |F'_c(nu)|^2 / sum_m Om^2 |pi(m)|^2 |(nu/m)_2|^2
 f = 1: the plain sum;  f = 3: the six unit classes separately (absorbs the unknown kappa_0(u));  f = 12: all ray
 characters mod 12 (mean over psi of |<pi, (nu/.)_2 psi>|^2 / random level).
 Random model: E Rat = 1.  Coherent: >> 1.  Structurally small (dual sum empty): << 1.
 Same for x = rho, and for exact calibration sequences (no noise).
"""
import numpy as np, sys, time
from a2lib import *
t0 = time.time()
Q = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
NUMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 4 * Q
R = Rows(Q); W = R.w; K = len(R.scales)
rows = np.where(R.sf)[0]; n = len(rows); Om = R.Om[rows]
print("Q=%d sf rows %d  scales %s  [%.1fs]" % (Q, n, R.Ys[R.scales].astype(int).tolist(), time.time() - t0), flush=True)

# --- columns nu: squarefree primary, coprime to 6 (= products of primary primes), N(nu) <= NUMAX
Pn = primes_upto(NUMAX)
nus = [(a, b, N) for (a, b, N, idx) in squarefree_products(Pn, NUMAX) if N >= 7]
nu_a = np.array([x[0] for x in nus], np.int64); nu_b = np.array([x[1] for x in nus], np.int64)
nu_N = np.array([x[2] for x in nus]); nnu = len(nus)
print("columns nu: %d (N(nu) in [7, %d])" % (nnu, NUMAX), flush=True)
used = R.used
E = np.zeros((len(R.P), nnu), np.int64); Z = np.zeros((len(R.P), nnu), np.int64)
for k in used:
    e = R.P[k].symbol(nu_a, nu_b); Z[k] = (e < 0); E[k] = np.where(e >= 0, e, 0)
Inc = R.Inc[rows].astype(np.int64)
par = (Inc @ (3 * E)) % 6            # exponent of the quadratic symbol (0 or 3)
zero = (Inc @ Z) > 0
TW = np.where(zero, 0.0, np.where(par == 0, 1.0, -1.0))      # (n, nnu)
print("twist matrix built [%.1fs]" % (time.time() - t0), flush=True)


def classes(f):
    if f == 1: return np.zeros(n, int)
    key = (R.ma[rows] % f) * f + (R.mb[rows] % f)
    return np.unique(key, return_inverse=True)[1]


fs = [1, 3, 12]
labs = {f: classes(f) for f in fs}
seqs_noisy = {"rho": R.V[:, rows], "pi": R.V[:, rows] * R.zetaZ[rows][None, :]}
rng = np.random.default_rng(5)
seqs_exact = {"1": np.ones(n), "conj(alpha)^2": np.conj(R.alpha[rows]) ** 2, "gamma2": R.G[2][rows],
              "gamma1": R.G[1][rows], "zeta": R.zetaZ[rows], "rand": np.exp(2j * np.pi * rng.random(n))}
TW2 = TW ** 2
out = {}
for name, X in seqs_noisy.items():
    # per-row |x|^2 estimate (cross-scale), denominator per nu
    d_row = cross_est((X * 1.0).T, W)                    # (n,)
    den = (Om ** 2 * d_row) @ TW2                          # (nnu,)
    for f in fs:
        lab = labs[f]; nc = lab.max() + 1
        num = np.zeros(nnu); reps = np.zeros((K, nnu))
        # T[nu, c, i]
        T = np.zeros((nnu, nc, K), complex)
        for c in range(nc):
            sel = lab == c
            T[:, c, :] = TW[sel].T @ (X[:, sel] * Om[sel]).T
        num = cross_est(T, W).sum(1)
        for i in range(K):
            keep = np.array([j for j in range(K) if j != i])
            dr = cross_est(X[keep].T, W[keep]); dn = (Om ** 2 * dr) @ TW2
            reps[i] = cross_est(T[:, :, keep], W[keep]).sum(1) / dn
        out[(name, f)] = (num / den, reps)
    print("%s done [%.1fs]" % (name, time.time() - t0), flush=True)
for name, x in seqs_exact.items():
    den = (Om ** 2 * np.abs(x) ** 2) @ TW2
    for f in fs:
        lab = labs[f]; nc = lab.max() + 1
        S = np.zeros((nnu, nc), complex)
        for c in range(nc):
            sel = lab == c
            S[:, c] = TW[sel].T @ (Om[sel] * x[sel])
        out[(name, f)] = ((np.abs(S) ** 2).sum(1) / den, None)

# --- report in dyadic bins of N(nu)
edges = [7, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
edges = [e for e in edges if e <= NUMAX] + [NUMAX + 1]
print("\nRat(nu) averaged over nu in bins of N(nu); noisy rows: mean (s.e. across nu) [bin-mean jackknife-over-scales s.e.]")
for f in fs:
    print("\n--- f = %d ---" % f)
    hdr = "%-15s" % "N(nu) bin" + "".join("%-22s" % s for s in ["#nu  N^2/Q", "rho", "pi"]) + "".join("%-10s" % s for s in seqs_exact)
    print(hdr)
    for b in range(len(edges) - 1):
        sel = (nu_N >= edges[b]) & (nu_N < edges[b + 1])
        if sel.sum() == 0: continue
        mid = np.sqrt(edges[b] * edges[b + 1])
        line = "[%5d,%5d)  " % (edges[b], edges[b + 1]) + "%-22s" % ("%4d  %8.3g" % (sel.sum(), mid ** 2 / Q))
        for name in ("rho", "pi"):
            r, reps = out[(name, f)]
            m = r[sel].mean(); se = r[sel].std() / np.sqrt(sel.sum())
            jk = reps[:, sel].mean(1); jse = np.sqrt((K - 1) / K * ((jk - jk.mean()) ** 2).sum())
            line += "%-22s" % ("%6.2f(%4.2f)[%4.2f]" % (m, se, jse))
        for name in seqs_exact:
            line += "%-10s" % ("%7.3f" % out[(name, f)][0][sel].mean())
        print(line)
np.savez("twist_Q%d.npz" % Q, nu_N=nu_N, nu_a=nu_a, nu_b=nu_b,
         **{"%s_f%d" % (k[0], k[1]): v[0] for k, v in out.items()})
print("[%.1fs]" % (time.time() - t0))
