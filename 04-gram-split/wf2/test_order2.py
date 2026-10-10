"""Second-order phase test at prime rows (motivated by Broeker-Hoffstein Conj. 5.7: tau_6(p)^2 = (6th root) x cubic Gauss
sum x a 'patternless' algebraic factor from a finite set).  If pi(p)^2 (or rho(p)^2 gamma_2^{-+1}) had a phase
concentrated (non-uniform) relative to alpha^k within classes mod f, the class sums of x(p)^2 h(p) would be coherent.

  Rat2(x; h, f) = sum_classes |sum_p Om x(p)^2 h(p)|^2 / sum_p Om^2 |x(p)|^4 |h|^2,
 x(p)^2 estimated by V_i V_j (i<j), |.|^2 and |x|^4 by products over 4 distinct scales (unbiased if noises independent).
 Null: 60 random +-1 twists h*eps;  z = (Rat2 - null mean)/null sd.  We flag |z| > 3 and check replication across Q.
"""
import numpy as np, sys, itertools, time
from a2lib import *
t0 = time.time()
rng = np.random.default_rng(11)
fs = [1, 3, 12]; ks = list(range(-9, 10))
allz = {}
for Q in (1500, 3000, 6000):
    R = Rows(Q); W = R.w; K = len(R.scales)
    rows = np.where(R.sf & (R.nfac == 1))[0]; Om = R.Om[rows]; n = len(rows)
    pairs = list(itertools.combinations(range(K), 2))
    wp = np.array([W[i] * W[j] for i, j in pairs])
    disj = [(a, b) for a in range(len(pairs)) for b in range(len(pairs)) if a != b and not set(pairs[a]) & set(pairs[b])]
    da = np.array([d[0] for d in disj]); db = np.array([d[1] for d in disj]); wd = wp[da] * wp[db]
    labs = {f: (np.zeros(n, int) if f == 1 else np.unique((R.ma[rows] % f) * f + R.mb[rows] % f, return_inverse=True)[1]) for f in fs}
    alpha = R.alpha_m[rows]; g2 = R.G[2][rows]
    for xname, zf in (("rho", np.ones(n)), ("pi", R.zetaZ[rows])):
        Vp = R.V[:, rows] * zf[None, :]
        PP = np.array([Vp[i] * Vp[j] for i, j in pairs])
        den0 = (Om ** 2 * (PP[da] * np.conj(PP[db])).real * wd[:, None]).sum() / wd.sum()

        def rat2(h, lab):
            nc = lab.max() + 1; Y = PP * (Om * h)[None, :]
            U = np.array([np.bincount(lab, y.real, nc) + 1j * np.bincount(lab, y.imag, nc) for y in Y])
            return ((U[da] * np.conj(U[db])).real * wd[:, None]).sum() / wd.sum() / (den0 * np.mean(np.abs(h) ** 2))
        null = {f: np.array([rat2(rng.choice([-1.0, 1.0], n), labs[f]) for _ in range(60)]) for f in fs}
        es = (0,) if xname == "pi" else (-1, 0, 1)
        for e in es:
            for k in ks:
                h = alpha ** k * g2 ** e
                for f in fs:
                    v = rat2(h, labs[f])
                    allz[(xname, e, k, f, Q)] = (v, (v - null[f].mean()) / null[f].std(), null[f].max())
        print("Q=%d x=%s  prime rows %d  null mean/sd: %s  [%.0fs]" % (Q, xname, n,
              "  ".join("f=%d %.2f/%.2f (max %.1f)" % (f, null[f].mean(), null[f].std(), null[f].max()) for f in fs), time.time() - t0), flush=True)

print("\nEntries with z > 3 (excess over null), and their values at the other Q:")
seen = set()
for key, (v, z, mx) in sorted(allz.items(), key=lambda kv: -kv[1][1]):
    if z < 3: break
    base = key[:4]
    if base in seen: continue
    seen.add(base)
    xname, e, k, f, Q = key
    reps = "  ".join("Q=%d: %6.2f (z=%5.1f)" % (q, allz[base + (q,)][0], allz[base + (q,)][1]) for q in (1500, 3000, 6000))
    print("  x=%-3s e=%+d k=%+d f=%-2d | %s" % (xname, e, k, f, reps))
nt = sum(1 for k in allz if k[4] == 6000)
print("\n(number of (x,e,k,f) cells per Q: %d;  expected # with z>3 under a Gaussian null ~ %.1f per Q, heavier tails likely)" % (nt, nt * 0.00135))
# replicated: z > 2 at all three Q
rep = [b for b in set(k[:4] for k in allz) if all(allz[b + (q,)][1] > 2 for q in (1500, 3000, 6000))]
print("cells with z > 2 at ALL three Q:", sorted(rep))
print("[%.0fs]" % (time.time() - t0))
