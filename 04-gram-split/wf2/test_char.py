"""Item 2: is pi(m) = rho(m) zeta_m 'simpler' than rho?  Coherence fingerprints against simple sequences.

Statistic (first order), for a row set S, a test function h and a class partition mod f:
   Rat(x; h, f) = sum_classes |sum_{m in class} Om x(m) h(m)|^2 / sum_m Om^2 |x(m) h(m)|^2
 = mean over the characters psi of (O/f)^x of |<x, h psi>|^2 / random level.
 random/independent phases -> 1;  x h constant on classes -> ~ (#rows per class);  0 = structurally cancelling.
For x = rho, pi (noisy) numerator and denominator are unbiased weighted cross-scale estimates (distinct scales);
s.e. = leave-one-scale-out jackknife.  Calibration sequences are exact.
Second order (prime rows): x = rho(p)^2 estimated by V_i V_j (i<j), |.|^2 by products over 4 distinct scales.
"""
import numpy as np, sys, time, itertools
from a2lib import *
t0 = time.time()
Q = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
R = Rows(Q)
W = R.w; K = len(R.scales)
print("Q=%d  scales Y=%s  sig2=%.2f  noise/scale=%s" % (Q, R.Ys[R.scales].astype(int).tolist(), R.sig2, np.round(R.noise, 1).tolist()))
rng = np.random.default_rng(0)


def classes(rows, f):
    if f == 1: return np.zeros(len(rows), int)
    key = (R.ma[rows] % f) * f + (R.mb[rows] % f)
    _, lab = np.unique(key, return_inverse=True); return lab


def rat_noisy(X, h, lab, Om):
    """X: (K, n) per-scale x; h: (n,); lab: class labels.  Returns ratio, jackknife se."""
    nc = lab.max() + 1
    Y = X * (Om * h)[None, :]
    T = np.zeros((nc, X.shape[0]), complex)
    for i in range(X.shape[0]): T[:, i] = np.bincount(lab, Y[i].real, nc) + 1j * np.bincount(lab, Y[i].imag, nc)
    D = (X * (Om * np.abs(h))[None, :]).T          # (n, K)
    def r(idx):
        w = W[idx]
        return cross_est(T[:, idx], w).sum() / cross_est(D[:, idx], w).sum()
    est = r(np.arange(K))
    reps = np.array([r(np.array([j for j in range(K) if j != i])) for i in range(K)])
    return est, np.sqrt((K - 1) / K * ((reps - reps.mean()) ** 2).sum())


def rat_exact(x, h, lab, Om):
    y = Om * x * h; nc = lab.max() + 1
    S = np.bincount(lab, y.real, nc) + 1j * np.bincount(lab, y.imag, nc)
    return (np.abs(S) ** 2).sum() / (np.abs(y) ** 2).sum()


sf = np.where(R.sf)[0]
pr = np.where(R.sf & (R.nfac == 1))[0]
alpha_m = R.alpha_m
ks = list(range(-6, 7)); fs = [1, 3, 12, 36]
for lab_rows, rows in (("ALL squarefree rows", sf), ("PRIME rows", pr)):
    Om = R.Om[rows]
    X_rho = R.V[:, rows]; X_pi = R.V[:, rows] * R.zetaZ[rows][None, :]
    cal = {"1": np.ones(len(rows)), "conj(alpha)^2": np.conj(R.alpha[rows]) ** 2, "gamma2": R.G[2][rows],
           "gamma1": R.G[1][rows], "zeta": R.zetaZ[rows],
           "rand": np.exp(2j * np.pi * rng.random(len(rows)))}
    print("\n==== %s (n=%d, Q=%d).  Rat(x; alpha(m)^k, classes mod f): estimate (jackknife se)" % (lab_rows, len(rows), Q))
    labs = {f: classes(rows, f) for f in fs}
    print("  rows per class: " + "  ".join("f=%d: %.1f" % (f, len(rows) / (labs[f].max() + 1)) for f in fs))
    for name, X in (("rho", X_rho), ("pi ", X_pi)):
        print("  x = %s" % name)
        for k in ks:
            h = alpha_m[rows] ** k
            line = []
            for f in fs:
                e, s = rat_noisy(X, h, labs[f], Om); line.append("%8.2f(%5.2f)" % (e, s))
            print("    k=%+d  " % k + "  ".join(line))
    # null: random +-1 twists (structure destroyed), same estimator
    for name, X in (("rho", X_rho), ("pi ", X_pi)):
        nul = {f: [] for f in fs}
        for sd in range(20):
            eps = rng.choice([-1.0, 1.0], len(rows))
            for f in fs: nul[f].append(rat_noisy(X, eps, labs[f], Om)[0])
        print("  NULL (x=%s, 20 random +-1 twists): " % name + "  ".join("f=%d: mean %.2f sd %.2f" % (f, np.mean(nul[f]), np.std(nul[f])) for f in fs))
    print("  exact calibration sequences: max over k in [-6,6] of Rat (argmax k), and Rat at k=0")
    for cname, x in cal.items():
        line = []
        for f in fs:
            v = [rat_exact(x, alpha_m[rows] ** k, labs[f], Om) for k in ks]
            line.append("f=%d: max %7.2f (k=%+d) k0 %7.2f" % (f, max(v), ks[int(np.argmax(v))], v[6]))
        print("    %-14s " % cname + " | ".join(line))
print("[%.1fs]" % (time.time() - t0), flush=True)

# ------------------------------------------------------------ second order on primes: rho(p)^2 * alpha^k * gamma_2^e
print("\n==== SECOND ORDER, prime rows: Rat2(h,f) = sum_cl |sum Om rho^2 h|^2 / sum Om^2 |rho|^4  (4 distinct scales)")
print("     note pi(p)^2 = -rho(p)^2 psi_1^2 alpha^{-3} gamma_2(p)  (from conj(gamma2)^8 = conj(alpha)^2 conj(gamma2)^2 = -conj(alpha)^3 gamma2 at primes)")
rows = pr; Om = R.Om[rows]; Vp = R.V[:, rows]
pairs = list(itertools.combinations(range(K), 2))
PP = np.array([Vp[i] * Vp[j] for i, j in pairs])          # (npairs, n) unbiased for rho^2
wp = np.array([W[i] * W[j] for i, j in pairs])
disj = [(a, b) for a in range(len(pairs)) for b in range(len(pairs)) if a != b and not set(pairs[a]) & set(pairs[b])]
da = np.array([d[0] for d in disj]); db = np.array([d[1] for d in disj]); wd = wp[da] * wp[db]
den = (Om ** 2 * (PP[da] * np.conj(PP[db])).real * wd[:, None]).sum() / wd.sum()


def rat2(h, lab):
    nc = lab.max() + 1; Y = PP * (Om * h)[None, :]
    U = np.array([np.bincount(lab, y.real, nc) + 1j * np.bincount(lab, y.imag, nc) for y in Y])   # (npairs, nc)
    num = ((U[da] * np.conj(U[db])).real * wd[:, None]).sum() / wd.sum()
    return num / (den * np.mean(np.abs(h) ** 2))


labs2 = {f: classes(rows, f) for f in (1, 3, 12)}
g2p = R.G[2][rows]; psi1sq = R.psi1[rows] ** 2
for e in (-1, 0, 1):
    print("  h = alpha^k gamma_2^%+d:" % e)
    for k in range(-9, 10):
        h = alpha_m[rows] ** k * g2p ** e
        print("    k=%+d  " % k + "  ".join("f=%d: %6.2f" % (f, rat2(h, labs2[f])) for f in (1, 3, 12)))
nul = {f: [] for f in (1, 3, 12)}
for sd in range(30):
    eps = rng.choice([-1.0, 1.0], len(rows))
    for f in (1, 3, 12): nul[f].append(rat2(eps, labs2[f]))
print("  NULL (30 random +-1 twists): " + "  ".join("f=%d: mean %.2f sd %.2f" % (f, np.mean(nul[f]), np.std(nul[f])) for f in (1, 3, 12)))
print("[%.1fs]" % (time.time() - t0))
