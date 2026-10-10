"""Attack E, (V'_LS) ground truth.

  F'_c(nu) = sum_{ideals R sf} Om(N R) pi(zeta6^c R) (nu/R)_2,      pi = rho * zeta  (zeta = Lemma Z, exact),
  nu squarefree primary, (nu,6)=1, 7 <= N(nu) <= numax.   Unit classes c = 0..5 kept separate (f = 3; absorbs kappa_0).

  num(nu) = sum_c |F'_c(nu)|^2       (cross-scale unbiased, scales Y >= ymin*Q only)
  den(nu) = sum_m Om^2 |pi(m)|^2 |(nu/m)_2|^2       (random-model expectation of num)
  Rat(nu) = num/den                  (random: 1;  coherent: >> 1)
  S(L)    = sum_{N nu <= L} num(nu)
  R_LS(L) = S(L) / (Q * sum_m Om |pi|^2)          [task normalisation; large sieve: R_LS <~ (L+Q)/Q]
  R_rand(L) = S(L) / sum_{N nu <= L} den(nu)       [1 for random]
Nulls: pi -> pi * eps(m), eps = random signs (keeps |pi| and the noise structure, destroys arithmetic structure).
Exact (noise-free) calibration sequences on the same rows: random phases, zeta alone, conj(alpha)^2 (degree 1), "1".
usage: python E_ls.py Q src(E|K12) ymin_over_Q numax_over_Q [nnull]
"""
import numpy as np, sys, time, json
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *

t0 = time.time()
Q = int(sys.argv[1]); src = sys.argv[2] if len(sys.argv) > 2 else "E"
ymin = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
numax = int(float(sys.argv[4]) * Q) if len(sys.argv) > 4 else Q
nnull = int(sys.argv[5]) if len(sys.argv) > 5 else 8
R = RowsE(Q, src=src, ymin_over_Q=ymin)
W = R.w; K = len(R.scales)
print("Q=%d src=%s rows %d (sf %d) ideals %d; scales %s; init %.0fs; zeta product-formula check: max rel err %.1e on %d rows"
      % (Q, src, R.n, R.sf.sum(), R.nideal, R.Ys[R.scales].astype(int).tolist(), R.tinit, R.zeta_check_err, R.zeta_check_n), flush=True)
print("||rho||^2 (cross-scale, Om weight, all rows) = %.3f ; per-scale noise/signal: %s" % (R.sig2, " ".join("%.2f" % x for x in R.noise / R.sig2)), flush=True)

# ---- ideal x unit-class layout (squarefree ideals with all six associates present)
sfi = np.unique(R.ideal[R.sf])
pos = {v: k for k, v in enumerate(sfi)}
nI = len(sfi)
IX = -np.ones((nI, 6), np.int64)
for t in np.where(R.sf)[0]:
    IX[pos[R.ideal[t]], R.ucls[t]] = t
assert (IX >= 0).all()
rep = IX[:, 0]
OmI = R.Om[rep]
zeta = R.zetaZ
Vp = R.V[:, IX] * zeta[IX][None]                      # (K, nI, 6)  pi estimates per scale
pi2 = cross_est(np.moveaxis(Vp, 0, -1), W)             # (nI, 6)  |pi|^2 per row (unbiased)
norm_Om2 = float((OmI[:, None] ** 2 * pi2).sum())
norm_Om = float((OmI[:, None] * pi2).sum())
print("ideals %d; sum Om^2|pi|^2 = %.3f, sum Om|pi|^2 = %.3f  [%.0fs]" % (nI, norm_Om2, norm_Om, time.time() - t0), flush=True)
Cpi = np.einsum("kic,lic->kl", Vp * OmI[None, :, None], Vp.conj()).real      # sum Om pi_k conj(pi_l)
offd = np.array([[Cpi[k, l] for l in range(K) if l != k] for k in range(K)])
print("bias diagnostic: off-diagonal scale Gram sum_m Om pi_k conj(pi_l) / sum Om|pi|^2 (should all be 1 if unbiased):")
for k in range(K):
    print("   Y=%8d: " % R.Ys[R.scales[k]] + " ".join("%6.3f" % (Cpi[k, l] / norm_Om) if l != k else "   -- " for l in range(K)))

# ---- columns and twist matrix
nu_a, nu_b, nu_N = nu_columns(numax)
nnu = len(nu_N)
TW = quad_twist(R, rep, nu_a, nu_b)                    # (nI, nnu)
print("columns nu: %d (N <= %d); twist matrix [%.0fs]" % (nnu, numax, time.time() - t0), flush=True)
TW2f = np.abs(TW).astype(np.float32)                  # |(nu/R)_2|^2 in {0,1}
den_w = (OmI[:, None] ** 2 * pi2).sum(1)               # (nI,)
den = den_w @ TW2f             # (nnu,)


def twisted_sums(Xc):
    """Xc: (K, nI, 6) complex -> T (nnu, 6, K) = sum_I TW[I,nu] Om_I Xc[k,I,c]."""
    Y = (Xc * OmI[None, :, None]).transpose(1, 2, 0).reshape(nI, -1)       # (nI, 6K)
    T = np.zeros((nnu, Y.shape[1]), complex)
    for c0 in range(0, nnu, 4000):
        T[c0:c0 + 4000] = TW[:, c0:c0 + 4000].T.astype(np.float64) @ Y
    return T.reshape(nnu, 6, K)


T = twisted_sums(Vp)
num = cross_est(T, W).sum(1)                            # (nnu,)
reps = np.zeros((K, nnu)); reps_den = np.zeros((K, nnu))
for i in range(K):
    keep = np.array([j for j in range(K) if j != i])
    reps[i] = cross_est(T[:, :, keep], W[keep] / W[keep].sum()).sum(1)
    p2 = cross_est(np.moveaxis(Vp[keep], 0, -1), W[keep] / W[keep].sum())
    reps_den[i] = ((OmI[:, None] ** 2 * p2).sum(1)) @ TW2f
Rat = num / den
print("twisted sums + jackknife [%.0fs]" % (time.time() - t0), flush=True)

# ---- nulls (random signs per row)
rng = np.random.default_rng(11)
null_num = []
for k in range(nnull):
    eps = rng.choice([-1.0, 1.0], size=(nI, 6))
    Tn = twisted_sums(Vp * eps[None])
    null_num.append(cross_est(Tn, W).sum(1))
null_num = np.array(null_num)
print("nulls [%.0fs]" % (time.time() - t0), flush=True)

# ---- exact calibration sequences (no noise)
exact = {"rand": np.exp(2j * np.pi * rng.random((nI, 6))),
         "zeta": zeta[IX],
         "conj(alpha)^2": np.conj(R.alpha_m[IX]) ** 2,
         "1": np.ones((nI, 6))}
ex_rat = {}
for name, x in exact.items():
    Y = (x * OmI[:, None])
    S = np.zeros((nnu, 6), complex)
    for c0 in range(0, nnu, 4000):
        S[c0:c0 + 4000] = TW[:, c0:c0 + 4000].T.astype(np.float64) @ Y
    d = ((OmI[:, None] ** 2 * np.abs(x) ** 2).sum(1)) @ TW2f
    ex_rat[name] = ((np.abs(S) ** 2).sum(1), d)

# ---- rho (zeta-free) for comparison, and Type I in m: rows restricted to d | R
Vr = R.V[:, IX]; rho2 = cross_est(np.moveaxis(Vr, 0, -1), W)
Tr = twisted_sums(Vr); num_r = cross_est(Tr, W).sum(1); den_r = ((OmI[:, None] ** 2 * rho2).sum(1)) @ TW2f
reps_r = np.zeros((K, nnu)); reps_rd = np.zeros((K, nnu))
for i in range(K):
    keep = np.array([j for j in range(K) if j != i])
    reps_r[i] = cross_est(Tr[:, :, keep], W[keep]).sum(1)
    reps_rd[i] = ((OmI[:, None] ** 2 * cross_est(np.moveaxis(Vr[keep], 0, -1), W[keep])).sum(1)) @ TW2f
typeIm = []
Pf = R.P
for kd in range(min(12, len(Pf))):
    d = Pf[kd]
    msk = np.array([kd in R.fac[r] for r in rep])
    if msk.sum() < 20: continue
    Y = (Vp[:, msk] * OmI[None, msk, None]).transpose(1, 2, 0).reshape(msk.sum(), -1)
    Td = (TW[msk].T.astype(np.float64) @ Y).reshape(nnu, 6, K)
    nd = cross_est(Td, W).sum(1); dd = ((OmI[msk, None] ** 2 * pi2[msk]).sum(1)) @ TW2f[msk]
    ok = dd > 0
    typeIm.append((d.N, int(msk.sum()), float(nd[ok].sum() / dd[ok].sum()), float(np.median(nd[ok] / dd[ok]))))

# ---- report
out = {"Q": Q, "src": src, "ymin": ymin, "scales": R.Ys[R.scales].tolist(), "nI": nI, "nrows_sf": int(R.sf.sum()),
       "norm_Om2": norm_Om2, "norm_Om": norm_Om, "nnu": nnu}
edges = [7, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144]
edges = [e for e in edges if e < numax] + [numax + 1]
print("\nper-nu ratio Rat = sum_c|F'_c|^2 / random expectation, dyadic bins of N(nu)   [jackknife s.e. of bin mean]")
print("%-17s %6s %9s %16s %14s   %s" % ("N(nu)", "#nu", "N/Q", "pi Rat [jk]", "null [sd]", "  ".join("%-12s" % k for k in exact)))
bins = []
for b in range(len(edges) - 1):
    sel = (nu_N >= edges[b]) & (nu_N < edges[b + 1])
    if not sel.any(): continue
    m = num[sel].sum() / den[sel].sum()
    jk = np.array([reps[i, sel].sum() / reps_den[i, sel].sum() for i in range(K)])
    se = np.sqrt((K - 1) / K * ((jk - jk.mean()) ** 2).sum())
    nm = null_num[:, sel].sum(1) / den[sel].sum()
    exs = {k: v[0][sel].sum() / v[1][sel].sum() for k, v in ex_rat.items()}
    bins.append(dict(lo=edges[b], hi=edges[b + 1], n=int(sel.sum()), rat=m, se=se, null_mean=nm.mean(), null_sd=nm.std(), **{"ex_" + k: v for k, v in exs.items()}))
    print("[%6d,%6d)  %6d %9.3g %8.3f [%5.3f] %7.3f [%5.3f]   %s" % (edges[b], edges[b + 1], sel.sum(), np.sqrt(edges[b] * edges[b + 1]) / Q, m, se,
          nm.mean(), nm.std(), "  ".join("%-12.3f" % exs[k] for k in exact)))
out["bins"] = bins
rb = []
for b in bins:
    sel = (nu_N >= b["lo"]) & (nu_N < b["hi"])
    jk = np.array([reps_r[i, sel].sum() / reps_rd[i, sel].sum() for i in range(K)])
    rb.append((num_r[sel].sum() / den_r[sel].sum(), np.sqrt((K - 1) / K * ((jk - jk.mean()) ** 2).sum())))
    b["rho"], b["rho_se"] = rb[-1]
print("rho (zeta-free) per bin [jk]: " + " ".join("%.3f[%.3f]" % x for x in rb))
print("Type I in m (rows with d | R): N(d), #ideals, sum-ratio over all nu, median per-nu ratio: " + "; ".join("%d %d %.3f %.3f" % t for t in typeIm))
out["typeI_m"] = typeIm

print("\ncumulative large-sieve ratios  (S(L) = sum_{N nu <= L} sum_c |F'_c|^2)")
print("%-10s %8s %8s %10s %12s %12s %12s %12s %12s" % ("L", "logL/logQ", "#nu(L)", "#nu/Q", "R_LS", "R_LS jk se", "R_rand", "R_rand jk", "null R_rand"))
cum = []
for ex in [0.5, 0.6, 0.7, 0.75, 0.79, 0.85, 0.9, 0.95, 1.0]:
    L = Q ** ex
    if L > numax: continue
    sel = nu_N <= L
    if sel.sum() < 3: continue
    S = num[sel].sum(); D = den[sel].sum()
    jkS = np.array([reps[i, sel].sum() for i in range(K)]); jkD = np.array([reps_den[i, sel].sum() for i in range(K)])
    jkR = jkS / jkD
    seLS = np.sqrt((K - 1) / K * ((jkS - jkS.mean()) ** 2).sum()) / (Q * norm_Om)
    seR = np.sqrt((K - 1) / K * ((jkR - jkR.mean()) ** 2).sum())
    nR = null_num[:, sel].sum(1) / D
    cum.append(dict(ex=ex, L=L, nnu=int(sel.sum()), R_LS=S / (Q * norm_Om), R_LS_se=seLS, R_rand=S / D, R_rand_se=seR,
                    null_mean=nR.mean(), null_sd=nR.std(), D_over_Qnorm=D / (Q * norm_Om)))
    print("Q^%-7.2f %8.2f %8d %10.4f %12.5f %12.5f %12.3f %12.3f %8.3f(%.3f)" % (ex, ex, sel.sum(), sel.sum() / Q, S / (Q * norm_Om), seLS,
          S / D, seR, nR.mean(), nR.std()))
out["cum"] = cum

# ---- extremes
o = np.argsort(-Rat)
nullR = null_num / den[None]
print("\nlargest Rat(nu): " + " ".join("%.1f(N=%d)" % (Rat[k], nu_N[k]) for k in o[:10]))
print("null: mean of max Rat over draws %.1f (sd %.1f); 99.9%% quantile of null Rat %.2f; of pi Rat %.2f" % (
    nullR.max(1).mean(), nullR.max(1).std(), np.quantile(nullR, 0.999), np.quantile(Rat, 0.999)))
print("fraction Rat > 5: pi %.4f  null %.4f;  Rat > 10: pi %.4f  null %.4f" % (np.mean(Rat > 5), np.mean(nullR > 5), np.mean(Rat > 10), np.mean(nullR > 10)))
out["top"] = [(float(Rat[k]), int(nu_N[k])) for k in o[:20]]
out["null_max_mean"] = float(nullR.max(1).mean())
out["frac5"] = [float(np.mean(Rat > 5)), float(np.mean(nullR > 5))]
np.savez(DATA + "ls_Q%d_%s_y%g.npz" % (Q, src, ymin), nu_N=nu_N, nu_a=nu_a, nu_b=nu_b, num=num, den=den, reps=reps, reps_den=reps_den,
         null_num=null_num, num_r=num_r, den_r=den_r, reps_r=reps_r, reps_rd=reps_rd, **{"ex_" + k: v[0] for k, v in ex_rat.items()}, **{"exden_" + k: v[1] for k, v in ex_rat.items()})
json.dump(out, open(DATA + "ls_Q%d_%s_y%g.json" % (Q, src, ymin), "w"), indent=1, default=float)
print("[%.0fs]" % (time.time() - t0))
