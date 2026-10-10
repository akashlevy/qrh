"""Attack E: collect ls_*.json and fit the large-sieve ratio R_LS(L,Q).
usage: python E_summary.py tag1 tag2 ...   (tags like 'Q12000_E_y4')"""
import json, sys, numpy as np
DATA = "/Users/akashlevy/Documents/riemann/gram_split/vls/data/"
tags = sys.argv[1:]
pts = []
print("%-18s %6s %8s %10s %10s %9s %8s %16s %14s" % ("run", "ex", "#nu", "R_LS", "se", "R_rand", "se", "null mean(sd)", "R_LS*Q/#nu"))
for t in tags:
    d = json.load(open(DATA + "ls_%s.json" % t))
    Q = d["Q"]
    for c in d["cum"]:
        pts.append((Q, c["L"], c["R_LS"], c["R_LS_se"], c["R_rand"], c["R_rand_se"], c["nnu"]))
        print("%-18s %6.2f %8d %10.5f %10.5f %9.3f %8.3f %8.3f(%.3f) %14.3f" % (t, c["ex"], c["nnu"], c["R_LS"], c["R_LS_se"], c["R_rand"], c["R_rand_se"],
              c["null_mean"], c["null_sd"], c["R_LS"] * Q / c["nnu"]))
pts = np.array(pts)
Q, L, RLS, se, Rr, ser, nn = pts.T
ok = (RLS > 0) & (se > 0)
x1 = np.log(L / Q)[ok]; x2 = np.log(Q)[ok]; y = np.log(RLS[ok]); w = (RLS[ok] / se[ok]) ** 2
X = np.stack([np.ones(ok.sum()), x1, x2], 1)
cov = np.linalg.inv(X.T @ (w[:, None] * X)); b = cov @ X.T @ (w * y)
print("\nweighted fit  log R_LS = a + b log(L/Q) + c log Q :  a = %.3f(%.3f)  b = %.3f(%.3f)  c = %.3f(%.3f)" %
      (b[0], np.sqrt(cov[0, 0]), b[1], np.sqrt(cov[1, 1]), b[2], np.sqrt(cov[2, 2])))
print("   random model: b = 1, c = 0;   saturation (Cauchy-Schwarz / large-sieve Q-term): b = 0")
# R_rand pooled
wr = 1 / ser[ok] ** 2
print("pooled R_rand (all runs, all L, weighted) = %.3f +- %.3f   [note: cumulative points are nested, not independent]" % ((Rr[ok] * wr).sum() / wr.sum(), 1 / np.sqrt(wr.sum())))

# ---- independent dyadic bins: Rat(bin) = A (N/Q)^{-gamma};  cumulative exponent b = 1 - gamma
rows = []
for t in tags:
    d = json.load(open(DATA + "ls_%s.json" % t))
    for b in d["bins"]:
        mid = np.sqrt(b["lo"] * b["hi"]) / d["Q"]
        err = np.sqrt(b["se"] ** 2 + b["null_sd"] ** 2)
        rows.append((d["Q"], mid, b["rat"], err, b.get("rho", np.nan), b.get("rho_se", np.nan), b["n"]))
rows = np.array(rows)
Qb, xb, rb, eb, rr, re, nb = rows.T
for lab, lo, hi in (("all bins", 0, 10), ("N/Q < Q^{-0.21} (L < Q^0.79)", 0, None)):
    sel = (xb > 0) & (eb > 0)
    if hi is None:
        sel &= xb < Qb ** -0.21
    y = np.log(rb[sel]); s = eb[sel] / rb[sel]; X = np.stack([np.ones(sel.sum()), np.log(xb[sel])], 1); w = 1 / s ** 2
    cov = np.linalg.inv(X.T @ (w[:, None] * X)); p = cov @ X.T @ (w * y)
    chi2 = np.sum(w * (y - X @ p) ** 2)
    print("bin fit (%s, %d bins): log Rat = %.3f(%.3f) + (-gamma) log(N/Q),  gamma = %.3f +- %.3f ;  chi2/dof = %.2f" %
          (lab, sel.sum(), p[0], np.sqrt(cov[0, 0]), -p[1], np.sqrt(cov[1, 1]), chi2 / max(sel.sum() - 2, 1)))
sel = (xb > 0) & (eb > 0)
wm = 1 / eb[sel] ** 2
print("weighted mean Rat over all independent bins: %.3f +- %.3f (chi2/dof %.2f)" % ((rb[sel] * wm).sum() / wm.sum(), 1 / np.sqrt(wm.sum()),
      np.sum(wm * (rb[sel] - (rb[sel] * wm).sum() / wm.sum()) ** 2) / (sel.sum() - 1)))
small = sel & (xb * Qb < 512)
print("pooled small-nu bins (N(nu) < 512): pi Rat weighted mean %.3f +- %.3f over %d bins;  rho (zeta-free) %.3f +- %.3f" % (
    (rb[small] / eb[small] ** 2).sum() / (1 / eb[small] ** 2).sum(), 1 / np.sqrt((1 / eb[small] ** 2).sum()), small.sum(),
    np.nansum(rr[small] / re[small] ** 2) / np.nansum(1 / re[small] ** 2), 1 / np.sqrt(np.nansum(1 / re[small] ** 2))))
