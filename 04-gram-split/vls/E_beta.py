"""Attack E, (V'_beta) and Type I sums with Patterson's cubic theta coefficients.

Columns (cubic_theta.py, from paper78 eq. reflection-infinity-coefficient / Dunn-Radziwill (5.7)):
  x = n b^3 primary, odd, n squarefree;  a(x) = |b| gt3(n);  mu = lambda^{-3} x,  lambda^4 mu = lambda x,
  beta_x = a(x) alpha(x) / sqrt(N x)       (constants 3^{5/2} sqrt27 alpha(lambda)^{-3} and theta(lambda^4 mu) dropped),
  psi_x(m) = (lambda x / m)_2  (exact, prod_{p | m} chi_p(lambda x)^3),  dyadic window w_L(x) = bump(N x / L, 1/2, 1).
Bilinear form (per unit class c, summed over c):
  TB_c(L) = sum_{R sf ideal} Om(R) pi(zeta6^c R) B0_R,   B0_R = sum_x beta_x w_L(x) psi_x(R)    (B0 exact, pi noisy)
  ratio_rand = sum_c |TB_c|^2 / sum_{R,c} Om^2 |pi|^2 |B0_R|^2      (random: 1)
  ratio_CS   = sum_c |TB_c|^2 / ( sum_R Om |B0_R|^2 * sum_{R,c} Om |pi|^2 )   (Cauchy-Schwarz: 1; random ~ 1/#ideals)
Type I (exact, no noise):  S_d(R) = sum_{d | x} beta_x w_L(x) psi_x(R),  d squarefree primary, N(d) <= DMAX.
usage: python E_beta.py Q src ymin
"""
import numpy as np, sys, time, json
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *
from cubic_theta import cubic_theta_coeffs

t0 = time.time()
Q = int(sys.argv[1]); src = sys.argv[2] if len(sys.argv) > 2 else "E"; ymin = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
DMAX = 400
R = RowsE(Q, src=src, ymin_over_Q=ymin)
W = R.w; K = len(R.scales)
IX, OmI, Vp, pi2 = pi_layout(R)
nI = len(OmI); rep = IX[:, 0]
Vr = R.V[:, IX]                                         # rho (zeta-free)
rho2 = cross_est(np.moveaxis(Vr, 0, -1), W)
print("Q=%d ideals %d scales %s [%.0fs]" % (Q, nI, R.Ys[R.scales].astype(int).tolist(), time.time() - t0), flush=True)
Lmax = Q
Lmin = int(np.sqrt(Q) / 2)
xa, xb, xN, ax, _, _ = cubic_theta_coeffs(Lmin, Lmax)
xa = xa.astype(np.int64); xb = xb.astype(np.int64); xN = xN.astype(np.int64)
alx = np.array([cval(a, b) for a, b in zip(xa, xb)]); alx /= np.abs(alx)
beta = ax * alx / np.sqrt(xN)
nua, nub = mul(np.int64(1), np.int64(2), xa, xb)        # lambda * x
TWb = quad_twist(R, rep, nua, nub).astype(np.float64)   # (nI, nx)
print("cubic theta columns: %d (N in (%d, %d]); |a(x)| range %.2f..%.2f [%.0fs]" % (len(xN), Lmin, Lmax, np.abs(ax).min(), np.abs(ax).max(), time.time() - t0), flush=True)
rng = np.random.default_rng(3)
exact = {"rand": np.exp(2j * np.pi * rng.random((nI, 6))), "zeta": R.zetaZ[IX]}
res = []
# classes of the primary ideal R mod 4 (12 classes; together with the 6 unit classes c this averages over 72 ray-type class functions)
clR = np.unique((R.Ra[rep] % 4) * 4 + (R.Rb[rep] % 4), return_inverse=True)[1]; ncl = clR.max() + 1
# additive column characters theta_h(x) = e(h x / 9) = exp(2 pi i (h x)_omega / 9): h = 1 plus 7 random residues mod 9
hs = [(1, 0), (0, 1), (1, 3), (2, 5), (4, 1), (7, 7), (5, 2), (3, 8)]
print("\n(V'_beta): bilinear form with cubic theta columns, dyadic windows N(x) in (L/2, L]; per (unit class c, R mod 4) and 8 additive column characters")
print("%-9s %6s %10s %13s %16s %16s %12s %12s %16s %10s %10s" % ("L", "#x", "sum|bw|^2", "<|B0|^2>/sum", "pi:rand(se_h)", "jk se(h=1)", "pi:CS", "1/#ideals", "rho:rand(se_h)", "rand(ex)", "zeta(ex)"))
for L in sorted(set([int(Q / 2 ** j) for j in range(0, 20) if Q / 2 ** j >= 2 * Lmin] + [int(Q ** 0.79)])):
    wL = bump(xN / L, 0.5, 1.0)
    sel = wL > 0
    out = {"L": L, "nx": int(sel.sum())}
    acc = {"pi": [], "rho": [], "pi_cs": [], "rand": [], "zeta": [], "jk": []}
    for (ha, hb) in hs:
        ya, yb = mul(np.int64(ha), np.int64(hb), xa[sel], xb[sel])
        th = np.exp(2j * np.pi * (yb % 9) / 9)
        bw = beta[sel] * wL[sel] * th
        B0 = TWb[:, sel] @ bw
        sb2 = float((np.abs(bw) ** 2).sum()); mB = float((OmI * np.abs(B0) ** 2).sum() / OmI.sum())
        for lab, Vx, x2 in (("pi", Vp, pi2), ("rho", Vr, rho2)):
            TB = np.zeros((ncl, 6, K), complex)
            for k in range(K):
                v = Vx[k] * (OmI * B0)[:, None]                         # (nI, 6)
                for c in range(6):
                    TB[:, c, k] = np.bincount(clR, weights=v[:, c].real, minlength=ncl) + 1j * np.bincount(clR, weights=v[:, c].imag, minlength=ncl)
            num = cross_est(TB, W).sum()
            den = float((OmI[:, None] ** 2 * x2 * np.abs(B0)[:, None] ** 2).sum())
            acc[lab].append(num / den)
            if lab == "pi":
                acc["pi_cs"].append(num / float((OmI * np.abs(B0) ** 2).sum() * (OmI[:, None] * x2).sum()))
                if (ha, hb) == (1, 0):
                    jk = []
                    for i in range(K):
                        keep = np.array([j for j in range(K) if j != i])
                        x2k = cross_est(np.moveaxis(Vx[keep], 0, -1), W[keep])
                        jk.append(cross_est(TB[:, :, keep], W[keep]).sum() / float((OmI[:, None] ** 2 * x2k * np.abs(B0)[:, None] ** 2).sum()))
                    jk = np.array(jk); acc["jk"].append(np.sqrt((K - 1) / K * ((jk - jk.mean()) ** 2).sum()))
        for k_, x in exact.items():
            v = x * (OmI * B0)[:, None]
            S_ = np.zeros((ncl, 6), complex)
            for c in range(6):
                S_[:, c] = np.bincount(clR, weights=v[:, c].real, minlength=ncl) + 1j * np.bincount(clR, weights=v[:, c].imag, minlength=ncl)
            acc[k_].append(float((np.abs(S_) ** 2).sum() / (OmI[:, None] ** 2 * np.abs(x) ** 2 * np.abs(B0)[:, None] ** 2).sum()))
    for k_ in acc: acc[k_] = np.array(acc[k_])
    out.update(sum_bw2=sb2, meanB2_over_sum=mB / sb2, pi=float(acc["pi"].mean()), pi_se_h=float(acc["pi"].std() / np.sqrt(len(hs))),
               pi_jk_h1=float(acc["jk"][0]), pi_cs=float(acc["pi_cs"].mean()), rho=float(acc["rho"].mean()), rho_se_h=float(acc["rho"].std() / np.sqrt(len(hs))),
               rand=float(acc["rand"].mean()), zeta=float(acc["zeta"].mean()), pi_all=acc["pi"].tolist())
    res.append(out)
    print("%-9d %6d %10.4f %13.3f %9.3f(%5.3f) %16.3f %12.2e %12.2e %9.3f(%5.3f) %10.3f %10.3f" % (L, sel.sum(), sb2, mB / sb2, out["pi"], out["pi_se_h"], out["pi_jk_h1"],
          out["pi_cs"], 1 / nI, out["rho"], out["rho_se_h"], out["rand"], out["zeta"]), flush=True)

# ---- Type I sums (exact): window L ~ Q^{0.79} if available else Lmax/2
L1 = int(Q ** 0.79)
wL = bump(xN / L1, 0.5, 1.0); sel = np.where(wL > 0)[0]
bw = beta[sel] * wL[sel]
Pd = primes_upto(DMAX)
ds = [(a, b, N) for (a, b, N, idx) in squarefree_products(Pd, DMAX) if N > 1]
print("\nType I sums S_d(R) = sum_{d|x} beta_x w(N x/L) (lambda x/R)_2  at L = %d (#x = %d), %d moduli d with N(d) <= %d" % (L1, len(sel), len(ds), DMAX))
print("%-14s %5s %8s %12s %12s %12s %12s" % ("d", "N(d)", "#x", "<|S_d|^2>/rnd", "max|S_d|^2/rnd", "|sum bw|^2/rnd", "(sum|bw|)^2/rnd"))
typeI = []
for (da, db, dN) in ds:
    msk = divides(np.int64(da), np.int64(db), dN, xa[sel], xb[sel])
    if msk.sum() < 2: continue
    Sd = TWb[:, sel[msk]] @ bw[msk]
    rnd = float((np.abs(bw[msk]) ** 2).sum())
    r1 = float((OmI * np.abs(Sd) ** 2).sum() / OmI.sum() / rnd)
    rmax = float(np.max(np.abs(Sd) ** 2) / rnd)
    coh = float(np.abs(bw[msk].sum()) ** 2 / rnd)
    triv = float(np.abs(bw[msk]).sum() ** 2 / rnd)
    typeI.append(dict(d=(da, db), N=dN, nx=int(msk.sum()), r1=r1, rmax=rmax, coh=coh, triv=triv))
typeI.sort(key=lambda t: t["N"])
for t in typeI[:12] + typeI[-6:]:
    print("%-14s %5d %8d %12.3f %12.2f %12.2f %12.1f" % (str(t["d"]), t["N"], t["nx"], t["r1"], t["rmax"], t["coh"], t["triv"]))
r1s = np.array([t["r1"] for t in typeI]); rmx = np.array([t["rmax"] for t in typeI]); cohs = np.array([t["coh"] for t in typeI])
nxs = np.array([t["nx"] for t in typeI])
print("over all %d d: <|S_d|^2>/rnd mean %.3f sd %.3f (min %.3f max %.3f);  max_R |S_d|^2/rnd: median %.1f (log #R = %.1f);  untwisted coherence median %.2f, max %.1f"
      % (len(typeI), r1s.mean(), r1s.std(), r1s.min(), r1s.max(), np.median(rmx), np.log(nI), np.median(cohs), cohs.max()))
# full untwisted theta sum (calibration: theta_3 main term)
print("untwisted, all x in window: |sum beta w|^2 / sum|beta w|^2 = %.2f  (#x = %d)" % (np.abs(bw.sum()) ** 2 / (np.abs(bw) ** 2).sum(), len(bw)))
json.dump({"Q": Q, "windows": res, "typeI": typeI, "L1": L1}, open(DATA + "beta_Q%d_%s.json" % (Q, src), "w"), indent=1, default=float)
print("[%.0fs]" % (time.time() - t0))
