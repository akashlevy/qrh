"""Referee: locate the rho excess in the (V'_beta) control at L=3000: cube-type columns (b != 1) vs squarefree x."""
import sys, numpy as np, time
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *
from cubic_theta import cubic_theta_coeffs
Q = int(sys.argv[1]); t0 = time.time()
R = RowsE(Q, ymin_over_Q=4.0, zeta_check=10)
W = R.w; K = len(R.scales)
IX, OmI, Vp, pi2 = pi_layout(R); nI = len(OmI); rep = IX[:, 0]
Vr = R.V[:, IX]; rho2 = cross_est(np.moveaxis(Vr, 0, -1), W)
xa, xb, xN, ax, _, _ = cubic_theta_coeffs(int(np.sqrt(Q) / 2), Q)
xa = xa.astype(np.int64); xb = xb.astype(np.int64); xN = xN.astype(np.int64)
alx = np.array([cval(a, b) for a, b in zip(xa, xb)]); alx /= np.abs(alx)
beta = ax * alx / np.sqrt(xN)
clR = np.unique((R.Ra[rep] % 4) * 4 + (R.Rb[rep] % 4), return_inverse=True)[1]; ncl = clR.max() + 1
def form(Vx, x2, B0):
    TB = np.zeros((ncl, 6, K), complex)
    for k in range(K):
        v = Vx[k] * (OmI * B0)[:, None]
        for c in range(6):
            TB[:, c, k] = np.bincount(clR, weights=v[:, c].real, minlength=ncl) + 1j * np.bincount(clR, weights=v[:, c].imag, minlength=ncl)
    num = cross_est(TB, W).sum(); den = float((OmI[:, None] ** 2 * x2 * np.abs(B0)[:, None] ** 2).sum())
    jk = []
    for i in range(K):
        keep = np.array([j for j in range(K) if j != i]); x2k = cross_est(np.moveaxis(Vx[keep], 0, -1), W[keep])
        jk.append(cross_est(TB[:, :, keep], W[keep]).sum() / float((OmI[:, None] ** 2 * x2k * np.abs(B0)[:, None] ** 2).sum()))
    jk = np.array(jk); return num / den, np.sqrt((K - 1) / K * ((jk - jk.mean()) ** 2).sum())
L = 3000
wL = bump(xN / L, 0.5, 1.0); sel = np.where(wL > 0)[0]
la, lb = mul(np.int64(1), np.int64(2), xa[sel], xb[sel])
TW = quad_twist(R, rep, la, lb).astype(np.float64)
bw = beta[sel] * wL[sel]
cube = np.abs(np.abs(ax[sel]) - 1) > 1e-6
print("Q=%d L=%d: #x %d, of which b!=1: %d carrying %.3f of sum|bw|^2; their N(x), |a|: %s" % (Q, L, len(sel), cube.sum(), (np.abs(bw[cube]) ** 2).sum() / (np.abs(bw) ** 2).sum(),
      sorted(set((int(xN[sel][k]), round(float(abs(ax[sel][k])), 2)) for k in np.where(cube)[0]))))
for lab_cols, msk in (("all x", np.ones(len(sel), bool)), ("b=1 only", ~cube), ("b!=1 only", cube)):
    B0 = TW[:, msk] @ bw[msk]
    for lab, Vx, x2 in (("rho", Vr, rho2), ("pi", Vp, pi2)):
        r, se = form(Vx, x2, B0)
        print("  %-10s %s: ratio %.3f +- %.3f" % (lab_cols, lab, r, se), flush=True)
# per-column with fine classes (R mod 4 x unit) for the b != 1 columns
for k in np.where(cube)[0]:
    B0 = TW[:, k] * bw[k]
    r, se = form(Vr, rho2, B0); r2, se2 = form(Vp, pi2, B0)
    print("  column N(x)=%d |a|=%.2f: rho ratio %.2f +- %.2f ; pi %.2f +- %.2f" % (xN[sel][k], abs(ax[sel][k]), r, se, r2, se2))
print("[%.0fs]" % (time.time() - t0))
