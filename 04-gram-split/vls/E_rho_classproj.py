"""Attack E check: is rho's small-conductor twist deficit (Rat ~ 0.75 at N(nu) < 512, Q >= 24000) its class-function
(constant-term, degree-1-like) component?  Remove per-scale class means of V_k mod F (F = 12, 36) and recompute the
pooled small-nu twist ratio for rho and pi.  usage: python E_rho_classproj.py Q ymin"""
import numpy as np, sys, time
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *
t0 = time.time()
Q = int(sys.argv[1]); ymin = float(sys.argv[2])
R = RowsE(Q, ymin_over_Q=ymin); W = R.w; K = len(R.scales)
V0 = R.V.copy()
nu_a, nu_b, nu_N = nu_columns(1024)
for F in (0, 12, 36):
    R.V = V0.copy()
    if F:
        key = (R.ma % F) * F + (R.mb % F); lab = np.unique(key, return_inverse=True)[1]
        sw = np.bincount(lab, weights=R.Om)
        for i in range(K):
            mr = np.bincount(lab, weights=R.Om * R.V[i].real) / sw; mi = np.bincount(lab, weights=R.Om * R.V[i].imag) / sw
            R.V[i] = R.V[i] - (mr + 1j * mi)[lab]
    IX, OmI, Vp, pi2 = pi_layout(R); rep = IX[:, 0]; nI = len(OmI)
    if F == 0:
        TW = quad_twist(R, rep, nu_a, nu_b).astype(np.float64)
    Vr = R.V[:, IX]; rho2 = cross_est(np.moveaxis(Vr, 0, -1), W)
    line = "F=%-3d" % F
    for lab_, X, x2 in (("rho", Vr, rho2), ("pi", Vp, pi2)):
        Y = (X * OmI[None, :, None]).transpose(1, 2, 0).reshape(nI, -1)
        T = (TW.T @ Y).reshape(len(nu_N), 6, K)
        num = cross_est(T, W).sum(1); den = ((OmI[:, None] ** 2 * x2).sum(1)) @ TW ** 2
        for lo, hi in ((7, 128), (128, 512), (512, 1025)):
            s = (nu_N >= lo) & (nu_N < hi)
            line += "  %s[%d,%d): %.3f" % (lab_, lo, hi, num[s].sum() / den[s].sum())
    print(line, flush=True)
print("[%.0fs]" % (time.time() - t0))
