"""Attack E bias diagnostic: the large-sieve sum S(L) = sum_{N nu <= L} sum_c |F'_c(nu)|^2 estimated from each single pair of
scales (k, l), k != l:  S_kl(L) = sum Re T_k conj(T_l).  If the noise of different scales is uncorrelated, every pair is an
unbiased estimate of the same S(L);  a deterministic small-Y component shows up as pairs involving the smallest Y deviating.
Also the untwisted per-scale class sums sum_m Om pi_k(m) (per unit class) as in the A2 referee report.
usage: python E_pairs.py Q ymin numax_over_Q
"""
import numpy as np, sys, time, json
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *

t0 = time.time()
Q = int(sys.argv[1]); ymin = float(sys.argv[2]); numax = int(float(sys.argv[3]) * Q)
R = RowsE(Q, ymin_over_Q=ymin)
W = R.w; K = len(R.scales); Ys = R.Ys[R.scales].astype(int)
IX, OmI, Vp, pi2 = pi_layout(R); nI = len(OmI); rep = IX[:, 0]
Vr = R.V[:, IX]
nu_a, nu_b, nu_N = nu_columns(numax)
TW = quad_twist(R, rep, nu_a, nu_b)
TW2f = np.abs(TW).astype(np.float32)
print("Q=%d scales %s, ideals %d, columns %d [%.0fs]" % (Q, Ys.tolist(), nI, len(nu_N), time.time() - t0), flush=True)
out = {}
for lab, X in (("pi", Vp), ("rho", Vr)):
    x2 = cross_est(np.moveaxis(X, 0, -1), W)
    den = ((OmI[:, None] ** 2 * x2).sum(1)) @ TW2f
    Y = (X * OmI[None, :, None]).transpose(1, 2, 0).reshape(nI, -1)
    T = np.zeros((len(nu_N), Y.shape[1]), complex)
    for c0 in range(0, len(nu_N), 4000):
        T[c0:c0 + 4000] = TW[:, c0:c0 + 4000].T.astype(np.float64) @ Y
    T = T.reshape(len(nu_N), 6, K)
    print("\n%s: pairwise-scale estimates of R_rand(L) = S(L) / sum den" % lab)
    for ex in (0.5, 0.6, 0.7, 0.79, 1.0):
        sel = nu_N <= Q ** ex
        D = den[sel].sum()
        M = np.array([[np.sum((T[sel, :, k] * T[sel, :, l].conj()).real) / D if k != l else np.nan for l in range(K)] for k in range(K)])
        out["%s_%.2f" % (lab, ex)] = M.tolist()
        print("  L = Q^%.2f (#nu %d):" % (ex, sel.sum()))
        for k in range(K):
            print("     Y=%8d: " % Ys[k] + " ".join("%7.3f" % M[k, l] if k != l else "    -- " for l in range(K)))
    # untwisted per-scale class sums (unit classes), normalised by sqrt of the random level
    cs = (X * OmI[None, :, None]).sum(1)                       # (K, 6)
    rl = np.sqrt((OmI[:, None] ** 2 * x2).sum(0))             # (6,)
    print("  untwisted per-scale unit-class sums / sqrt(random level):")
    for k in range(K):
        print("     Y=%8d: " % Ys[k] + " ".join("%7.2f%+7.2fi" % (z.real, z.imag) for z in cs[k] / rl))
json.dump(out, open(DATA + "pairs_Q%d_y%g.json" % (Q, ymin), "w"))
print("[%.0fs]" % (time.time() - t0))
