"""Numerical verdict for lemma (D): is ||A(Y) - Y^{-1/3} What rho||^2 = c Q/Y with no Y^{1/3} excess?
Two estimators of the residual, both unbiased if the scale-'rests' w_i = V_i - rho are mutually orthogonal (V_i = Y_i^{1/3} A_i / What):
 (T) triple estimator  <V_i - V_j, V_i - V_k>, averaged over ordered pairs j != k, both != i      (no rho needed)
 (R) rho_hat = mean of V_j over the top-4 scales J;  ||V_i - rho_hat||^2 - (1/|J|^2) sum_{j in J} ||w_j||^2_est   (i not in J)
"""
import numpy as np, itertools, pickle, sys
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi); uu = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (uu * (1 - uu)) + 4); return out
r = np.linspace(1, 2, 200001); What = np.trapezoid(bump(r, 1, 2) * r ** (-5 / 6), r)
K0 = pickle.load(open("k0.pkl", "rb"))["K0"].real; cdiag = pickle.load(open("exactlat.pkl", "rb"))["cdiag"]
rng = np.random.default_rng(7)
def estimators(C, K, top):
    T = np.array([np.mean([C[i, i] - C[i, j] - C[i, k] + C[j, k] for j, k in itertools.permutations([t for t in range(K) if t != i], 2)]) for i in range(K)])
    J = list(top)
    wJ = np.array([np.mean([C[j, j] - C[j, a] - C[j, b] + C[a, b] for a, b in itertools.permutations([t for t in J if t != j], 2)]) for j in J])
    bias = wJ.sum() / len(J) ** 2
    Rr = np.full(K, np.nan)
    for i in range(K):
        if i in J: continue
        # ||V_i - mean_J V_j||^2 = C_ii - 2 mean_j C_ij + mean_{j,k} C_jk
        Rr[i] = C[i, i] - 2 * np.mean([C[i, j] for j in J]) + np.mean([[C[j, k] for k in J] for j in J]) - bias
    return T, Rr
for Q in [1500, 3000, 6000]:
    A = np.load("data/A_Q%d_K12.npy" % Q); ma, mb, mN = np.load("data/m_Q%d_K12.npy" % Q)
    K = A.shape[0]; Ys = 256.0 * 2 ** np.arange(K); u = Ys ** (-1 / 3) * What; V = A / u[:, None]
    Om = bump(mN / Q, .5, 1); nm = len(Om)
    C = ((V * Om) @ V.conj().T).real
    top = range(K - 4, K)
    T, Rr = estimators(C, K, top)
    bt = []
    for b in range(300):
        w = np.bincount(rng.integers(0, nm, nm), minlength=nm) * Om
        bt.append(estimators(((V * w) @ V.conj().T).real, K, top))
    bt = np.array(bt)
    seT, seR = bt[:, 0].std(0), bt[:, 1].std(0)
    f = u ** 2 * Ys / Q
    rho2 = K0 * Q ** (5 / 6)
    print("\nQ=%d   rows=%d   K0*Q^{5/6}=%.3f   k=0 diagonal constant c_diag=%.4f" % (Q, nm, rho2, cdiag))
    print("    Y     P^1/6   G_ii*Y/Q   resid(T)*Y/Q       resid(R)*Y/Q       [full excess would add u^2 K0 Q^{5/6} Y/Q]")
    for i in range(K):
        sR = "   ---           " if np.isnan(Rr[i]) else "%.4f (%.4f)" % (Rr[i] * f[i], seR[i] * f[i])
        print("%8d  %5.2f    %.4f    %.4f (%.4f)    %s     %.4f" % (Ys[i], (Ys[i] ** 2 * 2.25 / Q) ** (1 / 6), C[i, i] * f[i], T[i] * f[i], seT[i] * f[i], sR, rho2 * f[i]))
    for lab, est, se in [("T", T, seT), ("R", Rr, seR)]:
        ok = ~np.isnan(est)
        y = est[ok] * f[ok]; s = se[ok] * f[ok]; x = What ** 2 * Ys[ok] ** (1 / 3) / Q
        Wt = 1 / s ** 2; X = np.stack([np.ones(ok.sum()), x], 1)
        cov = np.linalg.inv(X.T @ (Wt[:, None] * X)); beta = cov @ X.T @ (Wt * y)
        print("  fit resid(%s)*Y/Q = a + b What^2 Y^{1/3}/Q:  a = %.5f +- %.5f   b = %.3f +- %.3f   (b/(K0 Q^{5/6}) = %.3f +- %.3f)" %
              (lab, beta[0], np.sqrt(cov[0, 0]), beta[1], np.sqrt(cov[1, 1]), beta[1] / rho2, np.sqrt(cov[1, 1]) / rho2))
    off = C[np.triu_indices(K, 1)]
    big = [C[i, j] for i in range(K) for j in range(i + 1, K) if i >= K - 4]
    print("  off-diagonal <V_i,V_j>: mean all %.3f, mean top-4 block %.3f   (K0 Q^{5/6} = %.3f)" % (off.mean(), np.mean(big), rho2))
