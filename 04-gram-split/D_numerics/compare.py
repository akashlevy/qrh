import numpy as np, pickle
d = pickle.load(open("exactlat.pkl", "rb")); pred = d["pred"]
K0 = pickle.load(open("k0.pkl", "rb"))["K0"].real
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi); uu = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (uu * (1 - uu)) + 4); return out
r = np.linspace(1, 2, 200001); What = np.trapezoid(bump(r, 1, 2) * r ** (-5 / 6), r)
rng = np.random.default_rng(0)
for Q in [1500, 3000, 6000]:
    A = np.load("data/A_Q%d_K12.npy" % Q); ma, mb, mN = np.load("data/m_Q%d_K12.npy" % Q)
    K = A.shape[0]; Ys = 256.0 * 2 ** np.arange(K); u = Ys ** (-1 / 3) * What; V = A / u[:, None]
    Om = bump(mN / Q, .5, 1); nm = len(Om)
    C = ((V * Om) @ V.conj().T).real
    boots = []
    for b in range(300):
        w = np.bincount(rng.integers(0, nm, nm), minlength=nm) * Om
        boots.append(((V * w) @ V.conj().T).real)
    se = np.std(boots, axis=0)
    iu = np.triu_indices(K, 1)
    x = pred[Q][iu]; y = C[iu]; s = se[iu]
    for lo, lab in [(0, "all pairs"), (2, "pairs with i>=2"), (5, "pairs with i>=5")]:
        msk = iu[0] >= lo
        wts = 1 / s[msk] ** 2
        alpha = np.sum(wts * x[msk] * y[msk]) / np.sum(wts * x[msk] ** 2)
        chi2 = np.sum(wts * (y[msk] - alpha * x[msk]) ** 2) / (msk.sum() - 1)
        chi2c = np.sum(wts * (y[msk] - np.sum(wts * y[msk]) / np.sum(wts)) ** 2) / (msk.sum() - 1)
        print("Q=%d %-16s n=%2d  data mean %.2f  pred mean %.2f (asympt %.2f)  fit data=alpha*pred alpha=%.3f  chi2/dof %.2f  [const model chi2/dof %.2f]  corr(data,pred)=%.2f" %
              (Q, lab, msk.sum(), y[msk].mean(), x[msk].mean(), K0 * Q ** (5 / 6), alpha, chi2, chi2c, np.corrcoef(x[msk], y[msk])[0, 1]))
    print("   typical bootstrap se of kappa_ij (i>=5): %.2f" % s[iu[0] >= 5].mean())
