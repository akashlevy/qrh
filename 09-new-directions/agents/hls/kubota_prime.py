"""
hls: per-test-prime size of the theta residue at index p^j (sextic m = 1, 5; cubic m = 2).

Amplitude estimates R_r on disjoint X-ranges r (pole model T = R X^a, a = 1/2 + 1/n); the cross product
Q = Re(R_r conj R_s) (r != s) is an unbiased estimate of |R|^2 (independent noise).  Ratio to the j = 0 residue:
rho_j(p) = Q_j(p) / Q_0, fitted as c_j N(p)^{-gamma_j} by weighted least squares on binned primes.
Theta-poorness of the hybrid family needs gamma_3 > 0 (|tau_6(p^3)|^2 decaying in N(p)).
"""
import sys, json, numpy as np
from sk_read import read
D = read(sys.argv[1])
N = D["N"]; o = np.argsort(N, kind="stable")
N = N[o].astype(float); g = D["g"][o]; e = D["e"][o]; T = D["T"]
tp = D["tp"]; pN = np.array([a * a - a * b + b * b for a, b in tp], float)
X = N.max()
ranges = [(X / 64, X / 16), (X / 16, X / 4), (X / 4, X)]
order = {1: 6, 2: 3, 5: 6}
BINS = json.loads(sys.argv[3]) if len(sys.argv) > 3 else [(0, 40), (40, 110), (110, 330)]
out = {}
for m in (1, 5, 2):
    n = order[m]; a = 0.5 + 1.0 / n
    Q = np.zeros((6 if n == 6 else 3, T))
    for j in range(Q.shape[0]):
        for t in range(T):
            ok = e[:, t] >= 0
            ph = np.where(ok, np.exp(-2j * np.pi * (j * m) * e[:, t].astype(float) / 6), 0)
            cs = np.cumsum(g[:, m - 1] * ph)
            R = []
            for lo, hi in ranges:
                i0, i1 = np.searchsorted(N, lo, side="right"), np.searchsorted(N, hi, side="right")
                R.append((cs[i1 - 1] - cs[i0 - 1]) / (hi ** a - lo ** a))
            R = np.array(R)
            prs = [np.real(R[i] * np.conj(R[k])) for i in range(len(R)) for k in range(i + 1, len(R))]
            Q[j, t] = np.mean(prs)
    q0 = np.mean(Q[0])
    res = {}
    print(f"m = {m} (order {n}):  |R_0|^2 = {q0:.3e}")
    bins = BINS
    for j in range(1, Q.shape[0]):
        rho = Q[j] / q0
        bm = [float(np.mean(rho[(pN >= lo) & (pN < hi)])) for lo, hi in bins]
        bp = [float(np.exp(np.mean(np.log(pN[(pN >= lo) & (pN < hi)])))) for lo, hi in bins]
        # fit log(rho) ~ log c - gamma log N(p) on bins with positive mean
        ok = np.array(bm) > 0
        gam = float(-np.polyfit(np.log(np.array(bp)[ok]), np.log(np.array(bm)[ok]), 1)[0]) if ok.sum() >= 2 else float("nan")
        sem = [float(np.std(rho[(pN >= lo) & (pN < hi)]) / np.sqrt(((pN >= lo) & (pN < hi)).sum())) for lo, hi in bins]
        res[f"j{j}"] = dict(bin_geo_Np=bp, rho_bin=bm, sem=sem, gamma=gam, rho_all=rho.tolist())
        print(f"   j={j}: rho by N(p)-bin " + "  ".join(f"[{p:5.0f}] {r:+.4f}±{s:.4f}" for p, r, s in zip(bp, bm, sem))
              + f"   gamma_fit = {gam:.2f}")
    out[f"m{m}"] = dict(q0=float(q0), res=res)
out["pN"] = pN.tolist()
json.dump(out, open(sys.argv[2], "w"), indent=1)
