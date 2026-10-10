"""
hls: least-squares theta-pole amplitude for the twisted sextic/cubic Kubota sums.

For each (m, j) and test prime t: model T(X) = R X^{1/2+1/n} + random walk (variance ~ #{c <= X}).
Blocks: dyadic X-blocks B_i; unbiased block estimates r_i = B_i / (X_i^{a} - X_{i-1}^{a}), a = 1/2 + 1/n.
Statistic per (m,j): Z = sum_t |sum_i w_i r_i|^2 / sum_t sum_i w_i^2 |r_i - rbar|^2 ...
Simpler robust test: coherence between the two halves of the X-range,
   C = Re sum_t R1_t conj R2_t / sqrt(sum_t |R1_t|^2 sum_t |R2_t|^2),
R1 = amplitude fitted on [X/64, X/8], R2 on [X/8, X].  Pole: C ~ 1.  No pole: C ~ N(0, 1/sqrt(T)).
Index convention: gamma_m(c; pi0^j) = conj((pi0/c)_6)^{j m} gamma_m(c); n = order of chi^m.
"""
import sys, json, numpy as np
from sk_read import read
D = read(sys.argv[1])
N = D["N"]; o = np.argsort(N, kind="stable")
N = N[o].astype(float); g = D["g"][o]; e = D["e"][o]; T = D["T"]
Xmax = N.max()
order = {1: 6, 2: 3, 3: 2, 4: 3, 5: 6}
def amp(cs_ph, lo, hi, a):
    """least squares R for increments on (lo, hi]: sum (B - R dX^a)"""
    i0, i1 = np.searchsorted(N, lo, side="right"), np.searchsorted(N, hi, side="right")
    return (cs_ph[i1 - 1] - cs_ph[i0 - 1]) / (hi ** a - lo ** a)
out = {}
print(f"Xmax = {Xmax:.3g}, T = {T} test primes")
print(" m  j |  C(halves)   |R2|^2 mean  null-scale   ratio")
for m in (1, 5, 2):
    n = order[m]; a = 0.5 + 1.0 / n
    for j in range(6 if n == 6 else 3):
        R1, R2 = [], []
        for t in range(T):
            ok = e[:, t] >= 0
            ph = np.where(ok, np.exp(-2j * np.pi * (j * m) * e[:, t].astype(float) / 6), 0)
            cs = np.cumsum(g[:, m - 1] * ph)
            R1.append(amp(cs, Xmax / 64, Xmax / 8, a)); R2.append(amp(cs, Xmax / 8, Xmax, a))
        R1 = np.array(R1); R2 = np.array(R2)
        C = float(np.real(np.sum(R1 * np.conj(R2))) / np.sqrt(np.sum(abs(R1) ** 2) * np.sum(abs(R2) ** 2)))
        # null scale for |R2|^2: random walk variance 0.28 (Xmax - Xmax/8) / (Xmax^a - (Xmax/8)^a)^2
        dens = len(N) / Xmax
        null = dens * (Xmax - Xmax / 8) / (Xmax ** a - (Xmax / 8) ** a) ** 2
        r2 = float(np.mean(abs(R2) ** 2))
        out[f"m{m}j{j}"] = dict(C=C, R2=r2, null=null, ratio=r2 / null)
        print(f" {m}  {j} |  {C:+.3f}      {r2:.3e}   {null:.3e}   {r2/null:6.2f}")
json.dump(out, open(sys.argv[2], "w"), indent=1)
