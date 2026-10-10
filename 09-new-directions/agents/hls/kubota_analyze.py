"""
hls: which twisted Kubota sums have the theta main term?

T_{m,s,t}(X) = sum_{c sf odd primary, N(c) <= X, pi0_t not | c} gamma_m(c) * conj((pi0_t / c)_6)^s
  = Kubota sum of the order-(6/gcd(m,6)) Gauss sum gamma_m at the index pi0^j when s = j m (mod 6);
  for m in {2,4} (cubic) and s odd, or m = 3 and s not = 0 mod 3, it is a genuine character twist (no index).
Statistic: A(X) = mean_t |T(X)|^2 / #{c <= X}  (= 1 for square-root cancellation; ~ X^{2/n} if theta pole).
Fitted growth exponent of A over the last decades.  Expected pole growth: m=1,5: X^{1/3}; m=2,4: X^{2/3}.
"""
import sys, json, numpy as np
from sk_read import read
D = read(sys.argv[1])
N = D["N"]; o = np.argsort(N, kind="stable")
N = N[o]; g = D["g"][o]; e = D["e"][o]; T = D["T"]
Xs = [int(x) for x in np.geomspace(2e4, N.max(), 8)]
idx = [np.searchsorted(N, X, side="right") for X in Xs]
res = {}
print("X grid:", Xs)
for m in range(1, 6):
    for s in range(6):
        A = []
        for t in range(T):
            ok = e[:, t] >= 0
            ph = np.where(ok, np.exp(-2j * np.pi * s * e[:, t].astype(float) / 6), 0)
            cs = np.cumsum(g[:, m - 1] * ph)
            cnt = np.cumsum(ok)
            A.append([abs(cs[i - 1]) ** 2 / cnt[i - 1] for i in idx])
        A = np.array(A)            # T x len(Xs)
        mA = A.mean(axis=0)
        slope = np.polyfit(np.log(Xs[-4:]), np.log(mA[-4:]), 1)[0]
        res[f"m{m}s{s}"] = dict(A=mA.tolist(), slope=float(slope))
        print(f"m={m} s={s}:  A(X) = " + " ".join(f"{a:8.2f}" for a in mA) + f"   slope(last 4) = {slope:+.3f}")
json.dump(dict(X=Xs, res=res), open(sys.argv[2], "w"), indent=1)
