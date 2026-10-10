"""
hls: theta residue at composite cube indices (the hybrid rows k = p1 p2):
compare Q(p1^3 p2^3) with Q(p1^3) Q(p2^3) / Q(1)  (multiplicativity of |tau_6| at cubes), sextic m = 1, 5,
and also the index p1^3 p2^2 (row k = p1 with mark P = p2, m = 1) against Q(p1^3) Q(p2^2) / Q(1).
Q = unbiased cross-range estimate of |R|^2 as in kubota_prime.py.
"""
import sys, json, itertools, numpy as np
from sk_read import read
D = read(sys.argv[1])
N = D["N"]; o = np.argsort(N, kind="stable")
N = N[o].astype(float); g = D["g"][o]; e = D["e"][o].astype(np.int64); T = D["T"]
tp = D["tp"]; pN = np.array([a * a - a * b + b * b for a, b in tp], float)
X = N.max(); ranges = [(X / 64, X / 16), (X / 16, X / 4), (X / 4, X)]
a_ = 0.5 + 1 / 6
idx = [(np.searchsorted(N, lo, side="right"), np.searchsorted(N, hi, side="right")) for lo, hi in ranges]
def Q(m, expo):
    """expo: integer array (len N) of total sextic exponent of conj((h/c)_6)^m phase, -1 = excluded"""
    ok = expo >= 0
    ph = np.where(ok, np.exp(-2j * np.pi * m * expo / 6), 0)
    cs = np.cumsum(g[:, m - 1] * ph)
    R = np.array([(cs[i1 - 1] - cs[i0 - 1]) / (hi ** a_ - lo ** a_) for (i0, i1), (lo, hi) in zip(idx, ranges)])
    return float(np.mean([np.real(R[i] * np.conj(R[k])) for i in range(3) for k in range(i + 1, 3)]))
out = {}
for m in (1, 5):
    q1 = np.mean([Q(m, np.where(e[:, t] >= 0, 0, -1)) for t in range(T)])
    single3 = np.array([Q(m, np.where(e[:, t] >= 0, 3 * e[:, t] % 6, -1)) for t in range(T)])
    single2 = np.array([Q(m, np.where(e[:, t] >= 0, 2 * e[:, t] % 6, -1)) for t in range(T)])
    pq, pr, pq2, pr2 = [], [], [], []
    for t1, t2 in itertools.combinations(range(T), 2):
        ok = (e[:, t1] >= 0) & (e[:, t2] >= 0)
        ex = np.where(ok, (3 * e[:, t1] + 3 * e[:, t2]) % 6, -1)
        pq.append(Q(m, ex)); pr.append(single3[t1] * single3[t2] / q1)
        ex2 = np.where(ok, (3 * e[:, t1] + 2 * e[:, t2]) % 6, -1)
        pq2.append(Q(m, ex2)); pr2.append(single3[t1] * single2[t2] / q1)
    pq, pr, pq2, pr2 = map(np.array, (pq, pr, pq2, pr2))
    out[f"m{m}"] = dict(q1=float(q1), mean_Q_p3q3=float(pq.mean()), mean_pred=float(pr.mean()),
                        sem_Q_p3q3=float(pq.std() / np.sqrt(len(pq))), mean_Q_p3q2=float(pq2.mean()),
                        mean_pred2=float(pr2.mean()), sem_Q_p3q2=float(pq2.std() / np.sqrt(len(pq2))),
                        mean_single3=float(single3.mean()), mean_single2=float(single2.mean()))
    print(f"m={m}: Q(1)={q1:.3e}  mean Q(p^3)={single3.mean():.3e}  mean Q(p^2)={single2.mean():.3e}")
    print(f"   index p1^3 p2^3: mean Q = {pq.mean():.3e} +- {pq.std()/np.sqrt(len(pq)):.1e}   multiplicative prediction {pr.mean():.3e}")
    print(f"   index p1^3 p2^2: mean Q = {pq2.mean():.3e} +- {pq2.std()/np.sqrt(len(pq2)):.1e}   multiplicative prediction {pr2.mean():.3e}")
json.dump(out, open(sys.argv[2], "w"), indent=1)
