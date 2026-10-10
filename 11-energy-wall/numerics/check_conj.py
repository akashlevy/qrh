"""Is |S_k| = |S_kbar| (conjugate rows) for the theta coefficients?  Explains z rms ~ sqrt2 for symmetric designs."""
import numpy as np
from analyze import load
d = load("data/energy_N50000.npz")
A, B = d["rowA"].astype(int), d["rowB"].astype(int)
key = {(a, b): i for i, (a, b) in enumerate(zip(A, B))}
pairs = [(i, key[(a - b, -b)]) for i, (a, b) in enumerate(zip(A, B)) if (a - b, -b) in key and key[(a - b, -b)] > i]
I = np.array([p[0] for p in pairs]); J = np.array([p[1] for p in pairs])
S = d["S"]
for n in ["beta", "rand0", "bP97", "rP97", "X:s30:p78", "X:s300:p78", "X:s300:ones", "X:s300:rnd", "XR:s300:p78", "X:t31x61:p78x2", "g0"]:
    j = d["names"].index(n)
    a, b = np.abs(S[I, j]) ** 2, np.abs(S[J, j]) ** 2
    print(f"{n:16s} corr(|S_k|^2, |S_kbar|^2) = {np.corrcoef(a, b)[0,1]:+.3f}   max rel diff {np.max(np.abs(a-b)/(a+b+1e-9)):.2e}")
