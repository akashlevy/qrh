import numpy as np, sys, pickle, itertools
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi); uu = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (uu * (1 - uu)) + 4); return out
K = 12; Y0 = 256; Ymax = Y0 * 2 ** K
P = primes_upto(Ymax); S = squarefree_products(P, Ymax)
H0 = 0.6961675094177613
D = np.zeros(K)
for (sa, sb, sN, idx) in S:
    if sN < Y0: continue
    i = int(np.floor(np.log2(sN / Y0)))
    if i >= K: continue
    Yi = Y0 * 2 ** i; rr = sN / Yi
    w = bump(np.array([rr]), 1, 2)[0] * rr ** -0.5 / Yi
    phi = np.prod([1 - 1 / P[t].N for t in idx]) if idx else 1.0
    D[i] += w ** 2 * phi          # w^2 q_s^{-1} * phi(s) = w^2 * phi(s)/q_s
Ys = Y0 * 2.0 ** np.arange(K)
cd = H0 / 2 * D * Ys       # Diag*Y/Q
np.save("diag_exact.npy", cd)
r = np.linspace(1, 2, 200001); What = np.trapezoid(bump(r, 1, 2) * r ** (-5 / 6), r)
for Q in [1500, 3000, 6000]:
    A = np.load("data/A_Q%d_K12.npy" % Q); ma, mb, mN = np.load("data/m_Q%d_K12.npy" % Q)
    u = Ys ** (-1 / 3) * What; V = A / u[:, None]; Om = bump(mN / Q, .5, 1)
    C = ((V * Om) @ V.conj().T).real
    T = np.array([np.mean([C[i, i] - C[i, j] - C[i, k] + C[j, k] for j, k in itertools.permutations([t for t in range(K) if t != i], 2)]) for i in range(K)])
    res = T * u ** 2 * Ys / Q
    print("Q=%d  resid(T)*Y/Q / (Diag_exact*Y/Q):" % Q, np.round(res / cd, 3), " mean %.3f" % np.mean(res / cd))
print("Diag_exact*Y/Q per scale:", np.round(cd, 4))
