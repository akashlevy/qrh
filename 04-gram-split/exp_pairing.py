"""
Theta pairing test: T = sum_m Om rho(m) B_m with the (unreflected) completed-row model
  B_m^{(q)} = sum_{c sf primary, N c in [Z,2Z]} W(Nc/Z) gamma_2(c) conj(alpha(c)) chi_c(q m) / sqrt(N c)
for several fixed twists q.  rho enters only through cross-scale products of v_i = Y_i^{1/3} A(Y_i).
Ratio R = E|T|^2 / sum_m Om^2 |rho|^2 |B_m|^2 :  ~1 random,  >>1 coherent,  <<1 structurally cancelling.
"""
import numpy as np, sys, time
from zomega import *
Q = 3000; i0 = 4
Z = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
A = np.load("A_Q%d.npy" % Q); ma, mb, mN = np.load("m_Q%d.npy" % Q)
K = A.shape[0]; Ys = np.array([256 * 2 ** i for i in range(K)], float)
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (u * (1 - u)) + 4); return out
Om = bump(mN / Q, .5, 1)
V = A[i0:] * Ys[i0:, None] ** (1 / 3); nV = len(V)
t0 = time.time()
P = primes_upto(2 * Z)
# gamma_2 for primes
def gam2(Pr):
    if Pr.kind == "split":
        p = Pr.q; x = np.arange(1, p); jj = Pr.symbol(x, np.zeros_like(x))
        return np.sum(Z6[(2 * jj) % 6] * np.exp(-2j * np.pi * ((Pr.b * x) % p) / p)) / np.sqrt(Pr.N)
    q = Pr.q; c, d = np.meshgrid(np.arange(q), np.arange(q), indexing="ij"); c = c.ravel(); d = d.ravel()
    nz = (c != 0) | (d != 0); c = c[nz]; d = d[nz]; jj = Pr.symbol(c, d); sgn = -1 if Pr.a < 0 else 1
    return np.sum(Z6[(2 * jj) % 6] * np.exp(-2j * np.pi * sgn * d / q)) / np.sqrt(Pr.N)
g2 = np.array([gam2(Pr) for Pr in P])
print("gamma2 done %.1fs" % (time.time() - t0), flush=True)
twists = [(1, 0), (7, 3), (4, 3), (-2, 3), (13, 0), (10, 3), (-5, -3), (1, 3), (7, 0), (-11, -3)]
Jq = {}
for (qa, qb) in twists:
    xa, xb = mul(ma, mb, qa, qb)
    Jq[(qa, qb)] = np.stack([Pr.symbol(xa, xb) for Pr in P])
print("symbols done %.1fs" % (time.time() - t0), flush=True)
S = squarefree_products(P, 2 * Z)
Bs = {k: np.zeros(len(ma), complex) for k in twists}
for (sa, sb, sN, idx) in S:
    if not (Z <= sN <= 2 * Z): continue
    w = bump(np.array([sN / Z]), 1, 2)[0] / np.sqrt(sN)
    g = 1.0 + 0j
    for t, ii in enumerate(idx):
        g *= g2[ii]
        for jj in idx[t + 1:]:
            e1 = P[ii].symbol(np.array([P[jj].a]), np.array([P[jj].b]))[0]
            e2 = P[jj].symbol(np.array([P[ii].a]), np.array([P[ii].b]))[0]
            g *= Z6[(2 * (e1 + e2)) % 6]
    al = cval(sa, sb); al = np.conj(al / abs(al))
    for k in twists:
        rows = Jq[k][list(idx)]
        bad = (rows < 0).any(axis=0); e = rows.sum(axis=0) % 6
        val = Z6[e]; val[bad] = 0
        Bs[k] += w * g * al * val
print("B built %.1fs" % (time.time() - t0), flush=True)
pairs = [(i, j) for i in range(nV) for j in range(nV) if i != j]
rows_out = []
for k in twists:
    B = Bs[k]
    T = (V * Om * B).sum(axis=1)
    est = np.mean([(T[i] * np.conj(T[j])).real for i, j in pairs])
    rnd = np.mean([np.sum(Om ** 2 * V[i] * np.conj(V[j]) * np.abs(B) ** 2).real for i, j in pairs])
    cs = np.sqrt(np.mean([np.sum(Om ** 2 * V[i] * np.conj(V[j])).real for i, j in pairs]) * np.sum(Om ** 2 * np.abs(B) ** 2))
    rows_out.append((est / rnd, est / cs ** 2))
    print("twist %-10s  |T|^2/random = %7.3f   |T|^2/CS^2 = %.5f   (#m=%d)" % (k, est / rnd, est / cs ** 2, len(ma)))
r = np.array(rows_out)
print("mean |T|^2/random over twists = %.3f (se %.3f);  1/#m = %.5f" % (r[:, 0].mean(), r[:, 0].std() / np.sqrt(len(r)), 1 / len(ma)))
