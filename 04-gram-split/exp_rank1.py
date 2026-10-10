"""
N1/N2: is the P_a^{1/6} excess of the additive polynomial a rank-one (theta) component?

A_m(Y) = Y^{-1} sum_s W(N s/Y) (N s/Y)^{-1/2} gtilde(s,-m),  gtilde(s,h) = g(s,h)/sqrt(N s)
s: primary squarefree, (s,6)=1, N s in [Y, 2Y] (disjoint dyadic supports for Y_i = Y0 2^i).
Theta picture: A_m(Y) = Y^{-1/3} What(1/6) rho(m) + O(Y^{-1/2}),  What(1/6) = int W(r) r^{-1/2} r^{2/3} dr/r.
Prediction: for i != j (no common s, no diagonal):  G_ij (Y_i Y_j)^{1/3} -> kappa = What^2 sum_m Omega |rho|^2.
"""
import numpy as np, sys, time, json
from zomega import *

Q = float(sys.argv[1]) if len(sys.argv) > 1 else 1500.0
Y0 = int(sys.argv[2]) if len(sys.argv) > 2 else 256
K = int(sys.argv[3]) if len(sys.argv) > 3 else 8
Ymax = Y0 * 2 ** K
t0 = time.time()

def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t)
    inside = (t > lo) & (t < hi)
    u = (t[inside] - lo) / (hi - lo)
    out[inside] = np.exp(-1.0 / (u * (1 - u)) + 4.0)
    return out

ma, mb, mN = elements_in_annulus(Q / 2, Q)
Om = bump(mN / Q, 0.5, 1.0)
keep = Om > 1e-12
ma, mb, mN, Om = ma[keep], mb[keep], mN[keep], Om[keep]
nm = len(ma)
print("rows m:", nm, flush=True)

P = primes_upto(Ymax)
print("primes:", len(P), "t=%.1f" % (time.time() - t0), flush=True)
gt = np.array([Pr.gauss() / np.sqrt(Pr.N) for Pr in P])
print("prime gauss sums done t=%.1f" % (time.time() - t0), flush=True)
J = np.stack([Pr.symbol(-ma, -mb) for Pr in P])            # (nP, nm), exponent of (-m/pi)_6
print("symbol table done t=%.1f" % (time.time() - t0), flush=True)

S = squarefree_products(P, Ymax)
print("squarefree s:", len(S), flush=True)

Ys = [Y0 * 2 ** i for i in range(K)]
A = np.zeros((K, nm), complex)
r = np.linspace(1, 2, 20001)
What = np.trapezoid(bump(r, 1, 2) * r ** (-0.5 + 2 / 3) / r, r)
for (sa, sb, sN, idx) in S:
    if sN < Y0: continue
    i = int(np.floor(np.log2(sN / Y0)))
    if i >= K: continue
    # weight in each dyadic (supports [Y_i, 2 Y_i])
    w = 0.0
    Yi = Ys[i]
    rr = sN / Yi
    w = bump(np.array([rr]), 1, 2)[0] * rr ** -0.5 / Yi
    if w == 0: continue
    # gtilde(s) by twisted multiplicativity
    g = 1.0 + 0j
    for t, ii in enumerate(idx):
        g *= gt[ii]
        for jj in idx[t + 1:]:
            e1 = P[ii].symbol(np.array([P[jj].a]), np.array([P[jj].b]))[0]
            e2 = P[jj].symbol(np.array([P[ii].a]), np.array([P[ii].b]))[0]
            g *= Z6[(e1 + e2) % 6]
    rows = J[list(idx)] if idx else np.zeros((0, nm), np.int64)
    if len(idx):
        bad = (rows < 0).any(axis=0)
        e = rows.sum(axis=0) % 6
        val = np.conj(Z6[e])                                   # conj chi_s(-m)
        val[bad] = 0
    else:
        val = np.ones(nm)
    A[i] += w * g * val
print("A built t=%.1f" % (time.time() - t0), flush=True)

G = (A * Om) @ A.conj().T
np.save("A_Q%d.npy" % Q, A); np.save("m_Q%d.npy" % Q, np.stack([ma, mb, mN]))
Yv = np.array(Ys, float) * np.sqrt(2)   # geometric centre of [Y,2Y] for reporting
print("What(1/6) =", What)
print("diag G_ii * Y_i  (expect ~const + kappa Y^{1/3}):")
for i in range(K):
    print("  Y=%7d  G_ii*Y=%.4f   G_ii*Y^{2/3}=%.4f" % (Ys[i], G[i, i].real * Ys[i], G[i, i].real * Ys[i] ** (2 / 3)))
print("off-diagonal G_ij (Y_i Y_j)^{1/3}  (prediction: constant kappa, real>0):")
for i in range(K):
    print("  ", " ".join("%8.4f%+.3fi" % ((G[i, j] * (Ys[i] * Ys[j]) ** (1 / 3)).real, (G[i, j] * (Ys[i] * Ys[j]) ** (1 / 3)).imag) for j in range(K)))
# normalized correlation between Y^{1/3}A(Y) vectors
V = A * np.array(Ys)[:, None] ** (1 / 3)
C = (V * Om) @ V.conj().T
d = np.sqrt(np.real(np.diag(C)))
print("correlation of Y^{1/3}A(Y) across scales:")
for i in range(K):
    print("  ", " ".join("%6.3f" % (abs(C[i, j]) / d[i] / d[j]) for j in range(K)))
json.dump(dict(Q=Q, Ys=Ys, G_re=G.real.tolist(), G_im=G.imag.tolist(), What=What), open("rank1_Q%d.json" % Q, "w"))
print("done t=%.1f" % (time.time() - t0))
