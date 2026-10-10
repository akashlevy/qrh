"""
N3: coherence of twisted theta sums  F_j(mu, c) = sum_{m in class c} Omega(m) rho(m) (m/mu)_6^j.
rho is accessed only through cross-scale products (unbiased): v_i = Y_i^{1/3} A(Y_i),
T_i = sum_m Omega v_i(m) psi(m);  E[T_i conj T_j] = |What F|^2 for i != j (independent noise).
Coherence ratio = mean_{i<j} Re T_i conj T_j / (sum_m Omega^2 |What rho|^2 restricted to the class).
Random-sign twist => ratio ~ 1 on average over mu; coherent => ratio >> 1 (up to #m).
"""
import numpy as np, sys
from zomega import *
Q = int(sys.argv[1]); Lmu = int(sys.argv[2]); i0 = int(sys.argv[3]) if len(sys.argv) > 3 else 4
A = np.load("A_Q%d.npy" % Q); ma, mb, mN = np.load("m_Q%d.npy" % Q)
K = A.shape[0]; Y0 = 256; Ys = np.array([Y0 * 2 ** i for i in range(K)], float)
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); inside = (t > lo) & (t < hi)
    u = (t[inside] - lo) / (hi - lo); out[inside] = np.exp(-1.0 / (u * (1 - u)) + 4.0); return out
Om = bump(mN / Q, 0.5, 1.0)
V = A[i0:] * Ys[i0:, None] ** (1 / 3)
pairs = [(i, j) for i in range(len(V)) for j in range(i + 1, len(V))]
# class of m mod 12: use (a mod 12, b mod 12) up to nothing -- we use unit class mod 3 and residue mod 2 (F4) separately
def cls3(a, b):
    out = np.zeros(len(a), int)
    units = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]
    for k, (ua, ub) in enumerate(units):
        x, y = mul(a, b, ua, ub)
        out[(x % 3 == 1) & (y % 3 == 0)] = k        # m * unit_k is primary
    return out
c3 = cls3(ma, mb)
c2 = (ma % 2) * 2 + (mb % 2)                        # residue in F4 = O/2
c4 = ((ma % 4) * 4 + (mb % 4))                      # residue mod 4
classes = {"all": np.ones(len(ma), bool)}
for k in range(6): classes["u%d" % k] = c3 == k
def base_norm(mask):
    return np.mean([np.sum(Om[mask] ** 2 * V[i][mask] * np.conj(V[j][mask])).real for i, j in pairs])
P = [Pr for Pr in primes_upto(Lmu) if Pr.N > 50]
print("mu primes:", len(P), " rows:", len(ma))
for cname, mask in classes.items():
    bn = base_norm(mask)
    line = []
    for j in range(6):
        rat = []
        for Pr in P:
            e = Pr.symbol(ma[mask], mb[mask]); ok = e >= 0
            psi = np.where(ok, Z6[(j * e) % 6], 0)
            T = [np.sum(Om[mask] * V[i][mask] * psi) for i in range(len(V))]
            rat.append(np.mean([(T[i] * np.conj(T[k])).real for i, k in pairs]) / bn)
        rat = np.array(rat)
        line.append("j=%d: %7.2f±%5.2f" % (j, rat.mean(), rat.std() / np.sqrt(len(rat))))
    print("%4s |rho|^2-norm=%.3e  " % (cname, bn), "  ".join(line))
