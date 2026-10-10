"""Explicit constant K_0 of Lemma D1 for the numerical model (c(s)=mu^2(s) on primary s, (s,6)=1; rows m coprime to 6;
Omega(t)=bump(t,.5,1)), plus the 'exact-lattice' finite-scale version of the exceptional main term."""
import numpy as np, sys, pickle
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
from scipy.special import j0, gamma
from scipy.integrate import quad
PI = np.pi
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi); uu = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (uu * (1 - uu)) + 4); return out
# ---------- archimedean pieces
IOm = quad(lambda x: float(bump(np.array([x]), .5, 1)[0]) * x ** (-1 / 6), 0.5, 1)[0]
IH = (2 / np.sqrt(3)) ** (2 / 3) * PI ** (2 / 3) * gamma(1 / 6) / gamma(5 / 6) * IOm
ghat0 = (2 * PI / np.sqrt(3)) / 6 * IH
kK = PI / (3 * np.sqrt(3))
# ---------- class tables (L = 72 for unit/lambda/2 symbols, mod 4 for R)
P = primes_upto(30000)
A = np.array([Pr.a for Pr in P]); B = np.array([Pr.b for Pr in P])
xs = np.array([[1, 1], [1, -1], [2, 0], [-1, 0]])
V = np.array([Pr.symbol(xs[:, 0], xs[:, 1]) for Pr in P])
L = 72
key = (A % L) * L + (B % L)
cls, first = np.unique(key, return_index=True)          # the 432 classes of G_72
nG = len(cls); assert nG == 432
Vc = V[first]                                            # exponents for zeta6, lambda, 2, -1
k4 = (A % 4) * 4 + (B % 4)
n = 400
S = np.stack([P[j].symbol(A[:n], B[:n]) for j in range(n)], 1); Rm = (S - S.T) % 6
R4 = {}
for i in range(n):
    for j in range(n):
        if i != j: R4[(k4[i], k4[j])] = 1 - 2 * (Rm[i, j] == 3)
c4 = ((cls // L) % 4) * 4 + ((cls % L) % 4)              # class mod 4 of each G_72 class
RR = np.array([[R4[(c4[i], c4[j])] for j in range(nG)] for i in range(nG)], float)
def chi(e): return Z6[e % 6]
cz, cl, c2, cm = (chi(Vc[:, t]) for t in range(4))
# d | 6 : d in {1, lambda, 2, 2lambda}
dlist = [(1, 1, np.ones(nG)), (-1, 3, cl), (-1, 4, c2), (1, 12, cl * c2)]
def Rcol(Ckey4):  # R(x, C) for all x in G
    return np.array([R4[(c4[i], Ckey4)] for i in range(nG)], float)
# collapsed epsilon_0 sum  Etilde(C,d)
def Etilde(Ck4, dchi):
    rc = Rcol(Ck4)
    a1 = np.conj(dchi) * rc                     # x1 side: chi_{x1}(1/d) R(x1,C)
    a2 = np.conj(cm) * dchi * rc                # x2 side: conj chi_{x2}(-1) chi_{x2}(d) R(x2,C)
    same = (Vc[:, 0][:, None] % 6) == (Vc[:, 0][None, :] % 6)
    geo_l = (2 / 3) / (1 - 3 ** (-1 / 6) * cl[:, None] * np.conj(cl)[None, :])
    geo_2 = (3 / 4) / (1 - 4 ** (-1 / 6) * c2[:, None] * np.conj(c2)[None, :])
    M = same * geo_l * geo_2 * RR
    return (a1 @ M @ a2) / nG ** 2
# ---------- C sum
PC = primes_upto(200000)
Cs = squarefree_products(PC, 20000)
def euler_dens(Nmax=10 ** 6):
    pr = primes_upto(Nmax)
    qs = np.array([p.N for p in pr], float)
    f = np.prod(1 - 5 / qs ** 2 + 6 / qs ** 3 - 2 / qs ** 4)
    f0 = np.prod((1 + 2 / qs) * (1 - 1 / qs) ** 2)
    loc6 = (1 - 1 / 3) ** 2 * (1 - 1 / 4) ** 2
    return kK ** 2 * f * loc6, kK ** 2 * f0 * loc6
D1, D0 = euler_dens()
print("ghat0=%.6f  kappa_K=%.6f  Dens(1)=%.6f  Dens0(1)=%.6f" % (ghat0, kK, D1, D0))
Et_cache = {}
K0 = 0.0; terms = []
for (ca, cb, cN, idx) in Cs:
    qs = [PC[i].N for i in idx]
    phi = np.prod([q - 1 for q in qs]) if qs else 1.0
    h = np.prod([1 / (1 + 2 * (q - 1) / q ** 2) for q in qs]) if qs else 1.0
    Ck4 = (ca % 4) * 4 + (cb % 4)
    if Ck4 not in Et_cache:
        Et_cache[Ck4] = sum(mu * qd ** (-5 / 6) * Etilde(Ck4, dch) for (mu, qd, dch) in dlist)
    t = phi / cN ** 3 * h * Et_cache[Ck4]
    K0 += t; terms.append((cN, t))
K0 *= ghat0 * D1
print("K_0 = %.6e  (imag part should vanish)" % K0.real, K0.imag)
for (cN, t) in terms[:8]: print("   C norm %d  rel term %.4f" % (cN, (t * ghat0 * D1 / K0).real))
print("C=1 share: %.4f" % (terms[0][1] * ghat0 * D1 / K0).real)
pickle.dump(dict(K0=K0, ghat0=ghat0, D1=D1, D0=D0, IOm=IOm), open("k0.pkl", "wb"))
for Q in [400, 600, 1500, 3000, 6000]:
    print("Q=%5d  predicted ||rho||^2_Omega = K_0 Q^{5/6} = %.4f" % (Q, (K0 * Q ** (5 / 6)).real))
