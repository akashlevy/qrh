"""Finite-scale ('exact-lattice') exceptional main term and the k=0 diagonal constant, for comparison with the Gram data."""
import numpy as np, sys, pickle
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
from scipy.special import j0
from scipy.interpolate import interp1d
exec(open("k0.py").read().split("# ---------- C sum")[0])   # reuse class tables & pieces
kd = pickle.load(open("k0.pkl", "rb")); D0 = kd["D0"]; K0 = kd["K0"].real
# H on fine grid and lattice function Lambda(T) = (1/6) sum_{a != 0} H(N(a)^6/T)
rho = np.linspace(np.sqrt(0.5), 1, 20001); Omr = bump(rho ** 2, 0.5, 1.0) * rho
def H(t):
    k = 2 * PI * (2 / np.sqrt(3)) * np.sqrt(t)
    return (2 / np.sqrt(3)) * 2 * PI * np.trapezoid(Omr * j0(k * rho), rho)
tg = np.concatenate([np.linspace(0, 50, 5001), np.linspace(50, 3000, 6001)[1:]])
Hg = np.array([H(t) for t in tg]); Hf = interp1d(tg, Hg, bounds_error=False, fill_value=0.0)
Nmax = 200000
R_ = 600
aa, bb = np.meshgrid(np.arange(-R_, R_ + 1), np.arange(-R_, R_ + 1), indexing="ij")
NN = norm(aa, bb).ravel(); NN = NN[(NN > 0) & (NN <= Nmax)]
Nu, mult = np.unique(NN, return_counts=True)
Tgrid = np.exp(np.linspace(np.log(1e-4), np.log(1e14), 1500))
Lam = np.array([np.sum(mult * Hf(Nu.astype(float) ** 6 / T)) / 6 for T in Tgrid])
Lf = lambda T: np.interp(np.log(np.maximum(T, 1e-30)), np.log(Tgrid), Lam, left=0.0, right=np.nan)
print("check Lambda vs main at T=1e12:", Lf(1e12), kd["ghat0"] / 6 * 1e12 ** (1 / 6))
# primes of small norm for b-Mobius and for C
PS = primes_upto(5000); 
def sqfree(Nm): return squarefree_products(PS, Nm)
Bs = sqfree(400)    # b's: q_b^6 <= T/1e-3
def Lex(T, Cprimes):
    out = np.zeros_like(T)
    for (ba, bb_, bN, idx) in Bs:
        if any(i in Cprimes for i in idx): continue
        w = (-1) ** len(idx) * np.prod([2 / (PS[i].N + 2) for i in idx]) if idx else 1.0
        out += w * Lf(T / float(bN) ** 6)
    return out
# E(C, eps0, d) summed over units, as a function of (i,j) for eps0 = zeta^u lambda^i 2^j
def Eij(Ck4, dchi):
    rc = Rcol(Ck4)
    out = np.zeros((6, 6), complex)
    for uu in range(6):
        for i in range(6):
            for j in range(6):
                ce = cz ** uu * cl ** i * c2 ** j
                a1 = ce * np.conj(dchi) * rc
                a2 = np.conj(cm) * np.conj(ce) * dchi * rc
                out[i, j] += (a1 @ RR @ a2) / nG ** 2
    return out
r = np.linspace(1, 2, 241); W = bump(r, 1, 2)
r1, r2 = np.meshgrid(r, r, indexing="ij"); WW = (W[:, None] * W[None, :]) / (r1 * r2); dr = r[1] - r[0]
Cs = [c for c in squarefree_products(PS, 400)]
Ecache = {}
# precompute Phi(X) = int int W W/(r1 r2) Lambda_ex(X r1 r2) (C=1 version of Lambda_ex for all C)
Xg = np.exp(np.linspace(np.log(1e-4), np.log(1e13), 900))
rr = (r1 * r2).ravel(); wwf = WW.ravel() * dr * dr
Phi = np.array([np.sum(wwf * Lex(X * rr, set())) for X in Xg])
Phif = lambda X: np.interp(np.log(X), np.log(Xg), Phi, left=0.0, right=np.nan)
print("Phi done", flush=True)
def Gexc(Q, Y1, Y2):
    tot = 0.0
    for (ca, cb, cN, idx) in Cs:
        phi = np.prod([PS[i].N - 1 for i in idx]) if idx else 1.0
        Dc = D0 * np.prod([1 / (1 + 2 / PS[i].N) for i in idx]) if idx else D0
        Ck4 = (ca % 4) * 4 + (cb % 4)
        for (mu, qd, dch) in dlist:
            key = (Ck4, qd)
            if key not in Ecache: Ecache[key] = Eij(Ck4, dch)
            E = Ecache[key]
            for i in range(6):
                for j in range(6):
                    qe = 3.0 ** i * 4.0 ** j
                    X = Y1 * Y2 * qd / (Q * float(cN) ** 6 * qe)
                    if X < 1e-4: continue
                    val = Phif(X)
                    tot += mu / qd * phi / float(cN) ** 2 * Dc * E[i, j] * val
    return Q / np.sqrt(Y1 * Y2) * tot
r_ = np.linspace(1, 2, 200001); What = np.trapezoid(bump(r_, 1, 2) * r_ ** (-5 / 6), r_)
# k=0 diagonal constant: Diag*Y/Q = H(0) * (1/2) * dens_s * int W^2 / r dr
pr = primes_upto(10 ** 6); qs = np.array([p.N for p in pr], float)
dens_s = kK * (2 / 3) * (3 / 4) * np.prod(1 - 2 / qs ** 2 + 1 / qs ** 3)
cdiag = H(0.0) * 0.5 * dens_s * np.trapezoid(bump(r_, 1, 2) ** 2 / r_, r_)
print("diagonal (k=0) constant: Diag = c Q/Y with c = %.5f" % cdiag)
out = {}
for Q in [1500, 3000, 6000]:
    Ys = 256.0 * 2 ** np.arange(12)
    M = np.zeros((12, 12))
    for a in range(12):
        for b in range(a + 1, 12):
            g = Gexc(Q, Ys[a], Ys[b]).real
            M[a, b] = M[b, a] = g / (What ** 2 * (Ys[a] * Ys[b]) ** (-1 / 3))
    out[Q] = M
    np.set_printoptions(precision=2, suppress=True, linewidth=200)
    print("Q=%d  exact-lattice exceptional prediction for kappa_ij (asymptotic value K0 Q^{5/6} = %.3f):" % (Q, K0 * Q ** (5 / 6)))
    print(M)
pickle.dump(dict(pred=out, cdiag=cdiag), open("exactlat.pkl", "wb"))
