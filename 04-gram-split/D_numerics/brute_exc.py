"""Brute-force exceptional-frequency sum (C=1 part and C=prime parts) for the actual s-lists, vs data and vs density formula."""
import numpy as np, sys, pickle, math
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
from scipy.special import j0
from scipy.interpolate import interp1d
exec(open("k0.py").read().split("# ---------- C sum")[0])
Q = float(sys.argv[1]); pairs = [tuple(map(int, p.split(","))) for p in sys.argv[2:]]
Kmax = max(max(p) for p in pairs) + 1
Y0 = 256; Ymax = Y0 * 2 ** Kmax
rho_ = np.linspace(np.sqrt(0.5), 1, 20001); Omr = bump(rho_ ** 2, 0.5, 1.0) * rho_
def Hfun(t):
    k = 2 * PI * (2 / np.sqrt(3)) * np.sqrt(t)
    return (2 / np.sqrt(3)) * 2 * PI * np.trapezoid(Omr * j0(k * rho_), rho_)
tg = np.concatenate([np.linspace(0, 50, 5001), np.linspace(50, 1000, 4001)[1:]])
Hf = interp1d(tg, np.array([Hfun(t) for t in tg]), bounds_error=False, fill_value=0.0)
PP = primes_upto(Ymax)
S = squarefree_products(PP, Ymax)
cls_index = {c: i for i, c in enumerate(cls)}
sc = {}   # per scale: arrays
for (sa, sb, sN, idx) in S:
    if sN < Y0: continue
    i = int(np.floor(np.log2(sN / Y0)))
    if i >= Kmax: continue
    rr = sN / (Y0 * 2 ** i)
    w = bump(np.array([rr]), 1, 2)[0] * rr ** -0.5 / (Y0 * 2 ** i)
    if w == 0: continue
    sc.setdefault(i, []).append((sa, sb, sN, w, idx))
for i in sc:
    L_ = sc[i]
    sc[i] = dict(a=np.array([x[0] for x in L_]), b=np.array([x[1] for x in L_]), N=np.array([x[2] for x in L_], float),
                 w=np.array([x[3] for x in L_]), idx=[set(x[4]) for x in L_])
    k72 = (sc[i]["a"] % 72) * 72 + (sc[i]["b"] % 72)
    ci = np.array([cls_index[k] for k in k72])
    sc[i]["ci"] = ci
# small a's: elements of O (all), with norms; prime factor sets via small primes
smallP = [p for p in PP if p.N <= 400]
def afactors(a, b):
    out = set()
    for t, p in enumerate(PP[:len(smallP)]):
        # p | a  iff symbol == -1
        if p.symbol(np.array([a]), np.array([b]))[0] == -1: out.add(t)
    return out
R_ = 30
aa, bb = np.meshgrid(np.arange(-R_, R_ + 1), np.arange(-R_, R_ + 1), indexing="ij"); aa = aa.ravel(); bb = bb.ravel()
NN = norm(aa, bb); keep = (NN > 0) & (NN <= 400); aa, bb, NN = aa[keep], bb[keep], NN[keep]
# group by (norm, factor set restricted to primes outside S)
groups = {}
for a_, b_, n_ in zip(aa, bb, NN):
    fs = frozenset(t for t in afactors(a_, b_) if PP[t].N > 4 and PP[t].q != 3)
    groups[(int(n_), fs)] = groups.get((int(n_), fs), 0) + 1
print("a-groups:", len(groups), flush=True)
Hd = {}
def exc_pair(i, j, Conly1=True):
    A1, A2 = sc[i], sc[j]
    n1, n2 = len(A1["N"]), len(A2["N"])
    NN12 = A1["N"][:, None] * A2["N"][None, :]
    coprime = np.array([[len(x & y) == 0 for y in A2["idx"]] for x in A1["idx"]])
    # R(s1,s2), class characters
    Rs = RR[A1["ci"]][:, A2["ci"]]
    z1, z2 = Vc[A1["ci"], 0] % 6, Vc[A2["ci"], 0] % 6
    same = (z1[:, None] == z2[None, :])
    m2 = np.conj(cm[A2["ci"]])[None, :]
    tot = 0.0
    # divisibility masks of s by small primes
    for (mu, qd, dch) in dlist:
        dfac = np.conj(dch[A1["ci"]])[:, None] * dch[A2["ci"]][None, :]
        Qp = Q / qd
        for il in range(6):
            for j2 in range(6):
                qe = 3.0 ** il * 4.0 ** j2
                ce1 = (cl[A1["ci"]] ** il * c2[A1["ci"]] ** j2)[:, None]
                ce2 = np.conj(cl[A2["ci"]] ** il * c2[A2["ci"]] ** j2)[None, :]
                asum = np.zeros((n1, n2))
                for (n_, fs), mult in groups.items():
                    targ = Qp * qe * float(n_) ** 6 / NN12
                    if targ.min() > 1000: continue
                    hv = Hf(targ) * mult / 6
                    if fs:
                        ok1 = np.array([len(x & fs) == 0 for x in A1["idx"]]); ok2 = np.array([len(x & fs) == 0 for x in A2["idx"]])
                        hv = hv * ok1[:, None] * ok2[None, :]
                    asum += hv
                tot += mu * Qp * np.sum((A1["w"][:, None] * A2["w"][None, :]) / np.sqrt(NN12) * coprime * Rs * 6 * same * m2 * dfac * ce1 * ce2 * asum)
    return tot
r_ = np.linspace(1, 2, 200001); What = np.trapezoid(bump(r_, 1, 2) * r_ ** (-5 / 6), r_)
A = np.load("data/A_Q%d_K12.npy" % Q) if Q in (1500, 3000, 6000) else None
ma, mb, mN = np.load("data/m_Q%d_K12.npy" % Q); Om = bump(mN / Q, .5, 1)
G = ((A * Om) @ A.conj().T).real
pred = pickle.load(open("exactlat.pkl", "rb"))["pred"][int(Q)]
for (i, j) in pairs:
    g = exc_pair(i, j)
    u2 = What ** 2 * (256.0 ** 2 * 2 ** (i + j)) ** (-1 / 3)
    print("Q=%d (i,j)=(%d,%d): brute exceptional C=1: kappa=%.3f   density-formula (all C): %.3f   data kappa=%.3f" % (Q, i, j, g.real / u2, pred[i, j], G[i, j] / u2), flush=True)
