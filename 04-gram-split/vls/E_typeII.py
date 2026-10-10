"""Attack E: does pi(m) = rho(m) zeta_m have bilinear (Type II) structure?

For squarefree ideals R = a b (a prime, N(a) in a band [A1, A2), b = R/a), unit class c:
   M_c[a, b] = pi(zeta6^c a b),       twisted:  Mt = M * zeta6^{-j (e(a,b) + e(b,a))},   e(a,b) = exponent of chi_a(b).
Twisted multiplicativity  pi(ab) = f(a) g(b) zeta6^{j(e(a,b)+e(b,a))}  <=>  Mt is rank one.
(Calibration: zeta = psi1 conj(gamma2)^4 is of this form with j = 4, since gamma2(ab) = (chi_a(b)chi_b(a))^2 gamma2(a)gamma2(b).)

4-cycle statistic (unbiased: the four entries come from four distinct scales, averaged over the 24 orders of 4 scales):
   Gam_j = sum_c sum_{a != a', b != b'} Mt[a,b] conj(Mt[a,b']) conj(Mt[a',b]) Mt[a',b']
   Esym  = sum_c sum_{4-cycles} (P[a,b]P[a',b'] + P[a,b']P[a',b]) / 2,   P = |pi|^2 (cross-scale)
   Gam/Esym <= 1, = 1 iff rank one with aligned phases;  random: 0 (null s.d. from random-sign nulls).
Spectral test: top eigenvalue of the unbiased Gram  Mt^H Mt  (cross-scale)  vs the same with random signs.
usage: python E_typeII.py Q src ymin
"""
import numpy as np, sys, time, json, itertools
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *

t0 = time.time()
Q = int(sys.argv[1]); src = sys.argv[2] if len(sys.argv) > 2 else "E"; ymin = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
FPROJ = int(sys.argv[4]) if len(sys.argv) > 4 else 0      # if > 0: remove the class-function (mod FPROJ) component of rho first
NOCONJ = len(sys.argv) > 5 and sys.argv[5] == "noconj"    # keep one prime of each conjugate pair {a, abar} only
# (V(mbar) = conj V(m) exactly, so cycles (a, abar; b, bbar) contribute |pi(ab)|^2 |pi(a bbar)|^2 >= 0 for every j-twist that is trivial on them)
R = RowsE(Q, src=src, ymin_over_Q=ymin)
W = R.w; K = len(R.scales)
if FPROJ:
    key = (R.ma % FPROJ) * FPROJ + (R.mb % FPROJ)
    lab = np.unique(key, return_inverse=True)[1]
    for i in range(K):
        sw = np.bincount(lab, weights=R.Om, minlength=lab.max() + 1)
        mr = np.bincount(lab, weights=R.Om * R.V[i].real, minlength=lab.max() + 1) / sw
        mi = np.bincount(lab, weights=R.Om * R.V[i].imag, minlength=lab.max() + 1) / sw
        R.V[i] = R.V[i] - (mr + 1j * mi)[lab]
    print("removed class means mod %d (%d classes) from rho, per scale" % (FPROJ, lab.max() + 1))
IX, OmI, Vp, pi2 = pi_layout(R)
nI = len(OmI); rep = IX[:, 0]
Vr = R.V[:, IX]; rho2 = cross_est(np.moveaxis(Vr, 0, -1), W)
print("Q=%d ideals %d scales %s [%.0fs]" % (Q, nI, R.Ys[R.scales].astype(int).tolist(), time.time() - t0), flush=True)
# top 4 scales by weight (least noisy) for the 4-cycle statistic
top4 = np.argsort(-W)[:4]
orders = list(itertools.permutations(top4, 4))
rng = np.random.default_rng(7)


def build(A1, A2):
    ent = []                    # (ia, ib, I, e)
    akey, bkey = {}, {}
    for I in range(nI):
        f = R.fac[rep[I]]
        for k in f:
            Pa = R.P[k]
            if not (A1 <= Pa.N < A2): continue
            if NOCONJ and Pa.kind == "split" and Pa.b < 0: continue
            rest = [k2 for k2 in f if k2 != k]
            ba, bb = 1, 0
            for k2 in rest: ba, bb = mul(ba, bb, R.P[k2].a, R.P[k2].b)
            e1 = sym1(Pa, ba, bb)                                   # chi_a(b)
            e2 = sum(sym1(R.P[k2], Pa.a, Pa.b) for k2 in rest)      # chi_b(a)
            ia = akey.setdefault(k, len(akey)); ib = bkey.setdefault((ba, bb), len(bkey))
            ent.append((ia, ib, I, (e1 + e2) % 6))
    return np.array(ent, np.int64).reshape(-1, 4), len(akey), len(bkey)


def gam(Ms, Pm):
    """Ms: list of 4 (na, nb) complex matrices (one per scale, same twist); Pm: (na, nb) |pi|^2 estimate; returns (Gam, Esym)."""
    M1, M2, M3, M4 = Ms
    G13 = M1 @ M3.conj().T; H24 = M2.conj() @ M4.T
    full = np.sum(G13 * H24)
    aa = np.sum(np.diag(G13) * np.diag(H24))
    bb = np.sum((M1 * M2.conj()).sum(0) * (M3.conj() * M4).sum(0))
    ab = np.sum(M1 * M2.conj() * M3.conj() * M4)
    return full - aa - bb + ab


def esym(Pm, pres):
    X = Pm @ pres.T                       # X[a,a'] = sum_b P[a,b] pres[a',b]
    full = np.sum(X * X.T)                # sum P[a,b] pres[a',b] P[a',b'] pres[a,b']
    aa = np.sum(np.diag(X) ** 2)
    bb = np.sum((Pm * pres).sum(0) ** 2)
    ab = np.sum((Pm * pres) ** 2)
    return full - aa - bb + ab            # (the two pairings are equal by symmetry a <-> a')


bands = [(7, 60), (60, 400), (400, 3000), (3000, 20000)]
results = []
for (A1, A2) in bands:
    ent, na, nb = build(A1, A2)
    if len(ent) < 50 or na < 3 or nb < 3: continue
    ia, ib, II, ee = ent.T
    pres = np.zeros((na, nb)); pres[ia, ib] = 1
    line = {"band": (A1, A2), "na": na, "nb": nb, "entries": len(ent)}

    def stat(xs, x2, exact=False, signs=None):
        """xs: (K, nI, 6) per-scale values (or (nI,6) exact); returns Gam/Esym per j."""
        out = []
        for j in range(6):
            ph = np.conj(Z6[(j * ee) % 6])
            G = 0; E = 0
            for c in range(6):
                Pm = np.zeros((na, nb)); Pm[ia, ib] = (np.abs(xs[II, c]) ** 2 if exact else x2[II, c])
                if exact:
                    Mt = np.zeros((na, nb), complex); Mt[ia, ib] = xs[II, c] * ph * (1 if signs is None else signs[II, c])
                    G += gam([Mt, Mt, Mt, Mt], Pm)
                else:
                    Mk = []
                    for k in top4:
                        Mt = np.zeros((na, nb), complex); Mt[ia, ib] = xs[k, II, c] * ph * (1 if signs is None else signs[II, c]); Mk.append(Mt)
                    Md = dict(zip(top4, Mk))
                    G += np.mean([gam([Md[o[0]], Md[o[1]], Md[o[2]], Md[o[3]]], Pm) for o in orders])
                E += esym(Pm, pres)
            out.append((G / E).real)
        return np.array(out)

    def spec(xs, signs=None):
        """top-3 eigenvalues / trace of the unbiased cross-scale Gram Mt^H Mt (summed over c), per j."""
        out = []
        for j in range(6):
            ph = np.conj(Z6[(j * ee) % 6]); G = 0
            for c in range(6):
                Mk = np.zeros((K, na, nb), complex)
                Mk[:, ia, ib] = xs[:, II, c] * ph * (1 if signs is None else signs[II, c])
                S = np.einsum("k,kab->ab", W, Mk)
                G = G + (S.conj().T @ S - np.einsum("k,kab,kac->bc", W ** 2, Mk.conj(), Mk)) / (W.sum() ** 2 - (W ** 2).sum())
            G = (G + G.conj().T) / 2
            ev = np.linalg.eigvalsh(G)[::-1]
            out.append(ev[:3] / np.trace(G).real)
        return np.array(out)

    sp_pi = spec(Vp); sp_rho = spec(Vr)
    sp_null = np.array([spec(Vp, rng.choice([-1.0, 1.0], size=(nI, 6))) for _ in range(6)])
    st_pi = stat(Vp, pi2); st_rho = stat(Vr, rho2)
    nulls = np.array([stat(Vp, pi2, signs=rng.choice([-1.0, 1.0], size=(nI, 6))) for _ in range(6)])
    st_z = stat(R.zetaZ[IX], None, exact=True)
    st_rand = stat(np.exp(2j * np.pi * rng.random((nI, 6))), None, exact=True)
    st_g1 = stat(R.G[1][IX], None, exact=True)
    E4 = esym(np.ones((na, nb)) * pres, pres)
    line.update(pi=st_pi.tolist(), rho=st_rho.tolist(), null_sd=nulls.std(0).tolist(), null_mean=nulls.mean(0).tolist(),
                zeta=st_z.tolist(), rand=st_rand.tolist(), gamma1=st_g1.tolist(), cycles=float(E4),
                spec_pi=sp_pi[:, 0].tolist(), spec_rho=sp_rho[:, 0].tolist(), spec_null_mean=sp_null[:, :, 0].mean(0).tolist(), spec_null_sd=sp_null[:, :, 0].std(0).tolist())
    results.append(line)
    print("\nband N(a) in [%d,%d): #a=%d #b=%d entries=%d, 4-cycles=%d  [%.0fs]" % (A1, A2, na, nb, len(ent), E4, time.time() - t0), flush=True)
    print("   j:          " + " ".join("%8d" % j for j in range(6)))
    for lab, v in (("pi", st_pi), ("rho", st_rho), ("null mean", nulls.mean(0)), ("null sd", nulls.std(0)), ("zeta (exact)", st_z),
                   ("gamma1 (ex.)", st_g1), ("rand (ex.)", st_rand)):
        print("   %-12s" % lab + " ".join("%8.4f" % x for x in v))
    print("   spectral: top eigenvalue / trace of unbiased Gram Mt^H Mt (rank-r structure test)")
    for lab, v in (("pi", sp_pi[:, 0]), ("rho", sp_rho[:, 0]), ("null mean", sp_null[:, :, 0].mean(0)), ("null sd", sp_null[:, :, 0].std(0))):
        print("   %-12s" % lab + " ".join("%8.4f" % x for x in v))

# ---- spectral test on the largest band with enough entries: top eigenvalue of unbiased Gram, j = 0..5
json.dump({"Q": Q, "bands": results}, open(DATA + "typeII_Q%d_%s_p%d%s.json" % (Q, src, FPROJ, "_noconj" if NOCONJ else ""), "w"), indent=1, default=float)
print("[%.0fs]" % (time.time() - t0))
