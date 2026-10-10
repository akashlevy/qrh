"""Attack E side test (exact, noise-free, any Q): is the row scalar zeta_R coherent with the cubic-theta column sum B0?

  TBz(L; cl) = sum_{R in class cl} Om(N R/Q) zeta_R B0_R(L),   B0_R(L) = sum_x beta_x w_L(x) (lambda x / R)_2,
  R primary squarefree, (R,6)=1, N R in (Q/2, Q];  zeta_R = psi1(R) conj(gamma2(R))^4 (Lemma Z, kappa_0 dropped);
  classes cl = R mod F (F = 1, 4, 12, 36: absorbs ray factors and kappa_0).
  rat = sum_cl |TBz|^2 / sum_R Om^2 |B0_R|^2      (random: 1;  coherent / main term: grows with Q)
Also the same with zeta replaced by gamma2, conj(gamma2), random phases, and "1" (calibrations).
usage: python E_zeta_beta.py Q1 Q2 ...
"""
import numpy as np, sys, time, json
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *
from cubic_theta import cubic_theta_coeffs

t0 = time.time()
Qs = [int(x) for x in sys.argv[1:]] or [6000]
res = {}
for Q in Qs:
    P = primes_upto(Q)
    S = [t for t in squarefree_products(P, Q) if t[2] > Q / 2]
    Ra = np.array([t[0] for t in S], np.int64); Rb = np.array([t[1] for t in S], np.int64); RN = np.array([t[2] for t in S], np.int64)
    fac = [t[3] for t in S]
    used = sorted(set(k for f in fac for k in f))
    g2 = np.zeros(len(P), complex)
    e_lam = np.zeros(len(P), int); e_4 = np.zeros(len(P), int); e_m1 = np.zeros(len(P), int)
    for k in used:
        g2[k] = gamma_prime(P[k], 2)
        e_lam[k] = sym1(P[k], *LAM); e_4[k] = sym1(P[k], 4, 0); e_m1[k] = sym1(P[k], -1, 0)
    cross = {}
    G2 = np.zeros(len(S), complex); psi1 = np.zeros(len(S), complex)
    for t, f in enumerate(fac):
        g = 1 + 0j
        for x in range(len(f)):
            g *= g2[f[x]]
            for y in range(x + 1, len(f)):
                key = (f[x], f[y])
                if key not in cross:
                    cross[key] = sym1(P[f[x]], P[f[y]].a, P[f[y]].b) + sym1(P[f[y]], P[f[x]].a, P[f[x]].b)
                g *= Z6[(2 * cross[key]) % 6]
        G2[t] = g
        el = sum(e_lam[k] for k in f); e4 = sum(e_4[k] for k in f); em = sum(e_m1[k] for k in f)
        psi1[t] = np.conj(Z6[el % 6]) * np.conj(Z6[e4 % 6]) * Z6[em % 6]
    zeta = psi1 * np.conj(G2) ** 4
    Om = bump(RN / Q, .5, 1)
    alR = np.array([cval(a, b) for a, b in zip(Ra, Rb)]); alR /= np.abs(alR)
    print("Q=%d: %d primary sf ideals in (Q/2,Q] [%.0fs]" % (Q, len(S), time.time() - t0), flush=True)
    # twist matrix (lambda x / R)_2 via incidence over used primes
    xa, xb, xN, ax, _, _ = cubic_theta_coeffs(Q // 16, Q)
    xa = xa.astype(np.int64); xb = xb.astype(np.int64); xN = xN.astype(np.int64)
    alx = np.array([cval(a, b) for a, b in zip(xa, xb)]); alx /= np.abs(alx)
    beta = ax * alx / np.sqrt(xN)
    nua, nub = mul(np.int64(1), np.int64(2), xa, xb)
    import scipy.sparse as sps
    rows = np.repeat(np.arange(len(S)), [len(f) for f in fac]); cols = np.concatenate([np.array(f) for f in fac])
    Inc = sps.csr_matrix((np.ones(len(rows), np.int64), (rows, cols)), shape=(len(S), len(P)))
    Ep = np.zeros((len(P), len(xN)), np.int8); Zp = np.zeros((len(P), len(xN)), np.int8)
    for k in used:
        e = P[k].symbol(nua, nub); Zp[k] = e < 0; Ep[k] = np.where(e >= 0, e, 0) % 2
    TW = np.zeros((len(S), len(xN)), np.float64)
    for c0 in range(0, len(xN), 3000):
        par = (Inc @ Ep[:, c0:c0 + 3000].astype(np.int64)) % 2; zer = (Inc @ Zp[:, c0:c0 + 3000].astype(np.int64)) > 0
        TW[:, c0:c0 + 3000] = np.where(zer, 0, 1 - 2 * par)
    print("   twist matrix %s [%.0fs]" % (TW.shape, time.time() - t0), flush=True)
    rng = np.random.default_rng(5)
    seqs = {"zeta": zeta, "conj(g2)": np.conj(G2), "g2": G2, "conj(g2)*alpha": np.conj(G2) * alR, "rand": np.exp(2j * np.pi * rng.random(len(S))), "1": np.ones(len(S))}
    out = {}
    for F in (1, 12, 36):
        lab = np.unique((Ra % F) * F + (Rb % F), return_inverse=True)[1]; ncl = lab.max() + 1
        for L in [Q // 8, Q // 4, Q // 2, Q]:
            wL = bump(xN / L, 0.5, 1.0); sel = wL > 0
            B0 = TW[:, sel] @ (beta[sel] * wL[sel])
            den = (Om ** 2 * np.abs(B0) ** 2).sum()
            for name, z in seqs.items():
                v = Om * z * B0
                num = (np.abs(np.bincount(lab, weights=v.real, minlength=ncl) + 1j * np.bincount(lab, weights=v.imag, minlength=ncl)) ** 2).sum()
                out[(F, L, name)] = num / den
    res[Q] = out
    print("   ratio sum_cl|sum Om z B0|^2 / sum Om^2|B0|^2   (random = 1)")
    for F in (1, 12, 36):
        print("   F=%-3d " % F + "  ".join("%s: " % name + "/".join("%.2f" % out[(F, L, name)] for L in [Q // 8, Q // 4, Q // 2, Q]) for name in seqs))
json.dump({str(Q): {"%d_%d_%s" % k: v for k, v in o.items()} for Q, o in res.items()}, open(DATA + "zeta_beta.json", "w"), indent=1)
print("[%.0fs]" % (time.time() - t0))
