"""Step 1: exact row scalar zeta_m for the data rows; validation.
 (a) gamma_j(R) by twisted multiplicativity  vs  brute-force Gauss sum over O/R (incl. inert factors, conjugate pairs)
 (b) Lemma Z closed form psi_1(R) conj(gamma_2(R))^4  vs  the paper's product formula (reflection-row-phase), c_F=1
 (c) prime rows vs check_row_phase.py's w = conj(alpha)^2 tau^-_{p,1} tau^+_{p,3}:  w * gamma_2^4 constant on classes mod 72
"""
import numpy as np, sys, time
from a2lib import *
t0 = time.time()
Q = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
R = Rows(Q)
print("Q=%d rows %d, squarefree %d (%.1f%%), prime rows %d, max omega %d   [%.1fs]" %
      (Q, R.n, R.sf.sum(), 100 * R.sf.mean(), (R.sf & (R.nfac == 1)).sum(), R.nfac.max(), time.time() - t0))
# (a) brute force on a sample of composite squarefree rows of small norm
rng = np.random.default_rng(1)
cand = [t for t in range(R.n) if R.sf[t] and R.nfac[t] >= 2]
cand = sorted(cand, key=lambda t: R.mN[t])
sample = cand[:40] + list(rng.choice(cand, 20, replace=False))
# make sure inert factors and conjugate pairs are present
inert = [t for t in cand if any(R.P[k].kind == "inert" for k in R.fac[t])][:10]
conjp = [t for t in cand if len(set(R.P[k].q for k in R.fac[t])) < len(R.fac[t])][:10]
sample = list(dict.fromkeys(sample + inert + conjp))
err = {1: 0, 2: 0, 3: 0}
for t in sample:
    for j in (1, 2, 3):
        gb = gauss_bruteforce(int(R.Ra[t]), int(R.Rb[t]), [R.P[k] for k in R.fac[t]], j)
        err[j] = max(err[j], abs(gb - R.G[j][t]))
print("(a) %d composite rows (%d with inert factor, %d with pi*conj(pi)); max |gamma_j(R) twisted-mult - bruteforce|: j=1 %.1e  j=2 %.1e  j=3 %.1e"
      % (len(sample), len(inert), len(conjp), err[1], err[2], err[3]))
# signal identity gamma_1 gamma_2 = mu alpha G, G = conj(chi_R(4)) gamma_3, on all squarefree rows
sf = np.where(R.sf)[0]
mu = (-1.0) ** R.nfac
chi4 = np.array([R.chiR(t, R.e_4) for t in range(R.n)])
e_sig = np.abs(R.G[1][sf] * R.G[2][sf] - mu[sf] * R.alpha[sf] * np.conj(chi4[sf]) * R.G[3][sf]).max()
e_cub = np.abs(R.G[2][sf] ** 3 - mu[sf] * R.alpha[sf]).max()
print("    eq.signal on all %d sf rows: max|g1 g2 - mu alpha conj(chi4) g3| = %.1e ; max|g2^3 - mu alpha| = %.1e" % (len(sf), e_sig, e_cub))
# (b) Lemma Z vs product formula
ratio = R.zetaP[sf] / R.zetaZ[sf]
print("(b) zeta_P / zeta_Z over %d squarefree rows: max|ratio - 1| = %.2e   (|zeta|=1: max||zetaZ|-1| = %.1e)"
      % (len(sf), np.abs(ratio - 1).max(), np.abs(np.abs(R.zetaZ[sf]) - 1).max()))
# (c) prime rows vs check_row_phase.py convention
pr = [t for t in sf if R.nfac[t] == 1 and R.P[R.fac[t][0]].kind == "split"]
from collections import defaultdict
cls = defaultdict(list)
for t in pr:
    k = R.fac[t][0]; Pr = R.P[k]
    w = np.conj(R.alpha[t]) ** 2 * R.gm1[k] * R.g[3][k]
    cls[(Pr.a % 72, Pr.b % 72)].append((w * R.g[2][k] ** 4, R.zetaZ[t] * R.g[2][k] ** 4, w))
nc = bad = badZ = 0; spr = []
for key, L in cls.items():
    if len(L) < 2: continue
    nc += 1
    v = np.array([x[0] for x in L]); vz = np.array([x[1] for x in L]); W = np.array([x[2] for x in L])
    bad += np.max(np.abs(v - v[0])) > 1e-9; badZ += np.max(np.abs(vz - vz[0])) > 1e-9
    spr.append(np.max(np.abs(W - W[0])))
print("(c) %d distinct split prime rows R (primary); classes mod 72 with >=2 primes: %d" % (len(set(int(R.fac[t][0]) for t in pr)), nc))
print("    classes where w*gamma2^4 (check_row_phase) is non-constant: %d ; where zetaZ*gamma2^4 is non-constant: %d ; median spread of w itself: %.2f"
      % (bad, badZ, np.median(spr)))
np.savez("zeta_Q%d.npz" % Q, ma=R.ma, mb=R.mb, sf=R.sf, zetaZ=R.zetaZ, zetaP=R.zetaP, u=R.u, Ra=R.Ra, Rb=R.Rb)
print("saved zeta_Q%d.npz  [%.1fs]" % (Q, time.time() - t0))
