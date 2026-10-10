"""Systematic search of the rank-4 Weyl groupoid (t, A, m, s) of the theta pairing T
(wf2/A2_theory.md Sec. 4, wf2/weyl_groupoid.py) for functional-equation images in which the
double residue T = Res_t Res_s (...) is a short bilinear form.

Conventions (as in weyl_groupoid.py, Sawin-Whitehead 2507.08662 Sec 5.1, Heckenberger):
  n = 6, zeta = e^{2 pi i/6};  q_ii = -zeta^{M_ii} = zeta^{M_ii+3},  qt_ij = zeta^{M_ij}.
  Basis E: rows E_k = new simple roots written in the ORIGINAL simple roots (e_t, e_A, e_m, e_s).
  New Mellin variables u'_k = <E_k, x>  (centred), i.e. u' = E u.
  Lengths (log-sizes of the summation variables, in units of log Z):  l' = E^{-T} l.
  Residue-adapted chamber: e_t = +-E_{k_t} and e_s = +-E_{k_s}; then the double residue is a simple
  double residue in the new series, and the two remaining ('free') variables have lengths
  (E^{-T} l0)_k, l0 = (0, 1, M, 0), independent of the (infinite) residue lengths.
"""
import itertools, sys
from fractions import Fraction as Fr
import numpy as np

n = 6
NAMES = ['t', 'A', 'm', 's']
M0 = [[5, 4, 0, 0], [4, 0, 1, 0], [0, 1, 0, 5], [0, 0, 5, 4]]

def setup(M):
    r = len(M)
    B = [[0]*r for _ in range(r)]
    for i in range(r):
        B[i][i] = (M[i][i] + n // 2) % n
        for j in range(i + 1, r):
            B[i][j] = M[i][j] % n
    chi = lambda a, b: sum(a[i] * B[i][j] * b[j] for i in range(r) for j in range(r)) % n
    return r, chi

def enumerate_bases(M, cap=100000):
    r, chi = setup(M)
    def mij(E, i, j):
        qi = chi(E[i], E[i]); qt = (chi(E[i], E[j]) + chi(E[j], E[i])) % n
        for m in range(0, 4 * n + 1):
            if ((m + 1) * qi) % n == 0 or (m * qi + qt) % n == 0:
                return m
    E0 = tuple(tuple(int(i == k) for k in range(r)) for i in range(r))
    seen = {E0: None}; stack = [E0]; parent = {E0: (None, None)}
    while stack:
        E = stack.pop()
        for i in range(r):
            if chi(E[i], E[i]) == 0:
                continue
            newE = tuple(tuple(-x for x in E[i]) if j == i else
                         tuple(E[j][k] + mij(E, i, j) * E[i][k] for k in range(r)) for j in range(r))
            if newE not in seen:
                seen[newE] = None; parent[newE] = (E, i); stack.append(newE)
                if len(seen) > cap: raise RuntimeError("cap")
    return list(seen), parent, chi, mij

def Mof(E, chi):
    r = len(E)
    return tuple(tuple(((chi(E[i], E[i]) - n // 2) % n) if i == j else
                       ((chi(E[i], E[j]) + chi(E[j], E[i])) % n) for j in range(r)) for i in range(r))

def word(E, parent):
    w = []
    while parent[E][0] is not None:
        E, i = parent[E]; w.append(i)
    return list(reversed(w))

def order_of(q):          # order of zeta^q
    from math import gcd
    return n // gcd(n, q % n) if q % n else 1

if __name__ == "__main__":
    Ms = [Fr(80, 100), Fr(826, 1000), Fr(86, 100)]
    bases, parent, chi, mij = enumerate_bases(M0)
    print("bases:", len(bases))
    objs = {}
    for E in bases:
        objs.setdefault(Mof(E, chi), []).append(E)
    print("objects:", len(objs))
    for k, (Mm, Es) in enumerate(objs.items()):
        print(" object %d (%d chambers): M = %s" % (k, len(Es), [list(r) for r in Mm]))
    et = (1, 0, 0, 0); es = (0, 0, 0, 1)
    adapted = []
    for E in bases:
        kt = [k for k in range(4) if E[k] in (et, tuple(-x for x in et))]
        ks = [k for k in range(4) if E[k] in (es, tuple(-x for x in es))]
        if kt and ks:
            adapted.append((E, kt[0], ks[0]))
    print("residue-adapted chambers:", len(adapted))
    rows = {}
    for E, kt, ks in adapted:
        Em = np.array(E, dtype=float)
        EinvT = np.linalg.inv(Em).T
        Mm = Mof(E, chi)
        free = [k for k in range(4) if k not in (kt, ks)]
        # Heckenberger exponents in the target object (numerators of residue nodes)
        def m_in(i, j):
            return mij(E, i, j)
        lens = {}
        for Mv in Ms:
            l0 = np.array([0.0, 1.0, float(Mv), 0.0])
            lp = EinvT @ l0
            lens[float(Mv)] = tuple(round(lp[k], 6) for k in free)
        # symbolic: coefficient of (1, M) in each free length
        cA = EinvT @ np.array([0, 1.0, 0, 0]); cM = EinvT @ np.array([0, 0, 1.0, 0])
        sym = tuple((round(cA[k], 6), round(cM[k], 6)) for k in free)
        # which original variables the free/new variables are built from: row of E^{-T}?  rows of E give roots.
        desc = dict(
            word=''.join(NAMES[i] for i in word(E, parent)),
            kt=kt, ks=ks, free=free,
            M=[list(r) for r in Mm],
            ord_t=order_of(Mm[kt][kt] + 3), ord_s=order_of(Mm[ks][ks] + 3),
            numer_t={k: m_in(kt, k) for k in range(4) if k != kt},
            numer_s={k: m_in(ks, k) for k in range(4) if k != ks},
            free_types=[(Mm[k][k], order_of(Mm[k][k] + 3)) for k in free],
            free_edge=Mm[free[0]][free[1]],
            sym=sym, lens=lens, E=E)
        key = (sym, tuple(sorted(desc['numer_t'].items())), tuple(sorted(desc['numer_s'].items())), desc['free_edge'],
               tuple(desc['free_types']), Mm[kt][ks])
        rows.setdefault(key, []).append(desc)
    print("distinct residue-adapted representations (up to chamber multiplicity):", len(rows))
    for key, ds in sorted(rows.items(), key=lambda kv: max(kv[1][0]['lens'][0.826])):
        d = ds[0]
        print("-" * 100)
        print("word=%-12s mult=%d  kt=%d ks=%d free=%s  M_target=%s" % (d['word'] or 'id', len(ds), d['kt'], d['ks'], d['free'], d['M']))
        print("   residue t: order %d, numerator exps %s ;  residue s: order %d, numerator exps %s ; M[kt][ks]=%d"
              % (d['ord_t'], d['numer_t'], d['ord_s'], d['numer_s'], d['M'][d['kt']][d['ks']]))
        print("   free nodes (M_kk, order of Gauss char): %s ; free-free edge M=%d" % (d['free_types'], d['free_edge']))
        print("   free lengths = (a*1 + b*M):", d['sym'], " at M=.80/.826/.86:", d['lens'])
        print("   E =", d['E'])
