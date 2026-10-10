"""Test which encoding of the cubic-theta/column edge reproduces the paper's PROVEN reflection rule
(Prop. completed-reflection, paper78 l.1763ff): row twist chi_p^j on the column variable A  ==>  dual column
factor B_p = chi_p^{-j-2}, and the row acquires tau_{p,j}^- tau_{p,j+2}^+ (Gauss sums of chi^j, chi^{j+2}).
3-node model (t, A, m): t cubic Kubota [M_tt=5], edge t-A = e_tA, A Dirichlet [0], edge A-m = j, m: [M_mm].
We look for chambers with e_t simple, whose free variables have lengths (2*l_m - l_A, l_m) (the paper's
reflected lengths) and print the new edge and new M_mm."""
import sys
sys.path.insert(0, '/Users/akashlevy/Documents/riemann/gram_split/vls')
from groupoid_bfs import bfs
from groupoid_search import Mof
import sympy as sp
lA, lm = sp.symbols('lA lm')
for etA in (2, 4):
    for j in range(6):
        for Mmm in (0,):
            M = [[5, etA, 0], [etA, 0, j], [0, j, Mmm]]
            try:
                dist, chi, mij = bfs(M)
            except Exception as e:
                print(etA, j, 'fail', e); continue
            if len(dist) > 3000:
                print("e_tA=%d j=%d: groupoid large (%d)" % (etA, j, len(dist))); continue
            hits = []
            for E, w in dist.items():
                kt = [k for k in range(3) if E[k] in ((1, 0, 0), (-1, 0, 0))]
                if not kt: continue
                kt = kt[0]
                free = [k for k in range(3) if k != kt]
                lp = sp.Matrix(E).T.inv() * sp.Matrix([0, lA, lm])
                L = [sp.expand(lp[k]) for k in free]
                Mm = Mof(E, chi)
                hits.append((len(w), w, kt, free, L, Mm[free[0]][free[1]], [Mm[k][k] for k in free],
                             {k: mij(E, kt, k) for k in free}))
            hits.sort()
            refl = [h for h in hits if sp.expand(h[4][0] - (2*lm - lA)) == 0 or sp.expand(h[4][1] - (2*lm - lA)) == 0]
            print("e_tA=%d j=%d  (#chambers %d)  paper predicts dual edge %d" % (etA, j, len(dist), (-j - 2) % 6))
            for h in refl[:3]:
                print("     word=%s kt=%d free lengths=%s  new free edge=%d  free M_kk=%s  t-numerators=%s" % (h[1], h[2], h[4], h[5], h[6], h[7]))
            if not refl:
                print("     no chamber with reflected length 2lm-lA")
