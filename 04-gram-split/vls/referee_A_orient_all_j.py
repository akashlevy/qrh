"""Referee check: for 3-node (t, A, m) encodings with t-node (M_tt, e_tA) and A-m edge j, find chambers
with e_t simple and free lengths (2lm - lA, lm) [paper's reflected lengths], and compare the new
free edge and the new row-node exponent against the paper's Prop completed-reflection:
  B_p = chi^{-j-2} (j != 0,4)  -> dual edge (-j-2) mod 6
  row omega_{p,j} = tau_j tau_{j+2}  -> S-W exponent (j-3)+(j+2-3) = 2j-4 mod 6   (j != 0,4)
  j = 4: B_p = q^{-1/2}(-1 + q 1_{p|x}) (Ramanujan, 'principal', exponent 0 with Moebius-type), row tau_4 -> exponent 1
  j = 0 (active): B_p = chi^{-2} -> edge 4, row tau_2 -> exponent 5
Also the 3-node U-type direct model (A Kubota with M_AA) reflected at A."""
import sys
sys.path.insert(0, '/Users/akashlevy/Documents/riemann/gram_split/vls')
from groupoid_bfs import bfs
from groupoid_search import Mof
import sympy as sp
lA, lm = sp.symbols('lA lm')
def pred(j):
    if j == 4: return (0, 1)
    if j == 0: return (4, 5)
    return ((-j - 2) % 6, (2 * j - 4) % 6)
for (Mtt, etA, name) in [(5, 4, 'M0-type'), (1, 2, 'M1-type'), (5, 2, 'mixed52'), (1, 4, 'mixed14')]:
    print('====', name, 'M_tt=%d e_tA=%d' % (Mtt, etA))
    for j in range(6):
        M = [[Mtt, etA, 0], [etA, 0, j], [0, j, 0]]
        dist, chi, mij = bfs(M, cap=5000, maxcoef=30)
        res = set()
        for E, w in dist.items():
            kt = [k for k in range(3) if E[k] in ((1, 0, 0), (-1, 0, 0))]
            if not kt: continue
            kt = kt[0]; free = [k for k in range(3) if k != kt]
            lp = sp.Matrix(E).T.inv() * sp.Matrix([0, lA, lm])
            L = [sp.expand(lp[k]) for k in free]
            if L == [2*lm - lA, lm] or L == [lm, 2*lm - lA]:
                Mm = Mof(E, chi)
                # identify which free node is the row (length lm)
                ir = free[L.index(lm)]; ic = free[L.index(2*lm - lA)]
                res.add((Mm[ir][ic], Mm[ir][ir], Mm[ic][ic], tuple(sorted({k: mij(E, kt, k) for k in free}.items()))))
        print('  j=%d  #ch=%5d  paper (edge,row)=%s   found (edge,row,colnode,t-numerators)=%s' % (j, len(dist), pred(j), sorted(res)))
print('==== direct 3-node (m, A) with A Kubota M_AA, reflect at A')
from groupoid_search import setup, n
for MAA in (5, 1):
    for j in range(6):
        M = [[0, j], [j, MAA]]
        r, chi = setup(M)
        qi = chi((0, 1), (0, 1))
        if qi == 0: print('  M_AA=%d j=%d nonadmissible' % (MAA, j)); continue
        qt = (chi((0, 1), (1, 0)) + chi((1, 0), (0, 1))) % n
        for m in range(30):
            if ((m + 1) * qi) % n == 0 or (m * qi + qt) % n == 0: break
        E = ((1, m), (0, -1))
        Mm = Mof(E, chi)
        print('  M_AA=%d j=%d m_Am=%d -> edge=%d row=%d col=%d   paper=%s' % (MAA, j, m, Mm[0][1], Mm[0][0], Mm[1][1], pred(j)))
