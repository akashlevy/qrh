"""Complete chamber table of the exact 3-node model U = (s, m, A) [[4,5,0],[5,0,1],[0,1,5]]:
every chamber, its object, the lengths of ALL three summation variables (Lam = residue length -> infinity,
M = log Q/log Z), whether e_s is +-simple (residue-adapted), and Moebius (non-admissible) nodes in the object."""
import sys
sys.path.insert(0, '.')
from groupoid_bfs import bfs
from groupoid_search import Mof
import sympy as sp
Lam, M = sp.symbols('Lam M')
U = [[4, 5, 0], [5, 0, 1], [0, 1, 5]]
names = ['s', 'm', 'A']
dist, chi, mij = bfs(U, cap=10000, maxcoef=50)
objs = {}
for E in dist:
    objs.setdefault(Mof(E, chi), len(objs))
print("chambers:", len(dist), " objects:", len(objs))
for Mm, k in objs.items():
    print("  object %d: %s  Moebius nodes: %s" % (k, [list(r) for r in Mm], [i for i in range(3) if Mm[i][i] == 3]))
rows = []
for E, w in sorted(dist.items(), key=lambda kv: (len(kv[1]), kv[1])):
    lp = sp.Matrix(E).T.inv() * sp.Matrix([Lam, M, 1])
    adapted = any(E[k] in ((1, 0, 0), (-1, 0, 0)) for k in range(3))
    Mm = Mof(E, chi)
    rows.append((''.join(names[i] for i in w) or 'id', objs[Mm], E, [sp.expand(x) for x in lp], adapted,
                 [i for i in range(3) if Mm[i][i] == 3]))
print("%-10s %3s %-34s %-44s %-8s %s" % ('word', 'obj', 'E (rows = new simple roots in e_s,e_m,e_A)', 'lengths of the 3 variables', 'adapted', 'Moeb'))
for r in rows:
    print("%-10s %3d %-34s %-44s %-8s %s" % (r[0], r[1], str(r[2]), str(r[3]), r[4], r[5]))
