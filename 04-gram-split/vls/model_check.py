"""Which encodings reproduce the two PROVEN structural facts?
 (F-refl) paper Prop. completed-reflection: column twist chi^j (j=1) -> dual column factor chi^{-j-2}=chi^3,
          row node acquires tau_j tau_{j+2} = gamma1*gamma3 (S-W exponent (1-3)+(3-3) = 4);
 (F-sig)  V_gaps 1.6 [R]: m-Poisson/Voronoi of the unreflected pairing turns the column Gauss sum gamma2 into mu(A)
          (Moebius node, S-W exponent 3)."""
import sys
sys.path.insert(0, '/Users/akashlevy/Documents/riemann/gram_split/vls')
from groupoid_bfs import bfs
from groupoid_search import Mof
import sympy as sp

def reflect_once(M, i):
    from groupoid_search import setup, n
    r, chi = setup(M)
    def mij(E, i, j):
        qi = chi(E[i], E[i]); qt = (chi(E[i], E[j]) + chi(E[j], E[i])) % n
        for m in range(0, 4 * n + 1):
            if ((m + 1) * qi) % n == 0 or (m * qi + qt) % n == 0:
                return m
    E0 = tuple(tuple(int(a == b) for b in range(r)) for a in range(r))
    if chi(E0[i], E0[i]) == 0: return None, None
    E = tuple(tuple(-x for x in E0[i]) if j == i else tuple(E0[j][k] + mij(E0, i, j) * E0[i][k] for k in range(r)) for j in range(r))
    return E, Mof(E, chi)

print("3-node U = (s, m, A), A carries Patterson gamma_2 (M_AA=5), edge m-A = 1:")
U = [[4, 5, 0], [5, 0, 1], [0, 1, 5]]
E, Mn = reflect_once(U, 2); print("  reflect at A (cubic theta Voronoi):", [list(r) for r in Mn], " E=", E)
E, Mn = reflect_once(U, 1); print("  reflect at m (Poisson/Voronoi in m):", [list(r) for r in Mn], " E=", E)
print("  conjugate-theta variant M_AA=1:")
U1 = [[4, 5, 0], [5, 0, 1], [0, 1, 1]]
E, Mn = reflect_once(U1, 2); print("  reflect at A:", [list(r) for r in Mn])
E, Mn = reflect_once(U1, 1); print("  reflect at m:", [list(r) for r in Mn])

for name, M in [("M0 (A2_theory)", [[5, 4, 0, 0], [4, 0, 1, 0], [0, 1, 0, 5], [0, 0, 5, 4]]),
                ("M1 (t conj)   ", [[1, 2, 0, 0], [2, 0, 1, 0], [0, 1, 0, 5], [0, 0, 5, 4]])]:
    dist, chi, mij = bfs(M, cap=20000)
    objs = set(Mof(E, chi) for E in dist)
    print("\n4-node", name, ": chambers", len(dist), " objects", len(objs))
    # find the 'At' chamber
    for E, w in sorted(dist.items(), key=lambda kv: len(kv[1])):
        if w == [1, 0]:
            Mm = Mof(E, chi); print("   word At:", [list(r) for r in Mm], "E=", E)
        if w == [2, 3]:
            Mm = Mof(E, chi); print("   word ms:", [list(r) for r in Mm], "E=", E)
