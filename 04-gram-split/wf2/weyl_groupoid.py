"""Finiteness of the Heckenberger Weyl groupoid attached (Sawin-Whitehead, arXiv 2507.08662, Sec. 5.1)
to an axiomatic multiple Dirichlet series with data (n, M).
Dictionary (S-W Thm 5.2 / Sec 5.1): zeta = primitive n-th root of unity;
   q_ii = chi^(e_i,e_i) = -zeta^{M_ii},   qt_ij = chi^(e_i,e_j) chi^(e_j,e_i) = zeta^{M_ij}.
Node i carries a Gauss sum with character xi*chi^{M_ii} = chi^{M_ii+n/2}; M_ii = n/2 gives q_ii = 1
(non-admissible: 'Moebius node'); M_ii = 0 gives q_ii = -1 (Dirichlet-type node).
Bicharacter realised as chi^(a,b) = zeta^{a^T B b}, B upper triangular, B_ii = M_ii + n/2, B_ij = M_ij (i<j).
Reflection: m_ij = min{m>=0 : q_ii^{m+1}=1 or q_ii^m qt_ij = 1};  s_i(e_i) = -e_i, s_i(e_j) = e_j + m_ij e_i.
Returns (finite?, number of bases, number of objects=distinct M's, non-admissible encountered?)."""
import itertools, sys

def run(M, n, cap=4000, maxcoef=40, verbose=False):
    r = len(M)
    B = [[0]*r for _ in range(r)]
    for i in range(r):
        B[i][i] = (M[i][i] + n // 2) % n
        for j in range(i + 1, r):
            B[i][j] = M[i][j] % n
    def chi(a, b):  # exponent of zeta
        return sum(a[i] * B[i][j] * b[j] for i in range(r) for j in range(r)) % n
    def data(E):
        q = tuple(chi(E[i], E[i]) for i in range(r))
        qt = tuple((chi(E[i], E[j]) + chi(E[j], E[i])) % n for i in range(r) for j in range(i + 1, r))
        return q, qt
    def mij(E, i, j):
        qi = chi(E[i], E[i]); qt = (chi(E[i], E[j]) + chi(E[j], E[i])) % n
        for m in range(0, 4 * n + 1):
            if ((m + 1) * qi) % n == 0 or (m * qi + qt) % n == 0:
                return m
        return None
    E0 = tuple(tuple(int(i == k) for k in range(r)) for i in range(r))
    seen = {E0}; stack = [E0]; nonadm = False; objs = set()
    while stack:
        E = stack.pop()
        objs.add(data(E))
        for i in range(r):
            if chi(E[i], E[i]) == 0:   # q_ii = 1 : not admissible at i
                nonadm = True; continue
            newE = []
            for j in range(r):
                if j == i:
                    newE.append(tuple(-x for x in E[i]))
                else:
                    m = mij(E, i, j)
                    newE.append(tuple(E[j][k] + m * E[i][k] for k in range(r)))
            newE = tuple(newE)
            if max(abs(x) for v in newE for x in v) > maxcoef:
                return False, len(seen), len(objs), nonadm
            if newE not in seen:
                seen.add(newE); stack.append(newE)
                if len(seen) > cap:
                    return False, len(seen), len(objs), nonadm
    return True, len(seen), len(objs), nonadm

def show(name, M, n=6):
    fin, nb, no, na = run(M, n)
    print("%-62s finite=%-5s bases=%-6d objects=%-4d nonadmissible_seen=%s" % (name, fin, nb, no, na))
    return fin

if __name__ == "__main__":
    n = 6
    print("Sanity checks")
    show("A2 sextic WMDS  [[4,5],[5,4]]", [[4,5],[5,4]])
    show("A2 cubic WMDS (n=6 language, chi^2 nodes) [[5,4],[4,5]]", [[5,4],[4,5]])
    show("FHL n-th order twists [[0,1],[1,0]] (S-W example)", [[0,1],[1,0]])
    show("A3 sextic WMDS", [[4,5,0],[5,4,5],[0,5,4]])
    show("quadratic A2 (n=6 lang): [[0,3],[3,0]]", [[0,3],[3,0]])
    print("\nReflected theta pairing: nodes (s, m, mu); s: g6(s) [M=4]; m: gamma1*gamma3 [M=4];")
    print("edge s-m sextic^{-1} [5], edge m-mu quadratic [3]; mu: cubic theta coeff = cubic Gauss sum [1 or 5] or none [0]")
    for Msm in (1, 5):
        for Mmu in (1, 5, 0, 2, 4):
            show("(s,m,mu) M_sm=%d M_mumu=%d" % (Msm, Mmu), [[4, Msm, 0], [Msm, 4, 3], [0, 3, Mmu]])
    print("\nTwo-node pieces")
    for Mmu in (1, 5, 0):
        show("(m,mu) only: M_mm=4, edge 3, M_mumu=%d" % Mmu, [[4, 3], [3, Mmu]])
    show("(s,m) with m-node gamma1 only (A2)", [[4, 5], [5, 4]])
    print("\nUnreflected pairing: nodes (s, m, A); m: no Gauss sum [0]; edge m-A sextic [1 or 5]; A: cubic Gauss gamma_2 [5 or 1]")
    for MmA in (1, 5):
        for MAA in (5, 1, 3):
            show("(s,m,A) M_mA=%d M_AA=%d" % (MmA, MAA), [[4, 5, 0], [5, 0, MmA], [0, MmA, MAA]])
    print("\nFour-node (t, A, m, s): t cubic Kubota [5], t-A cubic edge [4 or 2], A: [0], A-m sextic [1], m: [0], m-s [5], s: [4]")
    for MtA in (2, 4):
        for MAm in (1, 5):
            show("(t,A,m,s) M_tA=%d M_Am=%d" % (MtA, MAm),
                 [[5, MtA, 0, 0], [MtA, 0, MAm, 0], [0, MAm, 0, 5], [0, 0, 5, 4]])

def objects(M, n=6, cap=4000):
    """List the distinct M-matrices (objects) reachable, with the basis that realises them and
    the nodes that are non-admissible (q_ii = 1, i.e. Moebius-type nodes, M_ii = n/2)."""
    r = len(M)
    B = [[0]*r for _ in range(r)]
    for i in range(r):
        B[i][i] = (M[i][i] + n // 2) % n
        for j in range(i + 1, r):
            B[i][j] = M[i][j] % n
    chi = lambda a, b: sum(a[i] * B[i][j] * b[j] for i in range(r) for j in range(r)) % n
    def Mof(E):
        return tuple(tuple(((chi(E[i], E[i]) - n // 2) % n) if i == j else
                           ((chi(E[i], E[j]) + chi(E[j], E[i])) % n) for j in range(r)) for i in range(r))
    def mij(E, i, j):
        qi = chi(E[i], E[i]); qt = (chi(E[i], E[j]) + chi(E[j], E[i])) % n
        for m in range(0, 4 * n + 1):
            if ((m + 1) * qi) % n == 0 or (m * qi + qt) % n == 0:
                return m
    E0 = tuple(tuple(int(i == k) for k in range(r)) for i in range(r))
    seen = {E0}; stack = [E0]; out = {}
    while stack and len(seen) < cap:
        E = stack.pop(); Mm = Mof(E)
        out.setdefault(Mm, E)
        for i in range(r):
            if chi(E[i], E[i]) == 0: continue
            newE = tuple(tuple(-x for x in E[i]) if j == i else
                         tuple(E[j][k] + mij(E, i, j) * E[i][k] for k in range(r)) for j in range(r))
            if newE not in seen:
                seen.add(newE); stack.append(newE)
    return out

def canon(M):
    """canonical form up to simultaneous permutation of nodes"""
    r = len(M)
    return min(tuple(tuple(M[p[i]][p[j]] for j in range(r)) for i in range(r)) for p in itertools.permutations(range(r)))
