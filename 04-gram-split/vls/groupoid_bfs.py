"""BFS version: shortest reflection words to each chamber; residue-adapted chambers with their shortest
words, the affine maps of the free variables on the double-residue plane P = {u_t = 1/3, u_s = 1/6},
and the target objects.  Centred coordinates: u_i = s_i - 1/2; Kubota pole of an order-N node at u = 1/N."""
import sys
from collections import deque
from fractions import Fraction as Fr
sys.path.insert(0, '/Users/akashlevy/Documents/riemann/gram_split/vls')
from groupoid_search import setup, Mof, order_of, M0, NAMES, n
import sympy as sp

def bfs(M, cap=5000, maxcoef=30):
    r, chi = setup(M)
    def mij(E, i, j):
        qi = chi(E[i], E[i]); qt = (chi(E[i], E[j]) + chi(E[j], E[i])) % n
        for m in range(0, 4 * n + 1):
            if ((m + 1) * qi) % n == 0 or (m * qi + qt) % n == 0:
                return m
    E0 = tuple(tuple(int(i == k) for k in range(r)) for i in range(r))
    dist = {E0: []}; dq = deque([E0])
    while dq:
        E = dq.popleft()
        for i in range(r):
            if chi(E[i], E[i]) == 0: continue
            newE = tuple(tuple(-x for x in E[i]) if j == i else
                         tuple(E[j][k] + mij(E, i, j) * E[i][k] for k in range(r)) for j in range(r))
            if max(abs(x) for v in newE for x in v) > maxcoef:
                continue
            if newE not in dist:
                dist[newE] = dist[E] + [i]; dq.append(newE)
                if len(dist) > cap:
                    return dist, chi, mij
    return dist, chi, mij

if __name__ == "__main__":
    dist, chi, mij = bfs(M0)
    print("chambers:", len(dist), " max word length:", max(len(w) for w in dist.values()))
    et = (1, 0, 0, 0); es = (0, 0, 0, 1)
    neg = lambda v: tuple(-x for x in v)
    uA, um = sp.symbols('uA um')
    pole = {'t': sp.Rational(1, 3), 's': sp.Rational(1, 6)}
    out = []
    for E, w in dist.items():
        kt = [k for k in range(4) if E[k] in (et, neg(et))]
        ks = [k for k in range(4) if E[k] in (es, neg(es))]
        if not (kt and ks): continue
        kt, ks = kt[0], ks[0]
        free = [k for k in range(4) if k not in (kt, ks)]
        Mm = Mof(E, chi)
        Em = sp.Matrix(E)
        u = sp.Matrix([pole['t'], uA, um, pole['s']])
        up = Em * u
        l0 = sp.Matrix([0, 1, sp.Symbol('M'), 0])
        lp = (Em.T).inv() * l0
        out.append((len(w), ''.join(NAMES[i] for i in w), E, kt, ks, free, Mm,
                    [sp.expand(up[k]) for k in free], [sp.simplify(lp[k]) for k in free],
                    'sgn_t=%+d sgn_s=%+d' % (sum(E[kt]), sum(E[ks]))))
    out.sort()
    for o in out:
        L, w, E, kt, ks, free, Mm, ufree, lfree, sg = o
        print("%2d %-28s kt=%d ks=%d free=%s %s" % (L, w or 'id', kt, ks, free, sg))
        print("     free vars u' =", ufree, "  lengths l' =", lfree)
        print("     target M =", [list(r) for r in Mm], " numerators: t->", {k: mij(E, kt, k) for k in range(4) if k != kt},
              " s->", {k: mij(E, ks, k) for k in range(4) if k != ks})
