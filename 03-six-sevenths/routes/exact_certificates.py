"""Exact (Fraction) LP-dual certificates for route R2 (bootstrap).  stdlib only.
Variables v = (lx, ly, ell, s, g, e).  A constraint is (coeffs, rhs, name) meaning coeffs.v >= rhs.
A certificate is a list of (name, y>=0); it is valid iff sum y*coeffs == (0,0,0,1,0,0) and then s >= sum y*rhs.
"""
from fractions import Fraction as F
from itertools import combinations

def low():      # s >= 5/6 + lx/6 - ell/6 + g/2 + e/2           (Prop probe-low with Gram and energy factors)
    return ([F(-1,6),0,F(1,6),1,F(-1,2),F(-1,2)], F(5,6), "low")
def g0():  return ([0,0,0,0,1,0], F(0), "g>=0")
def gb6(): return ([F(1,6),F(-1,6),0,0,1,0], F(0), "g>=b/6  (exceptional sixth-power frequencies, P_a^{1/6})")
def g2b(): return ([2,-1,0,0,1,0], F(0), "g>=2b-ly (nonexceptional, P_a^2/Y')")
def e0():  return ([0,0,0,0,0,1], F(0), "e>=0")
def eA():  return ([F(1,2),F(1,2),F(-3,4),0,0,1], F(1,4), "e>=(1+3ell-2M)/4 (case A)")
def eW(w=F(0)): return ([-1,-1,-1,0,0,1], -1-w, f"e>=M+ell-1-{w} (energy wall)")
def floor_bin(d0):
    # E(h) < 0 for the floor bin: s >= (1+d0)(1-ly)/2 - ell/2 + (d0/2) ell + h(1+d0/2-1/6),  h = 1-lx+ell
    c = 1 + d0/2 - F(1,6)
    return ([c, (1+d0)/2, F(1,2) - d0/2 - c, 1, 0, 0], (1+d0)/2 + c, f"floor bin (R=1) at delta0={d0}")

def verify(rows, cert):
    by = {r[2]: r for r in rows}
    co = [F(0)]*6; rhs = F(0)
    for name, y in cert:
        assert y >= 0, (name, y)
        c, b, _ = by[name]
        co = [a + y*t for a, t in zip(co, c)]; rhs += y*b
    assert co == [0,0,0,1,0,0], co
    return rhs

def best_vertex(rows):
    """max over basic dual-feasible y of sum y*rhs  (= LP optimum by duality; tiny LP, brute force)."""
    best = None
    for S in combinations(range(len(rows)), 6):
        A = [[F(rows[i][0][k]) for i in S] + [F(int(k == 3))] for k in range(6)]
        ok = True
        for col in range(6):
            piv = next((r for r in range(col, 6) if A[r][col] != 0), None)
            if piv is None: ok = False; break
            A[col], A[piv] = A[piv], A[col]
            for r in range(6):
                if r != col and A[r][col] != 0:
                    f = A[r][col]/A[col][col]
                    A[r] = [a - f*b for a, b in zip(A[r], A[col])]
        if not ok: continue
        y = [A[k][6]/A[k][k] for k in range(6)]
        if min(y) < 0: continue
        val = sum(yi*rows[i][1] for yi, i in zip(y, S))
        if best is None or val > best[0]:
            best = (val, [(rows[i][2], yi) for yi, i in zip(y, S) if yi != 0])
    return best

if __name__ == "__main__":
    base = [low(), g0(), gb6(), g2b(), e0(), eA(), eW()]
    # (1) GRH for all non-principal rows (floor at 1/2): 13/15
    rows = base + [floor_bin(F(0))]
    cert = [("low", F(4,5)), ("floor bin (R=1) at delta0=0", F(1,5)),
            ("g>=b/6  (exceptional sixth-power frequencies, P_a^{1/6})", F(1,5)), ("g>=0", F(1,5)),
            ("e>=M+ell-1-0 (energy wall)", F(1,15)), ("e>=0", F(1,3))]
    v = verify(rows, cert); print("(1) GRH rows, paper low side:           s >=", v, "=", float(v)); assert v == F(13,15)
    print("    LP optimum (vertex search):", best_vertex(rows)[0])
    # (2) actual detector floor 51/100 (delta0 = 1/50), any row counts: floor bin alone
    rows = base + [floor_bin(F(1,50))]
    bv = best_vertex(rows); print("(2) floor 51/100, rows otherwise perfect: s >=", bv[0], "=", float(bv[0]))
    for n, y in bv[1]: print("      ", y, "x", n)
    # (3) GRH rows AND Gram without the structural P_a^{1/6} term: exactly 6/7
    rows = [low(), g0(), g2b(), e0(), eA(), eW(), floor_bin(F(0))]
    cert = [("low", F(17,21)), ("floor bin (R=1) at delta0=0", F(4,21)),
            ("g>=2b-ly (nonexceptional, P_a^2/Y')", F(1,42)), ("g>=0", F(17,42)-F(1,42)),
            ("e>=M+ell-1-0 (energy wall)", F(1,14)), ("e>=0", F(17,42)-F(1,14))]
    v = verify(rows, cert); print("(3) GRH rows + Gram w/o P_a^{1/6}:      s >=", v, "=", float(v)); assert v == F(6,7)
    print("    LP optimum (vertex search):", best_vertex(rows)[0], " attained at (lx,ly,ell)=(2/7,4/7,1/7)")
    # check attainment of (3)
    lx, ly, ell = F(2,7), F(4,7), F(1,7); b = ly-lx; M = lx+ly; h = 1-lx+ell
    g = max(F(0), 2*b-ly); e = max(F(0), (1+3*ell-2*M)/4, M+ell-1)
    lowv = F(5,6)+lx/6-ell/6+g/2+e/2; hiv = (1-ly)/2-ell/2+5*h/6
    print("    at (2/7,4/7,1/7): low =", lowv, " GRH high =", hiv)
    # (4) GRH rows + energy wall relaxed by w (W(w)): 13/15 - w/... 
    for w in [F(1,10), F(1,5), F(1,3)]:
        rows = [low(), g0(), gb6(), g2b(), e0(), eA(), eW(w), floor_bin(F(0))]
        print(f"(4) GRH rows + wall relaxed by w={w}:     s >=", best_vertex(rows)[0], "=", float(best_vertex(rows)[0]))
