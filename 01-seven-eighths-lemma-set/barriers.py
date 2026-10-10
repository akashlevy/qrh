"""
Exact lower-bound certificates ("this method cannot go further") for the
OpenAI 7/8 quasi-Riemann argument (Part II of the 30 Sep 2026 preprint).

The argument needs, for a geometry (l_x, l_y, ell) and target sigma0:

  LOW  (Prop probe-low, Prop probe-gram, Lemma reflected-energy):
       sigma0 >= 5/6 + l_x/6 - ell/6 + g/2 + e/2
       g = max(0, b/6, 2b - l_y)           b = l_y - l_x   (additive Gram bound)
       e = max(0, (1+2ell-2M)/3, M+ell-1)  M = l_x + l_y
       NOTE: the true case-A term is the paper's (1+3ell-2M)/4 (powerful rows; see
       strengthening/energy_lp.py), which is >= the term used here near the optimum. Using the
       weaker constraint only lowers the LP value, so both certificates remain valid LOWER bounds.
  HIGH (eq. common-high-exponent at d = h, rows counted by Prop detector-counts):
       sigma0 >= (1+delta)(1-l_y)/2 - ell/2 + x*delta*ell + h*(R_*(delta,x) + delta/2 - 1/6)
       for every delta in (1/50, 5/6], x in [0, 1/2],   h = 1 - l_x + ell.

Each certificate below is a nonnegative rational combination of these linear
constraints that eliminates (l_x, l_y, ell, g, e) and leaves sigma0 >= const.
Because only *some* of the high-side constraints are used, the bound holds for
the full problem (a relaxation can only lower the optimum).

  (A) paper's ingredients as stated:                sigma0 >= 0.8749570...
  (B) same, but GRH for every row L-function
      (perfect row counting, detector floor -> 1/2):  sigma0 >= 13/15

Run: python3 barriers.py   (pure stdlib)
"""
from fractions import Fraction as F
from itertools import combinations

ALPHA, CP = F(5, 6), F(2, 9)

def Rstar(d, x):
    D = 3 - x - 4 * CP * x
    P = (1 - x) * (2 - 4 * CP * x)
    J = (ALPHA - d) * D + d * P
    return 1 - d + (ALPHA - d) * d * P / (2 * J)

# constraint rows:  a . (l_x, l_y, ell, sigma, g, e) >= rhs
LOW = [
    ([F(-1, 6), 0, F(1, 6), 1, F(-1, 2), F(-1, 2)], F(5, 6), "low side"),
    ([0, 0, 0, 0, 1, 0], F(0), "g >= 0"),
    ([F(1, 6), F(-1, 6), 0, 0, 1, 0], F(0), "g >= b/6"),
    ([2, -1, 0, 0, 1, 0], F(0), "g >= 2b - l_y"),
    ([0, 0, 0, 0, 0, 1], F(0), "e >= 0"),
    ([F(2, 3), F(2, 3), F(-2, 3), 0, 0, 1], F(1, 3), "e >= (1+2ell-2M)/3"),
    ([-1, -1, -1, 0, 0, 1], F(-1), "e >= M+ell-1"),
]

def high_row(d, x):
    c = Rstar(d, x) + d / 2 - F(1, 6)
    return ([c, (1 + d) / 2, F(1, 2) - x * d - c, 1, 0, 0], (1 + d) / 2 + c,
            f"high side at delta={d}, x={x}")

# GRH for rows: delta -> 0, R = 1, q = 0  =>  sigma >= (1-l_y)/2 - ell/2 + 5h/6
HIGH_GRH = ([F(5, 6), F(1, 2), F(1, 2) - F(5, 6), 1, 0, 0], F(1, 2) + F(5, 6),
            "high side, GRH for all rows")

def best_dual(rows):
    """max over basic dual-feasible y >= 0 with sum y_i a_i = e_sigma of sum y_i rhs_i."""
    best = None
    for S in combinations(range(len(rows)), 6):
        A = [[rows[i][0][k] for i in S] + [F(int(k == 3))] for k in range(6)]
        ok = True
        for col in range(6):
            piv = next((r for r in range(col, 6) if A[r][col] != 0), None)
            if piv is None:
                ok = False
                break
            A[col], A[piv] = A[piv], A[col]
            for r in range(6):
                if r != col and A[r][col] != 0:
                    f = A[r][col] / A[col][col]
                    A[r] = [a - f * b for a, b in zip(A[r], A[col])]
        if not ok:
            continue
        y = [A[k][6] / A[k][k] for k in range(6)]
        if min(y) < 0:
            continue
        val = sum(yi * rows[i][1] for yi, i in zip(y, S))
        if best is None or val > best[0]:
            best = (val, [(rows[i][2], yi) for yi, i in zip(y, S) if yi != 0])
    return best

def verify(cert_val, cert, rows):
    """independent check: combination reproduces exactly sigma >= cert_val."""
    byname = {r[2]: r for r in rows}
    coeff = [F(0)] * 6
    rhs = F(0)
    for name, y in cert:
        a, b, _ = byname[name]
        assert y >= 0
        coeff = [c + y * ai for c, ai in zip(coeff, a)]
        rhs += y * b
    assert coeff == [0, 0, 0, 1, 0, 0] and rhs == cert_val
    return True

if __name__ == "__main__":
    pts = [(F(387, 1000), F(1, 2)), (F(39, 100), F(1, 2))]
    rowsA = LOW + [high_row(d, x) for d, x in pts]
    vA, cA = best_dual(rowsA)
    verify(vA, cA, rowsA)
    print(f"(A) paper's ingredients:  sigma0 >= {vA} = {float(vA):.7f}")
    for n, y in cA:
        print(f"      {float(y):.5f} x [{n}]")
    rowsB = LOW + [HIGH_GRH]
    vB, cB = best_dual(rowsB)
    verify(vB, cB, rowsB)
    print(f"(B) GRH for all rows:     sigma0 >= {vB} = {float(vB):.7f}")
    for n, y in cB:
        print(f"      {y} x [{n}]")
    print("both certificates verified exactly")
