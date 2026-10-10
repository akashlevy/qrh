"""Pure-Fraction verification of the two exact 6/7 certificates in the GRH-for-rows regime of Part II.
Each certificate: sum_i y_i * (coef_i . v - rhs_i) >= 0 with y_i >= 0 and sum_i y_i coef_i = (0,0,0,1)
  ==> s >= sum_i y_i rhs_i = 6/7 for every geometry (lx, ly, ell)."""
from fractions import Fraction as F
def row(name):
    # returns (coef over (lx,ly,ell,s), rhs) for 'coef.v >= rhs'
    # low at d=0: s >= 5/6 + lx/6 - ell/6 + g/2 + e/2
    lowbase = ([-F(1,6), F(0), F(1,6), F(1)], F(5,6))
    def low(gc, gr, ec, er):  # g = gc.v + gr ; e = ec.v + er
        c = [lowbase[0][k] - gc[k]/2 - ec[k]/2 for k in range(4)]
        return (c, lowbase[1] + gr/2 + er/2)
    Z = [F(0)]*4
    g = {"g0": (Z, F(0)), "gb6": ([-F(1,6), F(1,6), F(0), F(0)], F(0)), "gq": ([F(-2), F(1), F(0), F(0)], F(0))}
    W = {"wall0": ([F(1), F(1), F(1), F(0)], F(-1)), "wall17": ([F(1), F(1), F(1), F(0)], F(-1) - F(1,7))}
    e = {"e0": (Z, F(0)), **W}
    if name == "high GRH":   # s >= 1/2 - ly/2 - ell/2 + 5/6 (1 - lx + ell)
        return ([F(5,6), F(1,2), F(1,2) - F(5,6), F(1)], F(1,2) + F(5,6))
    gn, en = name.split("+")
    return low(*g[gn], *e[en])
def check(cert):
    tot = [F(0)]*4; rhs = F(0)
    for name, y in cert:
        c, r = row(name); assert y >= 0
        tot = [t + y*ck for t, ck in zip(tot, c)]; rhs += y*r
    return tot, rhs
A = [("g0+e0", F(2,3)), ("g0+wall0", F(2,21)), ("gq+wall0", F(1,21)), ("high GRH", F(4,21))]
B = [("g0+e0", F(4,15)), ("g0+wall17", F(2,15)), ("gb6+e0", F(2,5)), ("high GRH", F(1,5))]
for nm, c in [("(a) GRH + Gram without P_a^{1/6} (no b/6 piece), wall M+ell<=1", A),
              ("(b) GRH + wall moved to M+ell<=8/7 (case A refined, inactive)", B)]:
    tot, rhs = check(c)
    print(nm, "\n    combination coefficients on (lx,ly,ell,s):", tot, " ==> s >=", rhs)
