"""Exponent-level consistency of the counterfactual configuration W_beta with every lemma conclusion
a bootstrap argument could use.  W_beta: one primitive target eta0 has a zero at Re s = beta, every other
primitive finite-order Hecke L-function over F satisfies GRH.  Only the principal sixth-power class
u = a^6 (count H^{1/6}) carries |A_u(D)|^2 ~ D^{2 beta}; all other rows are diagonal-size.
Exact (Fraction) checks; prints the largest beta for which each conclusion is consistent.
"""
from fractions import Fraction as F
from itertools import product
def thm_ms(theta):      # sum_{N(u)<=D^{1+theta}} |A_u|^2 << D^{2+theta}:  (1+theta)/6 + 2b <= 2+theta
    return (2 + theta - (1 + theta) / 6) / 2
def amplification(c):   # |M_1(D)|^2 = D^{2b-1} << D^{5(1+c)/6}
    return (1 + F(5, 6) * (1 + c)) / 2
def marked():
    """Lemma marked: sum_{q_u<<Z^m}|M_u(Z^r)Q_u|^2 << Z^m for r+2z<=m, 2r+8z<=3m.  Principal class with
    aligned slot coefficients: Z^{m/6} rows, |M|^2=Z^{r(2b-1)}, |Q|^2=Z^z.  Need r(2b-1)+z <= 5m/6 on the
    polytope; take m=1 (homogeneous).  Return sup beta allowed."""
    verts = [(F(1), F(0)), (F(1, 2), F(1, 4)), (F(0), F(3, 8)), (F(0), F(0))]
    best = F(10)
    for r, z in verts:
        if r > 0:   # r(2b-1) + z <= 5/6  ->  b <= (5/6 - z + r)/(2r)
            best = min(best, (F(5, 6) - z + r) / (2 * r))
        else:
            assert z <= F(5, 6)
    return best
def extraction(kappa):  # integer/bootstrapped extraction from M_kappa:  beta <= (6+5kappa)/12 is the consistent range
    return (6 + 5 * kappa) / 12
def dh_world_vs_Mkappa(beta, kappa):
    """DH-world: U^{2-2beta} sixth-power-free rows ~U with |A|^2=D^{2beta}.  Exponent (base D) of
    sum_{u<=D^kappa}|A_u(D)|^2 minus (1+kappa); >0 means M_kappa FAILS in the DH-world."""
    return max(kappa * (2 - 2 * beta) + 2 * beta, kappa / 6 + 2 * beta, 1 + kappa) - (1 + kappa)
print("thm:ms (paper2, H=D^{1+theta}) consistent with W_beta iff beta <=", thm_ms(F(0)), "+ 5theta/12")
print("Lemma inverse-amplification at U~1 consistent iff beta <=", amplification(F(0)), "+ 5c/12")
print("Lemma marked (all admissible m,r,z) consistent iff beta <=", marked())
for k in [F(6, 7), F(9, 10), F(1)]:
    print(f"extraction from M_kappa, kappa={k}: W_beta excluded only for beta > {extraction(k)}")
print("\nDH-world: excess exponent of S(D^kappa,D) over D^{1+kappa}  (>0  => M_kappa false there)")
for b in [F(6, 7) + F(1, 100), F(13, 15), F(7, 8)]:
    print("  beta=", b, [ (str(k), str(dh_world_vs_Mkappa(b, k))) for k in [F(6, 7), F(9, 10), F(1), F(11, 10)] ])
