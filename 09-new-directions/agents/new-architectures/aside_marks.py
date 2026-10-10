"""Variant (Idea 4): amplifier chi_p(m) on the Poisson variable, grouped with the Gauss-sum
polynomial in Cauchy-Schwarz:  |sum_m (A_m M_m) B_m| <= ||A M|| ||B||,  M_m = sum_p w_p chi_p(m).
Model: A M is a Gauss-sum polynomial of modulus s p (Y' = YP) on the fixed m-range Q, i.e.
X' = X/P; Gram  (Q/Y')(1 + (Y'/X')^{1/6} + (Y'/X')^2/Y');  B unmarked: wall M <= 1 only.
Signal: P^{2/3} (weights conj gamma_p), low: P^{1/2} x Gram'  => amplification P^{1/6} before Gram.
Rows as in P78 (range ZP/X, suppression P^{-1/2}).  GRH for rows."""
import numpy as np
from scipy.optimize import linprog

def solve(gram_slope=1/6, caseA=True, wallB=1.0, grh=True):
    A, B = [], []
    def ge(coef, const):
        A.append(list(coef) + [-1.0]); B.append(-const)
    def le(coef, const):
        A.append(list(coef) + [0.0]); B.append(-const)
    LX = np.array([1, 0, 0.]); LY = np.array([0, 1, 0.]); EL = np.array([0, 0, 1.])
    bp = LY - LX + 2 * EL                      # log(Y'/X')
    base = (LX / 6 - EL / 6, 5 / 6)
    gp = [np.zeros(3), gram_slope * bp, 2 * bp - LY - EL]
    M = LX + LY
    ep = [(np.zeros(3), 0.0), (M, -wallB), ((1/3) * M, -2/3)]
    if caseA:
        ep.append((-0.5 * M, 0.25))           # unmarked case A (1-2M)/4
    for g in gp:
        for e in ep:
            ge(base[0] + 0.5 * g + 0.5 * e[0], base[1] + 0.5 * e[1])
    ge(-0.5 * LY - 0.5 * EL + (5 / 6) * (-LX + EL), 0.5 + 5 / 6)
    le(-LX + EL, 0.0)                         # X' = X/P >= 1
    res = linprog([0, 0, 0, 1], A_ub=np.array(A), b_ub=np.array(B),
                  bounds=[(0, 2), (0, 3), (0, 2), (None, None)], method="highs")
    lx, ly, ell, sg = res.x
    return dict(theta=round(sg, 6), lx=round(lx, 4), ly=round(ly, 4), ell=round(ell, 4), M=round(lx + ly, 4))

if __name__ == "__main__":
    print("A-side amplifier, GRH:", solve())
    print("  ... with Gram slope halved:", solve(gram_slope=1/12))
    print("  ... no Gram at all:", solve(gram_slope=0))
