"""
Paper2-style route (11/12 paper architecture) transplanted to Q(i) with the DERIVED Kubota geometry.

Family A_u(D) = sum_n mu(n) nu(n) chi_n(u) W(N n / D), chi quartic.  Rows N(u) <= H = D^kappa.
 (E1) extraction (rows u = p^4, N p ~ Y = H^{1/4}):  |A_1|^2 << D H^{3/4} + D^2 H^{-1/2}
      [sextic: u = p^6, |A_1|^2 << D H^{5/6} + D^2 H^{-1/3}]
 (E2) Poisson in u: dual rows N(h) <~ calH = D^2/H, dual sums B_h(D) with coefficients Phi*conj(alpha)*gt4(h^3, n)
      target  E(calH, X) := sum_h |X^{-1/2} B_h(X)|^2 <~ X  at X = D.
 (E3) completed mean square: reflection (Kubota-Diaconu, conductor N(h)^3, dual length calH^3/X, quartic dual
      twist) + a quartic large sieve  Sieve(rows, cols); or the unreflected sieve on length X.
 (E4) zeta(4s-3) completion removal + transfer recursion (two Poissons, ratio calH/X preserved) as in paper2;
      the completed condition is monotone along the recursion (calH' = calH L/X, X' = L < X).
Reachable kappa = smallest kappa with E <~ X; theta = max over extraction terms.
"""
import numpy as np
from fractions import Fraction as F

SIEVES = {
    "deFaveri n=4 (proved, 2610.04045)": [(1, 0), (0, 1), (F(3, 4), F(1, 2)), (F(1, 2), F(3, 4))],
    "BGL (AB)^{2/3} (proved)": [(1, 0), (0, 1), (F(2, 3), F(2, 3))],
    "perfect quartic sieve A+B (conjectural, false in general by dFDH)": [(1, 0), (0, 1)],
    "LoA in cube-shift aspect (rows only)": [(1, 0)],
}

def completed_ok(eta, pieces):
    """eta = log_X calH.  Reflected: cols = 3 eta - 1 ; direct: cols = 1.  Need max piece <= 1 (normalized)."""
    refl = max(a * eta + b * (3 * eta - 1) for a, b in pieces)
    dirc = max(a * eta + b * 1 for a, b in pieces)
    return min(refl, dirc) <= 1 + 1e-15

def eta_max(pieces):
    lo, hi = F(0), F(1)
    # exact: each piece gives a linear condition; scan candidates
    cands = []
    for a, b in pieces:
        # reflected: a*eta + b*(3 eta - 1) <= 1  ->  eta <= (1+b)/(a+3b)
        if a + 3 * b > 0: cands.append(F(1 + b) / (a + 3 * b))
        # direct: a*eta + b <= 1 -> eta <= (1-b)/a
        if a > 0: cands.append(F(1 - b) / a)
    best = F(0)
    for c in sorted(set(cands)) + [F(1)]:
        if 0 <= c <= 1 and completed_ok(float(c), pieces):
            best = max(best, c)
    return min(best, F(1))

def theta_extraction(kappa, n=4):
    """|A_1|^2 << D^{1 + (1-1/n) kappa} + D^{2 - 2 kappa/n}."""
    return max((1 + (1 - F(1, n)) * kappa) / 2, 1 - kappa / n)

if __name__ == "__main__":
    print("PAPER2 ROUTE, quartic, derived geometry")
    for name, pcs in SIEVES.items():
        e = eta_max(pcs)
        kappa = 2 - e                       # calH = D^2/H = D^{2-kappa} <= X^{eta} = D^{eta}
        th = theta_extraction(kappa)
        print(f"  {name:66s}: dual rows calH <= X^{e} (= {float(e):.5f}) -> H >= D^{kappa} (= {float(kappa):.5f});"
              f"  extraction exponent theta = {th} = {float(th):.5f}" + ("  (trivial, no zero-free half-plane)" if th >= 1 else ""))
    # what kappa would 6/7 need
    k67 = F(20, 21)
    print(f"  6/7 would need H = D^{k67} (theta_ext({k67}) = {theta_extraction(k67)}), i.e. calH = D^{2-k67} > X = D:"
          " a below-diagonal (rows > columns) dual mean square, i.e. a quasi-RH-strength input (R5 obstruction O1).")
    print(f"  best possible from extraction alone: kappa = 4/5 -> theta = {theta_extraction(F(4,5))} (needs H = D^(4/5) < D)")
    # sextic check (paper2 itself): conductor N(h)^2, quadratic sieve
    pcs = [(1, 0), (0, 1)]
    # cubic reflected cols = 2 eta - 1 -> eta <= 1 ; kappa -> 1+, theta = 11/12
    print(f"  sextic check: theta_ext(1, n=6) = {theta_extraction(F(1), 6)} (= 11/12, paper2)")
    # the combined bound (no exact mean square): exponent of |A_1|^2 with dual bound inserted
    print("  combined bound |A_1|^2 << H^{-1/4}[D H + H F] + D^2 H^{-1/2}, F = min(direct, reflected) sieve, dF:")
    best = (9, None)
    for k in np.linspace(1.0, 2.0, 100001):
        eta = 2 - k
        pcs = SIEVES["deFaveri n=4 (proved, 2610.04045)"]
        Fx = min(max(float(a) * eta + float(b) * (3 * eta - 1) for a, b in pcs), max(float(a) * eta + float(b) for a, b in pcs))
        ex = max(1 + 0.75 * k, 2 - k / 2, 0.75 * k + Fx)
        if ex < best[0]: best = (ex, k)
    print(f"    min over kappa of the |A_1|^2 exponent = {best[0]:.6f} at kappa = {best[1]:.5f}  (>= 2 means trivial)")
