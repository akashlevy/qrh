"""
P2-shape (paper2 architecture) bookkeeping for an order-N Mobius family whose Poisson dual
is reflected by an automorphic functional equation of conductor exponent cB in the row h,
with the order-k large sieve (coefficient-blind) as the only base inequality.

Family (type O):  A_u = sum_{n~D} mu(n) nu(n) chi_n(u),  rows N(u) <= H = D^h.
Copies u = a^N  ->  |A_1|^2 <= (mean square) / H^{1/N}.
If  sum_{u<=H} |A_u|^2 << D^{1+eps} H  then  theta = 1/2 + (1 - 1/N) h / 2.

Moves on normalised exponents  r(F,L) = log(Energy/(F L))  (logs base D):
  Poisson     O(f,l) <-> D(2l-f, l):   rho_O(f,l) <= max(0, (l-f) + rho_D(2l-f,l)),  and symmetric.
  Reflection  D(f,l) -> D(f, cB f - l):  rho_D unchanged (exact identity).
Base:  order-k large sieve  F + L + (F L)^{1-1/k} ... for k=3 (cubic) F+L+(FL)^{2/3};
       for k=2 (quadratic, Heath-Brown) F+L.
Optionally a 'pole' term for duals whose Dirichlet series carries a normaliser 1/L(N s - c, chi^N):
       rho_D(f,l) <= max(0, (2 beta* - 1)/N * l)   (degree-one dual, all rows; see NOTES.md sec 3).
"""
from fractions import Fraction as Fr
import numpy as np


def ls_exp(f, l, k):
    """excess exponent of the order-k large sieve, family f, length l (normalised by F*L)."""
    if k == 2:
        return max(0.0, l - f)
    # Blomer-Goldmakher-Louvel: F + L + (FL)^{2/3} for every k >= 3 (sharp for k = 3, DFDH)
    return max(0.0, l - f, (2.0 / 3.0) * (f + l) - f)


def rho_D(f, l, cB, k, depth=6):
    """best bound for a type-D node (automorphic dual coefficients), using up to `depth` moves."""
    best = ls_exp(f, l, k)
    if depth == 0:
        return best
    # reflection then base / Poisson
    lr = cB * f - l
    if lr >= 0:
        best = min(best, ls_exp(f, lr, k))
        best = min(best, max(0.0, (lr - f) + rho_O(2 * lr - f, lr, cB, k, depth - 1)))
    # Poisson
    best = min(best, max(0.0, (l - f) + ls_exp(2 * l - f, l, k)))
    return best


def rho_O(f, l, cB, k, depth=6):
    best = ls_exp(f, l, k)
    if depth == 0:
        return best
    best = min(best, max(0.0, (l - f) + rho_D(2 * l - f, l, cB, k, depth - 1)))
    return best


def theta_P2(N, cB, k, grid=4001):
    """smallest h with rho_O(h,1) == 0 -> theta."""
    for h in np.linspace(1.0, 3.0, grid):
        if rho_O(h, 1.0, cB, k) <= 1e-12:
            return h, 0.5 + (1 - 1 / N) * h / 2
    return None, None


def theta_P2_pole(N, beta):
    """degree-one dual with harmful normaliser (A2-residue picture): bootstrap step."""
    # sum_h |B_h|^2 <~ Hd * D^{1 + (2 beta - 1)/N}; balance gives
    h = 1 + (2 * beta - 1) / N
    return h, 0.5 + (1 - 1 / N) * h / 2


if __name__ == "__main__":
    print("P2-shape, order-N family, dual reflected with conductor exponent cB, base = order-k LS")
    rows = [
        ("sextic -> cubic theta, quadratic residual (paper2)", 6, 2, 2),
        ("sextic, no reflection (BGL sextic LS only)", 6, 0, 6),
        ("cubic, no reflection (cubic LS only)", 3, 0, 3),
        ("cubic -> RS theta(x)theta, cB=2 (harmless normaliser)", 3, 2, 3),
        ("cubic -> RS, cB=3", 3, 3, 3),
        ("cubic -> sextic Kubota h^4, cB=4", 3, 4, 3),
        ("cubic, cB=2 with a PERFECT residual LS (hypothetical)", 3, 2, 2),
        ("cubic, cB=1 with perfect LS (hypothetical finite type)", 3, 1, 2),
    ]
    for lab, N, cB, k in rows:
        h, th = theta_P2(N, cB, k)
        print(f"  {lab:60s} h={h:.4f}  theta={th:.6f}  ~{Fr(th).limit_denominator(200)}")
    print()
    print("A2-residue (degree one, normaliser 1/L(3s-1,(alpha nu)^3)) bootstrap, N=3:")
    b = 1.0
    for it in range(12):
        h, th = theta_P2_pole(3, b)
        print(f"   beta*(chi^3)={b:.6f} -> beta(nu) <= {th:.6f}")
        b = th
    print("   fixed point 13/14 =", 13 / 14)
    print("   general N fixed point beta = 1/2 + N(N-1)/(2(N^2-N+1)):",
          {N: 0.5 + N * (N - 1) / (2 * (N * N - N + 1)) for N in (3, 4, 5, 6)})
