"""
Exact rational certificate: the OpenAI Part II argument (preprint
"The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane Re(s) > 7/8",
openai/math, 30 Sep 2026) closes at

        sigma0 = 7/8 - 1/40000 = 34999/40000

when its fixed geometry (old-eq:5.1) is replaced by

        b = 1/8,  ell = 1/6 + 1/10000,  M + ell = 1,
        l_x = (1-ell-b)/2,  l_y = l_x + b,  h = 1 - l_x + ell.

Every geometry-dependent inequality used in the paper's Part II
(Sections "compensated probe", "compensated low estimate",
"local compensation", "seven-eighths bound") is re-checked here in exact
rational arithmetic.  The single non-linear check, the analogue of the
paper's Lemma "Compensated endpoint certificate" (lem:balanced-endpoint),
is proved by Bernstein expansion with exact subdivision.

Run:  python3 certificate.py      (pure stdlib, no numpy/sympy needed)
"""
from fractions import Fraction as F
from itertools import product

# ---------------------------------------------------------------- geometry
b     = F(1, 8)
ell   = F(1, 6) + F(1, 10000)
M     = 1 - ell
lx    = (M - b) / 2
ly    = lx + b
h     = 1 - lx + ell
sigma0 = F(7, 8) - F(1, 40000)

ALPHA = F(5, 6)      # sixth-power amplification slope (Lemma inverse-amplification)
CP    = F(2, 9)      # plain capacity 1/(6 kappa) at kappa = 3/4 (Lemma plain)
Z0    = F(17, 50)    # z-contour (Euler region one)
D0    = F(1, 50)     # detector floor delta_0 = 2*(51/100) - 1
DMIN  = F(1, 100)    # small-row cutoff, as in the paper

checks = []
def check(name, cond, val=None):
    checks.append((name, bool(cond), val))

# ------------------------------------------------- Mellin exponent / low side
C = lambda s: s + lx / 2 - 1 + h / 6                    # Definition stage-high-data
low_exponent = lx / 2 + b / 12                          # Prop probe-low: X^{1/2} P_a^{1/12}
check("C(s) = s - 11/16 (unchanged)", C(0) == F(-11, 16), C(0))
check("low exponent == C(sigma0)", low_exponent == C(sigma0), low_exponent)
check("M + ell = 1 (energy: max(H, v+l_b) = H)", M + ell == 1)
# Prop probe-low length conditions, for all rescaled lengths 0 <= d <= ell (worst d = ell)
check("Q' = Z^{M-2d} >= 1        : M - 2 ell > 0", M - 2 * ell > 0, M - 2 * ell)
check("X' >= 1                    : l_x - ell > 0", lx - ell > 0, lx - ell)
check("Y' >= 1                    : l_y - ell > 0", ly - ell > 0, ly - ell)
check("Gram: P_a^2/Y' <= P_a^{1/6}: l_y - ell - 11b/6 > 0",
      ly - ell - 11 * b / 6 > 0, ly - ell - 11 * b / 6)
check("P_a = Z^b >= 1             : b > 0", b > 0)
# Lemma probe-row-norm: energy Z^{M' + ((1 + 3 l' - 2 M')/4)_+} with l' = ell - d,
# M' = M - 2d; tuple exponent f(d) = -d + (energy excess)/2 must be <= 0.
def f_tuple(d):
    excess = max(F(0), (1 + 3 * (ell - d) - 2 * (M - 2 * d)) / 4)
    return -d + excess / 2
grid = [ell * F(k, 1000) for k in range(1001)]
worst = max(f_tuple(d) for d in grid)
# f is piecewise linear in d with one kink, so the grid (which contains the kink
# region endpoints 0 and ell) gives the exact maximum; record the kink explicitly too.
kink = 1 - 5 * ell
check("rescaled tuples: max_d f(d) <= 0",
      worst <= 0 and (kink < 0 or kink > ell or f_tuple(kink) <= 0), worst)

# ----------------------------------------------------- principal / outer rows
check("principal margin m_w = l_y/20 > 0", ly / 20 > 0, ly / 20)
check("principal margin m_z = h/600 > 0", h / 600 > 0, h / 600)
small = h * (Z0 - F(1, 6)) - ly / 2 + F(63, 50) * DMIN
check("small rows: h(17/50-1/6) - l_y/2 + (63/50) d_min < 0", small < 0, small)
# Euler region two extended down to Re s >= sigma0: the four principal error
# exponents  -x, -6z, 4-5x-6z, 1-w-6z  at x>=sigma0, z>=33/200, w>=19/20
xr, zr, wr = sigma0, F(33, 200), F(19, 20)
errs = [-xr, -6 * zr, 4 - 5 * xr - 6 * zr, 1 - wr - 6 * zr]
check("principal local errors <= -sigma0 on extended region two",
      max(errs) <= -sigma0, max(errs))
check("region-two good-prime defect 4-6x-6z < -1", 4 - 6 * xr - 6 * zr < -1)
check("region-two good-prime defect 1-x-w-6z < -1", 1 - xr - wr - 6 * zr < -1)
check("region one at small-row line: beta*+1/2 > 1 (eps0>0)", sigma0 + F(1, 2) > 1)

# --------------------------------------------------------- high side pieces
def E_high(delta, q, R, d):
    """eq. (common-high-exponent), first line, relative to C(sigma0)."""
    a = (1 + delta) / 2
    return (a - sigma0 + h * (Z0 - F(1, 6)) - a * ly - (1 - a) * ell
            - (delta / 2 - q) * ell + d * (R + delta / 2 - Z0))

# supply: total prime length must exceed the maximal capacity 7/37
check("prime supply at d=h: ell/h > 7/37", ell / h > F(7, 37), ell / h - F(7, 37))
check("supply above h: 5*ell - h > 0 (room for zeta)", 5 * ell - h > 0, 5 * ell - h)
# floor bin: R = 1, q <= delta0/2
Efloor = E_high(D0, D0 / 2, F(1), h)
check("floor bin E(h) < 0", Efloor < 0, Efloor)
# intermediate rows d in [d_min, 1/2]: R = 76/75 - 2 delta/3, q = delta/2, linear in delta
Emid = [E_high(dl, dl / 2, F(76, 75) - 2 * dl / 3, F(1, 2)) for dl in (D0, ALPHA)]
check("intermediate rows E(1/2) < 0 at delta = 1/50 and 5/6", max(Emid) < 0, max(Emid))
check("intermediate slope 101/150 - delta/6 > 0", F(101, 150) - ALPHA / 6 > 0)

# ------------------------------------- the endpoint certificate (non-linear)
# R_*(delta, x) from Prop detector-counts after balancing t:
#   D_x = 3 - x - 4 cP x,  P_x = (1-x)(2 - 4 cP x),  J = (alpha-delta) D_x + delta P_x
#   R_* = 1 - delta + (alpha-delta) delta P_x / (2 J)
# -E(h) = sigma0 - [(1+delta)(1-l_y)/2 - ell/2 + x delta ell + h (R_* + delta/2 - 1/6)]
# Multiply by J > 0:  N(delta,x) = J*(-E) - m*J  must be >= 0 on the box.
MARGIN = F(9, 100000)   # claimed uniform margin -E(h) >= 9e-5

# bivariate polynomials as dict {(i,j): coeff}  for delta^i x^j
def padd(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}
def pmul(p, q):
    out = {}
    for (i1, j1), a in p.items():
        for (i2, j2), c in q.items():
            k = (i1 + i2, j1 + j2)
            out[k] = out.get(k, 0) + a * c
    return {k: v for k, v in out.items() if v != 0}
def const(c): return {(0, 0): F(c)}
def scal(c, p): return {k: F(c) * v for k, v in p.items()}
DEL = {(1, 0): F(1)}
X = {(0, 1): F(1)}

Dx = padd(const(3), scal(-(1 + 4 * CP), X))
Px = pmul(padd(const(1), scal(-1, X)), padd(const(2), scal(-4 * CP, X)))
J  = padd(pmul(padd(const(ALPHA), scal(-1, DEL)), Dx), pmul(DEL, Px))
# bracket = sigma0 - (1+d)(1-ly)/2 + ell/2 - x d ell - h (1 - d/2 - 1/6)   [R_* = 1-d+...]
bracket = padd(const(sigma0 - (1 - ly) / 2 + ell / 2 - h * F(5, 6)),
               scal(-(1 - ly) / 2 + h / 2, DEL),
               scal(-ell, pmul(DEL, X)))
N = padd(pmul(J, bracket),
         scal(-h / 2, pmul(pmul(padd(const(ALPHA), scal(-1, DEL)), DEL), Px)),
         scal(-MARGIN, J))

from math import comb
def bernstein_min(poly, d0, d1, x0, x1):
    """min Bernstein coefficient of poly on [d0,d1]x[x0,x1] (exact)."""
    n = max(i for i, _ in poly) if poly else 0
    m = max(j for _, j in poly) if poly else 0
    # substitute delta = d0 + (d1-d0) s, x = x0 + (x1-x0) t, expand in s,t
    c = {}
    for (i, j), a in poly.items():
        for k in range(i + 1):
            for l in range(j + 1):
                term = (a * comb(i, k) * d0 ** (i - k) * (d1 - d0) ** k
                          * comb(j, l) * x0 ** (j - l) * (x1 - x0) ** l)
                c[(k, l)] = c.get((k, l), 0) + term
    best = None
    for p_, q_ in product(range(n + 1), range(m + 1)):
        val = sum(c.get((k, l), 0) * F(comb(p_, k), comb(n, k)) * F(comb(q_, l), comb(m, l))
                  for k in range(p_ + 1) for l in range(q_ + 1))
        best = val if best is None or val < best else best
    return best

def certify(poly, box, depth=0, maxdepth=14, stats=[0]):
    d0, d1, x0, x1 = box
    stats[0] += 1
    if bernstein_min(poly, d0, d1, x0, x1) >= 0:
        return True
    if depth >= maxdepth:
        return False
    if d1 - d0 >= x1 - x0:
        mid = (d0 + d1) / 2
        return (certify(poly, (d0, mid, x0, x1), depth + 1, maxdepth, stats)
                and certify(poly, (mid, d1, x0, x1), depth + 1, maxdepth, stats))
    mid = (x0 + x1) / 2
    return (certify(poly, (d0, d1, x0, mid), depth + 1, maxdepth, stats)
            and certify(poly, (d0, d1, mid, x1), depth + 1, maxdepth, stats))

stats = [0]
ok_endpoint = certify(N, (F(0), ALPHA, F(0), F(1, 2)), stats=stats)
# J > 0 on the box (needed to divide): J is linear in delta, check corners/edges exactly
Jpos = bernstein_min(J, F(0), ALPHA, F(0), F(1, 2)) > 0
check("J(delta,x) > 0 on [0,5/6]x[0,1/2]", Jpos)
check(f"endpoint certificate: -E(h) >= {MARGIN} on [0,5/6]x[0,1/2]  "
      f"({stats[0]} boxes)", ok_endpoint)
# high-side margin vs Delta: E_actual - Delta <= -m - (1 - h/4) Delta  needs 1 - h/4 > 0
check("Delta absorption: 1 - h/4 > 0", 1 - h / 4 > 0, 1 - h / 4)

# ------------------------------------------------------------------ report
if __name__ == "__main__":
    print(f"sigma0 = {sigma0} = {float(sigma0):.7f}")
    print(f"b={b}  ell={ell}  M={M}  l_x={lx}  l_y={ly}  h={h}")
    print()
    allok = True
    for name, ok, val in checks:
        allok &= ok
        v = "" if val is None else f"   [{val} ~ {float(val):.6g}]"
        print(("PASS " if ok else "FAIL ") + name + v)
    print()
    print("ALL CHECKS PASS" if allok else "SOME CHECK FAILED")
