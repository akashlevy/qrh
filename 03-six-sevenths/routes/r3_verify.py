"""
R3 (zero-density route) -- exact-arithmetic checks.  Pure stdlib (fractions).

Sections
  1. Dictionary theta <-> multiplicity-counted density exponent A_eff:  theta = 1 - 1/(6 A_eff).
  2. M_kappa as a multiplicity-counted large-values estimate: A_kappa(sigma) = (1+kappa-2sigma)/(kappa(1-sigma)).
  3. Part II's *pure* zero-density estimate for the distinct (sixth-power-free) row family
     (Prop detector-counts with no prime-amplitude information, x = 0, t optimised):
        R(delta) = 1 - delta + delta(5-6delta)/(3(5-2delta)),  A(sigma) = R/(1-sigma), delta = 2sigma-1.
     Compared with: DH (A=2), BGL Cor 1.6, Halasz-Montgomery variants, optimal-sextic-LS (dFDH conj.).
  4. Lemma `marked' used in the simple sixth-power extraction (slots as amplifier):
        theta(z) = 1/2 + (5m/6 - z)/2 with m = max(1+2z,(2+8z)/3)  -> minimum 11/12 at z = 0.
  5. Part II density saturation:
     (a) LP certificate: low side + the FLOOR-BIN high side (count R=1, any floor delta0>=0, any q>=0)
         => sigma0 >= 13/15 for every geometry.  The floor bin needs no zero witness and every row
         lies in it or above it, so NO row-count hypothesis (DH, beyond-DH, GRH) can remove it.
     (b) DH for rows (R = 1-delta for all delta in (0,5/6], any amplitude x) with floor -> 1/2
         attains exactly 13/15 at l_x = l_y = 2/5, ell = 1/5 (h = 4/5): every model constraint checked.
"""
from fractions import Fraction as F
from itertools import product

def hdr(s): print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)

# ---------------------------------------------------------------- 1
hdr("1. theta <-> A_eff  (target multiplicity H^{1/6};  H^{1/6} <= H^{A(1-beta)}  =>  beta <= 1-1/(6A))")
def theta_of_A(A): return 1 - 1 / (6 * A)
def A_of_theta(t): return 1 / (6 * (1 - t))
for name, t in [("11/12 (paper2 / Part I)", F(11, 12)), ("7/8 (Part II)", F(7, 8)),
                ("7/8 - 1/40000 (QRH-analysis)", F(34999, 40000)), ("13/15 (GRH-for-rows barrier)", F(13, 15)),
                ("6/7 (target)", F(6, 7))]:
    A = A_of_theta(t)
    assert theta_of_A(A) == t
    print(f"  theta = {str(t):>12s} = {float(t):.6f}   <->   A_eff = {str(A):>12s} = {float(A):.6f}   [{name}]")

# ---------------------------------------------------------------- 2
hdr("2. M_kappa => #{u<=H=D^kappa : |A_u(D)| >= D^sigma} << D^{1+kappa-2sigma} = H^{A_kappa(sigma)(1-sigma)}")
def A_kappa(k, s): return (1 + k - 2 * s) / (k * (1 - s))
for k in [F(1), F(9, 10), F(6, 7), F(3, 5)]:
    vals = ", ".join(f"A({s})={str(A_kappa(k, s))}" for s in [F(6, 7), F(7, 8), F(11, 12)])
    beta = (6 + 5 * k) / 12
    assert A_kappa(k, beta) * (1 - beta) == F(1, 6)       # beats multiplicity exactly above (6+5k)/12
    print(f"  kappa={str(k):>5s}: {vals};  multiplicity beaten iff sigma > (6+5kappa)/12 = {beta}")
assert A_kappa(F(6, 7), F(6, 7)) == F(7, 6)
print("  check: kappa = 6/7 gives A = 7/6 at sigma = 6/7  (the number quoted for route R3)  OK")
# zero-free consequence of M_kappa for each individual row:
print("  single-row consequence of M_kappa: |A_u(D)|^2 <= D^{1+kappa} => every row zero-free for Re s > (1+kappa)/2;"
      f" kappa=6/7 -> {(1 + F(6, 7)) / 2}")

# ---------------------------------------------------------------- 3
hdr("3. Part II's best PURE density estimate for distinct rows (Prop detector-counts, x = 0)")
ALPHA = F(5, 6)
def R_short(t, d, x=F(0)):
    Dx = 3 - F(17, 9) * x; Px = (2 - F(8, 9) * x) * (1 - x)
    return 1 - d + d * Px / Dx * (F(3, 2) - t)
def L_long(t, d): return 1 - d + (ALPHA - d) * (t - 1)
def R_pure(d):
    tstar = 5 / (5 - 2 * d)                       # crossing of R_short(.,d,0) and L(.,d)
    assert 1 <= tstar <= F(3, 2)
    assert R_short(tstar, d) == L_long(tstar, d)
    # brute-force min over a rational t-grid never beats the crossing value
    grid = [1 + F(j, 2000) for j in range(1001)]
    best = min(max(R_short(t, d), L_long(t, d)) for t in grid)
    assert best >= L_long(tstar, d)
    closed = 1 - d + d * (5 - 6 * d) / (3 * (5 - 2 * d))
    assert closed == L_long(tstar, d)
    return closed
def bgl_g(s):     # Blomer-Goldmakher-Louvel IMRN 2014, Cor 1.6 (5/6 < s <= 1)
    return 2 * (10 * s - 7) * (1 - s) / (24 * s - 12 * s * s - 11)
print("   sigma     delta   R_PartII(delta)  A_PartII   A_DH  A_BGL(Cor1.6)  A_HM(vanish)  A_HM(convex)")
for s in [F(6, 7), F(13, 15), F(7, 8), F(9, 10), F(11, 12)]:
    d = 2 * s - 1
    R = R_pure(d)
    A = R / (1 - s)
    Abgl = bgl_g(s) / (1 - s)
    Ahm_v = F(4)                       # length N >= H^{2+eps} kills all non-duplicate off-diagonals: R <= N^{2-2s}
    Ahm_c = 2 / (4 * s - 3)            # convexity for off-diagonal L(1/2, psi_u psibar_u'), conductor <= H^2
    print(f"  {str(s):>6s}  {str(d):>6s}   {str(R):>10s}     {float(A):6.4f}   2.0    {float(Abgl):6.4f}"
          f"        {float(Ahm_v):5.3f}       {float(Ahm_c):6.4f}")
assert R_pure(F(5, 7)) / (1 - F(6, 7)) == F(7, 3)
assert R_pure(F(3, 4)) / (1 - F(7, 8)) == F(16, 7)
assert R_pure(F(5, 6)) == 1 - F(5, 6)          # DH exactly at sigma = 11/12
print("  exact: A_PartII(6/7) = 7/3,  A_PartII(7/8) = 16/7,  A_PartII(11/12) = 2 (DH);  needed A_eff = 7/6")

# optimal sextic large sieve (de Faveri-Dunn-Hoffstein conj. (1.8), n=6), coefficient-blind, distinct rows:
#   Delta(A,B) = A + B + A^{5/6} B^{1/3} + A^{1/3} B^{5/6};  rows A = U, columns B = U^r;
#   count exponent  min_r max(1, r, 5/6 + r/3, 1/3 + 5r/6) + r(1-2s)
def A_optLS(s, steps=4000):
    best = None
    for j in range(1, steps + 1):
        r = F(3 * j, steps)
        e = max(F(1), r, F(5, 6) + r / 3, F(1, 3) + 5 * r / 6) + r * (1 - 2 * s)
        best = e if best is None or e < best else best
    return best / (1 - s)
print(f"  optimal-sextic-LS (dFDH conj. n=6) coefficient-blind density: A(6/7) ~ {float(A_optLS(F(6,7))):.4f},"
      f" A(7/8) ~ {float(A_optLS(F(7,8))):.4f}")

# ---------------------------------------------------------------- 4
hdr("4. Lemma marked as an AMPLIFIED simple extraction (rows u = a^6 included, slots a_i(p)=conj nu(p))")
# |M_1|^2 <= Z^{5m/6 - z}, columns Z^r with r = 1; constraints r+2z <= m, 2r+8z <= 3m (margins -> 0)
def theta_marked(z):
    m = max(1 + 2 * z, (2 + 8 * z) / 3)
    return F(1, 2) + (F(5, 6) * m - z) / 2
for z in [F(0), F(1, 20), F(1, 6), F(1, 4), F(1, 2), F(1)]:
    print(f"  slot length z = {str(z):>5s}:  theta = {str(theta_marked(z)):>8s} = {float(theta_marked(z)):.5f}")
assert all(theta_marked(F(j, 100)) >= F(11, 12) for j in range(0, 301))
assert theta_marked(F(0)) == F(11, 12)
print("  => slots cost 2z rows but amplify by only Z^z: theta(z) = 11/12 + z/3 (z<=1/2); never below 11/12.")

# ---------------------------------------------------------------- 5
hdr("5a. Part II: low side + FLOOR-BIN high side  =>  sigma0 >= 13/15 (any geometry, any row counts)")
# variables (l_x, l_y, ell, sigma, g, e); constraint rows a.v >= rhs
low      = ([F(-1, 6), 0, F(1, 6), 1, F(-1, 2), F(-1, 2)], F(5, 6))   # sigma >= 5/6 + lx/6 - ell/6 + g/2 + e/2
g_pos    = ([0, 0, 0, 0, 1, 0], F(0))                                  # g >= 0
g_b6     = ([F(1, 6), F(-1, 6), 0, 0, 1, 0], F(0))                     # g >= b/6,  b = l_y - l_x
e_pos    = ([0, 0, 0, 0, 0, 1], F(0))                                  # e >= 0
e_wall   = ([-1, -1, -1, 0, 0, 1], F(-1))                              # e >= M + ell - 1
def floor_row(d0, q):
    # floor bin: count R = 1 (no witness), amplitude q in [0, d0/2]:
    # sigma >= (1+d0)(1-l_y)/2 - ell/2 + q*ell + h(1 + d0/2 - 1/6),  h = 1 - l_x + ell
    c = 1 + d0 / 2 - F(1, 6)
    return ([c, (1 + d0) / 2, F(1, 2) - q - c, 1, 0, 0], (1 + d0) / 2 + c)
def combine(rows_w):
    coef = [F(0)] * 6; rhs = F(0)
    for (a, b), w in rows_w:
        assert w >= 0
        coef = [ci + w * ai for ci, ai in zip(coef, a)]; rhs += w * b
    return coef, rhs
W = [F(4, 5), F(1, 5), F(1, 5), F(1, 3), F(1, 15), F(1, 5)]
for d0, q in [(F(0), F(0)), (F(1, 50), F(0)), (F(1, 50), F(1, 100)), (F(1, 10), F(1, 20))]:
    coef, rhs = combine(list(zip([low, g_pos, g_b6, e_pos, e_wall, floor_row(d0, q)], W)))
    # leftover coefficients multiply l_x, l_y, ell (must be killed or have the right sign).  Write the
    # combination as  sigma + c.(lx,ly,ell) >= rhs ; since floor row at d0>0 differs from d0=0 only by
    # nonnegative multiples of (1-l_y) and h  (both >= 0 on the admissible region), we report exactly.
    extra = rhs - F(13, 15)
    print(f"  floor delta0={str(d0):>5s}, q={str(q):>5s}: combination  sigma + ({', '.join(str(c) for c in coef[:3])}).(lx,ly,ell)"
          f" + ({coef[4]})g + ({coef[5]})e >= {rhs}")
coef, rhs = combine(list(zip([low, g_pos, g_b6, e_pos, e_wall, floor_row(F(0), F(0))], W)))
assert coef == [0, 0, 0, 1, 0, 0] and rhs == F(13, 15)
print("  delta0 = 0: all geometric coefficients vanish exactly -> sigma0 >= 13/15 for EVERY geometry.")
print("  delta0 > 0, q >= 0: floor row = (delta0=0 row) + delta0[(1-l_y)/2 + h/2] + q*ell  >= (delta0=0 row)")
print("  on the admissible region (l_y <= 1, h >= 0, ell >= 0), so the bound 13/15 persists.")

hdr("5b. DH for rows (R = 1-delta), floor -> 1/2: the point l_x=l_y=2/5, ell=1/5 attains exactly 13/15")
lx, ly, ell = F(2, 5), F(2, 5), F(1, 5)
h = 1 - lx + ell; b = ly - lx; M = lx + ly
g = max(F(0), b / 6, 2 * b - ly)
e = max(F(0), (1 + 3 * ell - 2 * M) / 4, M + ell - 1)        # paper's case-A term (1+3ell-2M)/4
s_low = F(5, 6) + lx / 6 - ell / 6 + g / 2 + e / 2
def s_high(d, x): return (1 + d) * (1 - ly) / 2 - ell / 2 + x * d * ell + h * ((1 - d) + d / 2 - F(1, 6))
# affine in x for fixed d and affine in d for fixed x: check the four vertices of [0,5/6]x[0,1/2] + a grid
vals = [s_high(F(i, 600) * F(5, 6) * 600 / 600, F(j, 50)) for i in range(0, 601, 50) for j in range(0, 26)]
mx = max(s_high(F(5, 6) * F(i, 600), F(j, 50)) for i in range(601) for j in range(26))
print(f"  h = {h}, b = {b}, M+ell = {M+ell}, Gram g = {g}, energy e = {e}")
print(f"  low side = {s_low};  max high side over delta in [0,5/6], x in [0,1/2] = {mx}")
print(f"  high side identity: s_high(delta,x) = 13/15 + delta*(x/5 - 1/10)  ->  check at (5/6,1/2): {s_high(F(5,6),F(1,2))}")
assert s_low == F(13, 15) and mx == F(13, 15)
for i in range(0, 601):
    for j in range(0, 26):
        d, x = F(5, 6) * F(i, 600), F(j, 50)
        assert s_high(d, x) == F(13, 15) + d * (x / 5 - F(1, 10))
side = {"l_y/2 - h(17/50-1/6) > 0 (principal margin)": ly / 2 - h * (F(17, 50) - F(1, 6)),
        "l_y - ell >= 0": ly - ell, "M - 2 ell >= 0": M - 2 * ell, "5 ell - h > 0 (room above h)": 5 * ell - h}
for k, v in side.items():
    assert v > 0, k
    print(f"  side constraint {k:45s}: {v}")
print("  => in Part II, the density hypothesis for the row family already reaches the GRH-for-rows barrier 13/15;")
print("     any zero-density estimate stronger than DH (A < 2), even GRH, gains nothing more there.")
print("\nALL CHECKS PASSED")
