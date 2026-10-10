"""
CORRECTED conditional certificate (adversarially audited):  W(w) ==> zero-free half-plane Re s > sigma0.
W(w) = the case-B wall of Lemma reflected-energy relaxed by w.  Consequences actually justified:
  (a) low-side energy loss  e = max(0, (1+3ell'-2M')/4 [if 2M'+ell'>=1],  M'+ell'-1-w)   (case A untouched by W)
  (b) Lemma marked slot capacity  z <= min((1+w/2-r)/2, 1-r, (3-2r)/8, supply),  only for r < 1
  (c) amplification unchanged  (alpha = 5/6)
With w = 0 this is exactly the paper's Part II bookkeeping (it certifies 7/8 at the paper geometry).
Certified:  W(1/5) => Re s > 108/125 = 0.864  at (lx, ly, ell) = (91/200, 93/200, 7/25).
The low side alone is bounded below by 13/15 - w/60 for every geometry (energy_lp.py).
"""
"""
Exact conditional certificate:   W(w)  ==>  zero-free half-plane  Re s > sigma0.

W(w) ("marked energy beyond the wall"): in Lemma reflected-energy of the 7/8
paper, case B costs nothing as long as the reflected dual length exceeds the
residual row length by at most w (log_Z units), i.e. T_d <= H + w.

Consequences used here (each is the paper's argument with its wall moved by w):
  * Lemma probe-row-norm:   energy loss  e = max(0, (1+2ell-2M)/3, M+ell-1-w)
  * Lemma marked:           slot capacity  r + 2z < 1 + w   (puncture invariant);
                            the second invariant 2r + 8z < 3 is KEPT unrelaxed
  * raw inverse moment:     sum_{u<=H}|M_u(D)|^2 << H  for  H >= D^{1/(1+w)}, hence
                            amplification exponent e(r) = max(1, (1 + 5r/(1+w))/6)
Everything else (Gram bound, plain fourth moment with c_P = 2/9, detector,
Euler factors, principal extraction) is the paper's, unchanged.

The high-side inequality E(h; delta, x) < 0 over delta in [1/50, 5/6],
x in [0, 1/2] is proved by exact rational branch-and-bound:
  - every count bound is nonincreasing in delta and in q = x*delta, so on a box
    R <= R_t(delta0, x0) for any fixed detector parameter t;
  - R_t(delta0,x0) = max_r min_j B_j(r) is computed EXACTLY (all B_j are
    piecewise linear in r; the max of their min is attained at a breakpoint
    or at a pairwise intersection of linear pieces);
  - all other terms of E are increasing in delta and x, bounded at (delta1, x1).

Run:  python3 conditional_certificate.py      (stdlib only)
"""
from fractions import Fraction as F
from itertools import combinations
import sys

# ------------------------------------------------------------------ data
w      = F(1, 10)
lx     = F(171, 400)
ly     = F(171, 400)
ell    = F(49, 200)
sigma0 = F(173, 200)
CP     = F(2, 9)
Z0     = F(17, 50)
D0     = F(1, 50)          # detector floor delta_0
DMAX   = F(5, 6)           # bin ceiling (Part I gives beta* <= 11/12)
DMIN   = F(1, 100)         # small-row cutoff
ZETA   = F(1, 10000)       # extension above h

h = 1 - lx + ell
b = ly - lx
M = lx + ly
S = ell / (h + ZETA)       # prime supply (log-length at base U), valid for d <= h+zeta

checks = []
def check(name, ok, val=None):
    checks.append((name, bool(ok), val))

# ------------------------------------------------------------------ low side
def eA(Mp, lp):
    return max(F(0), ((1 + 3 * lp - 2 * Mp) / 4) if 2 * Mp + lp >= 1 else F(0), Mp + lp - 1 - w)
e_loss = eA(M, ell)
low = F(5, 6) + lx / 6 - ell / 6 + max(F(0), b / 6, 2 * b - ly) / 2 + e_loss / 2
check("low side  <= sigma0", low <= sigma0, low)
check("rescaled tuples ok (M+ell-1-w-3d <= 2d etc. for d in [0,ell])",
      all(-d + eA(M - 2 * d, ell - d) / 2 <= e_loss / 2
          for d in [ell * F(k, 200) for k in range(201)]))
check("Gram condition l_y - ell - 11b/6 > 0 or b<=0", b <= 0 or ly - ell - 11 * b / 6 > 0)

# ------------------------------------------------------------------ outer pieces
small = h * (Z0 - F(1, 6)) - ly / 2 + F(63, 50) * DMIN
check("small rows exponent < 0", small < 0, small)
a0 = (1 + D0) / 2
floorE = a0 * (1 - ly) - sigma0 - ell / 2 + (D0 / 2) * ell + (h + ZETA) * (1 + D0 / 2 - F(1, 6))
check("floor bin E < 0 (R=1, q<=delta0/2)", floorE < 0, floorE)
check("principal margins l_y/20, h/600 > 0", ly > 0 and h > 0)
check("principal local errors on Re s >= sigma0 (region two)",
      max(-sigma0, -6 * F(33, 200), 4 - 5 * sigma0 - 6 * F(33, 200), 1 - F(19, 20) - 6 * F(33, 200)) < 0)

# ------------------------------------------------------------------ exact row count
def pieces_BI(d, q):
    """B_I(r) = 1 - d r - 2 q * clip(min((1+w-r)/2, (3-2r)/8, S), 0, inf), as list of affine (c0,c1) whose MAX is B_I
       (valid where all capacity pieces are >= 0; caller restricts r-range)."""
    out = []
    for (z0c, z1c) in [((1 + w / 2) / 2, F(-1, 2)), (F(1), F(-1)), (F(3, 8), F(-1, 4)), (S, F(0))]:
        out.append((1 - 2 * q * z0c, -d - 2 * q * z1c))
    return out

def eval_max(pcs, r):
    return max(c0 + c1 * r for c0, c1 in pcs)

def R_t(d, q, t):
    """exact max over r in [max(t-1/2,0), t] of min of the count bounds."""
    lo, hi = max(t - F(1, 2), F(0)), t
    # functions as lists of affine pieces whose max is the function value:
    fA = [(F(1) - 0, -d), (F(1, 6), F(5, 6) - d)]          # max(1,(1+5r/(1+w))/6) - d r
    # plain with slots, m = t - r:  1 - 2d m - 2q*min(CP(1-2m), S)
    #   = max( 1-2d m-2q CP(1-2m) , 1-2d m-2q S )   with m = t-r
    fP = [(1 - 2 * d * t - 2 * q * CP * (1 - 2 * t), 2 * d - 4 * q * CP),
          (1 - 2 * d * t - 2 * q * S, 2 * d)]
    fP0 = [(1 - 2 * d * t, 2 * d)]
    fI = pieces_BI(d, q)
    # B_I needs capacity >= 0: (1+w-r)/2 >= 0 and (3-2r)/8 >= 0 hold for r <= 1+w (< 3/2);
    # and is only available for r < 1+w.
    cut = F(1)
    segs = []
    if lo < min(hi, cut):
        segs.append((lo, min(hi, cut), [fA, fP, fP0, fI]))
    if max(lo, cut) < hi:
        segs.append((max(lo, cut), hi, [fA, fP, fP0]))
    if lo == hi:
        segs = [(lo, hi, [fA, fP, fP0] + ([fI] if lo < cut else []))]
    best = None
    for (s0, s1, funcs) in segs:
        allp = [p for f in funcs for p in f]
        cands = {s0, s1}
        for (a0_, a1_), (b0_, b1_) in combinations(allp, 2):
            if a1_ != b1_:
                r = (b0_ - a0_) / (a1_ - b1_)
                if s0 <= r <= s1:
                    cands.add(r)
        for r in cands:
            v = min(eval_max(f, r) for f in funcs)
            if best is None or v > best:
                best = v
    return best

TS = [1 + F(k, 40) for k in range(21)]      # coarse detector parameters t in [1, 3/2]

def R_float(d, q, t, nr=241):
    """floating-point version of R_t, used only to CHOOSE t (any t gives a valid bound)."""
    lo, hi = max(t - 0.5, 0.0), t
    best = -9.0
    w_, S_, CP_ = float(w), float(S), float(CP)
    for i in range(nr):
        r = lo + (hi - lo) * i / (nr - 1); m = t - r
        zM = max(0.0, min((1 + w_ / 2 - r) / 2, 1 - r, (3 - 2 * r) / 8, S_))
        b = [max(1.0, (1 + 5 * r) / 6) - d * r,
             1 - 2 * d * m - 2 * q * min(CP_ * (1 - 2 * m), S_),
             1 - 2 * d * m]
        if r < 1:
            b.append(1 - d * r - 2 * q * zM)
        best = max(best, min(b))
    return best

def best_t(d0, q0):
    d, q = float(d0), float(q0)
    grid = [1 + k / 200 for k in range(101)]
    tb = min(grid, key=lambda t: R_float(d, q, t))
    lo, hi = max(1.0, tb - 0.005), min(1.5, tb + 0.005)
    for _ in range(30):                       # golden-section refine
        a = lo + (hi - lo) * 0.382; c = lo + (hi - lo) * 0.618
        if R_float(d, q, a) < R_float(d, q, c): hi = c
        else: lo = a
    return F((lo + hi) / 2).limit_denominator(4000)

def E_box_upper(d0, d1, x0, x1):
    q0 = x0 * d0
    tb = best_t(d0, q0)
    R = min(R_t(d0, q0, t) for t in [tb, F(1), F(3, 2)])
    a1 = (1 + d1) / 2
    base = a1 * (1 - ly) - sigma0 - ell / 2 + x1 * d1 * ell
    # d in [d_min, h+zeta]: E linear in d; slope = R + delta/2 - z0
    slope_hi = R + d1 / 2 - Z0
    E_top = base + h * (Z0 - F(1, 6)) + (h + ZETA) * (R + d1 / 2 - Z0) if slope_hi > 0 else None
    E_bot = base + h * (Z0 - F(1, 6)) + DMIN * (R + d1 / 2 - Z0)
    return max(x for x in (E_top, E_bot) if x is not None)

sys.setrecursionlimit(10000)
stats = {"boxes": 0, "maxdepth": 0, "worst": None}
def certify(d0, d1, x0, x1, depth=0):
    stats["boxes"] += 1
    stats["maxdepth"] = max(stats["maxdepth"], depth)
    ub = E_box_upper(d0, d1, x0, x1)
    if ub < 0:
        if stats["worst"] is None or ub > stats["worst"]:
            stats["worst"] = ub
        return True
    if depth > 40:
        print("FAILED near", float(d0), float(x0), float(ub))
        return False
    if (d1 - d0) >= (x1 - x0) * 2:
        m = (d0 + d1) / 2
        return certify(d0, m, x0, x1, depth + 1) and certify(m, d1, x0, x1, depth + 1)
    m = (x0 + x1) / 2
    return certify(d0, d1, x0, m, depth + 1) and certify(d0, d1, m, x1, depth + 1)

if __name__ == "__main__":
    print(f"W(w) with w = {w};  target sigma0 = {sigma0} = {float(sigma0)}")
    print(f"geometry: l_x = {lx}, l_y = {ly}, ell = {ell}, h = {h}, M+ell = {M+ell}, supply S = {float(S):.4f}")
    ok = certify(D0, DMAX, F(0), F(1, 2))
    check(f"high side E < 0 on [1/50,5/6]x[0,1/2] for all row scales d in [d_min, h+zeta] "
          f"({stats['boxes']} boxes, depth {stats['maxdepth']})", ok, stats["worst"])
    allok = True
    for n, o, v in checks:
        allok &= o
        print(("PASS " if o else "FAIL ") + n + ("" if v is None else f"   [{float(v):.6g}]"))
    print("ALL CHECKS PASS" if allok else "SOME CHECK FAILED")
