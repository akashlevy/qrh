"""
R1 (direct mean square) -- exact exponent bookkeeping for the duality calculus of paper2.

All exponents are log_D of the quantity.  A "form" is a positive mean square
    Q = sum_{rows, N <= D^r} | sum_{cols, N ~ D^c} c_n psi_n(row) |^2
normalised so that the column coefficients are unimodular (sum |c_n|^2 = D^c).
t is the exponent of the TARGET bound needed for Q (or for its off-diagonal, see below).

Natural sizes (with an *optimal* large sieve granted at every level -- generous):
    diagonal  = r + c,        column term = 2c,        LS bound = max(r,c) + c.
    delta = t - (r + c)   (target minus diagonal / main term)
    eps   = t - 2c        (target minus column term)
The LS alone proves the target iff min(delta, eps) >= 0.

Moves (both are exact identities in paper2):
  P  (Poisson summation in the row variable; sextic or quadratic):
       off-diagonal of Q  =  D^{r-c} * Q',   Q' : rows 2c - r, cols c  (coefficients mu <-> Gauss)
       so (r, c, t) -> (2c - r, c, t - r + c).   [paper2 Lemma poisson / eq. initial-column-output]
       The diagonal of Q' is the 'fake diagonal' produced by removing coprimality with mu(v).
  R  (theta reflection in the column variable; only for cubic-Gauss/theta coefficients):
       column sum of length D^c  =  D^{c/2} * (sum of length D^{2r-c} with coefficients l^{-1/2}),
       i.e. in unimodular normalisation  Q = D^{c - (2r - c)} Q'',
       (r, c, t) -> (r, 2r - c, t + 2r - 2c).   [paper2 Prop. reflection, outline (intro-dual-theta-sum)]
Positivity moves (row enlargement r -> r' >= r with t fixed; Cauchy-Schwarz; separating weights)
can only decrease delta and eps.
"""
from fractions import Fraction as F

def state(r, c, t):
    return dict(r=r, c=c, t=t, delta=t - (r + c), eps=t - 2 * c)

def P(s):
    r, c, t = s["r"], s["c"], s["t"]
    return state(2 * c - r, c, t - r + c)

def R(s):
    r, c, t = s["r"], s["c"], s["t"]
    return state(r, 2 * r - c, t + 2 * r - 2 * c)

# coefficient types: 'mu' (A-level), 'gauss' (B-level: conj(alpha) gamma_2, sextic twist),
# 'theta_q' (C-level: theta coefficients, quadratic twist).  R is defined only on 'gauss' and 'theta_q'.
NEXT_TYPE = {("mu", "P"): "gauss", ("gauss", "P"): "mu",
             ("gauss", "R"): "theta_q", ("theta_q", "R"): "gauss",
             ("theta_q", "P"): "theta_q"}   # quadratic Poisson: quadratic Gauss sums are characters

def forced_path(kappa, steps=8):
    """From A = (kappa, 1, 1 + kappa): never undo the previous move (undoing returns to the parent)."""
    s = state(kappa, F(1), 1 + kappa); typ = "mu"; last = None
    out = [("A", typ, s)]
    for _ in range(steps):
        opts = [m for m in ("P", "R") if (typ, m) in NEXT_TYPE and m != last]
        assert len(opts) <= 1, (typ, last, opts)
        if not opts:          # only the undo move remains: the reachable graph is a finite path
            break
        m = opts[0]
        s = P(s) if m == "P" else R(s)
        typ = NEXT_TYPE[(typ, m)]
        last = m
        out.append((m, typ, s))
    return out

def show(kappa, steps=8):
    print(f"\nkappa = {kappa} (= {float(kappa):.4f}),  rho = 1 - kappa = {1 - kappa}")
    print(f"  {'move':4s} {'type':8s} {'rows r':>10s} {'cols c':>10s} {'target t':>10s} {'delta':>8s} {'eps':>8s}  min")
    for m, typ, s in forced_path(kappa, steps):
        mn = min(s["delta"], s["eps"])
        print(f"  {m:4s} {typ:8s} {str(s['r']):>10s} {str(s['c']):>10s} {str(s['t']):>10s} "
              f"{str(s['delta']):>8s} {str(s['eps']):>8s}  {str(mn):>6s}{'  <-- LS suffices' if mn >= 0 else ''}")

def invariants(kappa, steps=40):
    rho = 1 - kappa
    path = forced_path(kappa, steps)
    for m, typ, s in path:
        assert abs(s["eps"] - s["delta"]) == abs(rho), "|eps - delta| = |1 - kappa| must be invariant"
    mins = [min(s["delta"], s["eps"]) for _, _, s in path]
    return mins

if __name__ == "__main__":
    # 1. transformation rules on (delta, eps)
    import random
    for _ in range(1000):
        r, c, t = (F(random.randint(1, 60), random.randint(1, 12)) for _ in range(3))
        s = state(r, c, t)
        p, q = P(s), R(s)
        assert (p["delta"], p["eps"]) == (s["eps"], s["delta"])            # P swaps
        assert (q["delta"], q["eps"]) == (s["delta"], 2 * s["delta"] - s["eps"])   # R: eps -> 2 delta - eps
        assert P(P(s)) == s and R(R(s)) == s                                 # involutions
    print("checked: P swaps (delta,eps); R maps (delta,eps) -> (delta, 2 delta - eps); P,R involutions")

    # 2. paper2's own regime (kappa > 1) and the R1 regime (kappa < 1)
    for k in (F(11, 10), F(1), F(9, 10), F(6, 7), F(1, 2)):
        show(k, steps=6)

    # 3. along the forced path, min(delta, eps) never becomes >= 0 for kappa < 1 (exact, 40 steps)
    for k in (F(1, 100), F(1, 2), F(6, 7), F(99, 100)):
        mins = invariants(k)
        assert all(x < 0 for x in mins), k
        print(f"kappa={k}: path has {len(mins)} nodes (A..A'); max of min(delta,eps) = {max(mins)}  (= -(1-kappa) = {-(1-k)})")
    mins = invariants(F(11, 10))
    print(f"kappa=11/10: min(delta,eps) along path = {[str(x) for x in mins[:4]]} ... (>=0 from node C on)")

    # 4. what WOULD work: a column functional equation for mu (an R-move at level A) -- not available
    k = F(6, 7); s = state(k, F(1), 1 + k)
    sR = R(s)
    print(f"\nhypothetical R at level A (kappa=6/7): (delta,eps) {s['delta'],s['eps']} -> {sR['delta'],sR['eps']}"
          f"  => min = {min(sR['delta'], sR['eps'])} >= 0.  This is exactly the missing ingredient.")

    # 5. what the paper2 arrangement actually gives for kappa < 1
    #    (a) positivity at eq. (weighted): sum_u |A_u|^2 << D^{1+eps}(H + D)          -> exponent max(1+k, 2)
    #    (b) un-separated + reflected + quadratic LS (Prop R shape):
    #        sum_u|A_u|^2 << (H/D) * X * (calH + calH^2/X), X = D, calH = D^2/H  -> max(2, 3-k)
    for k in (F(1, 2), F(6, 7), F(19, 20), F(21, 20)):
        a = max(1 + k, F(2)); b = max(F(2), 3 - k, 1 + k)
        print(f"kappa={k}: target 1+kappa={1+k};  positivity route gives {a};  reflected-quadratic-LS route gives {b}")

    # 6. extraction formulas (paper2 Sec. 3 prime extraction; 7/8 paper Lemma inverse-amplification at U~1)
    def theta_prime(k): return max((6 + 5 * k) / 12, 1 - k / 6)
    def theta_int(k): return (6 + 5 * k) / 12
    for k in (F(6, 7), F(9, 10), F(1), F(11, 10), F(3, 5)):
        print(f"kappa={k}: prime extraction theta={theta_prime(k)}, integer extraction theta={theta_int(k)}")
    # minimum of the prime extraction over kappa:
    kstar = F(6, 7)
    assert (6 + 5 * kstar) / 12 == 1 - kstar / 6 == F(6, 7)
    print("prime extraction: (6+5k)/12 = 1 - k/6 exactly at k = 6/7, value 6/7 (its absolute floor).")
