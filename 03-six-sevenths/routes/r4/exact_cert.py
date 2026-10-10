"""
Exact (stdlib Fraction) certificate for the Part II bookkeeping with switchable lemma inputs.
Generalizes /Users/akashlevy/Documents/riemann/strengthening/conditional_certificate.py:
  wE     : energy wall relaxation (case B of Lemma reflected-energy: loss (M'+l'-1-wE)_+ and the (M'+4l'-2)/3-wE piece)
  wM     : relaxation of Lemma marked's first width condition r + 2z < 1 + wM (second kept: 2r+8z < 3)
  caseA  : 'paper'  -> (1+3l'-2M')/4  ;  'refined' -> (1+2l'-2M')/3   (refined count of powerful row parts)
  gram_b : 1 (paper, g >= b/6) or 0 (no P_a^{1/6} term)
  kA     : amplification slope; e(r) = max(1, (1+5 kA r)/6)   (kA=1 paper; kA<1 = raw moment for H >= D^kA)
  D0     : detector floor delta_0 (paper 1/50)
Checks (all exact):
  LOW   sigma0 >= 5/6 + lx/6 - ell/6 - d + g(d)/2 + e(d)/2   at d = 0 and d = ell (affine in d between)
  FLOOR, SMALL ROWS, and HIGH: E(d) < 0 for delta in [D0, DMAX], x in [0,1/2], all row scales d in [d_min, h+zeta]
  (branch and bound; R is nonincreasing in delta and in q = x delta, every other term of E is increasing).
Usage: python3 exact_cert.py key=value ...   (fractions as a/b)
"""
from fractions import Fraction as F
from itertools import combinations
import sys, json, time

P = dict(lx=F(4, 7), ly=F(4, 7), ell=F(3, 7), sigma0=F(6, 7), wE=F(4, 7), wM=F(0), caseA="paper",
         gram_b=F(1), kA=F(1), D0=F(1, 50), DMAX=F(5, 6), CP=F(2, 9), Z0=F(17, 50), DMIN=F(1, 100), ZETA=F(1, 10000))
for a in sys.argv[1:]:
    k, v = a.split("=")
    if k in ("caseA",):
        P[k] = v
    else:
        n, _, d = v.partition("/"); P[k] = F(int(n), int(d or 1))
lx, ly, ell, sigma0 = P["lx"], P["ly"], P["ell"], P["sigma0"]
wE, wM, kA, D0, DMAX, CP, Z0, DMIN, ZETA = (P[k] for k in ("wE", "wM", "kA", "D0", "DMAX", "CP", "Z0", "DMIN", "ZETA"))
h = 1 - lx + ell; b = ly - lx; M = lx + ly
checks = []
def check(name, ok, val=None):
    checks.append((name, bool(ok), val))

# ---------------------------------------------------------------- low side
def energy(Mp, lp):
    ca = (1 + 3 * lp - 2 * Mp) / 4 if P["caseA"] == "paper" else (1 + 2 * lp - 2 * Mp) / 3
    return max(F(0), ca, Mp + lp - 1 - wE, (Mp + 4 * lp - 2) / 3 - wE)
def gram(d):
    return max(F(0), P["gram_b"] * b / 6, 2 * b - ly + d) if b > 0 else F(0)
lows = []
for d in (F(0), ell):
    lows.append(F(5, 6) + lx / 6 - ell / 6 - d + gram(d) / 2 + energy(M - 2 * d, ell - d) / 2)
low = max(lows)
check("low side <= sigma0", low <= sigma0, low)
check("geometry: M-2ell>0, ly-ell>0, lx-ell>0", M - 2 * ell > 0 and ly - ell > 0 and lx - ell > 0)
small = h * (Z0 - F(1, 6)) - ly / 2 + F(63, 50) * DMIN
check("small rows exponent < 0", small < 0, small)
a0 = (1 + D0) / 2
floorE = a0 * (1 - ly) - sigma0 - ell / 2 + (D0 / 2) * ell + (h + ZETA) * (1 + D0 / 2 - F(1, 6))
check("floor bin E(h+zeta) < 0", floorE < 0, floorE)
S = ell / (h + ZETA)          # prime supply at base U, valid for every d <= h + zeta

# ---------------------------------------------------------------- exact row count
def R_t(d, q, t):
    """exact max over r in [max(t-1/2,0), t] of min{fA, fI (r<1+wM), fP, fP0}."""
    lo, hi = max(t - F(1, 2), F(0)), t
    # each function = list of (affine pieces) whose MIN or MAX gives it; we evaluate exactly at candidates
    pieces = [(F(1), -d), (F(1, 6), 5 * kA / 6 - d),                                   # fA
              (1 - 2 * q * (1 + wM) / 2, -d + q), (1 - 2 * q * F(3, 8), -d + q / 2), (1 - 2 * q * S, -d), (F(1), -d),  # fI
              (1 - 2 * d * t - 2 * q * CP * (1 - 2 * t), 2 * d - 4 * q * CP), (1 - 2 * d * t - 2 * q * S, 2 * d),  # fP
              (1 - 2 * d * t, 2 * d)]                                                   # fP0
    cands = {lo, hi}
    for bp in (1 / kA, 1 + wM, t - F(1, 2), t - (1 - S / CP) / 2, 2 * ((1 + wM) / 2 - S), 4 * (F(3, 8) - S),
               1 + wM - 4 * F(3, 8) + F(3, 2) - 1):
        cands.add(bp)
    for (a1, b1), (a2, b2) in combinations(pieces, 2):
        if b1 != b2:
            cands.add((a2 - a1) / (b1 - b2))
    best = None
    for r in cands:
        if not (lo <= r <= hi):
            continue
        for side in (0, 1):     # evaluate at r and at the right limit of the marked range
            if side == 1 and r != 1 + wM:
                continue
            fA = max(F(1), (1 + 5 * kA * r) / 6) - d * r
            vals = [fA]
            if r < 1 + wM and side == 0:
                zI = min((1 + wM - r) / 2, (3 - 2 * r) / 8, S)
                zI = max(zI, F(0))
                vals.append(1 - d * r - 2 * q * zI)
            m = t - r
            if m < F(1, 2):
                zP = min(CP * (1 - 2 * m), S)
                vals.append(1 - 2 * d * m - 2 * q * zP)
            vals.append(1 - 2 * d * m)
            v = min(vals)
            if best is None or v > best:
                best = v
    return best

def _Rf(d, q, ts):
    import numpy as np
    ts = np.asarray(ts, float)[:, None]
    lo = np.maximum(ts - 0.5, 0.0); r = lo + (ts - lo) * np.linspace(0, 1, 401)[None, :]; m = ts - r
    fA = np.maximum(1.0, (1 + 5 * float(kA) * r) / 6) - d * r
    zI = np.clip(np.minimum(np.minimum((1 + float(wM) - r) / 2, (3 - 2 * r) / 8), float(S)), 0, None)
    fI = np.where(r < 1 + float(wM), 1 - d * r - 2 * q * zI, np.inf)
    fP = np.where(m < 0.5, 1 - 2 * d * m - 2 * q * np.minimum(float(CP) * (1 - 2 * m), float(S)), np.inf)
    fP0 = 1 - 2 * d * m
    return np.minimum(np.minimum(fA, fI), np.minimum(fP, fP0)).max(axis=1)

def R_float_best_t(d, q):
    import numpy as np
    from model import R_batch, DEFAULT
    cfg = dict(DEFAULT); cfg.update(kA=float(kA), w1=float(wM), w2=0.0, cP=float(CP))
    d, q, Sf = float(d), float(q), float(S)
    ts = np.linspace(1.0, 1.5, 201)
    v = R_batch(np.full_like(ts, d), np.full_like(ts, q), Sf, ts, cfg)
    k = int(v.argmin()); a, c = ts[max(k - 1, 0)], ts[min(k + 1, 200)]
    for _ in range(50):
        m1 = a + (c - a) * 0.382; m2 = a + (c - a) * 0.618
        f = R_batch(np.array([d, d]), np.array([q, q]), Sf, np.array([m1, m2]), cfg)
        if f[0] < f[1]: c = m2
        else: a = m1
    return F((a + c) / 2).limit_denominator(10**9)

def E_box_upper(d0, d1, x0, x1):
    q0 = x0 * d0
    tb = R_float_best_t(d0, q0)
    R = min(R_t(d0, q0, t) for t in {tb, F(1), F(3, 2)})
    a1 = (1 + d1) / 2
    base = a1 * (1 - ly) - sigma0 - ell / 2 + x1 * d1 * ell
    slope_hi = R + d1 / 2 - Z0
    vals = [base + h * (Z0 - F(1, 6)) + DMIN * (R + d1 / 2 - Z0)]
    vals.append(base + h * (Z0 - F(1, 6)) + (h + ZETA) * slope_hi)
    return max(vals)

stats = dict(boxes=0, worst=None)
def certify(d0, d1, x0, x1, depth=0):
    stats["boxes"] += 1
    ub = E_box_upper(d0, d1, x0, x1)
    if ub < 0:
        stats["worst"] = ub if stats["worst"] is None else max(stats["worst"], ub)
        return True
    if depth > 30:
        print("FAILED near delta=%.5f x=%.4f ub=%.3g" % (float(d0), float(x0), float(ub)))
        return False
    if (d1 - d0) >= (x1 - x0) * 2:
        m = (d0 + d1) / 2
        return certify(d0, m, x0, x1, depth + 1) and certify(m, d1, x0, x1, depth + 1)
    m = (x0 + x1) / 2
    return certify(d0, d1, x0, m, depth + 1) and certify(d0, d1, m, x1, depth + 1)

if __name__ == "__main__":
    t0 = time.time()
    ok = certify(D0, DMAX, F(0), F(1, 2))
    check(f"high side E<0 on [D0,DMAX]x[0,1/2], all d in [d_min,h+zeta] ({stats['boxes']} boxes)", ok, stats["worst"])
    allok = all(o for _, o, _ in checks)
    print(json.dumps({k: (str(v) if isinstance(v, F) else v) for k, v in P.items()}))
    print(f"h = {h}, supply S = {S} ~ {float(S):.4f}, low = {low} ~ {float(low):.6f}")
    for n, o, v in checks:
        print(("PASS " if o else "FAIL ") + n + ("" if v is None else f"   [{float(v):.6g}]"))
    print("ALL CHECKS PASS" if allok else "SOME CHECK FAILED", f"({time.time()-t0:.0f}s)")
