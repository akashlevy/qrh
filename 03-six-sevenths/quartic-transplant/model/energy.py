"""
Generalized reflected-energy LP (analogue of Lemma reflected-energy + Lemma probe-row-norm of the 7/8 paper).

Block data (log_Z lengths, all >= 0):
  O  : powerful part of the row (count Z^{O/2}: every powerful ideal is x^2 y^3)
  Cp : conductor contribution of the frozen powerful part (Cp <= cpow*O)
  H  : residual squarefree row length
  za : active mark (slot) length, za <= l'
  v  : dual squarefree dyad,  lb : dual completion dyad (cube part b in the cubic case)
Reflection geometry:
  T_d = crow*H + Cp + cmark*za - N_*,       N_* = 1 + l'      (completed length incl. marked slots)
  retained dual:  y = v + kb*lb <= T_d
Energy (relative to the row count M'):
  cnt*O + Sieve(H, v[+lb]) - lb + Marks(v, za) - v - ks*(T_d - y) - M'
  Sieve = max over affine pieces (rows x dual large sieve), Marks = max over pieces (mark sieve) - v.
The maximum over block data of a max of affine pieces = max over piece choices of an LP.
Configurations:
  cubic (7/8 paper): crow=cmark=2, Sieve = quadratic (H, v+lb) [Goldmakher-Louvel/HB], Marks = cubic HB
                     (v, za, 2(v+za)/3), ks=1/2 (kernel m_ker = Z^{-D/4}), kb=3, lb inside the row sieve.
  quartic (derived): crow=cmark=3 (row shift h^3; mark p|c -> shift p^2 on cofactor + length/N(p)),
                     Sieve = de Faveri n=4 (H, v, 3H/4+v/2, H/2+3v/4) or perfect (H, v) or LoA (H),
                     Marks = quadratic (v, za) [mark twist chi_m(p)^2], ks=3/2 (U^natural(x)~x^{1/4}, dyad
                     normalization gives Z^{-3(T_d-y)/4}) or the conservative 1/2, d^4 completion outside sieve.
"""
import numpy as np, itertools
from scipy.optimize import linprog

names = ["O", "Cp", "H", "za", "v", "lb"]
IX = {k: i for i, k in enumerate(names)}
NV = len(names)

def vec(**kw):
    r = np.zeros(NV)
    for k, c in kw.items():
        r[IX[k]] += c
    return r

CUBIC = dict(crow=2.0, cmark=2.0, cpow=1.0, cnt=0.5, ks=0.5, kb=3.0, lb_in_sieve=True,
             sieve=[(1, 0), (0, 1)], marks=[(1, 0), (0, 1), (2/3, 2/3)], wall=0.0)
QUARTIC_DF = dict(crow=3.0, cmark=3.0, cpow=1.0, cnt=0.5, ks=1.5, kb=4.0, lb_in_sieve=False,
                  sieve=[(1, 0), (0, 1), (0.75, 0.5), (0.5, 0.75)], marks=[(1, 0), (0, 1)], wall=0.0)
QUARTIC_PERFECT = dict(QUARTIC_DF, sieve=[(1, 0), (0, 1)])
QUARTIC_LOA = dict(QUARTIC_DF, sieve=[(1, 0)])
QUARTIC_BGL = dict(QUARTIC_DF, sieve=[(1, 0), (0, 1), (2/3, 2/3)])

def energy_loss(Mp, lp, cfg, ret=False):
    """max over block data of (energy - M'), floored at 0 (an empty dual range costs nothing)."""
    Ns = 1.0 + lp
    Td = vec(H=cfg["crow"], Cp=1.0, za=cfg["cmark"])
    y = vec(v=1.0, lb=cfg["kb"])
    A = [vec(H=1, O=1), vec(Cp=1, O=-cfg["cpow"]), vec(za=1), y - Td]
    b = [Mp, 0.0, lp, -Ns]
    if not cfg["lb_in_sieve"]:
        A.append(vec(lb=1)); b.append(0.0)        # completion part outside the sieve only shortens v: worst at lb=0
    A = np.array(A); b = np.array(b)
    best = (0.0, None)
    for sp, mp in itertools.product(cfg["sieve"], cfg["marks"]):
        c = vec(O=cfg["cnt"])
        if cfg["lb_in_sieve"]:
            c += sp[0] * vec(H=1) + sp[1] * vec(v=1, lb=1) - vec(lb=1)
        else:
            c += sp[0] * vec(H=1) + sp[1] * vec(v=1)
        c += mp[0] * vec(v=1) + mp[1] * vec(za=1) - vec(v=1)
        c -= cfg["ks"] * (Td - y)
        c0 = cfg["ks"] * Ns - Mp
        # optional wall relaxation (only meaningful for the column pieces): subtract w when the dual piece is used
        if sp[1] > 0:
            c0 -= cfg.get("wall", 0.0) * (1.0 if sp == (0, 1) else 0.0)
        res = linprog(-c, A_ub=A, b_ub=b, bounds=[(0, None)] * NV, method="highs")
        if res.status == 0:
            val = -res.fun + c0
            if val > best[0] + 1e-12:
                best = (val, (sp, mp, dict(zip(names, np.round(res.x, 5)))))
        elif res.status == 3:
            return (np.inf, None) if ret else np.inf
    return best if ret else best[0]

def cubic_closed(Mp, lp, w=0.0):
    return max(0.0, (1 + 3 * lp - 2 * Mp) / 4, Mp + lp - 1 - w, (Mp + 4 * lp - 2) / 3 - w)

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    worst = 0; arg = None
    for _ in range(3000):
        Mp = rng.uniform(0.2, 1.4); lp = rng.uniform(0, 0.6)
        e = energy_loss(Mp, lp, CUBIC); c = cubic_closed(Mp, lp)
        if abs(e - c) > worst:
            worst = abs(e - c); arg = (round(Mp, 4), round(lp, 4), e, c, energy_loss(Mp, lp, CUBIC, True)[1])
    print("cubic: max |LP - paper closed form max(0,(1+3l'-2M')/4, M'+l'-1, (M'+4l'-2)/3)| =", worst, arg)

# ---------------------------------------------------------------------------------------------
# Refined count of the frozen powerful row part, by valuation type (quartic):
#   a prime with row valuation e >= 2 and log-norm a contributes  O += e*a,  count += a,
#   conductor (class of k^3 mod 4th powers) Cp += c_e*a  with c_e = (3e mod 4):  e=2:2, 3:1, 4:0, 5:3, 6:2, 7:1
# (the crude count O/2 with Cp=0 is only attained by 4th powers, whose true count is O/4).
QTYPES = [(2, 2), (3, 1), (4, 0), (5, 3), (6, 2), (7, 1)]
names2 = ["H", "za", "v"] + ["a%d" % e for e, _ in QTYPES]
IX2 = {k: i for i, k in enumerate(names2)}
NV2 = len(names2)
def vec2(**kw):
    r = np.zeros(NV2)
    for k, c in kw.items():
        r[IX2[k]] += c
    return r

def energy_loss_refined(Mp, lp, cfg, ret=False, Ns=None):
    """quartic energy with the refined powerful count (no completion part: it lies outside the sieve)."""
    Ns = 1.0 + lp if Ns is None else Ns
    Ovec = sum(e * vec2(**{"a%d" % e: 1}) for e, _ in QTYPES)
    Cvec = sum(c * vec2(**{"a%d" % e: 1}) for e, c in QTYPES)
    cnt = sum(vec2(**{"a%d" % e: 1}) for e, _ in QTYPES)
    Td = cfg["crow"] * vec2(H=1) + Cvec + cfg["cmark"] * vec2(za=1)
    y = vec2(v=1)
    A = np.array([vec2(H=1) + Ovec, vec2(za=1), y - Td]); b = np.array([Mp, lp, -Ns])
    best = (0.0, None)
    for sp, mp in itertools.product(cfg["sieve"], cfg["marks"]):
        c = cnt + sp[0] * vec2(H=1) + sp[1] * vec2(v=1) + mp[0] * vec2(v=1) + mp[1] * vec2(za=1) - vec2(v=1)
        c = c - cfg["ks"] * (Td - y)
        c0 = cfg["ks"] * Ns - Mp
        res = linprog(-c, A_ub=A, b_ub=b, bounds=[(0, None)] * NV2, method="highs")
        if res.status == 0:
            val = -res.fun + c0
            if val > best[0] + 1e-12:
                best = (val, (sp, mp, {k: round(float(x), 5) for k, x in zip(names2, res.x) if abs(x) > 1e-9}))
        elif res.status == 3:
            return (np.inf, None) if ret else np.inf
    return best if ret else best[0]
