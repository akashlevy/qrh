"""
Parametrized GRH-for-rows template of the P78 Part II architecture (new-architectures agent).

Exponents are in units of log Z.  Geometry: X = Z^lx, Y = Z^ly, P = Z^ell, M = lx+ly, b = ly-lx.

Physical meaning of the knobs (base values reproduce the r4 model with grh=True, i.e. 13/15):
  c     copy exponent: the principal frequency has copies H = u a^N, N = 1/c (sextic: c = 1/6).
        Signal  X^{1/2-c} Z^{beta-1+c} P^{a};  low-side base  sigma >= (1-c) + c lx - a ell.
  a     slot amplification per unit ell (P78: a = c, it is the growth F^c of the number of copies
        when the frequency range F = ZP/X grows by P).
  rho   relative slot suppression of a non-principal row (P78: 1/2, square-root cancellation
        of sum_p chibar_p(u) p^{z-1}).
  kU    growth of the row range per unit ell: h = 1 - lx + kU*ell  (P78: 1).
  gb    Gram slope: ||A||^2 << (Q/Y)(1 + P_a^{gb} + P_a^2/Y), P_a = Y/X  (sextic: 1/6 = c).
  gq    weight of the P_a^2/Y Gram piece (1 = on).
  r     degree of the dual object used in the theta reflection: the reflected dual length of
        a row of length ZP twisted with modulus QP is (QP)^r/(ZP); diagonal energy needs
        (r-1)(M'+l') <= 1 + wE.  (P78: r = 2.)
  caseA include the paper's case-A energy piece (1+3l'-2M')/4 (pure-cube dual).  'paper',
        'refined' ((1+2l'-2M')/3) or 'none'.
  third include the (M'+4l'-2)/3 energy piece (only meaningful for r = 2).
  wE    wall relaxation.
  margin side-condition margin (paper 0.05).
HIGH (GRH for rows):  sigma >= (1-ly)/2 - rho*ell + (1-c)*h.
"""
import json, sys
import numpy as np
from scipy.optimize import linprog

BASE = dict(c=1/6, a=1/6, rho=0.5, kU=1.0, gb=1/6, gq=1.0, r=2, caseA="paper", third=True,
            wE=0.0, margin=0.05, extra=None)


def solve(cfg_over=None, verbose=False):
    cfg = dict(BASE); cfg.update(cfg_over or {})
    c, a, rho, kU, gb, gq, r = (cfg[k] for k in ("c", "a", "rho", "kU", "gb", "gq", "r"))
    wE = cfg["wE"]
    # variables v = (lx, ly, ell, sigma); constraints A v <= b
    A, B = [], []

    def ge(coef, const):            # sigma >= coef.(lx,ly,ell) + const
        A.append([coef[0], coef[1], coef[2], -1.0]); B.append(-const)

    def le(coef, const):            # coef.(lx,ly,ell) + const <= 0
        A.append([coef[0], coef[1], coef[2], 0.0]); B.append(-const)

    LX = np.array([1, 0, 0.]); LY = np.array([0, 1, 0.]); EL = np.array([0, 0, 1.])
    for lam in (0.0, 1.0):
        d = lam * EL
        Mp = LX + LY - 2 * d; lp = EL - d
        base = (c * LX - a * EL - d, 1 - c)
        gp = [(np.zeros(3), 0.0)]
        if gb > 0:
            gp.append((gb * (LY - LX), 0.0))
        if gq > 0:
            gp.append((gq * (2 * (LY - LX) - LY + d), 0.0))
        ep = [(np.zeros(3), 0.0)]
        if cfg["caseA"] == "paper":
            ep.append((0.75 * lp - 0.5 * Mp, 0.25))
        elif cfg["caseA"] == "refined":
            ep.append(((2 / 3) * lp - (2 / 3) * Mp, 1 / 3))
        ep.append(((r - 1) * (Mp + lp), -1 - wE))
        if cfg["third"] and r == 2:
            ep.append(((1 / 3) * Mp + (4 / 3) * lp, -2 / 3 - wE))
        for g in gp:
            for e in ep:
                ge(base[0] + 0.5 * g[0] + 0.5 * e[0], base[1] + 0.5 * g[1] + 0.5 * e[1])
    # high side, GRH for rows
    h = (-LX + kU * EL, 1.0)
    ge(-0.5 * LY - rho * EL + (1 - c) * h[0], 0.5 + (1 - c) * h[1])
    if cfg.get("extra"):
        for coef, const, kind in cfg["extra"]:
            (ge if kind == "ge" else le)(np.array(coef, float), const)
    m = cfg["margin"]
    le(-(LX + LY) + 2 * EL, m)        # M - 2 ell >= margin
    le(-LY + EL, m)                   # ly - ell >= margin
    le(-LX + EL, 0.0)                 # lx >= ell
    res = linprog([0, 0, 0, 1], A_ub=np.array(A), b_ub=np.array(B),
                  bounds=[(0, 2), (0, 3), (0, 2), (None, None)], method="highs")
    if res.status != 0:
        return None
    lx, ly, ell, sg = res.x
    return dict(theta=sg, lx=lx, ly=ly, ell=ell, M=lx + ly, h=1 - lx + kU * ell)


if __name__ == "__main__":
    over = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    print(json.dumps(solve(over)))
