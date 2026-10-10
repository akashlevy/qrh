"""
P78 Part II template (GRH for rows) with a pluggable low side, to test the cubic
'product of two Gauss sums' architecture:

  twist family  chi_A(m) = (m/A)_3   (cubic, copies u = a^3, c = 1/3),
  Gauss side    A_m = Y^{-1} sum_s c(s) gamma(psi_s) psi_s(m)bar   (cubic Gauss sums),
  theta side    B_m = sum_A T(A)^2 eta(A) psi_A(m) [marks]   (T = cubic theta coefficients,
                T(A)^2 = Hecke * mu(A) * conj-Gauss: the Rankin-Selberg square theta(x)theta),
  absorption    T(A)^2 * gamma(psi_A) = mu(A) * Hecke(A)  (Gauss: g^3 = p J).

Low side (Cauchy-Schwarz in m):
  Gram  ||A||^2 << (Q/Y)(1 + P_a^{1/3})        (cubic large sieve; DFDH: sharp)
  Energy ||B||^2 / (Q * diag) = Z^{E}, where E is the best of four 'orbit nodes'
        (direct cubic LS / Poisson / RS-reflection with conductor (QP)^{cB} / reflection+Poisson),
        each a max of linear pieces; the node is enumerated (min of max is not an LP).
High side (GRH for rows):  sigma >= (1-ly)/2 - rho*ell + (1-c)(1 - lx + kU*ell)   [template_lp.py]

Exponents in units of log Z; X=Z^lx, Y=Z^ly, P=Z^ell, M=lx+ly, b=ly-lx.
Validation: mode='sextic' reproduces template_lp.py (13/15 under GRH for rows).
"""
import itertools, json, sys
import numpy as np
from scipy.optimize import linprog

LX = np.array([1, 0, 0.]); LY = np.array([0, 1, 0.]); EL = np.array([0, 0, 1.])
Z3 = np.zeros(3)


def energy_nodes(Mp, lp, cB, k):
    """list of nodes; each node = list of (coef, const) linear pieces of the excess exponent.
    f = M' (family), l = 1 + l' (row length), lR = cB(M'+l') - (1+l') (reflected length)."""
    f = (Mp, 0.0)
    l = (lp, 1.0)
    lR = (cB * (Mp + lp) - lp, -1.0)

    def lin(*terms):  # sum of (scalar, (coef,const))
        c = Z3.copy(); k0 = 0.0
        for s, (cc, kk) in terms:
            c = c + s * cc; k0 += s * kk
        return (c, k0)
    zero = (Z3, 0.0)
    nodes = []
    if k == 2:   # quadratic residual: perfect large sieve
        nodes.append([zero, lin((1, lR), (-1, f))])                           # reflect + quadratic LS
    else:
        nodes.append([zero, lin((1, l), (-1, f)), lin((2 / 3, l), (-1 / 3, f))])        # direct LS
        nodes.append([zero, lin((1, l), (-2 / 3, f))])                                  # Poisson
        nodes.append([zero, lin((1, lR), (-1, f)), lin((2 / 3, lR), (-1 / 3, f))])      # reflect
        nodes.append([zero, lin((1, lR), (-2 / 3, f))])                                 # reflect+Poisson
    return nodes


def solve(cfg, verbose=False):
    c, a, rho, kU = cfg["c"], cfg["a"], cfg["rho"], cfg["kU"]
    gb, gq = cfg["gb"], cfg["gq"]
    cB, k = cfg["cB"], cfg["k"]
    best = None
    lams = (0.0, 1.0)
    # node choice per lambda-case
    n_nodes = len(energy_nodes(LX, EL, cB, k))
    for choice in itertools.product(range(n_nodes), repeat=len(lams)):
        A, B = [], []

        def ge(coef, const):
            A.append([coef[0], coef[1], coef[2], -1.0]); B.append(-const)

        def le(coef, const):
            A.append([coef[0], coef[1], coef[2], 0.0]); B.append(-const)
        for lam, ch in zip(lams, choice):
            d = lam * EL
            Mp = LX + LY - 2 * d; lp = EL - d
            base = (c * LX - a * EL - d, 1 - c)
            gp = [(Z3, 0.0)]
            if gb > 0: gp.append((gb * (LY - LX), 0.0))
            if gq > 0: gp.append((gq * (2 * (LY - LX) - LY + d), 0.0))
            ep = list(energy_nodes(Mp, lp, cB, k)[ch])
            if cfg.get("caseA") == "paper":
                ep.append((0.75 * lp - 0.5 * Mp, 0.25))
            if cfg.get("third"):
                ep.append(((1 / 3) * Mp + (4 / 3) * lp, -2 / 3))
            if cfg.get("wall_extra"):        # optional extra explicit wall piece
                for (cc, kk) in cfg["wall_extra"](Mp, lp):
                    ep.append((cc, kk))
            for g in gp:
                for e in ep:
                    ge(base[0] + 0.5 * g[0] + 0.5 * e[0], base[1] + 0.5 * g[1] + 0.5 * e[1])
        h = (-LX + kU * EL, 1.0)
        ge(-0.5 * LY - rho * EL + (1 - c) * h[0], 0.5 + (1 - c) * h[1])
        m = cfg["margin"]
        le(-(LX + LY) + 2 * EL, m)
        le(-LY + EL, m)
        le(-LX + EL, 0.0)
        res = linprog([0, 0, 0, 1], A_ub=np.array(A), b_ub=np.array(B),
                      bounds=[(0, 2), (0, 3), (0, 2), (None, None)], method="highs")
        if res.status != 0:
            continue
        if best is None or res.x[3] < best["theta"] - 1e-12:
            lx, ly, ell, sg = res.x
            best = dict(theta=sg, lx=lx, ly=ly, ell=ell, M=lx + ly, nodes=choice)
    return best


SEXTIC = dict(c=1 / 6, a=1 / 6, rho=0.5, kU=1.0, gb=1 / 6, gq=1.0, cB=2, k=2,
              caseA="paper", third=True, margin=0.05)
CUBIC = dict(c=1 / 3, a=1 / 3, rho=0.5, kU=1.0, gb=1 / 3, gq=0.0, cB=2, k=3,
             caseA=None, third=False, margin=0.05)

if __name__ == "__main__":
    from fractions import Fraction as Fr

    def show(lab, cfg):
        r = solve(cfg)
        if r is None:
            print(f"{lab:66s} infeasible"); return
        print(f"{lab:66s} theta={r['theta']:.6f} ~{Fr(r['theta']).limit_denominator(300)}  "
              f"lx={r['lx']:.3f} ly={r['ly']:.3f} ell={r['ell']:.3f} M={r['M']:.3f} nodes={r['nodes']}")
    print("== validation ==")
    show("sextic P78 (quadratic residual, caseA, third): expect 13/15", SEXTIC)
    print("== cubic RS architecture (GRH for rows) ==")
    for cB in (1, 2, 3, 4):
        show(f"cubic, cB={cB}, cubic LS base, no caseA/third", dict(CUBIC, cB=cB))
    show("cubic, cB=2, + caseA (paper form)", dict(CUBIC, caseA="paper"))
    show("cubic, cB=2, + third piece", dict(CUBIC, third=True))
    show("cubic, cB=2, + P78-type gq Gram piece", dict(CUBIC, gq=1.0))
    show("cubic, cB=2, hypothetical perfect residual LS (k=2)", dict(CUBIC, k=2))
    show("cubic, cB=2, margin 0", dict(CUBIC, margin=0.0))
    show("cubic, cB=2, no marks (ell=0 forced via a=0,rho=0,kU=0)", dict(CUBIC, a=0.0, rho=0.0, kU=0.0))
