"""Apples-to-apples: add the trivial 'Poisson node' (energy diagonal once Q >= (ZP)^2,
BGL (FL)^{2/3} in between) to BOTH the sextic template and the cubic model, and restrict
the geometry box to see where the GRH floor comes from."""
import numpy as np
from fractions import Fraction as Fr
import p78_cubic_lp as P

LX, LY, EL, Z3 = P.LX, P.LY, P.EL, P.Z3

def nodes_sextic_plus(Mp, lp, cB, k):
    f = (Mp, 0.0); l = (lp, 1.0); lR = (cB * (Mp + lp) - lp, -1.0)
    def lin(*terms):
        c = Z3.copy(); k0 = 0.0
        for s, (cc, kk) in terms:
            c = c + s * cc; k0 += s * kk
        return (c, k0)
    zero = (Z3, 0.0)
    return [[zero, lin((1, lR), (-1, f))],              # theta reflection + quadratic LS (paper)
            [zero, lin((1, l), (-2 / 3, f))]]           # Poisson in m + BGL on Mobius dual

def run(lab, cfg, nodes_fn=None, box=None):
    old = P.energy_nodes
    if nodes_fn: P.energy_nodes = nodes_fn
    # crude box restriction through linprog bounds: monkeypatch
    import scipy.optimize as so
    orig = so.linprog
    if box:
        def lp2(*a, **kw):
            kw["bounds"] = box
            return orig(*a, **kw)
        P.linprog = lp2
    r = P.solve(cfg)
    P.energy_nodes = old; P.linprog = orig
    if r is None:
        print(f"{lab:70s} infeasible"); return
    print(f"{lab:70s} theta={r['theta']:.6f} ~{Fr(r['theta']).limit_denominator(300)} "
          f"lx={r['lx']:.3f} ly={r['ly']:.3f} ell={r['ell']:.3f} M={r['M']:.3f} nodes={r['nodes']}")

big = [(0, 2), (0, 3), (0, 2), (None, None)]
small = [(0, 1), (0, 1), (0, 1), (None, None)]
print("== sextic, paper low side + Poisson node ==")
run("sextic paper (box big)", P.SEXTIC, box=big)
run("sextic + Poisson node (box big)", P.SEXTIC, nodes_sextic_plus, box=big)
run("sextic + Poisson node, no caseA (box big)", dict(P.SEXTIC, caseA=None), nodes_sextic_plus, box=big)
run("sextic + Poisson node, no caseA, no third (box big)", dict(P.SEXTIC, caseA=None, third=False), nodes_sextic_plus, box=big)
print("== cubic, box restricted lx,ly,ell <= 1 ==")
for cB in (1, 2, 3, 4):
    run(f"cubic cB={cB} (box small)", dict(P.CUBIC, cB=cB), box=small)
run("cubic cB=2 + third (box small)", dict(P.CUBIC, third=True), box=small)
