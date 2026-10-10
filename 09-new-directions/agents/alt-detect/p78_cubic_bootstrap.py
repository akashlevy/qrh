"""P78-shape cubic architecture where the B-row energy carries the normaliser term of
1/L(3s-1,(alpha nu)^3): every row B_m of length ZP' contains (ZP')^{(1+beta*)/3} coherent terms,
so the normalised energy excess is >= (2 beta* - 1)(1 + l')/3  (in addition to the node bound).
GRH for the rows 1/L(x, eta psi_u) is assumed (floor); beta* = zero-free abscissa for L(s,(alpha nu)^3).
We report theta(beta*) and the bootstrap fixed point theta = theta(theta)."""
import numpy as np
from fractions import Fraction as Fr
import p78_cubic_lp as P
Z3 = P.Z3

def make_nodes(beta, cB):
    def nodes(Mp, lp, cB_, k):
        base = P.energy_nodes.__wrapped__(Mp, lp, cB, k) if hasattr(P.energy_nodes, "__wrapped__") else ORIG(Mp, lp, cB, k)
        pole = ((2 * beta - 1) / 3 * lp, (2 * beta - 1) / 3)
        return [n + [pole] for n in base]
    return nodes

ORIG = P.energy_nodes
def theta_of(beta, cB=2, box=[(0, 1), (0, 1), (0, 1), (None, None)]):
    P.energy_nodes = make_nodes(beta, cB)
    orig_lp = P.linprog
    def lp2(*a, **kw):
        kw["bounds"] = box; return orig_lp(*a, **kw)
    P.linprog = lp2
    r = P.solve(dict(P.CUBIC, cB=cB))
    P.energy_nodes = ORIG; P.linprog = orig_lp
    return r

for cB in (1, 2):
    print(f"cB = {cB}")
    for beta in (0.5, 0.6, 0.7, 0.8, 0.875, 0.9, 1.0):
        r = theta_of(beta, cB)
        print(f"   beta*={beta:.3f}: theta_GRH-rows = {r['theta']:.5f}  (lx={r['lx']:.3f}, ly={r['ly']:.3f}, ell={r['ell']:.3f})")
    b = 1.0
    for it in range(60):
        b = theta_of(b, cB)["theta"]
    print(f"   bootstrap fixed point (GRH for rows, normaliser zeros bootstrapped): {b:.6f}")
