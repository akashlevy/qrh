"""Variant scans for the k-Poisson probes (MODEL).  Output: runs/variants.out"""
import itertools
from fractions import Fraction as Fr
import kvar_lp as K

def run(lab, kind, **over):
    c = K.default_cfg(kind); c.update(over)
    s = K.solve(c)
    if s is None:
        print(f"{lab:70s} infeasible"); return None
    print(f"{lab:70s} theta={s['theta']:.6f} ~{Fr(s['theta']).limit_denominator(400)}  "
          f"l={[round(v,3) for v in s['l']]} y={[round(v,3) for v in s['y']]} ell={s['ell']:.3f} M={s['M']:.3f} nodes={s['nodes']}")
    return s

print("== reference: one Poisson variable (k=1, cubic copies, alt-detect template) ==")
for cB in (2, 3, 4):
    run(f"k=1 cubic, row theta(x)theta, cB={cB}" + ("  [ACTUAL, cB-theory]" if cB == 4 else "  [hypothetical]"), "cubic1", cB=cB, ymax=3)
run("k=1 cubic, cB=2, no marks", "cubic1", marks=False, ymax=3)
run("k=1 cubic, perfect energy (no wall/excess), marks", "cubic1", Eperfect=True, third=False, ymax=3)
run("k=1 cubic, perfect energy, no marks", "cubic1", Eperfect=True, third=False, marks=False, ymax=3)

print("\n== k=2 (Poisson, Poisson, theta_3): single cubic theta row, cB=2 ==")
run("k=2 default (third, marks, comp, partial rows)", "cubic2")
run("k=2 reflection node only (energy below the wall M+2l/3<=2/3)", "cubic2", node_filter=[2])
run("k=2 reflect / reflect+Poisson nodes only", "cubic2", node_filter=[2, 3])
run("k=2 no third piece", "cubic2", third=False)
run("k=2 + gq Gram piece", "cubic2", gq=1.0)
run("k=2 no marks", "cubic2", marks=False)
run("k=2 marks forced ell=0.1", "cubic2", fix={"e": 0.1})
run("k=2 marks forced ell=0.2", "cubic2", fix={"e": 0.2})
run("k=2 symmetric geometry", "cubic2", sym=True)
run("k=2 perfect energy", "cubic2", Eperfect=True, third=False)
run("k=2 perfect energy, no marks", "cubic2", Eperfect=True, third=False, marks=False)
run("k=2 Gram exceptional x0.5", "cubic2", gb=[1/6, 1/6])
run("k=2 Gram exceptional removed (hypothetical; it is the signal's shadow)", "cubic2", gb=[0, 0])
run("k=2 Gram exc. removed on var 2 only", "cubic2", gb=[1/3, 0])
print("-- (a) unequal Y1 != Y2 / (b) short second variable --")
for y2 in (0.0, 0.1, 0.25, 0.5):
    run(f"k=2 y2 fixed = {y2}", "cubic2", fix={"y2": y2}, margin_y=0.0)
for l2 in (0.0, 0.25, 0.5, 0.75, 1.0):
    run(f"k=2 l2 fixed = {l2}", "cubic2", fix={"l2": l2})
run("k=2 Q2 tiny (l2=y2=0.02)", "cubic2", fix={"l2": 0.02, "y2": 0.02}, margin_y=0.0)
print("-- (c) m2 smooth, Poisson-transformed, no s2-polynomial (Y2=1) --")
run("k=2 y2 = 0 (A' = smooth weight)", "cubic2", fix={"y2": 0.0}, margin_y=0.0)
print("-- (d) one Poisson + Gauss-sum orthogonality (m2 complete: F2=1, s2=A, Y2=ZP) --")
run("k=2 l2 = y2 = 1 (F2 = 1, orthogonality), ell=0", "cubic2", fix={"l2": 1.0, "y2": 1.0, "e": 0.0})
run("   same, no third piece", "cubic2", fix={"l2": 1.0, "y2": 1.0, "e": 0.0}, third=False)
run("   same, perfect energy", "cubic2", fix={"l2": 1.0, "y2": 1.0, "e": 0.0}, Eperfect=True, third=False)

print("\n== k=3 (Poisson, Poisson, Poisson): Hecke row lambda-bar eta, degree 1, cB=1 ==")
run("k=3 default", "cubic3")
run("k=3 no marks", "cubic3", marks=False)
run("k=3 perfect energy", "cubic3", Eperfect=True)
run("k=3 Gram removed + perfect energy (= low-side base)", "cubic3", Eperfect=True, gb=[0, 0, 0])

print("\n== sextic analogues (copy exponent 1/ord(chi^e) per variable) ==")
# k=3 Hecke rows, triples a+b+c = 0 mod 6
run("sextic (1,1,4) three Poisson, Hecke row", "cubic3", c=[1/6, 1/6, 1/3])
run("sextic (2,2,2) three Poisson, Hecke row (= cubic k=3)", "cubic3", c=[1/3, 1/3, 1/3])
run("sextic (1,2,3) three Poisson, Hecke row", "cubic3", c=[1/6, 1/3, 1/2])
run("sextic (1,2)+ray: two Poisson, Hecke*ray row", "cubic2", c=[1/6, 1/3], cB=1, third=False)
# k=2 theta rows with a quadratic residual (thin): (1,3,2)
run("sextic (1,3 | 2): Poisson e=1, Poisson e=3, cubic theta, quadr. residual", "cubic2",
    c=[1/6, 1/2], ls=2, caseA="paper", gq=1.0, third=True)
run("   same, var-2 forced trivial (l2=1, y2=0)", "cubic2",
    c=[1/6, 1/2], ls=2, caseA="paper", gq=1.0, third=True, fix={"l2": 1.0, "y2": 0.0}, margin_y=0.0)
