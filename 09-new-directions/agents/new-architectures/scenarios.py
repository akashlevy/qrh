"""Scenario table for the architecture template (GRH for rows).  Run: python3 scenarios.py"""
from fractions import Fraction as Fr
from template_lp import solve

def show(name, over):
    r = solve(over)
    if r is None:
        print(f"{name:62s}  infeasible (no sigma < 2)")
        return
    th = r["theta"]
    print(f"{name:62s}  theta={th:.6f} ~ {Fr(th).limit_denominator(400)}  "
          f"lx={r['lx']:.3f} ly={r['ly']:.3f} ell={r['ell']:.3f} M={r['M']:.3f} h={r['h']:.3f}")

if __name__ == "__main__":
    print("== sensitivity of the GRH floor to each scaling (others at P78 values) ==")
    show("P78 base (c=a=gb=1/6, rho=1/2, kU=1, r=2)", {})
    for a in (0.2, 0.25, 1/3):
        show(f"amplification a={a:.3f}", dict(a=a))
    for rho in (0.6, 2/3, 0.75, 1.0):
        show(f"slot suppression rho={rho:.3f}", dict(rho=rho))
    for kU in (0.9, 0.75, 0.5, 0.0):
        show(f"row-range growth kU={kU:.2f}", dict(kU=kU))
    for c in (0.2, 0.25, 1/3):
        show(f"copy exponent c={c:.3f} (a=gb=c, GL2 wall)", dict(c=c, a=c, gb=c))
    for c in (0.2, 0.25):
        show(f"copy exponent c={c:.3f} only in signal/rows (a=gb=1/6)", dict(c=c))
    show("no wall (energy diagonal for all M,ell), caseA kept", dict(wE=10))
    show("no wall, no case A", dict(wE=10, caseA="none"))
    show("no wall, no case A, no Gram", dict(wE=10, caseA="none", gb=0, gq=0))
    print()
    print("== Idea 1: higher-rank dual objects (degree r reflection), optimistic: no caseA/third, margin 0 ==")
    for (lab, c, r) in [("sextic twists, GL3^(3)-type dual (r=3)", 1/6, 3),
                        ("sextic twists, r=4", 1/6, 4),
                        ("cubic twists -> sextic Gauss sums, GL5^(6) (r=5)", 1/3, 5),
                        ("cubic twists -> sextic Gauss sums, GL6^(6) (r=6)", 1/3, 6),
                        ("order-10 twists -> quintic, GL4^(5) (r=4)", 1/10, 4),
                        ("order-10 twists -> quintic, GL5^(5) (r=5)", 1/10, 5),
                        ("order-10 twists, hypothetical GL2 quintic theta (r=2)", 1/10, 2),
                        ("cubic twists, hypothetical GL2 dual (r=2)", 1/3, 2),
                        ("cubic twists -> sextic Kubota series, shift h^4 (cond N(h)^4, r=4)", 1/3, 4),
                        ("quartic twists -> conj quartic Kubota (cond N(h)^3, r=3), optimistic LS", 1/4, 3)]:
        show(lab, dict(c=c, a=c, gb=c, r=r, caseA="none", third=False, margin=0.0))
    print()
    print("== Idea 4: other placements of the slot amplifier (GRH) ==")
    show("quadratic m-twist (m/p)_2, weights conj gamma_3(p), grouped with B: P78 wall, no case A",
         dict(caseA="none"))
    show("sextic m-twist chi_p(m), grouped with B: reflected dual needs BGL, wall 3M/2+ell<=1",
         dict(wE=10, caseA="none", third=False, extra=[([1.5, 1.5, 1], -1.0, "le")]))
    import aside_marks
    print(f"{'m-twist grouped with A (Gram (YP^2/X)^(1/6) eats P^(1/6)); see aside_marks.py':62s}  ",
          aside_marks.solve())
    show("s-marks (p|s, Y->YP) = reparametrization Y->YP: base", {})
    show("unmarked probe (ell = 0 forced), for reference",
         dict(margin=0, extra=[([0, 0, 1], 0.0, "le")]))
    print()
    print("== P2-type mean square (paper2 architecture), order-N twists, degree-r dual ==")
    print("   need sum_{h<=Hd}|B_h|^2 << D^2 with Hd = D^2/H and dual length Hd^r/D (optimistic")
    print("   diagonal large sieve for the reflected dual): Hd <= D^{2/r}, i.e. H >= D^{2-2/r};")
    print("   extraction with H^{1/N} copies:  theta = 1/2 + (1-1/N)(1-1/r)")
    for (lab, N, r) in [("sextic -> cubic theta GL2^(3) [paper2]", 6, 2),
                        ("sextic -> GL3^(3) first-row (if Gauss sums)", 6, 3),
                        ("quartic -> conj. quartic Gauss (Kubota, cond N(h)^3) [refuted transplant]", 4, 3),
                        ("cubic -> sextic Gauss, shift h^4 (Kubota, cond N(h)^4)", 3, 4),
                        ("cubic -> sextic Gauss on GL5^(6) theta", 3, 5),
                        ("order 10 -> quintic Gauss on GL4^(5) theta", 10, 4),
                        ("order 10 -> hypothetical GL2 quintic theta", 10, 2),
                        ("cubic -> hypothetical GL2 object with sextic Gauss coefficients", 3, 2),
                        ("quartic -> hypothetical GL2 theta", 4, 2)]:
        th = Fr(1, 2) + (1 - Fr(1, N)) * (1 - Fr(1, r))
        print(f"   {lab:74s} theta = {th} = {float(th):.4f}")

    print()
    print("== generic wall-free device: low gain gL per unit length, row cost gR (device.py) ==")
    import device
    for k in (2, 3, 4, 14):
        x = device.solve_dev(1/12, k/12)
        print(f"   gR/gL = {k:3d}:  theta = {x[-1]:.5f}  (device length used {x[3]:.3f})")
