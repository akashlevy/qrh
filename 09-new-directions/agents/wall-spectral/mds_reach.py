"""
Affine MDS bookkeeping for E (Route 1).  Coordinates x_i = s_i - 1/2, y = w - 1/2 of
Z(s1,s2,w) = sum_k q_k^{-w} D_{P1}(s1,k) conj D_{P2}(s2,k),  D_P(s,k) = sum_n tau(n) psi_k(n) chi_P(n)^{-2} q_n^{-s}.
FEs: s1 (theta twist, conductor (q_k q_P)^2, dual supported on P|mu => level cost L^{-x_1}),
     s2 likewise, s3 (k-sum = L(w, chi_{n1 n2}), conductor n1 n2).
(1) the group preserves I = x1 + x2 + y  (so Bochner hulls never lower inf I);
(2) orbit of the diagonal polar divisor {y = 1/2}: linear forms  y o g  = a*I + (finite root);
(3) per-(P1,P2) contour bound: error exponent at p=(x,x,y) after optimal reflections, compared with K*N.
"""
import numpy as np, itertools
s1 = np.array([[-1,0,0],[0,1,0],[2,0,1]]); s2 = np.array([[1,0,0],[0,-1,0],[0,2,1]])
s3 = np.array([[1,0,1],[0,1,1],[0,0,-1]]); gens = [s1, s2, s3]
I = np.array([1,1,1])
print("(1) I o g == I for generators:", all((I @ g == I).all() for g in gens))
# (2) orbit of the form y under right multiplication (phi -> phi o g)
forms = {(0,0,1)}; frontier = [np.array([0,0,1])]
for depth in range(8):
    new = []
    for f in frontier:
        for g in gens:
            h = tuple(f @ g)
            if h not in forms: forms.add(h); new.append(np.array(h))
    frontier = new
dec = {}
for f in sorted(forms, key=lambda t: sum(map(abs, t))):
    f = np.array(f)
    # write f = a*I + r with r in the finite root set {+-y, +-(x1-x2)..., +-x1+..}: choose a = min coordinate
    a = int(min(f)); r = tuple(f - a * I)
    dec.setdefault(r, []).append(a)
print("(2) polar forms y o g = a*I + r  (r: finite part):")
for r, al in sorted(dec.items()):
    print(f"     r = {r}:  a in {sorted(set(al))[:12]}")
print("     divisor {a*I + r(v) = 1/2}: for |a| -> inf these accumulate at I = 0 (natural boundary).")

# (3) exponent bookkeeping at the target geometry (energy-wall s.1, wE = 0.064)
M, ell = 0.3662 + 0.4936, 0.2042
kap, lam = M, ell; nu = 2 * M + ell - 1; w = nu - kap
print(f"(3) geometry: kappa=M={kap:.4f}, nu=2M+l-1={nu:.4f}, w={w:.4f}, l={lam:.4f}")
def excess(x, y, cprime, pair):
    """log_Z(error / (K N)) for a contour point (x, x, y) with level-cost c' at the reflected point."""
    c = 2 * max(0.0, -x) + cprime           # each reflection of a negative x_i costs L^{|x_i|}
    return kap * (0.5 + y) + nu * (1 + 2 * x) + lam * ((1 if pair else 0) + c) - (kap + nu)
for pair in (False, True):
    best = min(((excess(x, y, 0.0, pair), x, y) for x in np.linspace(-1, 1, 401) for y in np.linspace(-1, 1, 401)
                if 2 * x + y >= 0.5 - 1e-12))      # large-sieve / Lindelof-on-average region I >= 1/2
    print(f"    {'cross pairs P1!=P2' if pair else 'diagonal P1=P2   '}: min excess over I>=1/2 with c'=0: {best[0]:+.4f} at x={best[1]:+.3f}, y={best[2]:+.3f}")
    Ineed = 0.5 - lam * (1 if pair else 0) / kap
    print(f"      needs a bound at (0,0,I) with I <= {Ineed:.3f}   (I-invariance: unreachable by FE + convexity from I >= 1/2)")
