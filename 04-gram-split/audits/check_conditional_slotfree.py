"""Pure rational certificate for the note's CONDITIONAL 0.8665 geometry.

This verifies exponent arithmetic only. In particular it does not establish
the theta-pairing estimate V' used on the low side. The old Voronoi proof
of that estimate used the wrong coefficient sequence.

Instead of floating-point optimization and rectangle subdivision, use the
exact zero-slot row count R(delta)=(15-16*delta)/(15-6*delta). Its denominator
is positive on the entire allowed interval, so each high-side inequality
reduces to the maximum of one concave quadratic with rational coefficients.
"""

from fractions import Fraction as F

lx, ly, ell, sigma = F(989, 5000), F(3141, 5000), F(0), F(1733, 2000)
h, M, b = 1 - lx, lx + ly, ly - lx
delta0, delta1, dmin, zeta, z0 = F(1, 50), F(5, 6), F(1, 100), F(1, 10000), F(17, 50)

def report(label, value, relation="negative"):
    ok = value < 0 if relation == "negative" else value <= 0
    print(f"{'PASS' if ok else 'FAIL'} {label}: {value} = {float(value):.12g}")
    assert ok

# Assumed clean Gram split plus the UNPROVED large-sieve-shaped V' estimate.
g = max(F(0), b / 6 + max(M - 1, 7 * M / 3 - 2, -M / 6))
e = max(F(0), (1 - 2 * M) / 4, M - 1, (M - 2) / 3)
low = F(5, 6) + lx / 6 + (g + e) / 2
report("conditional low exponent minus sigma", low - sigma, "nonpositive")
assert lx > ell and ly > ell

small = h * (z0 - F(1, 6)) - ly / 2 + F(63, 50) * dmin
report("small rows", small)
floor = (1 + delta0) * (1 - ly) / 2 - sigma + (h + zeta) * (1 + delta0 / 2 - F(1, 6))
report("conservative floor bin", floor)

for d in (dmin, h + zeta):
    # E(delta,d)=C+A*delta+d*(15-16*delta)/(15-6*delta).
    C = (1 - ly) / 2 - sigma + h * (z0 - F(1, 6)) - d * z0
    A = (1 - ly + d) / 2
    c0, c1, c2 = 15 * (C + d), 15 * A - 6 * C - 16 * d, -6 * A
    assert c2 < 0 and 15 - 6 * delta1 > 0
    vertex = -c1 / (2 * c2)
    points = [delta0, delta1]
    if delta0 <= vertex <= delta1:
        points.append(vertex)
    peak = max(c0 + c1 * x + c2 * x * x for x in points)
    report(f"high-side numerator at row endpoint d={d}", peak)

print("PASS continuous high-side range: affine in row scale, no amplitude dependence at ell=0")
print("CONDITIONAL arithmetic certificate complete; no analytic V' estimate proved.")
