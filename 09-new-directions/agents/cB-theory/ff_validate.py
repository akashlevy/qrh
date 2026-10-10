"""
Function-field sanity check of the Gauss-sum formula used by ff_kubota.c (cB-theory).
F_7[T], sextic symbol (x/c)_6 = iota(Res(c,x)) with iota(3^k) = exp(2 pi i k/6)  (q = 7, (q-1)/6 = 1),
additive character psi_c(x) = e_7(coefficient of T^{d-1} of x mod c)   (= e_7(res_infty-coefficient of x/c)).
Claim (HD + Stickelberger, derived in NOTES.md):  for monic squarefree c of degree d with omega(c) prime factors
      g(1,c) = (-1)^(d - omega(c)) * G^d * iota(Res(c, c')),   G = sum_{y in F_7^x} iota(y) e_7(y).
Brute force over all x mod c for d <= 4.
"""
import itertools, numpy as np
from sympy import Poly, symbols, GF, factor_list
T = symbols('T')
p = 7
LOG = {pow(3, k, p): k for k in range(6)}
Z6 = np.exp(2j * np.pi / 6)
iota = lambda y: Z6 ** LOG[y % p]
G = sum(iota(y) * np.exp(2j * np.pi * y / p) for y in range(1, p))


def res(c, x):
    """Res(c, x) = prod_{c(b)=0} x(b) for monic c (sympy resultant convention matches for monic c)."""
    return int(Poly(c, T, modulus=p).resultant(Poly(x, T, modulus=p))) % p


def brute(ccoef):
    d = len(ccoef) - 1
    c = Poly(ccoef, T, modulus=p)
    tot = 0j
    for xs in itertools.product(range(p), repeat=d):
        x = Poly(list(xs), T, modulus=p)          # degree < d, leading coeff xs[0] is T^{d-1}
        r = int(c.resultant(x)) % p if not x.is_zero else 0
        if r == 0: continue
        tot += iota(r) * np.exp(2j * np.pi * xs[0] / p)
    return tot


rng = np.random.default_rng(1)
bad = 0; n = 0
for d in (1, 2, 3, 4):
    for trial in range(6 if d < 4 else 3):
        while True:
            cc = [1] + list(rng.integers(0, p, d))
            c = Poly(cc, T, modulus=p)
            if c.discriminant() % p != 0: break
        om = len(factor_list(c.as_expr(), modulus=p)[1])
        R = res(c.as_expr(), c.diff(T).as_expr())
        pred = (-1) ** (d - om) * G ** d * iota(R)
        bf = brute(cc)
        n += 1; bad += abs(bf - pred) > 1e-6
        print(f"d={d} c={cc} omega={om}  brute={bf:.4f}  formula={pred:.4f}")
print(f"{n - bad}/{n} agree")
