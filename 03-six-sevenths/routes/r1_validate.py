import sys, math, cmath, random
import numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/strengthening")
from eisenstein import eisenstein_primes, ResidueField, gauss_prime, primary, norm, cx, mul, Z6, e_add
from r1_fake_diagonal import sextic_exponents
random.seed(1)
primes = [p for p in eisenstein_primes(400) if p[2] != 2]
fields = {pr[1]: ResidueField(*pr) for pr in primes}
# 1. vectorised sextic exponents vs reference (split and inert)
ra = np.array([random.randint(-300, 300) for _ in range(3000)], dtype=np.int64)
rb = np.array([random.randint(-300, 300) for _ in range(3000)], dtype=np.int64)
bad = 0
for p, F in fields.items():
    v = sextic_exponents(F, ra, rb)
    for i in range(0, 3000, 37):
        k = F.sextic_exp((int(ra[i]), int(rb[i])))
        bad += (v[i] != (-1 if k is None else k))
print("sextic exponent mismatches:", bad, " (inert primes tested:", sum(F.kind=='inert' for F in fields.values()), ")")
# 2. gamma_2 of composite n = p1 p2 (split, distinct norms): brute force over Z/N(n) vs twisted multiplicativity
split = [p for p, F in fields.items() if F.kind == "split"]
def chi_exp_n(fac, x):
    tot = 0
    for p in fac:
        k = fields[p].sextic_exp(x)
        if k is None: return None
        tot += k
    return tot % 6
errs = []
for trial in range(12):
    p1, p2 = random.sample(split, 2)
    if fields[p1].q == fields[p2].q: continue
    n = mul(p1, p2); Nn = norm(n); nc = cx(n)
    s = 0
    for x in range(Nn):
        k = chi_exp_n((p1, p2), (x, 0))
        if k is None: continue
        s += Z6[(2 * k) % 6] * e_add(x / nc)
    brute = s / math.sqrt(Nn)
    g1 = gauss_prime(fields[p1], 2); g2 = gauss_prime(fields[p2], 2)
    k = fields[p2].sextic_exp(p1)          # chi_{p2}(p1)
    tm = g1 * g2 * Z6[(4 * k) % 6]
    errs.append(abs(brute - tm))
print("max |gamma_2(p1 p2) brute - twisted-mult| =", max(errs), " |brute| ~", abs(brute))
# 3. paper identity mu(n) gamma_{-1}(n) = chi_n(-1) G(n)^{-1} conj(alpha(n)) gamma_2(n): check |.| and gamma_1 gamma_2 = -alpha G at primes
