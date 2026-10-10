"""
Route 1 check: are quadratic twists of the cubic theta coefficients Fourier coefficients of the
sextic Kubota Eisenstein series (up to a Groessencharacter)?  We verify numerically, for squarefree
n built from split primes with distinct norms (so Z[w]/n = Z/N(n)), the identity

   beta(n) psi_k(n) chi_P(n)^{-2} = mu(n) * G(n) * r(n,k) * conj( gamma_1(n; k^3 P^4) )      (*)

where beta(n) = gamma_2(n) conj(alpha(n))  (P78 theta coefficient, no Moebius),
      psi_k(n) = chi_k(n)^3, gamma_j(n;h) = N(n)^{-1/2} sum_{x mod n} chi_n(x)^j e(h x/n),
      G(n) = conj(chi_n(4)) Gamma(n), r(a,b) = Gamma(ab)/(Gamma(a)Gamma(b))  (finite ray data mod 12).
gamma_1(n; h) is the sextic Gauss sum = numerator of the h-th Fourier coefficient of the
sextic Kubota Eisenstein series.  The factor mu(n) is NOT a Groessencharacter.
"""
import sys, math, cmath, itertools
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/energy-wall")
from eisenstein import eisenstein_primes, ResidueField, mul, cx, norm, Z6, primary

W = complex(-0.5, math.sqrt(3) / 2)
def e_add(z): return cmath.exp(4j * math.pi * z.imag / math.sqrt(3))
def Gamma(c):
    a, b = c
    return (1 + 1j ** (-b % 4) + 1j ** (a % 4) + 1j ** ((b - a) % 4)) / 2
def chi_exp(F, e):           # sextic exponent of (e / p)_6, None if p | e
    return F.sextic_exp((e[0] % F.q, e[1] % F.q)) if F.kind == "split" else F.sextic_exp(e)
def chi(fields, den, num):   # (num / den)_6 for den = list of prime factors (fields)
    v = 1 + 0j
    for F in den:
        k = chi_exp(F, num)
        if k is None: return 0
        v *= Z6[k]
    return v
def epow(e, k):
    r = (1, 0)
    for _ in range(k): r = mul(r, e)
    return r

pr = [x for x in eisenstein_primes(400) if x[0] == "split"]
# one prime per rational norm (distinct norms => Z[w]/n cyclic of order N(n))
seen = set(); P1 = []
for kind, p, q in pr:
    if q not in seen: seen.add(q); P1.append(ResidueField(kind, p, q))
def gam(Fs, j, h):
    """gamma_j(n; h) for n = prod of the primes in Fs (split, distinct norms)."""
    n = (1, 0)
    for F in Fs: n = mul(n, F.p)
    Nn = norm(n); nc = cx(n)
    s = 0
    for t in range(Nn):
        c = chi(None, Fs, (t, 0))
        if c == 0: continue
        s += c ** j * e_add(complex(*(0, 0)) + (cx(mul(h, (t, 0)))) / nc)
    return s / math.sqrt(Nn), n

fails = 0; tests = 0; fails_nomu = 0
rows = []
for r in (1, 2, 3):
    for Fs in itertools.combinations(P1[:7], r):
        if math.prod(F.q for F in Fs) > 60000: continue
        g1, n = gam(Fs, 1, (1, 0)); g2, _ = gam(Fs, 2, (1, 0))
        mu = (-1) ** r; al = cx(n) / abs(cx(n))
        Gn = (1 / chi(None, Fs, (4, 0))) * Gamma(n)
        # (a) paper identity gamma_1 gamma_2 = mu alpha G
        ok_a = abs(g1 * g2 - mu * al * Gn) < 1e-6
        beta = g2 * al.conjugate()
        # (b) twisted identity for a few k, P coprime to n
        others = [F for F in P1[7:14]]
        for kF in others[:3]:
            for PF in others[3:6]:
                k, P = kF.p, PF.p
                psi = chi(None, [kF], n) ** 3
                chiP = chi(None, [PF], n)
                lhs = beta * psi * chiP ** (-2)
                h = mul(epow(k, 3), epow(P, 4))
                g1h, _ = gam(Fs, 1, h)
                rr = Gamma(mul(n, k)) / (Gamma(n) * Gamma(k))
                rhs = mu * Gn * rr * g1h.conjugate()
                tests += 1
                fails += abs(lhs - rhs) > 1e-4   # float error grows with |h|
                fails_nomu += abs(lhs - rhs / mu) > 1e-4
        rows.append((r, norm(n), ok_a))
print("paper identity gamma1*gamma2 = mu*alpha*G at", len(rows), "squarefree n:", all(x[2] for x in rows))
print(f"twisted identity (*) : {tests} tests, failures {fails}")
print(f"same identity WITHOUT mu(n): failures {fails_nomu} (all tests with omega(n) odd fail)")
