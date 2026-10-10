"""Lemma C0 check: for squarefree primary A and any m,
   c_A chi_A(m) = mu(A) * q_A^{-1/2} * sum_{y mod A} conj(chi_A(y)) e(-m y / A),
   c_A = gamma_2(A) conj(alpha(A)) conj(G(A)),  G = conj(chi_A(4)) gamma_3(A)   (ray factor nu omitted).
i.e. the unreflected completed row B_m = sum_A c_A chi_A(m) V/sqrt(q_A) is a Moebius-twisted conjugate
sextic Kubota sum with numerator -m.  Brute force over O/A (zomega conventions)."""
import sys, numpy as np, itertools
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split"); sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
from a2lib import *
P = [p for p in primes_upto(400) if p.kind == "split"][:8]
def G(Aa, Ab, prs, j, ha, hb):
    N = norm(Aa, Ab); g = np.gcd(Aa, Ab); g = abs(int(g)) if g else 1
    x, y = np.meshgrid(np.arange(N // g, dtype=np.int64), np.arange(g, dtype=np.int64), indexing="ij"); x = x.ravel(); y = y.ravel()
    tot = np.zeros(len(x), np.int64); ok = np.ones(len(x), bool)
    for p in prs:
        s = p.symbol(x, y); ok &= s >= 0; tot += np.where(s >= 0, s, 0)
    va, vb = mul(x, y, ha, hb)                      # h*v
    ca, cb = mul(va, vb, Aa - Ab, -Ab)              # (h v) conj(A)
    ph = np.exp(2j * np.pi * (cb % N) / N)
    return np.sum(np.where(ok, Z6[(j * tot) % 6], 0) * ph) / np.sqrt(N)
def chi(prs, a, b):
    e = 0
    for p in prs:
        s = sym1(p, a, b)
        if s < 0: return 0
        e += s
    return Z6[e % 6]
worst = 0; n = 0
for r in (1, 2, 3):
    for prs in itertools.combinations(P, r):
        Aa, Ab = 1, 0
        for p in prs: Aa, Ab = mul(Aa, Ab, p.a, p.b)
        Aa, Ab = int(Aa), int(Ab)
        mu = (-1) ** r
        g1 = G(Aa, Ab, prs, 1, 1, 0); g2 = G(Aa, Ab, prs, 2, 1, 0); g3 = G(Aa, Ab, prs, 3, 1, 0)
        al = cval(Aa, Ab); al = al / abs(al)
        Gc = np.conj(chi(prs, 4, 0)) * g3
        cA = g2 * np.conj(al) * np.conj(Gc)
        for (ma, mb) in [(1, 0), (2, 1), (5, -3), (7, 2), (-4, 9), (prs[0].a, prs[0].b), (3 * prs[0].a, 3 * prs[0].b)]:
            lhs = cA * chi(prs, ma, mb)
            rhs = mu * G(Aa, Ab, prs, 5, -ma, -mb)
            worst = max(worst, abs(lhs - rhs)); n += 1
print("Lemma C0: %d (A,m) pairs (A = 1, 2 or 3 split primes, incl. m divisible by a prime of A); max |lhs-rhs| = %.2e" % (n, worst))
