"""Fast normalized sextic Gauss sums gtilde(pi,1) for split primary primes, via primitive-root power tables.
Checked against zomega.Prime.gauss on small primes."""
import numpy as np, sys
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *

def prim_root(p):
    fac = []; n = p - 1; d = 2
    while d * d <= n:
        if n % d == 0:
            fac.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: fac.append(n)
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in fac): return g

def powtable(g, p):
    pw = np.ones(1, np.int64); n = 1
    while n < p - 1:
        gn = pow(g, n, p)
        pw = np.concatenate([pw, (pw * gn) % p]); n *= 2
    return pw[:p - 1]

def gauss_split(args):
    a, b, p, r = args
    g = prim_root(p)
    z = (1 + r) % p
    roots = [pow(z, j, p) for j in range(6)]
    jg = roots.index(pow(g, (p - 1) // 6, p))
    pw = powtable(g, p)
    k = np.arange(p - 1, dtype=np.int64)
    ph = np.exp(-2j * np.pi * ((b * pw) % p) / p)
    return np.sum(Z6[(k * jg) % 6] * ph) / np.sqrt(p)

if __name__ == "__main__":
    P = primes_upto(5000)
    err = 0
    for Pr in P[:200]:
        if Pr.kind != "split": continue
        g1 = Pr.gauss() / np.sqrt(Pr.N); g2 = gauss_split((Pr.a, Pr.b, Pr.q, Pr.r))
        err = max(err, abs(g1 - g2))
    print("max err vs zomega:", err)
    # conjugate relation check
    d = {}
    for Pr in P[:200]:
        if Pr.kind == "split": d.setdefault(Pr.q, []).append(Pr)
    e2 = 0
    for p, L in d.items():
        if len(L) == 2:
            e2 = max(e2, abs(L[0].gauss() - np.conj(L[1].gauss())))
    print("max |g(pi) - conj g(pibar)|:", e2)
