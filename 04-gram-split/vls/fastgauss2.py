"""Chunked, cache-friendly version of D_numerics/fastgauss.gauss_split (same output; checked in E_primes2.py).
g(pi,1) = sum_{k=0}^{p-2} zeta6^{k jg} exp(-2 pi i b g^k / p),  chi_pi(g) = zeta6^{jg}.
Sum split by k mod 6 (Gaussian periods), chunks of C = 6*2^13 consecutive k."""
import numpy as np, sys
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/D_numerics")
from fastgauss import prim_root

C = 6 * 8192


def gauss_split2(args):
    a, b, p, r = args
    g = prim_root(p)
    z = (1 + r) % p
    roots = [pow(z, j, p) for j in range(6)]
    jg = roots.index(pow(g, (p - 1) // 6, p))
    n = p - 1
    # base table g^0..g^{C-1}
    m = min(C, n)
    base = np.ones(m, np.int64); x = 1
    # build by doubling within the chunk
    L = 1
    while L < m:
        gl = pow(g, L, p)
        take = min(L, m - L)
        base[L:L + take] = (base[:take] * gl) % p
        L += take
    bb = (b * base) % p                              # b g^k for k in chunk 0
    gC = pow(g, m, p)
    S = np.zeros(6, complex)
    two_pi_over_p = 2 * np.pi / p
    k0 = 0; cur = bb
    while k0 < n:
        cnt = min(m, n - k0)
        ang = cur[:cnt] * two_pi_over_p
        c = np.cos(ang); s = np.sin(ang)
        full = (cnt // 6) * 6
        # k0 is a multiple of m, m a multiple of 6 (or the whole range) -> class of index i is (k0 + i) mod 6 = i mod 6
        S += (c[:full].reshape(-1, 6).sum(0) - 1j * s[:full].reshape(-1, 6).sum(0))
        for i in range(full, cnt):
            S[(k0 + i) % 6] += c[i] - 1j * s[i]
        k0 += cnt
        cur = (cur * gC) % p
    jj = (np.arange(6) * jg) % 6
    return np.sum(np.exp(1j * np.pi / 3 * jj) * S) / np.sqrt(p)
