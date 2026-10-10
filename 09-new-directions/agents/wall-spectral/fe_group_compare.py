"""Compare FE groups: (a) our energy MDS (quadratic twists of the cubic theta, NOTES s.16);
(b) second moment of cubic Hecke L-functions sum_c L(s1,chi_c) L(s2,conj chi_c) N(c)^{-w} (DdFDS / DIPP family);
(c) second moment of quadratic twists of a GL2 form (Soundararajan-Young / Li)."""
import numpy as np
def closure(gens, maxdepth=14):
    seen = {tuple(np.eye(3, dtype=int).ravel())}; frontier = [np.eye(3, dtype=int)]; sizes = []
    for d in range(maxdepth):
        new = []
        for M in frontier:
            for g in gens:
                P = g @ M; t = tuple(P.ravel())
                if t not in seen: seen.add(t); new.append(P)
        frontier = new; sizes.append(len(seen))
        if not new: break
    return sizes
def refl(i, coeffs):   # x_i -> -x_i, y -> y + c x_i ;  coordinates (x1, x2, y)
    g = np.eye(3, dtype=int); g[i, i] = -1; g[2, i] = coeffs; return g
def refl_w(a1, a2):    # x1 -> x1 + a1 y, x2 -> x2 + a2 y, y -> -y
    g = np.eye(3, dtype=int); g[0, 2] = a1; g[1, 2] = a2; g[2, 2] = -1; return g
cases = {
 "(a) theta x quadratic twists (ours)": [refl(0, 2), refl(1, 2), refl_w(1, 1)],
 "(c) GL2 x quadratic twists (SY/Li)": [refl(0, 2), refl(1, 2), refl_w(1, 1)],
 "(b) cubic Hecke L second moment":    [refl(0, 1), refl(1, 1), refl_w(1, 1)],
 "(b2) same, conductor n1 n2^2 variant": [refl(0, 1), refl(1, 1), refl_w(1, 2)],
}
for name, gens in cases.items():
    s = closure(gens)
    print(f"{name:40s} group sizes by word length: {s}  -> {'FINITE, order %d' % s[-1] if s[-1] == s[-2] else 'infinite'}")
