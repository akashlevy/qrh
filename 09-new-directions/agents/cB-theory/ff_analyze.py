"""
Analyse ff_kubota output (cB-theory).  q = 7, sextic.  Normalised coefficients a(d) = S(d)/q^{d/2}
(the coefficient of X^d, X = q^{-u}, of sum g~ N^{-u}).  Prints:
  (1) the coprime parts T_j at h = T (index T^j) and the full series psi_l at index T^l (built by the exact
      p-part algebra), each multiplied by the normaliser 1 - q^3 X^6 (local factor of zeta(6u-2) at all
      places is 1/(1 - q^{3} X^6)?  we just test which multiplier kills high coefficients);
  (2) the same for h = T^2+1, T(T-1), T^2  (index h^j).
"""
import sys, numpy as np
q = 7
Z6 = np.exp(2j * np.pi / 6)
LOG = {pow(3, k, q): k for k in range(6)}
G = sum(Z6 ** LOG[y] * np.exp(2j * np.pi * y / q) for y in range(1, q))
fn = sys.argv[1] if len(sys.argv) > 1 else 'ff_D9.txt'
rows = [list(map(int, l.split())) for l in open(fn)]
D = max(r[1] for r in rows)
N = {}
for r in rows: N[(r[0], r[1])] = np.array(r[2:]).reshape(6, 6)


def S(h, j, d):
    if d < 0: return 0
    M = N[(h, d)]
    E, Sx = np.meshgrid(range(6), range(6), indexing='ij')
    return G ** d * np.sum(M * Z6 ** ((E - j * Sx) % 6))


def norm_series(h, j):
    return np.array([S(h, j, d) / q ** (d / 2) for d in range(D + 1)])


def Gq(j):
    j %= 6
    return -1 if j == 0 else sum(Z6 ** ((j * LOG[y]) % 6) * np.exp(2j * np.pi * y / q) for y in range(1, q))


def full_T(l):
    """full series at index T^l (all monic c incl. powers of T), normalised coefficients."""
    out = np.zeros(D + 1, complex)
    for d in range(D + 1):
        tot = 0
        for k in range(0, min(d, l + 1) + 1):
            if k <= l and k % 6 == 0:
                gk = 1 if k == 0 else q ** (k - 1) * (q - 1)
            elif k == l + 1:
                gk = q ** l * Gq(l + 1)
            else:
                continue
            tot += gk * (-1) ** ((d - k) * k) * S(1, l - 2 * k, d - k)
        out[d] = tot / q ** (d / 2)
    return out


def show(lab, a):
    print(f"  {lab:34s}" + " ".join(f"{abs(v):9.3g}" for v in a))


def times(a, poly):   # poly: dict power -> coeff
    b = np.zeros_like(a)
    for k, c in poly.items():
        b[k:] += c * a[:len(a) - k]
    return b


print(f"q = {q}, sextic, degrees d = 0..{D};  |normalised coefficient| of X^d")
print("  " + " " * 34 + " ".join(f"{d:9d}" for d in range(D + 1)))
Znorm = {0: 1, 6: -q ** 3}       # candidate global normaliser 1 - q^3 X^6
print("== h = 1 (index 1, all sf c)")
show("psi(1)", norm_series(0, 0))
show("psi(1) * (1 - q^3 X^6)", times(norm_series(0, 0), Znorm))
print("== full series at index T^l (built from coprime parts by exact p-part algebra), times (1 - q^3 X^6)")
for l in range(6):
    show(f"psi(T^{l}) (1-q^3X^6)", times(full_T(l), Znorm))
print("== coprime parts at h = T, index T^j:  T_j * (1 - q^3 X^6)")
for j in range(6):
    show(f"T_{j} (1-q^3X^6)", times(norm_series(1, j), Znorm))
for h, lab in [(2, "T^2+1"), (3, "T(T-1)"), (4, "T^2")]:
    print(f"== coprime parts at h = {lab}, index h^j,  times (1 - q^3 X^6)")
    for j in range(6):
        show(f"[{lab}]^{j} coprime (1-q^3X^6)", times(norm_series(h, j), Znorm))
