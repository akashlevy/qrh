"""
(b) "Moments of different lengths + near-orthogonality" (OpenAI cubic first moment, dual.tex,
Prop. prop:dual + Lemma lem:character-gram) transplanted to the n-th order family.

Normalisation (log_U units): family of U rows (n-th-power-free u, q_u <= U), a polynomial F of
length Z = U^z with divisor-bounded coefficients, normalised |F~|^2 = |F|^2 / Z.
On a row set of size U^rho the values satisfy |F~| = U^W (W = normalised excess).

Inputs available to the technique:
  (LS)  for every product of the factors with multiplicities: rho + 2 W(z) <= L_n(z),
        L_n(z) = max(1, z, 1 - 1/n + 2z/n, 2/n + (1 - 1/n) z)   [de Faveri Thm 1.1, n-th power free]
  (BHM) Bombieri-Halasz-Montgomery on the U^rho rows with Gram row sums
        sum_{u'} |T_{uu'}|^2 = U^{g(z)}, where truth is g = 1 + z (square-root sums, z < 2):
        rho + 2W <= max(z, (rho + g)/2).
The theta profile rho = 1 - 1/n, W(z) = z/n saturates the theta term of (LS) for every z and
every factorisation (it is linear in z).  It violates the diagonal (rho + 2W > 1) iff z > 1/2.

This script checks, for n = 3, 4, 6, the z-range where BHM (with the TRUE Gram size g = 1 + z,
i.e. the most optimistic possible near-orthogonality) excludes the theta profile.
Result: n = 3 -> 1 < z < 2 (OpenAI use z in (4/3, 2) with their weaker Gram bound);
        n = 4 -> empty; n = 6 -> empty.  Hence the technique cannot even recover the sextic
large sieve's sharpness for structured coefficients, let alone DIAG.
"""
import numpy as np

def L(n, z):
    return np.maximum.reduce([np.ones_like(z), z, 1 - 1 / n + 2 * z / n, 2 / n + (1 - 1 / n) * z])

zs = np.linspace(0.0, 3.0, 30001)
for n in (3, 4, 6):
    s = 1.0 / n
    rho = 1 - s
    lhs = rho + 2 * s * zs                     # theta profile
    ls_ok = np.all(lhs <= L(n, zs) + 1e-12)
    bhm_rhs = np.maximum(zs, (rho + 1 + zs) / 2)
    excl = zs[lhs > bhm_rhs + 1e-12]
    rng = (round(float(excl.min()), 4), round(float(excl.max()), 4)) if excl.size else None
    # also with an R-restricted square-root Gram (sum_{u' in R} |T| = R sqrt(Z)):
    bhm2 = np.maximum(zs, rho + zs / 2)
    excl2 = zs[lhs > bhm2 + 1e-12]
    rng2 = (round(float(excl2.min()), 4), round(float(excl2.max()), 4)) if excl2.size else None
    print(f"n={n}: theta profile consistent with all LS constraints: {ls_ok}; "
          f"violates diagonal for z > 1/2: {bool(np.all(lhs[zs > 0.5001] > 1))}; "
          f"BHM(true Gram) excludes it for z in {rng}; BHM(R-restricted sqrt Gram) excludes for z in {rng2}")
print("General: exclusion needs slope s = 1/n > 1/4 (BHM with g = 1 + z needs z(2s - 1/2) > s/2).")
