"""Referee check: e_lib.quad_twist (sparse, prod chi_p(nu)^3) vs an independent Jacobi-symbol computation.
For a primary ideal R = a + b w whose norm n is a squarefree rational integer (all prime factors split, distinct
rational primes), O/R = Z/n via w -> r = -a b^{-1} mod n, and (nu/R)_2 = Jacobi((nu_a + nu_b r) mod n | n)."""
import sys, numpy as np, time
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *
from sympy import jacobi_symbol, factorint
Q = int(sys.argv[1])
t0 = time.time()
R = RowsE(Q, ymin_over_Q=4.0, zeta_check=10)
IX, OmI, Vp, pi2 = pi_layout(R)
rep = IX[:, 0]
rng = np.random.default_rng(7)
# pick rows with squarefree rational norm
cand = [i for i in range(len(rep)) if all(e == 1 for e in factorint(int(R.mN[rep[i]])).values()) and int(R.mN[rep[i]]) % 3 != 0]
sub = rng.choice(cand, size=min(400, len(cand)), replace=False)
nu_a, nu_b, nu_N = nu_columns(Q)
cols = np.concatenate([np.arange(40), rng.choice(len(nu_N), 160, replace=False)])
TW = quad_twist(R, rep[sub], nu_a[cols], nu_b[cols])
mism = 0; tot = 0
for ii, i in enumerate(sub):
    a, b, n = int(R.Ra[rep[i]]), int(R.Rb[rep[i]]), int(R.mN[rep[i]])
    r = (-a * pow(b, -1, n)) % n
    assert (r * r + r + 1) % n == 0
    for jj, j in enumerate(cols):
        y = (int(nu_a[j]) + int(nu_b[j]) * r) % n
        J = jacobi_symbol(y, n) if np.gcd(y, n) == 1 else 0
        tot += 1; mism += (J != int(TW[ii, jj]))
print("Q=%d: quad_twist vs Jacobi: %d mismatches out of %d (rows %d, cols %d) [%.0fs]" % (Q, mism, tot, len(sub), len(cols), time.time() - t0))
# beta columns: (lambda x / R)_2
from cubic_theta import cubic_theta_coeffs
xa, xb, xN, ax, _, _ = cubic_theta_coeffs(int(np.sqrt(Q) / 2), Q)
xa = xa.astype(np.int64); xb = xb.astype(np.int64)
cx = rng.choice(len(xa), 150, replace=False)
la, lb = mul(np.int64(1), np.int64(2), xa[cx], xb[cx])
TWb = quad_twist(R, rep[sub], la, lb)
mism = 0; tot = 0
for ii, i in enumerate(sub):
    a, b, n = int(R.Ra[rep[i]]), int(R.Rb[rep[i]]), int(R.mN[rep[i]])
    r = (-a * pow(b, -1, n)) % n
    for jj in range(len(cx)):
        y = (int(la[jj]) + int(lb[jj]) * r) % n
        J = jacobi_symbol(y, n) if np.gcd(y, n) == 1 else 0
        tot += 1; mism += (J != int(TWb[ii, jj]))
print("Q=%d: beta-column twist (lambda x/R)_2 vs Jacobi: %d mismatches out of %d [%.0fs]" % (Q, mism, tot, time.time() - t0))
# |zeta| = 1 and zeta depends only on the ideal
z = R.zetaZ[IX]
print("max | |zeta|-1 | = %.1e ; max spread of zeta over the 6 associates = %.1e" % (np.abs(np.abs(z) - 1).max(), np.abs(z - z[:, :1]).max()))
