"""Sanity checks of the Gauss-sum conventions used in conductor_test.py:
 (1) |G~(pi)| = 1;  (2) J := G^3/N(pi) is an Eisenstein integer of norm N(pi) (Gauss/Jacobi);
 (3) twisted multiplicativity G(AB) = G(A)G(B)(A/B)_3^2 against a direct sum over O/AB;
 (4) mu-absorption for cubic twists:  mu(A) G~(A) = alpha-type Hecke * conj(G~(A))^2, i.e.
     G~^3 * conj(alpha_J) is a ray-class quantity: we print G~(pi)^3 / (J/|J|)."""
import numpy as np
import conductor_test as C
pr = C.make_primes(400)
for q in pr: q.G = C.gauss(q)
print("max ||G~|-1| =", max(abs(abs(q.G) - 1) for q in pr))
for q in pr[:8]:
    J = q.G ** 3 * np.sqrt(q.N) ** 3 / q.N
    a = J.real + J.imag / np.sqrt(3); b = 2 * J.imag / np.sqrt(3)   # J = a + b w
    print(f"N={q.N:4d} pi=({q.a},{q.b}) J~({a:+.6f},{b:+.6f}) N(J)={C.norm(round(a),round(b))} "
          f" J mod 3 = ({round(a)%3},{round(b)%3})  J/pi or J/conj(pi): "
          f"{complex(round(a)-round(b)/2, round(b)*np.sqrt(3)/2)/C.cval(q.a,q.b):.3f} "
          f"{complex(round(a)-round(b)/2, round(b)*np.sqrt(3)/2)/C.cval(*C.conj(q.a,q.b)):.3f}")
# (3) direct check for A = pi1 pi2 (split, distinct norms)
sp = [q for q in pr if q.kind == "split"]
q1, q2 = sp[0], sp[2]
a, b = C.mul(q1.a, q1.b, q2.a, q2.b); NA = q1.N * q2.N
# O/A = Z/NA (omega -> r mod NA via CRT)
from sympy.ntheory.modular import crt
r = int(crt([q1.p, q2.p], [q1.r, q2.r])[0])
x = np.arange(1, NA)
def sym(q, x):
    v = C.vpow(x % q.p, (q.p - 1) // 3, q.p)
    j = np.where(v == 1, 0, np.where(v == q.r % q.p, 1, 2)); j = np.where(x % q.p == 0, -1, j)
    return j
j1, j2 = sym(q1, x), sym(q2, x)
okm = (j1 >= 0) & (j2 >= 0)
# Tr(x/A) = x(2a' - b')/NA where A = a + b w ; (x/A) = x conj(A)/NA, Re(conj A) = a - b/2
tr = (2 * a - b) % NA
Gdir = np.sum((C.W3 ** (j1 + j2))[okm] * np.exp(2j * np.pi * (x[okm] * tr % NA) / NA)) / np.sqrt(NA)
j = C.symk(q2, q1.a, q1.b, 3)   # (pi1/pi2)_3
Gtm = q1.G * q2.G * C.W3 ** (2 * j)
print("twisted multiplicativity: direct", np.round(Gdir, 6), " formula", np.round(Gtm, 6))
