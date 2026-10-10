"""Attack D, check (N1): the Gauss-class calculus behind mu-absorption.

For split primary primes p (norm <= NMAX) compute gamma_j(p), j = 0..5 (gamma_0 := 1), and for every
pair (j, i) the product f = gamma_j(p) gamma_i(p).  We ask:

  (a) is f * alpha(p)^{-k} a function of the class of p mod 72 for some k in -3..3?  ("class 0")
      (alpha(p) = p/|p|; we use alpha = -gamma_2^3 at primes, paper Lemma prime-gauss-identities)
  (b) if so, is that class function a CHARACTER of the ray group (sign +1: f extends to a Hecke x ray
      character on squarefree A) or MINUS a character (sign -1: f extends to mu(A) x Hecke x ray)?
      Test: c(x) c(y) / c(xy) over pairs of classes x, y (all three populated).

Prediction (paper eq. signal + Jacobi sums / Stickelberger):
  * class 0 iff i + j = 0 mod 3;
  * sign -1 (Moebius) iff {i, j} = {1, 2} or {4, 5}  (i + j = 3 mod 6, Jacobi sign);
  * sign +1 iff i + j = 0 mod 6 or (i, j) in {(0,3),(3,0),(3,3)...};
  * single gamma_j (pair (j, 0)) has nontrivial class for j in {1,2,4,5}.
"""
import sys, numpy as np
from collections import defaultdict
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
from zomega import *          # noqa
from a2lib import gamma_prime

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
MOD = int(sys.argv[2]) if len(sys.argv) > 2 else 72
P = [Q for Q in primes_upto(NMAX) if Q.kind == "split"]
G = np.array([[1.0 + 0j] + [gamma_prime(Q, j) for j in range(1, 6)] for Q in P])   # (np, 6)
alpha = np.array([cval(Q.a, Q.b) / abs(cval(Q.a, Q.b)) for Q in P])
# the paper's ray function G(p) = conj(chi_p(4)) gamma_3(p): a QUADRATIC REFINEMENT of the reciprocity
# bicharacter R (G(vw) = G(v)G(w)R(v,w)), not a character.  Products with i+j = 3 mod 6 contain exactly
# one gamma_3, so we divide by G before testing the character property.
from a2lib import sym1
Gp = np.array([np.conj(Z6[sym1(Q, 4, 0)]) * gamma_prime(Q, 3) for Q in P])
assert np.max(np.abs(G[:, 2] ** 3 + alpha)) < 1e-9, "gamma_2^3 = -alpha check failed"
cls = [(Q.a % MOD, Q.b % MOD) for Q in P]

def class_function(vals):
    d = defaultdict(list)
    for c, v in zip(cls, vals):
        d[c].append(v)
    spread = [np.ptp(np.angle(np.array(L) * np.conj(L[0]))) for L in d.values() if len(L) >= 3]
    spread = max(spread) if spread else np.nan
    return spread, {c: np.mean(L) for c, L in d.items()}

def sign_test(cf):
    keys = list(cf.keys()); out = []
    for x in keys[:60]:
        for y in keys[:60]:
            xy = mul(x[0], x[1], y[0], y[1]); xy = (xy[0] % MOD, xy[1] % MOD)
            if xy in cf:
                out.append(cf[x] * cf[y] / cf[xy])
    out = np.array(out)
    return np.mean(out.real), np.max(np.abs(out - np.mean(out)))

print("split primes: %d (norm <= %d), classes mod %d" % (len(P), NMAX, MOD))
print("(for i+j = 3 mod 6 the product is divided by G = conj(chi(4)) gamma_3 before the sign test)")
print(" j  i  i+j mod6 | best k  max within-class spread | sign c(x)c(y)/c(xy)   verdict")
for j in range(0, 6):
    for i in range(j, 6):
        f = G[:, j] * G[:, i]
        if (i + j) % 6 == 3:
            f = f * np.conj(Gp)          # strip the quadratic refinement G (|G| = 1)
        best = None
        for k in range(-3, 4):
            sp, cf = class_function(f * np.conj(alpha) ** k)
            if best is None or sp < best[1]:
                best = (k, sp, cf)
        k, sp, cf = best
        if sp < 1e-6:
            sg, dev = sign_test(cf)
            verdict = ("MOEBIUS x Hecke x ray" if sg < -0.99 else "Hecke x ray (no mu)") + "  (dev %.1e)" % dev
            print("%2d %2d   %d       | %+d     %.2e              | %+.3f   %s" % (j, i, (i + j) % 6, k, sp, sg, verdict))
        else:
            print("%2d %2d   %d       | (none) %.2f rad              |          nontrivial Gauss class" % (j, i, (i + j) % 6, sp))
