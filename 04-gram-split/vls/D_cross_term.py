"""Attack D, check (N2): the s--A cross term of the Poisson coefficient forces j = a (mod 6).

Probe family: s-slot Gauss sum g_a(s,h) = sum_{d mod s} chi_s(d)^a e(h d/s); A-twist chi_A(m)^j.
For split primes r (s-modulus) and p (A-modulus) of distinct norms, brute-force
   T(H) = sum_{m mod rp} chi_p(m)^j g_a(r,-m) e(H m/(rp))
and compare with the closed form (paper eq. local-F, t = k = 1, H a unit at r and p):
   T(H) = q_r * chi_r(H)^a conj(chi_p(H))^j G_j(p) * X(r,p),   X(r,p) = conj(chi_r(p))^a chi_p(r)^j.
Then test whether X(r,p) is a function of the classes of (r,p) mod 12 (i.e. a ray-class phase that
the paper's R(A,s)-insertion can cancel): predicted iff j = a mod 6, because
   X = chi_p(r)^{j-a} * R(p,r)^{-a}   (sextic reciprocity chi_r(p) = R(p,r) chi_p(r)).
"""
import sys, numpy as np
from collections import defaultdict
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
from zomega import *          # noqa
from a2lib import sym1

def efrac(xa, xb, N):
    """e((xa + xb w)/N) = exp(2 pi i xb / N) (paper's e, Im-based)."""
    return np.exp(2j * np.pi * (np.asarray(xb) % N) / N)

def Tbrute(R, Pp, a, j, H):
    qr, qp = R.N, Pp.N; N = qr * qp
    pa, pb = mul(R.a, R.b, Pp.a, Pp.b)             # pi = r p
    m = np.arange(N, dtype=np.int64)
    sp = Pp.symbol(m, np.zeros_like(m)); chip = np.where(sp >= 0, Z6[(j * sp) % 6], 0)
    if j % 6 == 0: chip = np.where(sp >= 0, 1.0, 0.0)
    d = np.arange(qr, dtype=np.int64); sr = R.symbol(d, np.zeros_like(d))
    chir = np.where(sr >= 0, Z6[(a * sr) % 6], 0)
    # e(-m d / r) = e(-m d conj(r) / q_r)
    ca, cb = R.a - R.b, -R.b                         # conj(r)
    md = np.outer(m, d)
    g = (np.exp(2j * np.pi * ((-md * cb) % qr) / qr) * chir[None, :]).sum(1)
    # e(H m / pi) = e(H m conj(pi) / N)
    hc = mul(H[0], H[1], pa - pb, -pb)
    ph = np.exp(2j * np.pi * ((m * hc[1]) % N) / N)
    return np.sum(chip * g * ph)

def Gj(Pp, j):
    x = np.arange(1, Pp.N, dtype=np.int64); s = Pp.symbol(x, np.zeros_like(x))
    ca, cb = Pp.a - Pp.b, -Pp.b
    return np.sum(Z6[(j * s) % 6] * np.exp(2j * np.pi * ((x * cb) % Pp.N) / Pp.N))

P = [Q for Q in primes_upto(400) if Q.kind == "split"]
# --- 1. formula check on a few pairs
err = 0.0; ntest = 0
for (a, j) in [(1, 1), (1, 5), (2, 1), (2, 2), (5, 5), (4, 1)]:
    for R in P[:5]:
        for Pp in P[:7]:
            if R.N == Pp.N: continue
            for H in [(1, 0), (2, 1), (-1, 3)]:
                if sym1(R, *H) < 0 or sym1(Pp, *H) < 0: continue
                X = np.conj(Z6[sym1(R, Pp.a, Pp.b)]) ** a * Z6[sym1(Pp, R.a, R.b)] ** j
                pred = R.N * Z6[sym1(R, *H)] ** a * np.conj(Z6[sym1(Pp, *H)]) ** j * Gj(Pp, j) * X
                err = max(err, abs(Tbrute(R, Pp, a, j, H) - pred) / abs(pred)); ntest += 1
print("closed form of the Poisson coefficient: %d cases, max rel. error %.2e" % (ntest, err))

# --- 2. is X(r,p) a ray-class phase?
P = [Q for Q in primes_upto(6000) if Q.kind == "split"]
print("cross term X(r,p) over %d split primes: max within-(class r, class p mod 12) phase spread" % len(P))
for (a, j) in [(1, 1), (5, 5), (2, 2), (1, 5), (5, 1), (2, 1), (1, 2), (4, 1)]:
    d = defaultdict(list)
    for R in P[:150]:
        for Pp in P[:150]:
            if R.N == Pp.N: continue
            X = np.conj(Z6[sym1(R, Pp.a, Pp.b)]) ** a * Z6[sym1(Pp, R.a, R.b)] ** j
            d[(R.a % 12, R.b % 12, Pp.a % 12, Pp.b % 12)].append(X)
    sp = max(np.ptp(np.angle(np.array(L) * np.conj(L[0]))) for L in d.values() if len(L) >= 2)
    print("  (a, j) = (%d, %d), j-a = %+d mod 6: spread %.2e  -> %s" %
          (a, j, (j - a) % 6, sp, "ray phase (cancellable)" if sp < 1e-9 else "genuine sextic symbol chi_p(r)^(j-a): NOT cancellable"))
