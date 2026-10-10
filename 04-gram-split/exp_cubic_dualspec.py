"""
Cubic analogue: direct look at the Voronoi dual side through additive twists.
For a split prime mu (O/mu = F_p) and a in F_p^*, A(a) = sum_x a(x) W(N x/Q) e(x a/mu).
Voronoi (multiplier chi_mu(a)^4, see exp_cubic_voronoi.py):
    A(a) chi_mu(a)^{-4} = M + sum_{n in Lambda} c_n e(-n abar/mu),     abar = a^{-1} mod mu,
with c_n negligible unless N(n) <~ C N(mu)^2/Q.  So f(b) := A(1/b) chi_mu(1/b)^{-4}, as a function
of b in F_p, has its Fourier transform concentrated on the residues r(n) of SHORT dual vectors n.
We measure: energy fraction of f-hat captured by the residues of n in (1/c)O, N(n) <= X, for several c,
and the radius X needed to capture 90%: Voronoi predicts X_90 ~ N(mu)^2/Q (times a constant).
"""
import numpy as np, sys, os
from zomega import *
from cubic_theta import squarefree_with_gauss, primary_odd_upto, cubic_gauss_echeck

Q = int(float(sys.argv[1])) if len(sys.argv) > 1 else 200000
cache = sys.argv[2] if len(sys.argv) > 2 else "cubic_gauss_cache_400000.npz"
Qc = int(cache.split("_")[-1].split(".")[0])
primesC = primes_upto(Qc); gpC = np.load(cache)["gp"]
keep = [i for i, P in enumerate(primesC) if P.N <= Q]
primes = [primesC[i] for i in keep]; gp = gpC[keep]
A_, B_, N_, G_ = squarefree_with_gauss(Q, primes, gp)
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (u * (1 - u)) + 4); return out
xs = []
for (ba, bb, bn) in primary_odd_upto(int(round(Q ** (1 / 3))) + 1):
    if bn ** 3 > Q: continue
    k = (N_ * bn ** 3 > Q / 2) & (N_ * bn ** 3 <= Q)
    ca, cb = mul(ba, bb, ba, bb); ca, cb = mul(ca, cb, ba, bb)
    xa, xb = mul(A_[k], B_[k], ca, cb); xs.append((xa, xb, N_[k] * bn ** 3, np.sqrt(bn) * G_[k]))
xa = np.concatenate([t[0] for t in xs]); xb = np.concatenate([t[1] for t in xs])
xN = np.concatenate([t[2] for t in xs]); w = bump(xN / Q, .5, 1) * np.concatenate([t[3] for t in xs])
den = np.sum(np.abs(w) ** 2)
print("Q=%d #x=%d" % (Q, len(xN)))

R = 40
ga, gb = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
ga = ga.ravel(); gb = gb.ravel(); gN = norm(ga, gb); o = np.argsort(gN, kind="stable")
ga, gb, gN = ga[o][1:], gb[o][1:], gN[o][1:]          # nonzero lattice points by norm

cs = {"1": (1, 0), "2": (2, 0), "lam": (1, 2), "3": (3, 0), "6": (6, 0), "3lam": (3, 6), "9": (9, 0), "18": (18, 0)}
for P in [p for p in primes if p.kind == "split" and 60 <= p.N <= 2000][::6][:14]:
    p = P.q; pa, pb = P.a - P.b, -P.b
    X, Yc = mul(xa, xb, pa, pb); Yr = Yc % p
    A = np.array([np.sum(w * np.exp(2j * np.pi * ((a * Yr) % p) / p)) for a in range(p)])
    a_s = np.arange(1, p); ea = P.symbol(a_s, np.zeros_like(a_s))
    D = A[1:] * np.conj(Z6[(4 * ea) % 6])                # A(a) chi^{-4}(a)
    f = np.zeros(p, complex); binv = np.array([pow(int(a), -1, p) for a in a_s])
    f[binv] = D                                          # f(b) = D(1/b)
    fh = np.fft.fft(f)                                   # fh[r] = sum_b f(b) exp(-2 pi i r b/p)
    E = np.abs(fh) ** 2; E0 = E[0]; Etot = E.sum()
    # residue of a dual vector n (in O) for the character b -> e(-n b/mu): e(-n b/mu) = exp(-2 pi i b Y(n mubar)/p)
    _, Yn = mul(ga, gb, pa, pb); rn = Yn % p               # fh[r] picks exp(+2 pi i r b/p)... sign absorbed: use both
    line = []
    for cname, (c1, c2) in cs.items():
        # n = n'/c : residue r = Y(n' * cinv * mubar);  cinv = inverse of c mod mu as an integer (c's residue mod mu)
        cres = (c1 + c2 * P.r) % p
        cinv = pow(int(cres), -1, p)
        r = (rn * cinv) % p
        seen = np.zeros(p, bool); seen[0] = True; cum = E0; X90 = None
        for idx in range(len(r)):
            for rr in (r[idx], (-r[idx]) % p):
                if not seen[rr]: seen[rr] = True; cum += E[rr]
            if X90 is None and cum >= 0.9 * Etot: X90 = gN[idx]
        frac_at = lambda Xc: None
        line.append("%s:%s" % (cname, X90))
    print("N(mu)=%5d  N^2/Q=%7.3f  E0/Etot=%.3f  X_90 by lattice scale c -> %s   (random: X_90 ~ %d)"
          % (P.N, P.N ** 2 / Q, E0 / Etot, "  ".join(line), int(0.9 * P.N * 0.55)))
