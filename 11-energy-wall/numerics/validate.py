"""
numerics / step 1: independent re-validation of the edge-theta library.

(a) beta(n) from edge_lib (twisted multiplicativity from gamma_2 at primes) against a DIRECT
    evaluation gamma_2(n) conj(alpha(n)), gamma_j(n;h) = N(n)^{-1/2} sum_{x mod n} chi_n(x)^j e(hx/n),
    summed over x = a + b w, 0 <= a,b < N(n) (each class N(n) times), for squarefree primary n
    with omega(n) = 1..3, including inert prime factors (norm q^2).
(b) the identity (*) of WALL_SPECTRAL s.1.1:
       beta(n) psi_k(n) chi_P(n)^{-2} = mu(n) G(n) r(n,k) conj gamma_1(n; k^3 P^4),
    G(n) = conj chi_n(4) Gamma(n), r(n,k) = Gamma(nk)/(Gamma(n)Gamma(k)), for random primary k, P coprime to n
    (k squarefree with 1..2 prime factors, P a split prime), all by direct sums.
(c) the Moebius partner a(n) = mu(n) conj G(n) = conj(gamma_1(n) beta(n)) (P78 eq. signal).
"""
import sys, math, itertools, time
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/energy-wall")
from edge_lib import setup, Symb, Columns, Gamma, W
from eisenstein import mul, norm, cx

fields, pnorm, plist = setup(3000)
symb = {}
def S(p):
    if p not in symb: symb[p] = Symb(fields[p])
    return symb[p]

def chi_exp(den_primes, A, B):
    """sextic exponent of (x/n)_6, n = prod den_primes; -1 where not a unit."""
    e = np.zeros(np.shape(A), np.int64); z = np.zeros(np.shape(A), bool)
    for p in den_primes:
        k = S(p)(A, B).astype(np.int64)
        z |= k < 0; e += np.where(k < 0, 0, k)
    return np.where(z, -1, e % 6)

def red_mod(h, n):
    """exact reduction of h modulo n in Z[w] (Python ints): h - round(h nbar / N(n)) n."""
    Nn = norm(n); nb = (n[0] - n[1], -n[1])
    t = mul(h, nb)
    qa = (2 * t[0] + Nn) // (2 * Nn); qb = (2 * t[1] + Nn) // (2 * Nn)
    qn = mul((qa, qb), n)
    return (h[0] - qn[0], h[1] - qn[1])

def gam_direct(fac, j, h=(1, 0)):
    n = (1, 0)
    for p in fac: n = mul(n, p)
    h = red_mod(h, n)
    Nn = norm(n); nc = cx(n)
    a = np.arange(Nn)
    A, B = np.meshgrid(a, a, indexing="ij"); A = A.ravel(); B = B.ravel()
    e = chi_exp(fac, A, B)
    hx_a = A * h[0] - B * h[1]; hx_b = A * h[1] + B * h[0] - B * h[1]        # h * x
    z = (hx_a + hx_b * W) / nc
    val = np.where(e < 0, 0, np.exp(2j * np.pi * j * np.where(e < 0, 0, e) / 6)) * np.exp(4j * np.pi * z.imag / math.sqrt(3))
    return val.sum() / Nn / math.sqrt(Nn), n

def epow(x, k):
    r = (1, 0)
    for _ in range(k): r = mul(r, x)
    return r

t0 = time.time()
good = [p for p in plist if pnorm[p] not in (3, 4)]
split = [p for p in good if fields[p].kind == "split"]
inert = [p for p in good if fields[p].kind == "inert"]
rng = np.random.default_rng(1)
cases = []
# omega = 1: 6 split + 2 inert ; omega = 2: some with an inert factor ; omega = 3 small
for p in split[:6] + inert[:2]: cases.append([p])
for _ in range(10):
    f = list(rng.choice(len(split[:25]), 2, replace=False)); cases.append([split[i] for i in f])
cases.append([inert[0], split[0]]); cases.append([inert[0], split[3]]); cases.append([inert[1], split[1]])
for _ in range(4):
    f = list(rng.choice(len(split[:8]), 3, replace=False)); cases.append([split[i] for i in f])
cases = [c for c in cases if math.prod(pnorm[p] for p in c) <= 1200]
# library beta via Columns over a norm window containing each n
res_a = []; res_b = []; res_c = []
for fac in cases:
    g2, n = gam_direct(fac, 2); g1, _ = gam_direct(fac, 1)
    Nn = norm(n); al = cx(n) / abs(cx(n)); mu = (-1) ** len(fac)
    beta_dir = g2 * np.conj(al)
    C = Columns(Nn - 1, Nn, fields, pnorm, plist, None)
    idx = [i for i in range(C.n) if (int(C.A[i]), int(C.B[i])) == n]
    if not idx:      # n may be excluded by the (n,6)=1 filter (inert 2) -- skip
        continue
    i = idx[0]
    res_a.append(abs(C.beta[i] - beta_dir))
    # (c) Moebius partner
    res_c.append(abs(C.amob[i] - np.conj(g1 * beta_dir)))
    # (b) identity (*), 3 random (k, P) per n
    chi4 = np.exp(2j * np.pi * chi_exp(fac, np.array([4]), np.array([0]))[0] / 6)
    Gn = np.conj(chi4) * Gamma(n)
    others = [p for p in split[6:60] if p not in fac]
    for _ in range(3):
        kf = [others[i] for i in rng.choice(len(others), int(rng.integers(1, 3)), replace=False)]
        P = others[int(rng.integers(len(others)))]
        if P in kf: continue
        k = (1, 0)
        for p in kf: k = mul(k, p)
        nA, nB = np.array([n[0]]), np.array([n[1]])
        psi = 1.0
        for p in kf:                       # psi_k(n) = prod_{p|k} (n/p)_6^3
            e = S(p)(nA, nB)[0]; psi *= (-1) ** int(e)
        eP = S(P)(nA, nB)[0]
        lhs = beta_dir * psi * np.exp(-2j * np.pi * 2 * int(eP) / 6)
        h = mul(epow(k, 3), epow(P, 4))
        g1h, _ = gam_direct(fac, 1, h)
        rr = Gamma(mul(n, k)) / (Gamma(n) * Gamma(k))
        rhs = mu * Gn * rr * np.conj(g1h)
        res_b.append((len(fac), any(fields[p].kind == "inert" for p in fac), abs(lhs - rhs), abs(lhs + rhs)))
print(f"(a) library beta vs direct gamma_2 conj(alpha): {len(res_a)} squarefree n (omega 1..3, incl. inert factors), "
      f"max |diff| = {max(res_a):.2e}")
print(f"(c) Moebius partner a = mu conj G vs conj(gamma_1 beta) (direct): max |diff| = {max(res_c):.2e}")
fails = sum(r[2] > 1e-6 for r in res_b)
print(f"(b) identity (*): {len(res_b)} tests (with inert factor: {sum(r[1] for r in res_b)}), failures {fails}, "
      f"max |lhs-rhs| = {max(r[2] for r in res_b):.2e}; without mu(n) the omega-odd cases would fail: "
      f"{sum((r[0] % 2 == 1) and r[3] > 1e-6 for r in res_b)} of {sum(r[0] % 2 == 1 for r in res_b)} odd cases differ in sign")
print(f"time {time.time()-t0:.1f}s")
