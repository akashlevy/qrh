"""N2.  Is the paper's completed-row moment sharp for the actual cubic-theta column weights?
   R(beta) := sum_m w(N m/Q) |B0_m|^2 / ( sum_m w(N m/Q) * sum_mu |beta_mu|^2 ),   B0_m = sum_mu beta_mu (mu/m)_2,
with beta_mu = gamma_2(mu) alpha(mu) N(mu)^{-1/2} v(N mu / L) on squarefree primary mu (the Patterson / reflected
column at squarefree indices; the ray factors theta, B_p-normalisations and the lambda^4 shift are dropped:
(lambda^4 mu/m)_2 = (mu/m)_2 since lambda^4 is a square).  m runs over primary m, (N m, 6) = 1.
R ~ 1 means the diagonal Q*||beta||^2 is the whole story: Cauchy-Schwarz in m is sharp, so no manipulation
of the mu-side alone (Type I/II, DR correction, ...) can give a saving (Lemma D of B_bilinear.md).
Controls: beta with random phases; beta = |beta| (coherent); the conjugate orientation conj(gamma_2)."""
import sys, numpy as np, time
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
from a2lib import gamma_prime, sym1, bump

def primary_elements(lo, hi):
    a, b, N = elements_in_annulus(lo, hi)
    keep = (a % 3 == 1) & (b % 3 == 0)
    return a[keep], b[keep], N[keep]

def factor_primary(a, b, N, byq, P):
    """list of prime indices (with multiplicity) of the primary element a + b w."""
    out = []; rem = int(N); ca, cb = int(a), int(b); p = 2
    n = rem; fac = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            fac.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: fac.append(n)
    for p in fac:
        for k in byq.get(p, []):
            Pr = P[k]
            while True:
                if Pr.kind == "split":
                    if (ca + cb * Pr.r) % p: break
                else:
                    if ca % p or cb % p: break
                xa, xb = mul(ca, cb, Pr.a - Pr.b, -Pr.b)
                ca, cb = xa // Pr.N, xb // Pr.N
                out.append(k)
    assert norm(ca, cb) == 1
    return out

def run(Q, ell, seed=0, verbose=True):
    t0 = time.time()
    L = Q ** ell
    P = primes_upto(int(2 * Q) + 1)
    byq = {}
    for k, Pr in enumerate(P): byq.setdefault(Pr.q, []).append(k)
    # columns: squarefree primary mu, N in (L, 2L]
    Pc = [Pr for Pr in P if Pr.N <= 2 * L]
    cols = [c for c in squarefree_products(Pc, int(2 * L)) if c[2] > L]
    g2p = {}
    for c in cols:
        for i in c[3]:
            if i not in g2p: g2p[i] = gamma_prime(Pc[i], 2)
    beta = []; mua = []; mub = []
    for (a, b, N, idx) in cols:
        g = 1.0 + 0j
        for x in range(len(idx)):
            g *= g2p[idx[x]]
            for y in range(x + 1, len(idx)):
                e1 = sym1(Pc[idx[x]], Pc[idx[y]].a, Pc[idx[y]].b); e2 = sym1(Pc[idx[y]], Pc[idx[x]].a, Pc[idx[x]].b)
                g *= Z6[(2 * (e1 + e2)) % 6]
        al = cval(a, b) / abs(cval(a, b))
        beta.append(g * al / np.sqrt(N) * bump(N / L, 1, 2) ** 0.5)
        mua.append(a); mub.append(b)
    beta = np.array(beta); mua = np.array(mua, np.int64); mub = np.array(mub, np.int64)
    # rows
    ra, rb, rN = primary_elements(Q, 2 * Q)
    w = bump(rN / Q, 1, 2)
    # quadratic symbols (mu/P)_2 for all primes P up to 2Q, as +-1 (0 if P | mu)
    nP = len(P)
    S = np.zeros((nP, len(beta)), np.int8)
    for k, Pr in enumerate(P):
        s = Pr.symbol(mua, mub)
        S[k] = np.where(s < 0, 0, np.where(s % 2 == 0, 1, -1))
    # B0_m
    B0 = {}
    norms = {}
    variants = {"true": beta, "conj_gamma2": np.conj(beta) * (np.array([cval(a, b) / abs(cval(a, b)) for a, b in zip(mua, mub)]) ** 2),
                "coherent |beta|": np.abs(beta)}
    rng = np.random.default_rng(seed)
    for r in range(3):
        variants["random phase %d" % r] = np.abs(beta) * np.exp(2j * np.pi * rng.random(len(beta)))
    acc = {k: 0.0 for k in variants}
    for t in range(len(ra)):
        fl = factor_primary(ra[t], rb[t], rN[t], byq, P)
        chi = np.ones(len(beta))
        for k in fl: chi = chi * S[k]
        for kname, bv in variants.items():
            acc[kname] += w[t] * abs(np.dot(bv, chi)) ** 2
    out = {}
    for kname, bv in variants.items():
        out[kname] = acc[kname] / (w.sum() * np.sum(np.abs(bv) ** 2))
    if verbose:
        print("Q=%d L=%.0f (ell=%.2f): rows=%d cols=%d  [%.0fs]" % (Q, L, ell, len(ra), len(beta), time.time() - t0))
        for k, v in out.items(): print("   R[%s] = %.3f" % (k, v))
    return out

if __name__ == "__main__":
    for Q, ell in [(20000, 0.79), (20000, 0.6), (20000, 0.95), (60000, 0.79)]:
        run(Q, ell)

def l1_ratio(Q=20000, ell=0.79):
    """E|B0| / sqrt(E|B0|^2) and E|B0|^4/(E|B0|^2)^2 for the true weights (Gaussian: 0.886 and 2)."""
    L = Q ** ell
    P = primes_upto(int(2 * Q) + 1); byq = {}
    for k, Pr in enumerate(P): byq.setdefault(Pr.q, []).append(k)
    Pc = [Pr for Pr in P if Pr.N <= 2 * L]
    cols = [c for c in squarefree_products(Pc, int(2 * L)) if c[2] > L]
    beta = []; mua = []; mub = []
    for (a, b, N, idx) in cols:
        g = 1.0 + 0j
        for x in range(len(idx)):
            g *= gamma_prime(Pc[idx[x]], 2)
            for y in range(x + 1, len(idx)):
                e1 = sym1(Pc[idx[x]], Pc[idx[y]].a, Pc[idx[y]].b); e2 = sym1(Pc[idx[y]], Pc[idx[x]].a, Pc[idx[x]].b)
                g *= Z6[(2 * (e1 + e2)) % 6]
        beta.append(g * cval(a, b) / abs(cval(a, b)) / np.sqrt(N) * bump(N / L, 1, 2) ** 0.5); mua.append(a); mub.append(b)
    beta = np.array(beta); mua = np.array(mua, np.int64); mub = np.array(mub, np.int64)
    ra, rb, rN = primary_elements(Q, 2 * Q)
    S = np.zeros((len(P), len(beta)), np.int8)
    for k, Pr in enumerate(P):
        s = Pr.symbol(mua, mub); S[k] = np.where(s < 0, 0, np.where(s % 2 == 0, 1, -1))
    vals = []
    for t in range(len(ra)):
        chi = np.ones(len(beta))
        for k in factor_primary(ra[t], rb[t], rN[t], byq, P): chi = chi * S[k]
        vals.append(np.dot(beta, chi))
    v = np.abs(np.array(vals)); m2 = np.mean(v ** 2)
    print("Q=%d ell=%.2f: E|B0|/sqrt(E|B0|^2) = %.3f (Gaussian 0.886);  E|B0|^4/(E|B0|^2)^2 = %.3f (Gaussian 2)"
          % (Q, ell, np.mean(v) / np.sqrt(m2), np.mean(v ** 4) / m2 ** 2))
