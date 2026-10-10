"""Mean square over quadratic twists of the TRUE row sequence a(m) = rho(m) conj(gamma_2(m))^4 (squarefree m),
versus the zeta-free sequence rho(m), at Q = 3000, with the cross-scale estimator of exp_sextic_jtable.py.
R(mu) = |sum_m Om(m) a(m) psi_mu(m)|^2 / sum_m Om^2 |a(m)|^2   (each estimated by products over distinct scales).
Random-like sequence: E R = 1 for every mu.  Large-sieve saturation would need R ~ Q/#mu >> 1 on many mu.
psi_mu(m) = (m/mu)_6^j with mu squarefree primary (j=3: quadratic; j=1: coherent control for rho).
Ray factors (units, reciprocity signs, kappa0 psi1) are omitted: they do not change these statistics."""
import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
Q = 3000; i0 = int(sys.argv[1]) if len(sys.argv) > 1 else 4
A = np.load("/Users/akashlevy/Documents/riemann/gram_split/A_Q%d.npy" % Q)
ma, mb, mN = np.load("/Users/akashlevy/Documents/riemann/gram_split/m_Q%d.npy" % Q)
K = A.shape[0]; Ys = np.array([256 * 2 ** i for i in range(K)], float)
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (u * (1 - u)) + 4); return out
Om = bump(mN / Q, .5, 1)
V = A[i0:] * Ys[i0:, None] ** (1 / 3)
pairs = [(i, j) for i in range(len(V)) for j in range(i + 1, len(V))]
PR = primes_upto(Q)
def gam(P, j):
    if P.kind == "split":
        p = P.q; x = np.arange(1, p, dtype=np.int64)
        s = P.symbol(x, np.zeros_like(x)); ph = np.exp(-2j * np.pi * ((P.b * x) % p) / p)
        return np.sum(Z6[(j * s) % 6] * ph) / np.sqrt(p)
    q = P.q; c, d = np.meshgrid(np.arange(q), np.arange(q), indexing="ij"); c = c.ravel(); d = d.ravel()
    nz = (c != 0) | (d != 0); c = c[nz]; d = d[nz]; s = P.symbol(c, d)
    sgn = -1 if P.a < 0 else 1
    return np.sum(Z6[(j * s) % 6] * np.exp(-2j * np.pi * sgn * d / q)) / q
G2 = {}
# factor each m: list of prime indices dividing m, multiplicities via norms
M = len(ma); g2 = np.zeros(M, complex); sqf = np.zeros(M, bool)
divs = [[] for _ in range(M)]
for k, P in enumerate(PR):
    e = P.symbol(ma, mb)
    for idx in np.nonzero(e < 0)[0]:
        divs[idx].append(k)
units = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]
for idx in range(M):
    ks = divs[idx]
    if int(np.prod([PR[k].N for k in ks])) != int(mN[idx]):
        continue                                  # not squarefree (or not coprime to 6)
    sqf[idx] = True
    val = 1.0 + 0j
    for k in ks:
        if k not in G2: G2[k] = gam(PR[k], 2)
        val *= G2[k]
    for x in range(len(ks)):                       # twisted multiplicativity cross terms (j = 2)
        for y in range(x + 1, len(ks)):
            P1, P2 = PR[ks[x]], PR[ks[y]]
            e12 = int(P1.symbol(np.array([P2.a]), np.array([P2.b]))[0]); e21 = int(P2.symbol(np.array([P1.a]), np.array([P1.b]))[0])
            val *= Z6[(2 * (e12 + e21)) % 6]
    # unit: m = u * m0 with m0 primary; gamma_2(m) = chi_{m0}(u)^2 gamma_2(m0)
    a0, b0 = primary(int(ma[idx]), int(mb[idx]))
    for (ua, ub) in units:
        if mul(a0, b0, ua, ub) == (int(ma[idx]), int(mb[idx])):
            eu = sum(int(PR[k].symbol(np.array([ua]), np.array([ub]))[0]) for k in ks)
            val *= Z6[(2 * eu) % 6]; break
    g2[idx] = val
print("m in annulus:", M, " squarefree:", sqf.sum(), " max||gamma2|-1| =", np.max(np.abs(np.abs(g2[sqf]) - 1)))
zeta = np.where(sqf, np.conj(g2) ** 4, 0)
# mu: squarefree primary products of primes, 1 < N(mu) <= 600
mus = [t for t in squarefree_products(PR, 600) if t[2] > 1]
def stat(weight, j):
    W = Om * weight
    bn = np.mean([np.sum(Om ** 2 * np.abs(weight) ** 2 * V[i] * np.conj(V[k])).real for i, k in pairs])
    res = []
    for (a, b, N, idxs) in mus:
        e = np.zeros(M, np.int64); bad = np.zeros(M, bool)
        for k in idxs:
            s = PR[k].symbol(ma, mb); bad |= s < 0; e += np.where(s < 0, 0, s)
        psi = np.where(bad, 0, Z6[(j * e) % 6])
        T = [np.sum(W * V[i] * psi) for i in range(len(V))]
        res.append((N, np.mean([(T[i] * np.conj(T[k])).real for i, k in pairs]) / bn))
    return np.array(res)
edges = [1, 60, 200, 601]
for label, wt, j in [("TRUE a = rho*conj(g2)^4, quadratic twist j=3", zeta, 3),
                     ("zeta-free rho (sqfree m), quadratic twist j=3", sqf.astype(float), 3),
                     ("TRUE a, sextic twist j=1", zeta, 1),
                     ("control: zeta-free rho, coherent sextic twist j=1", sqf.astype(float), 1)]:
    r = stat(wt, j)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (r[:, 0] >= lo) & (r[:, 0] < hi); x = r[sel, 1]
        out.append("N(mu) in [%d,%d): %6.2f (%.2f) [%d]" % (lo, hi, x.mean(), x.std() / np.sqrt(len(x)), len(x)))
    print(label); print("    " + "   ".join(out))

if __name__ == "__main__" and len(sys.argv) > 2:
    r = stat(zeta, 3); x = r[:, 1]
    print("TRUE a, j=3: per-mu R: max %.2f, 95%% quantile %.2f, fraction R>5: %.3f, sum_mu R / #mu = %.2f, Q/#mu = %.1f"
          % (x.max(), np.quantile(x, .95), np.mean(x > 5), x.mean(), Q / len(x)))
    # bilinear form with model column weights beta_mu = conj(gamma_2(mu)) N(mu)^{-1/2} (cubic-theta-like), all mu
    print("(model bilinear form uses the same cross-scale estimator; ratio to Cauchy-Schwarz scale ||beta|| (Q ||a||^2)^{1/2})")

if __name__ == "__main__" and len(sys.argv) > 3:
    for lab, wt in [("zeta-free rho, j=3 (baseline)", sqf.astype(float))]:
        x = stat(wt, 3)[:, 1]
        print(lab + ": per-mu R: max %.2f, 95%% quantile %.2f, fraction R>5: %.3f, mean %.2f" % (x.max(), np.quantile(x, .95), np.mean(x > 5), x.mean()))
