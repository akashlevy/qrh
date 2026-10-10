"""
Numerical test of the effective conductor exponent (in the row m) of the dual of the CUBIC Mobius family.

Objects, for primary primes m and smooth dyadic weight W(N(A)/X), sums over primary squarefree A
coprime to 6m (O = Z[omega], primary = 1 mod 3):
  T_mu+ (m,X) = sum mu(A) G~(A) (m/A)_3^{-1} W     [dual of sum_n mu(n) (u/n)_3 ]
  T_mu- (m,X) = sum mu(A) G~(A) (m/A)_3     W     [other orientation]
  K1    (m,X) = sum      G~(A) (m/A)_3^{-1} alpha(A) W   control: Kubota series, numerator m,
                                                  angular twist (no theta pole); conductor N(m)
  P2    (m,X) = sum      G~(A) (m/A)_6     alpha(A)^{-1} W  control: cubic-theta coefficients with
                                                  a SEXTIC twist (paper2's B_h): dual length N(m)^2/X
  RND   (m,X) = sum mu(A) e(random) W                 control: random phases
G~(A) = cubic Gauss sum / sqrt(N A),  alpha(A) = A/|A|.
Statistic: E(X) = mean_m |T(m,X)|^2 / sum_A mu^2(A) W^2.   Random-size ~ 1.
An object with effective conductor N(m)^c and no poles in Re s > 1/2 should drop below ~1 once
X >~ N(m)^c (beyond that only harmless normaliser residues remain, of size <~ X^{5/12}).
"""
import numpy as np, sys, json, time
from math import isqrt

W3 = np.exp(2j * np.pi / 3)
W6 = np.exp(2j * np.pi / 6)
SQ3 = np.sqrt(3.0)


def norm(a, b): return a * a - a * b + b * b
def mul(a1, b1, a2, b2): return a1 * a2 - b1 * b2, a1 * b2 + a2 * b1 - b1 * b2
def cval(a, b): return complex(a - b / 2.0, b * SQ3 / 2.0)


def primary(a, b):
    for (ua, ub) in [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]:
        x, y = mul(a, b, ua, ub)
        if x % 3 == 1 and y % 3 == 0:
            return x, y
    raise ValueError


def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, isqrt(n) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]


def vpow(x, e, p):
    x = x.astype(np.int64) % p; r = np.ones_like(x)
    while e:
        if e & 1: r = (r * x) % p
        x = (x * x) % p; e >>= 1
    return r


class P:
    pass


def make_primes(Nmax):
    """primary primes coprime to 6 with norm <= Nmax: split (norm p) and inert (norm q^2)."""
    out = []
    for p in sieve(Nmax):
        p = int(p)
        if p % 3 == 1:
            # find a,b with norm p
            for b in range(1, isqrt(4 * p // 3) + 2):
                d = 4 * p - 3 * b * b
                if d < 0: break
                r = isqrt(d)
                if r * r == d:
                    a = (b + r) // 2
                    if norm(a, b) == p: break
                    a = (b - r) // 2
                    if norm(a, b) == p: break
            a, b = primary(a, b)
            for (x, y) in {(a, b), primary(*conj(a, b))}:
                q = P(); q.a, q.b, q.N, q.kind, q.p = x, y, p, "split", p
                q.r = (-x * pow(y, -1, p)) % p            # omega mod pi
                out.append(q)
        elif p % 3 == 2 and p != 2 and p * p <= Nmax:
            q = P(); q.a, q.b = primary(p, 0); q.N, q.kind, q.p = p * p, "inert", p
            out.append(q)
    return out


def conj(a, b):  # conj(a + b w) = a + b w^2 = (a - b) - b w
    return a - b, -b


def sym3_split(q, x_mod_p):
    """(x/pi)_3 as power j of omega (x given mod p, nonzero)."""
    e = (q.p - 1) // 3
    v = pow(int(x_mod_p) % q.p, e, q.p)
    if v == 1: return 0
    if v == q.r % q.p: return 1
    return 2


def ff_mul(u, v, q):  # F_q[w]/(w^2+w+1)
    a1, b1 = u; a2, b2 = v
    return ((a1 * a2 - b1 * b2) % q, (a1 * b2 + a2 * b1 - b1 * b2) % q)


def ff_pow(u, e, q):
    r = (1, 0)
    while e:
        if e & 1: r = ff_mul(r, u, q)
        u = ff_mul(u, u, q); e >>= 1
    return r


def symk_inert(q, a, b, k):
    """(x/q)_k for x = a + b w mod q, as power j of exp(2 pi i/k) (k = 3 or 6), q inert prime."""
    v = ff_pow((a % q.p, b % q.p), (q.N - 1) // k, q.p)
    z = (1, 0);
    # zeta_k: for k=3 omega=(0,1); k=6: zeta6 = 1 + omega = (1,1)
    g = (0, 1) if k == 3 else (1, 1)
    for j in range(k):
        if z == v: return j
        z = ff_mul(z, g, q.p)
    raise ValueError


def symk(q, a, b, k):
    if q.kind == "split":
        p = q.p
        x = (a + b * q.r) % p
        if x == 0: return None
        v = pow(x, (p - 1) // k, p)
        # omega ≡ r; zeta6 = 1 + omega
        g = q.r if k == 3 else (1 + q.r) % p
        z = 1
        for j in range(k):
            if z == v: return j
            z = (z * g) % p
        raise ValueError
    else:
        if a % q.p == 0 and b % q.p == 0: return None
        return symk_inert(q, a, b, k)


def gauss(q):
    """normalised cubic Gauss sum G(pi)/sqrt(N pi) with e(Tr(x/pi))."""
    if q.kind == "split":
        p = q.p
        x = np.arange(1, p, dtype=np.int64)
        v = vpow(x, (p - 1) // 3, p)
        j = np.where(v == 1, 0, np.where(v == q.r % p, 1, 2))
        tr = (2 * q.a - q.b) % p
        return complex(np.sum(W3 ** j * np.exp(2j * np.pi * (x * tr % p) / p)) / np.sqrt(p))
    else:
        qq = q.p
        s = 0j
        for a in range(qq):
            for b in range(qq):
                if a == 0 and b == 0: continue
                j = symk_inert(q, a, b, 3)
                # x/(-q): Tr = -(2a - b)/q
                s += W3 ** j * np.exp(-2j * np.pi * ((2 * a - b) % qq) / qq)
        return s / qq


def run(Nmax, mlist_norms, Xs, seed=1):
    t0 = time.time()
    primes = make_primes(Nmax)
    primes.sort(key=lambda q: q.N)
    for q in primes:
        q.G = gauss(q)
    print(f"# {len(primes)} primes, Gauss sums in {time.time()-t0:.1f}s", flush=True)
    # test rows m: primes with given norms
    ms = [q for q in primes if q.N in mlist_norms]
    # per prime: symbols (m/pi)_3, (m/pi)_6
    for q in primes:
        q.s3 = [symk(q, m.a, m.b, 3) for m in ms]
        q.s6 = [symk(q, m.a, m.b, 6) for m in ms]
    # enumerate squarefree products by DFS: state (a,b,N, mu, G, list s3, list s6)
    rng = np.random.default_rng(seed)
    recs = []   # (N, mu, Gt, alpha, s3 vector, s6 vector)
    nm = len(ms)

    def dfs(start, a, b, N, mu, G, s3, s6, fac):
        recs.append((N, mu, G, cval(a, b) / np.sqrt(N), s3.copy(), s6.copy(), fac))
        for i in range(start, len(primes)):
            q = primes[i]
            if N * q.N > Nmax: break
            # twisted multiplicativity: G(A pi) = G(A) G(pi) (A/pi)_3^2
            j = symk(q, a, b, 3)
            if j is None: continue
            Gn = G * q.G * W3 ** (2 * j)
            an, bn = mul(a, b, q.a, q.b)
            s3n = [(x + (y if y is not None else 99)) for x, y in zip(s3, q.s3)]
            s6n = [(x + (y if y is not None else 99)) for x, y in zip(s6, q.s6)]
            dfs(i + 1, an, bn, N * q.N, -mu, Gn, s3n, s6n, fac + [i])
    dfs(0, 1, 0, 1, 1, 1 + 0j, [0] * nm, [0] * nm, [])
    print(f"# {len(recs)} squarefree A, enumeration {time.time()-t0:.1f}s", flush=True)
    Ns = np.array([r[0] for r in recs], float)
    mu = np.array([r[1] for r in recs], float)
    Gt = np.array([r[2] for r in recs])
    al = np.array([r[3] for r in recs])
    S3 = np.array([r[4] for r in recs])   # (nA, nm), >=99 means m | A
    S6 = np.array([r[5] for r in recs])
    ok = S3 < 99
    ch3 = np.where(ok, W3 ** (S3 % 3), 0)
    ch6 = np.where(ok, W6 ** (S6 % 6), 0)
    rnd = np.exp(2j * np.pi * rng.random(len(recs)))
    res = {}

    def Wt(x):  # smooth bump on (1, 2)
        y = np.zeros_like(x)
        m = (x > 1) & (x < 2)
        t = (x[m] - 1.5) / 0.5
        y[m] = np.exp(-1 / (1 - t * t))
        return y
    for X in Xs:
        w = Wt(Ns / X)
        diag = np.sum(w ** 2 * mu ** 2)
        if diag == 0: continue
        out = {}
        for lab, coef in [("mu+", mu[:, None] * Gt[:, None] * np.conj(ch3)),
                          ("mu-", mu[:, None] * Gt[:, None] * ch3),
                          ("K1", Gt[:, None] * np.conj(ch3) * al[:, None]),
                          ("P2", Gt[:, None] * ch6 * np.conj(al)[:, None]),
                          ("RND", (mu * rnd)[:, None] * np.abs(ch3))]:
            T = np.sum(coef * w[:, None], axis=0)          # per m
            out[lab] = [float(abs(t) ** 2 / diag) for t in T]
        res[X] = out
    return ms, res


if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
    norms = [7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 109, 127, 139, 151, 157, 163, 181, 193, 199]
    Xs = [x for x in [20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960] if 2 * x <= Nmax]
    ms, res = run(Nmax, norms, Xs)
    json.dump(dict(m=[(m.a, m.b, m.N) for m in ms], res={str(k): v for k, v in res.items()}),
              open(f"runs/conductor_test_{Nmax}.json", "w"))
    # summary: mean over m grouped by norm bands
    bands = [(7, 13), (19, 43), (61, 109), (127, 199)]
    for lab in ["mu+", "mu-", "K1", "P2", "RND"]:
        print(f"== {lab}")
        for lo, hi in bands:
            idx = [i for i, m in enumerate(ms) if lo <= m.N <= hi]
            row = []
            for X in Xs:
                v = np.mean([res[X][lab][i] for i in idx])
                row.append(f"{v:6.2f}")
            print(f"  N(m) in [{lo},{hi}] (n={len(idx)}):  " + " ".join(row))
    print("  X = " + " ".join(f"{x:6d}" for x in Xs))
