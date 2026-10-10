"""
Numerical sanity check of the "adversarial slot" obstruction (diag agent).

Claim: for the n-th order (n = 3 or 6) residue-symbol family over Z[omega],
    S(u) = P^{-1/2} * sum_{pi ~ P} a(pi) * (u/pi)_n ,   rows u in Z[omega], N(u) <= U,
with a(pi) = conj(normalized Gauss sum g_n(pi)) has mean square
    (1/#u) sum_u |S(u)|^2  ~  1 + c * (P^2/U)^{1/n} ... (theta excess),
i.e. the slot alone violates the diagonal bound once P^2 >> U (P >> U^{1/2}),
whereas Hecke-type coefficients a(pi) = 1 or random phases stay diagonal.

Mechanism: Poisson in u turns sum_u g(pi1,u) conj(g(pi2,u)) into a sum over dual
frequencies h with characters chi_{pi1} conj(chi_{pi2})(h); with Gauss-phase coefficients
the Gauss sums cancel and the dual rows h = w^n (n-th powers) are fully coherent.

Usage: python3 theta_slot_test.py n U P [mode]
"""
import sys, math, cmath
import numpy as np

def primes_upto(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

def eis_norm(a, b):            # N(a + b w), w = e^{2 pi i/3}
    return a * a - a * b + b * b

def cplx(a, b):
    return a + b * complex(-0.5, math.sqrt(3) / 2)

def split_primes(P, n):
    """primary primes pi = a + b w with N(pi) = p, p = 1 mod 3 (and p = 1 mod n), p <= P.
    primary: a = 2 mod 3, b = 0 mod 3 (pi = 2 mod 3)."""
    out = []
    for p in primes_upto(P):
        if p % 3 != 1 or (p - 1) % n != 0: continue
        found = set()
        R = int(math.isqrt(4 * p // 3)) + 2
        for a in range(-R, R + 1):
            D = 4 * p - 3 * a * a
            if D < 0: continue
            s_ = math.isqrt(D)
            if s_ * s_ != D: continue
            for b2 in (a + s_, a - s_):
                if b2 % 2: continue
                b = b2 // 2
                if eis_norm(a, b) == p and a % 3 == 2 and b % 3 == 0:
                    found.add((a, b))
        for (a, b) in found:
            out.append((int(p), a, b))
    return out

def modpow(x, e, p):
    x = x % p
    res = np.ones_like(x)
    while e:
        if e & 1: res = (res * x) % p
        x = (x * x) % p
        e >>= 1
    return res

def residue_data(p, a, b, n):
    """image r of w in F_p with a + b r = 0, and table of n-th roots of unity."""
    # roots of x^2 + x + 1 mod p
    rs = [x for x in range(p) if (x * x + x + 1) % p == 0] if p < 2000 else None
    if rs is None:
        # Tonelli-free: find via random cube root of unity
        g = 2
        while True:
            r = pow(g, (p - 1) // 3, p)
            if r != 1: break
            g += 1
        rs = [r, (r * r) % p]
    r = [x for x in rs if (a + b * x) % p == 0][0]
    z6 = (1 + r) % p                      # image of zeta6 = 1 + w = e^{i pi/3}
    table = {pow(z6, j, p): j for j in range(6)}
    return r, table

def symbols(us_a, us_b, p, a, b, n):
    """(u/pi)_n as complex array for u = us_a + us_b w."""
    r, table = residue_data(p, a, b, n)
    red = (us_a + us_b * r) % p
    e = (p - 1) // n
    vals = modpow(red.astype(np.int64), e, p)
    out = np.zeros(len(vals), complex)
    for v, j in table.items():
        out[vals == v] = cmath.exp(1j * math.pi * j / 3)
    out[red == 0] = 0
    return out

def gauss_sum(p, a, b, n, addchar):
    """g = sum_{x mod pi} (x/pi)_n e(Tr(x / (pi * d))) with d = 1 or sqrt(-3)."""
    xs = np.arange(p)                     # integers 0..p-1 represent Z[w]/pi
    chi = symbols(xs, np.zeros_like(xs), p, a, b, n)
    pic = cplx(a, b)
    d = 1.0 if addchar == 1 else complex(0, math.sqrt(3))
    z = xs / (pic * d)
    return np.sum(chi * np.exp(2j * math.pi * 2 * z.real))

def rows(U):
    R = int(math.isqrt(4 * U)) + 2
    A, B = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
    A = A.ravel(); B = B.ravel()
    Nn = eis_norm(A, B)
    m = (Nn > 0) & (Nn <= U)
    return A[m], B[m]

def powerfree_mask(ua, ub, U, n):
    """True where u is not divisible by pi^n for any prime pi (u in Z[w])."""
    keep = np.ones(len(ua), bool)
    lim = int(U ** (1.0 / n)) + 1
    # prime elements of small norm: lambda = 1 - w (norm 3), 2 (norm 4), split primes
    plist = [(1, -1), (2, 0)]
    for p in primes_upto(max(lim, 7)):
        if p % 3 == 1 and p <= lim:
            for t in split_primes(p, 1):
                if t[0] == p: plist.append((t[1], t[2]))
    for (a, b) in plist:
        if eis_norm(a, b) ** n > U: continue
        # pi^n
        ca, cb = 1, 0
        for _ in range(n):
            ca, cb = ca * a - cb * b, ca * b + cb * a - cb * b
        Nq = eis_norm(ca, cb)
        # u * conj(q): conj(a + b w) = (a - b) - b w
        qa, qb = ca - cb, -cb
        ra = ua * qa - ub * qb
        rb = ua * qb + ub * qa - ub * qb
        keep &= ~((ra % Nq == 0) & (rb % Nq == 0))
    return keep

def run(n, U, P, Plo_frac=0.5, addchar=2, seed=0, filt=1, pfree=False):
    ua, ub = rows(U)
    if pfree:
        m = powerfree_mask(ua, ub, U, n); ua, ub = ua[m], ub[m]
    prs = [t for t in split_primes(P, n) if t[0] > Plo_frac * P and t[0] % filt == 1 % filt]
    rng = np.random.default_rng(seed)
    M = np.zeros((len(ua), len(prs)), complex)
    g = np.zeros(len(prs), complex)
    cls = []
    for j, (p, a, b) in enumerate(prs):
        M[:, j] = symbols(ua, ub, p, a, b, n)
        g[j] = gauss_sum(p, a, b, n, addchar) / math.sqrt(p)
        cls.append(((a % 9), (b % 9)))
    res = {}
    coeffs = {
        "gauss": np.conj(g),
        "one": np.ones(len(prs)),
        "random": np.exp(2j * math.pi * rng.random(len(prs))),
    }
    K = len(prs)
    for name, c in coeffs.items():
        S = M @ c / math.sqrt(K)
        res[name] = float(np.mean(np.abs(S) ** 2))
    return dict(n=n, U=U, P=P, rows=len(ua), primes=K, **res)

if __name__ == "__main__":
    n = int(sys.argv[1]); U = int(sys.argv[2]); P = int(sys.argv[3])
    addchar = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    filt = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    pfree = bool(int(sys.argv[6])) if len(sys.argv) > 6 else False
    print(run(n, U, P, addchar=addchar, filt=filt, pfree=pfree))
