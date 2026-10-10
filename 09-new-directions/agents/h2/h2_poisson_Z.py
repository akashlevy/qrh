"""
h2, idea 1 (b): Poisson involution of the k-moment with EXACT smooth weights (rational model, Jacobi symbols).
W(x) = sin^2(2 pi (x - 1/2)) on (1/2, 1]  (C^1, hat W(xi) = O(xi^-3)),  hat W(0) = 1/4.

Check A (per pair, exactness):  for m = n1 n2 odd squarefree,
   sum_k W(k/K) (k/m) = (K/m) sum_{h in Z} tau_m(h) hatW(hK/m),   tau_m(h) = sum_{x mod m} (x/m) e(hx/m).
Check B ("dual diagonal"): for n1 = n2 = n, the h != 0 part of the Poisson formula equals
   E(n) = sum_k W(k/K) 1_{(k,n)=1} - K hatW(0) phi(n)/n = sum_{d|n} mu(d) (K/d) sum_{l != 0} hatW(lK/d)
        = -hatW(0) sum_{d|n, d>K} mu(d) K/d + (Schwartz tails from d < K)  =  O(tau(n)),
   while crude bookkeeping (positive dual weight on |h| <= H = N^2/K, diagonal H per column, times K/N) assigns N per column.
Consequence (stated in H2.md): with exact weights the involution maps off-diagonal to off-diagonal (both diagonals of the
mixed moment are the same small character sum); Phi o Phi = id, so no contraction / fixed point is available.
"""
import math, cmath, json, sys
import numpy as np


def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0


def W(x):
    return np.where((x > 0.5) & (x <= 1), np.sin(2 * np.pi * (x - 0.5)) ** 2, 0.0)


def What(xi):
    """int W(x) e(-x xi) dx, closed form for the sin^2 bump on (1/2,1]"""
    xi = np.asarray(xi, float)
    out = np.empty(xi.shape, complex)
    # W = 1/2 - cos(4 pi (x-1/2))/2 = 1/2 - (e(2(x-1/2)) + e(-2(x-1/2)))/4
    def I(c, xi):   # int_{1/2}^{1} e(c (x - 1/2)) e(-x xi) dx
        f = c - xi
        with np.errstate(divide="ignore", invalid="ignore"):
            val = np.where(np.abs(f) < 1e-12, 0.5 + 0j,
                           (np.exp(2j * np.pi * f * 0.5) - 1) / (2j * np.pi * f + 1e-300))
        return val * np.exp(-2j * np.pi * xi * 0.5)
    out = 0.5 * I(0, xi) - 0.25 * (I(2, xi) + I(-2, xi))
    return out


def sqfree_odd(lo, hi):
    out = []
    for n in range(lo + 1, hi + 1):
        if n % 2 == 0: continue
        m, ok, p = n, True, 3
        while p * p <= m:
            if m % (p * p) == 0: ok = False; break
            if m % p == 0: m //= p
            p += 2
        if ok: out.append(n)
    return out


def factor(n):
    f, p = [], 2
    while p * p <= n:
        if n % p == 0:
            f.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: f.append(n)
    return f


def run(N=600, K=150, npairs=6):
    res = {}
    ks = np.arange(1, K + 1)
    Wk = W(ks / K)
    # Check A
    rng = np.random.default_rng(1)
    cols = sqfree_odd(N, 2 * N)
    errs = []
    for _ in range(npairs):
        while True:
            n1, n2 = rng.choice(cols, 2, replace=False)
            if math.gcd(int(n1), int(n2)) == 1: break
        m = int(n1) * int(n2)
        lhs = float(np.sum(Wk * np.array([jacobi(int(k), m) for k in ks])))
        chi = np.array([jacobi(x, m) for x in range(m)], float)
        tau = np.fft.fft(chi)            # tau[h] = sum_x chi(x) e(-hx/m)
        tau = np.conj(tau) if False else tau
        H = int(m / K * 400)
        hs = np.arange(-H, H + 1)
        tauh = tau[(-hs) % m]            # sum_x chi(x) e(+hx/m)
        rhs = (K / m) * np.sum(tauh * What(hs * K / m))
        errs.append(dict(n1=int(n1), n2=int(n2), direct=lhs, poisson=float(rhs.real), imag=float(rhs.imag)))
    res["checkA"] = errs
    print("Check A (per-pair Poisson, exact weights):")
    for e in errs:
        print(f"  n1={e['n1']} n2={e['n2']}: direct {e['direct']:+.6f}  Poisson {e['poisson']:+.6f}  (imag {e['imag']:.1e})")
    # Check B
    rows = []
    for n in cols[:: max(1, len(cols) // 25)]:
        f = factor(n)
        exact = float(np.sum(Wk * np.array([1.0 if math.gcd(int(k), n) == 1 else 0.0 for k in ks])))
        phi_over_n = np.prod([1 - 1 / p for p in f])
        E = exact - K * 0.25 * phi_over_n
        divs = [1]
        for p in f: divs += [d * p for d in divs]
        mu = lambda d: (-1) ** len(factor(d)) if d > 1 else 1
        pred = -0.25 * sum(mu(d) * K / d for d in divs if d > K)
        rows.append(dict(n=n, omega=len(f), E=E, pred_big_d=pred, crude=float(N)))
    res["checkB"] = rows
    Es = np.array([r["E"] for r in rows]); Ps = np.array([r["pred_big_d"] for r in rows])
    print(f"Check B (dual diagonal per column, N={N}, K={K}, crude bookkeeping would give ~{N}):")
    for r in rows[:8]:
        print(f"  n={r['n']:<6d} omega={r['omega']}  E(n) = {r['E']:+.4f}   -hatW(0) sum_(d|n,d>K) mu(d)K/d = {r['pred_big_d']:+.4f}")
    print(f"  max |E(n)| over {len(rows)} columns: {np.abs(Es).max():.4f};  max |E - pred| {np.abs(Es-Ps).max():.4f}")
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    r = run(**args)
    json.dump(r, open("/home/user/qrh/09-new-directions/agents/h2/runs/poisson_Z.json", "w"), indent=1)
