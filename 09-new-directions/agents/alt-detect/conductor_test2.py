"""Scaling-collapse version of conductor_test.py (see docstring there).
For every primary prime m with N(m) <= Mmax and every X on a log grid, compute
E = |T(m,X)|^2 / diag and bin by u = log X / log N(m).  A family with effective conductor
N(m)^c (and no poles right of 1/2) leaves the random level ~1 at u ~ c."""
import numpy as np, sys, json, time
import scipy.sparse as sps
import conductor_test as C

def main(Nmax, Mmax, nX=48, seed=3):
    t0 = time.time()
    primes = C.make_primes(Nmax); primes.sort(key=lambda q: q.N)
    for q in primes: q.G = C.gauss(q)
    print(f"# {len(primes)} primes, Gauss sums {time.time()-t0:.1f}s", flush=True)
    ms = [q for q in primes if q.kind == "split" and q.N <= Mmax]
    nm = len(ms)
    # per-prime symbol tables (m/pi)_3, (m/pi)_6  (None -> m == pi)
    S3p = np.zeros((len(primes), nm), np.int64); S6p = np.zeros((len(primes), nm), np.int64)
    for i, q in enumerate(primes):
        for k, m in enumerate(ms):
            j3 = C.symk(q, m.a, m.b, 3); j6 = C.symk(q, m.a, m.b, 6)
            S3p[i, k] = 10 ** 6 if j3 is None else j3
            S6p[i, k] = 10 ** 6 if j6 is None else j6
    print(f"# symbols for {nm} rows {time.time()-t0:.1f}s", flush=True)
    recs = []
    def dfs(start, a, b, N, mu, G, fac):
        recs.append((N, mu, G, C.cval(a, b) / np.sqrt(N), fac))
        for i in range(start, len(primes)):
            q = primes[i]
            if N * q.N > Nmax: break
            j = C.symk(q, a, b, 3)
            if j is None: continue
            an, bn = C.mul(a, b, q.a, q.b)
            dfs(i + 1, an, bn, N * q.N, -mu, G * q.G * C.W3 ** (2 * j), fac + [i])
    sys.setrecursionlimit(10000)
    dfs(0, 1, 0, 1, 1, 1 + 0j, [])
    nA = len(recs)
    print(f"# {nA} squarefree A {time.time()-t0:.1f}s", flush=True)
    rows, cols = [], []
    for r_, rec in enumerate(recs):
        for i in rec[4]:
            rows.append(r_); cols.append(i)
    F = sps.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(nA, len(primes)))
    S3 = F @ S3p; S6 = F @ S6p                       # (nA, nm) sums; >=1e6 means m | A
    ok = S3 < 10 ** 6
    ch3 = np.where(ok, C.W3 ** (S3 % 3), 0); ch6 = np.where(ok, C.W6 ** (S6 % 6), 0)
    Ns = np.array([r[0] for r in recs], float); mu = np.array([r[1] for r in recs], float)
    Gt = np.array([r[2] for r in recs]); al = np.array([r[3] for r in recs])
    rnd = np.exp(2j * np.pi * np.random.default_rng(seed).random(nA))
    objs = {"mu+": mu[:, None] * Gt[:, None] * np.conj(ch3),
            "mu-": mu[:, None] * Gt[:, None] * ch3,
            "K1": Gt[:, None] * np.conj(ch3) * al[:, None],
            "P2": Gt[:, None] * ch6 * np.conj(al)[:, None],
            "RND": (mu * rnd)[:, None] * np.abs(ch3)}
    Xs = np.exp(np.linspace(np.log(20), np.log(Nmax / 2.05), nX))
    def Wt(x):
        y = np.zeros_like(x); m = (x > 1) & (x < 2); t = (x[m] - 1.5) / 0.5
        y[m] = np.exp(-1 / (1 - t * t)); return y
    logNm = np.log([m.N for m in ms])
    out = {k: [] for k in objs}
    for X in Xs:
        w = Wt(Ns / X); diag = np.sum(w ** 2 * mu ** 2)
        u = np.log(X) / logNm
        for k, coef in objs.items():
            T = (w @ coef)
            out[k].append(np.c_[u, np.abs(T) ** 2 / diag, np.log(np.array([m.N for m in ms]))])
    print(f"# sums done {time.time()-t0:.1f}s", flush=True)
    bins = np.arange(0.0, 7.01, 0.5)
    summ = {}
    for k in objs:
        D = np.vstack(out[k])
        row = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            s = (D[:, 0] >= lo) & (D[:, 0] < hi)
            row.append((float(D[s, 1].mean()) if s.sum() > 5 else float("nan"), int(s.sum())))
        summ[k] = row
    print("u-bin:  " + " ".join(f"{lo:5.1f}" for lo in bins[:-1]))
    for k, row in summ.items():
        print(f"{k:5s}:  " + " ".join(f"{v:5.2f}" for v, n in row))
    print("count:  " + " ".join(f"{n:5d}" for v, n in summ['mu+']))
    np.savez(f"runs/conductor_test2_{Nmax}.npz", **{k: np.vstack(v) for k, v in out.items()})
    return summ

if __name__ == "__main__":
    Nmax = int(sys.argv[1]); Mmax = int(sys.argv[2])
    main(Nmax, Mmax)
