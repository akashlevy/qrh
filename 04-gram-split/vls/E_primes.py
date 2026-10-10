"""Attack E, stage 1: all primary primes of O = Z[omega] (outside 2, lambda) with norm <= X = 2^22,
with r (omega = r mod pi, split), sixth roots table, and normalised sextic Gauss sums gtilde(pi) = g(pi,1)/sqrt(N pi).

Conventions exactly those of zomega.py / D_numerics/fastgauss.py / expA.py (validated in E_validate.py):
  split pi = a + b w (primary: a = 1, b = 0 mod 3), p = N pi, r = -a b^{-1} mod p, zeta6 = 1 + w = 1 + r mod pi,
  g(pi,1) = sum_{x mod p} (x/pi)_6 exp(-2 pi i b x / p);  inert pi = -q, q = 2 mod 3, q^2 <= X (zomega.Prime).
Only one of each conjugate pair is summed: g(pibar) = conj(g(pi)) (pibar = (a-b) - b w), checked on the first 2000.
Output: vls/data/primes_X{X}.npz
"""
import numpy as np, sys, os, time
from math import isqrt
GS = "/Users/akashlevy/Documents/riemann/gram_split"
sys.path.insert(0, GS); sys.path.insert(0, GS + "/D_numerics")
from zomega import sieve, Prime, Z6, norm
from fastgauss import gauss_split
from multiprocessing import Pool

DATA = GS + "/vls/data/"


def powmod_vec(base, e, mod):
    """elementwise base^e mod mod (int64 arrays, mod < 2^31)."""
    base = np.asarray(base, np.int64) % mod; e = np.asarray(e, np.int64)
    res = np.ones_like(base)
    nb = int(e.max()).bit_length() if e.size else 0
    for bit in range(nb):
        m = ((e >> bit) & 1).astype(bool)
        res = np.where(m, (res * base) % mod, res)
        base = (base * base) % mod
    return res


def gen_primes(X):
    t0 = time.time()
    isp = np.zeros(X + 1, bool); isp[sieve(X)] = True
    R = isqrt(4 * X // 3) + 3
    av = np.arange(-R, R + 1); av = av[av % 3 == 1]
    bv = np.arange(-R, R + 1); bv = bv[bv % 3 == 0]
    A, B = np.meshgrid(av, bv, indexing="ij"); A = A.ravel().astype(np.int64); B = B.ravel().astype(np.int64)
    N = A * A - A * B + B * B
    keep = (N <= X) & (N > 3)
    A, B, N = A[keep], B[keep], N[keep]
    sel = isp[N] & (N % 3 == 1)
    sa, sb, sN = A[sel], B[sel], N[sel]
    o = np.lexsort((sa, sb, sN)); sa, sb, sN = sa[o], sb[o], sN[o]
    assert np.all(np.bincount(sN)[np.unique(sN)] == 2)
    # r = -a b^{-1} mod p
    binv = powmod_vec(sb % sN, sN - 2, sN)
    r = (-sa % sN) * binv % sN
    assert np.all((r * r + r + 1) % sN == 0)
    z = (1 + r) % sN
    roots = np.ones((len(sN), 6), np.int64)
    for j in range(1, 6):
        roots[:, j] = roots[:, j - 1] * z % sN
    # inert
    qs = [int(q) for q in sieve(isqrt(X)) if q % 3 == 2 and q != 2 and q * q <= X]
    ia = np.array([-q for q in qs], np.int64); ib = np.zeros(len(qs), np.int64); iN = np.array([q * q for q in qs], np.int64)
    a = np.concatenate([sa, ia]); b = np.concatenate([sb, ib]); Nn = np.concatenate([sN, iN])
    q = np.concatenate([sN, np.array(qs, np.int64)])
    kind = np.concatenate([np.zeros(len(sN), np.int8), np.ones(len(qs), np.int8)])     # 0 split, 1 inert
    rr = np.concatenate([r, np.zeros(len(qs), np.int64)])
    roots = np.concatenate([roots, np.zeros((len(qs), 6), np.int64)])
    o = np.argsort(Nn, kind="stable")
    out = dict(a=a[o], b=b[o], N=Nn[o], q=q[o], kind=kind[o], r=rr[o], roots=roots[o])
    print("primes: %d (split %d, inert %d)  [%.1fs]" % (len(o), len(sN), len(qs), time.time() - t0), flush=True)
    return out


def conj_index(P):
    """index of the conjugate prime (split); itself for inert."""
    key = {}
    for i, (a, b) in enumerate(zip(P["a"].tolist(), P["b"].tolist())):
        key[(a, b)] = i
    ci = np.array([key[(a - b, -b)] if k == 0 else i for i, (a, b, k) in
                   enumerate(zip(P["a"].tolist(), P["b"].tolist(), P["kind"].tolist()))])
    return ci


if __name__ == "__main__":
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 2 ** 22
    t0 = time.time()
    P = gen_primes(X)
    nP = len(P["N"])
    ci = conj_index(P)
    gt = np.zeros(nP, complex)
    # representatives: split with b > 0 (exactly one of each pair since conj(b) = -b)
    rep = np.where((P["kind"] == 0) & (P["b"] > 0))[0]
    assert len(rep) * 2 == (P["kind"] == 0).sum()
    # check conj relation on the first 2000 split primes (both computed)
    chk = np.where(P["kind"] == 0)[0][:2000]
    gchk = np.array([gauss_split((int(P["a"][i]), int(P["b"][i]), int(P["q"][i]), int(P["r"][i]))) for i in chk])
    gt[chk] = gchk
    print("max |g(pi) - conj g(pibar)| on first 2000 split:", np.max(np.abs(gt[chk] - np.conj(gt[ci[chk]]))) if np.all(np.isin(ci[chk], chk)) else "n/a", flush=True)
    todo = rep[~np.isin(rep, chk)]
    # big primes first for load balance
    todo = todo[np.argsort(-P["N"][todo])]
    args = [(int(P["a"][i]), int(P["b"][i]), int(P["q"][i]), int(P["r"][i])) for i in todo]
    with Pool(11) as pool:
        gs = pool.map(gauss_split, args, chunksize=16)
    gt[todo] = np.array(gs)
    sp = np.where(P["kind"] == 0)[0]
    missing = sp[gt[sp] == 0]
    gt[missing] = np.conj(gt[ci[missing]])
    for i in np.where(P["kind"] == 1)[0]:
        Pr = Prime(int(P["a"][i]), int(P["b"][i]), int(P["N"][i]), "inert", int(P["q"][i]))
        gt[i] = Pr.gauss() / np.sqrt(Pr.N)
    print("gauss sums done [%.1fs]  max||g|-1| = %.2e" % (time.time() - t0, np.max(np.abs(np.abs(gt) - 1))), flush=True)
    np.savez(DATA + "primes_X%d.npz" % X, gt=gt, ci=ci, **P)
    print("saved", DATA + "primes_X%d.npz" % X)
