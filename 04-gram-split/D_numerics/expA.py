"""Extended version of gram_split/exp_rank1.py (same conventions): Ymax up to 2^20, several Q at once.
A_i(m) = Y_i^{-1} sum_{s in [Y_i,2Y_i]} W(Ns/Y_i)(Ns/Y_i)^{-1/2} gtilde(s,-m), s squarefree primary, (s,6)=1, c(s)=1, v=0.
usage: python expA.py K Q1 Q2 ...   (Y_i = 256*2^i, i<K)"""
import numpy as np, sys, time, os
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
from fastgauss import gauss_split
from multiprocessing import Pool
OUT = os.path.dirname(os.path.abspath(__file__)) + "/data/"
os.makedirs(OUT, exist_ok=True)
K = int(sys.argv[1]); Qs = [float(x) for x in sys.argv[2:]]
Y0 = 256; Ymax = Y0 * 2 ** K
t0 = time.time()
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi)
    u = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1.0 / (u * (1 - u)) + 4.0); return out

def symblock(args):
    """sextic symbol exponents of x=(xa,xb) for a block of split primes (p, r, roots)."""
    ps, rs, roots, xa, xb = args
    p = ps[:, None]; v = (xa[None, :] + xb[None, :] * rs[:, None]) % p
    e = ((ps - 1) // 6)[:, None]
    res = np.ones_like(v); b = v.copy()
    for bit in range(22):
        m = ((e >> bit) & 1).astype(bool)
        res = np.where(m, (res * b) % p, res); b = (b * b) % p
    out = np.full(v.shape, -1, np.int8)
    for j in range(6):
        out[res == roots[:, j:j + 1]] = j
    out[v == 0] = -1
    return out

if __name__ == "__main__":
    P = primes_upto(Ymax)
    nP = len(P); print("primes", nP, "t=%.1f" % (time.time() - t0), flush=True)
    gfile = OUT + "gt_%d.npy" % Ymax
    if os.path.exists(gfile):
        gt = np.load(gfile)
    else:
        split = [i for i, Pr in enumerate(P) if Pr.kind == "split"]
        with Pool(11) as pool:
            gs = pool.map(gauss_split, [(P[i].a, P[i].b, P[i].q, P[i].r) for i in split], chunksize=64)
        gt = np.zeros(nP, complex); gt[split] = gs
        for i, Pr in enumerate(P):
            if Pr.kind == "inert": gt[i] = Pr.gauss() / np.sqrt(Pr.N)
        np.save(gfile, gt)
    print("gauss sums done t=%.1f  max||g|-1|=%.2e" % (time.time() - t0, np.max(np.abs(np.abs(gt) - 1))), flush=True)
    S = squarefree_products(P, Ymax)
    print("squarefree s:", len(S), "t=%.1f" % (time.time() - t0), flush=True)
    Ys = [Y0 * 2 ** i for i in range(K)]
    # gtilde(s,1) and scale weights, once
    cache = {}
    def cross(ii, jj):
        e1 = P[ii].symbol(np.array([P[jj].a]), np.array([P[jj].b]))[0]
        e2 = P[jj].symbol(np.array([P[ii].a]), np.array([P[ii].b]))[0]
        return Z6[(e1 + e2) % 6]
    slist = []
    for (sa, sb, sN, idx) in S:
        if sN < Y0: continue
        i = int(np.floor(np.log2(sN / Y0)))
        if i >= K: continue
        rr = sN / Ys[i]
        w = bump(np.array([rr]), 1, 2)[0] * rr ** -0.5 / Ys[i]
        if w == 0: continue
        g = 1.0 + 0j
        for t, ii in enumerate(idx):
            g *= gt[ii]
            for jj in idx[t + 1:]:
                g *= cross(ii, jj)
        slist.append((i, w * g, np.array(idx, np.int64)))
    print("s weights done:", len(slist), "t=%.1f" % (time.time() - t0), flush=True)
    splitidx = np.array([i for i, Pr in enumerate(P) if Pr.kind == "split"])
    for Q in Qs:
        ma, mb, mN = elements_in_annulus(Q / 2, Q)
        Om = bump(mN / Q, 0.5, 1.0); keep = Om > 1e-12
        ma, mb, mN = ma[keep], mb[keep], mN[keep]; nm = len(ma)
        J = np.zeros((nP, nm), np.int8)
        blocks = []
        B = 2000
        for s0 in range(0, len(splitidx), B):
            ids = splitidx[s0:s0 + B]
            ps = np.array([P[i].q for i in ids], np.int64); rs = np.array([P[i].r for i in ids], np.int64)
            roots = np.array([P[i].roots for i in ids], np.int64)
            blocks.append((ps, rs, roots, (-ma) % 1, None))
        args = []
        for s0 in range(0, len(splitidx), B):
            ids = splitidx[s0:s0 + B]
            ps = np.array([P[i].q for i in ids], np.int64); rs = np.array([P[i].r for i in ids], np.int64)
            roots = np.array([P[i].roots for i in ids], np.int64)
            args.append((ps, rs, roots, (-ma).astype(np.int64), (-mb).astype(np.int64)))
        with Pool(11) as pool:
            outs = pool.map(symblock, args)
        for s0, o in zip(range(0, len(splitidx), B), outs):
            J[splitidx[s0:s0 + B]] = o
        for i, Pr in enumerate(P):
            if Pr.kind == "inert": J[i] = Pr.symbol(-ma, -mb).astype(np.int8)
        print("Q=%d rows %d, symbol table t=%.1f" % (Q, nm, time.time() - t0), flush=True)
        A = np.zeros((K, nm), complex)
        for (i, wg, idx) in slist:
            rows = J[idx]
            bad = (rows < 0).any(axis=0)
            e = rows.astype(np.int64).sum(axis=0) % 6
            val = np.conj(Z6[e]); val[bad] = 0
            A[i] += wg * val
        np.save(OUT + "A_Q%d_K%d.npy" % (Q, K), A); np.save(OUT + "m_Q%d_K%d.npy" % (Q, K), np.stack([ma, mb, mN]))
        print("Q=%d A built t=%.1f" % (Q, time.time() - t0), flush=True)
