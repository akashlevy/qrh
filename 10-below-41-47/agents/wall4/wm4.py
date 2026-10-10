"""wall4 model: weighted high side (same formulas as model/gmodel.py, vectorized as in marked-cap/vm.py)
with a general LOW side

    B_low(M, ell, b) = 5/6 + M/12 - ell/6 + exc(M, ell, b)/2 ,   w := M + ell - 1,

where exc is the low-side excess (in Z-units) of the separated probe beyond the energy wall:
  'paper'      exc = max(w,0)                         (P78 lem:reflected-energy + probe-gram)
  'W'          exc = max(w - wW, 0)                   (lossless wall relaxation W(wW))
  'slope'      exc = f*max(w,0)                       (excess grows at rate f past the wall)
  'theta'      exc = max(w - wT, w - b/6, 0)          (Corollary D of 04 + rank-one wall W_theta(wT);
                                                       the rest part tolerates w <= b/6, see NOTES.md sec. 5)
Gram domination (P78 needs l_y - ell - 11b/6 >= 0 for P_a^2/Y << P_a^{1/6}) is reported/enforced
optionally ('gram=True'); Corollary D removes it.

High side: E_old/E_wt of the BRIEF, count R(delta,x) as in gmodel.count (paper knobs) on a fine grid;
optionally the marked-capacity transfer of a shared wall saving (r + 2z <= 1 + wM).
"""
import numpy as np
zc = 17/50

def count_v(dl, x, L, cM=2.0, wM=0.0, cP=2/9, alpha=5/6, cap82=True, nt=101, nr=201):
    q = dl*x
    t = np.linspace(1, 1.5, nt)[:, None]; s = np.linspace(0, 1, nr)[None, :]
    r = np.maximum(0, t-0.5) + (t-np.maximum(0, t-0.5))*s; m = t-r
    zp = np.clip((1-2*m)*cP, 0, L); bb = 1-dl*2*m-2*q*zp
    zm = np.clip((1+wM-r)/cM, 0, L)
    if cap82: zm = np.minimum(zm, np.maximum(0, (3-2*r)/8))
    a = 1-dl*r-2*q*zm
    c = 1-dl+(alpha-dl)*(r-1)
    w = np.where(r < 1, a, c)
    v = np.minimum(w, bb)
    return max(float(np.min(np.max(v, axis=1))), 1-dl)

GD = np.linspace(0.05, 0.75, 71)
GX = (0, 0.2, 0.3, 0.4, 0.45, 0.5)

def high_worst(ell, M, b, B, kappa=0.75, gd=GD, gx=GX, **kw):
    lx = (M-b)/2; ly = (M+b)/2; h = 1-lx+ell; d = h; L = ell/h
    worst = -9; arg = None
    for dl in gd:
        for x in gx:
            q = dl*x; a = (1+dl)/2; Rv = count_v(dl, x, L, **kw)
            base = a-B+h*(zc-1/6)-ell/2+q*ell
            Eold = base-a*ly+d*(Rv+dl/2-zc)
            n0 = 1-ly/d; zN = min(L, max(0, (1-2*n0)/(6*kappa)))
            Ewt = base-ly/2+d*(Rv+(1-Rv)/4-q*zN/2-zc)
            e = min(Eold, Ewt)
            if e > worst: worst = e; arg = (round(float(dl), 3), x, round(Rv, 4))
    return worst, arg

def exc(M, ell, b, mode='paper', wW=0.0, f=1.0, wT=0.0):
    w = M+ell-1
    if mode == 'paper': return max(w, 0.0)
    if mode == 'W': return max(w-wW, 0.0)
    if mode == 'slope': return f*max(w, 0.0)
    if mode == 'theta': return max(w-wT, w-b/6, 0.0)
    if mode == 'thetaW': return max(w-wT, w-wW-b/6, 0.0)   # W_theta(wT) plus full W(wW) for the rest part
    raise ValueError(mode)

def Mmax(ell, b, B, **ekw):
    """largest M with B_low(M) <= B (B_low increasing in M)."""
    lo, hi = 0.3, 1.5
    if 5/6+lo/12-ell/6+exc(lo, ell, b, **ekw)/2 > B: return None
    for _ in range(60):
        mid = (lo+hi)/2
        if 5/6+mid/12-ell/6+exc(mid, ell, b, **ekw)/2 <= B: lo = mid
        else: hi = mid
    return lo

def scan(B=41/47, ells=np.linspace(0.165, 0.215, 11), bs=np.linspace(0, 0.15, 7), Ms=None,
         gram=False, hkw=None, **ekw):
    """best (min over ell, b, M <= Mmax) of the worst high exponent at target B."""
    hkw = hkw or {}
    best = (9, None)
    for ell in ells:
        for b in bs:
            Mm = Mmax(ell, b, B, **ekw)
            if Mm is None: continue
            cands = [Mm] if Ms is None else [Mm-dM for dM in Ms if Mm-dM > 0.5]
            for M in cands:
                ly = (M+b)/2
                if gram and ly-ell-11*b/6 < 0: continue
                wv, arg = high_worst(ell, M, b, B, **hkw)
                if wv < best[0]: best = (wv, dict(ell=round(float(ell), 4), b=round(float(b), 3),
                                                   M=round(M, 4), w=round(M+ell-1, 4), arg=arg))
    return best


# ---------- faster vectorized high side (all delta at once) ----------
def count_all(gd, x, L, cM=2.0, wM=0.0, cP=2/9, alpha=5/6, cap82=True, nt=101, nr=201):
    dl = np.asarray(gd)[:, None, None]; q = dl*x
    t = np.linspace(1, 1.5, nt)[None, :, None]; s = np.linspace(0, 1, nr)[None, None, :]
    r = np.maximum(0, t-0.5) + (t-np.maximum(0, t-0.5))*s; m = t-r
    zp = np.clip((1-2*m)*cP, 0, L); bb = 1-dl*2*m-2*q*zp
    zm = np.clip((1+wM-r)/cM, 0, L)
    if cap82: zm = np.minimum(zm, np.maximum(0, (3-2*r)/8))
    a = 1-dl*r-2*q*zm
    c = 1-dl+(alpha-dl)*(r-1)
    v = np.minimum(np.where(r < 1, a, c), bb)
    R = np.min(np.max(v, axis=2), axis=1)
    return np.maximum(R, 1-np.asarray(gd))

def high_worst_fast(ell, M, b, B, kappa=0.75, gd=GD, gx=GX, **kw):
    lx = (M-b)/2; ly = (M+b)/2; h = 1-lx+ell; d = h; L = ell/h
    worst = -9; arg = None
    for x in gx:
        Rv = count_all(gd, x, L, **kw); dl = np.asarray(gd)
        q = dl*x; a = (1+dl)/2
        base = a-B+h*(zc-1/6)-ell/2+q*ell
        Eold = base-a*ly+d*(Rv+dl/2-zc)
        n0 = 1-ly/d; zN = min(L, max(0, (1-2*n0)/(6*kappa)))
        Ewt = base-ly/2+d*(Rv+(1-Rv)/4-q*zN/2-zc)
        e = np.minimum(Eold, Ewt); i = int(np.argmax(e))
        if e[i] > worst: worst = float(e[i]); arg = (round(float(dl[i]), 3), x, round(float(Rv[i]), 4))
    return worst, arg

def scanf(B=41/47, ells=np.linspace(0.165, 0.215, 11), bs=np.linspace(0, 0.15, 7), dMs=(0.0,),
          gram=False, hkw=None, **ekw):
    hkw = hkw or {}
    best = (9, None)
    for ell in ells:
        for b in bs:
            Mm = Mmax(ell, b, B, **ekw)
            if Mm is None: continue
            for dM in dMs:
                M = Mm-dM
                ly = (M+b)/2
                if gram and ly-ell-11*b/6 < 0: continue
                wv, arg = high_worst_fast(ell, M, b, B, **hkw)
                if wv < best[0]: best = (wv, dict(ell=round(float(ell), 4), b=round(float(b), 3),
                                                   M=round(M, 4), w=round(M+ell-1, 4), arg=arg))
    return best

if __name__ == "__main__":
    import sys
    B = 41/47
    print("paper low side:", scan(B))
