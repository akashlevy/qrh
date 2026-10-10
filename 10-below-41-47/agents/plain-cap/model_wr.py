"""Weighted model with a conditional slot input WR(lam):
    WR(lam):  sum_k |S(n1)S(n2)Q_z|^2 << Z^{M + lam z}  whenever n1+n2 <= M   (lam = kappa: PROVED)
Moment excess available for a symmetric plain product of length 2n (n=m plain witness, n=n0 numerator):
    E(n,z) = min( (kappa z - (1-2n)/6)_+          [P78 induction + pointwise beyond capacity],
                  ((1+5 lam) z/6 - (1-2n)/6)_+     [descend to width 2n+(1-lam)z, then WR terminal],
                  lam z )                          [WR at the top]
Plain witness count  b = min_z [1 - 2 dl m - 2 q z + E(m,z)],  z in [0,L]
Numerator            sig = (1-R)/4 + min_z [E(n0,z)/4 - q z/2]."""
import numpy as np

zc = 17 / 50
ZG = np.linspace(0, 1, 41)


def Emin(n, zs, kappa, lam):
    a = np.maximum(0, kappa * zs - (1 - 2 * n) / 6)
    b = np.maximum(0, (1 + 5 * lam) * zs / 6 - (1 - 2 * n) / 6)
    c = lam * zs
    return np.minimum(np.minimum(a, b), c)


def count(dl, x, L, kappa, lam, mk=0.5, cap82=True, alpha=5/6, ts=np.linspace(1, 1.5, 26)):
    q = dl * x; best = 9
    zs = ZG * L
    for t in ts:
        worst = -9
        for r in np.linspace(max(0, t - 0.5), t, 61):
            m = t - r
            b = np.min(1 - 2 * dl * m - 2 * q * zs + Emin(m, zs, kappa, lam))
            if r < 1:
                zm = min(L, max(0, (1 - r) * mk))
                if cap82: zm = min(zm, max(0, (3 - 2 * r) / 8))
                a = 1 - dl * r - 2 * q * zm; v = min(a, b)
            else:
                c = 1 - dl + (alpha - dl) * (r - 1); v = min(c, b)
            worst = max(worst, v)
        best = min(best, worst)
    return max(best, 1 - dl)


def E(ell, b, B, dl, x, kappa=0.75, lam=0.75, **kw):
    M = 1 - ell; lx = (M - b) / 2; ly = (M + b) / 2; h = 1 - lx + ell; d = h; L = ell / h
    q = dl * x; a = (1 + dl) / 2; Rv = count(dl, x, L, kappa, lam, **kw)
    base = a - B + h * (zc - 1 / 6) - ell / 2 + q * ell
    Eold = base - a * ly + d * (Rv + dl / 2 - zc)
    n0 = 1 - ly / d
    zs = ZG * L
    sig = (1 - Rv) / 4 + np.min(Emin(n0, zs, kappa, lam) / 4 - q * zs / 2)
    Ewt = base - ly / 2 + d * (Rv + sig - zc)
    return min(Eold, Ewt), Rv


def best(ell=0.1773, B=41/47, fine=False, **kw):
    gd = np.linspace(0.05, 0.75, 29 if fine else 15)
    gx = (0, 0.125, 0.25, 0.375, 0.5) if fine else (0, 0.25, 0.5)
    res = []
    for b in np.linspace(0, 0.15, 13 if fine else 7):
        w = max(E(ell, b, B, dl, x, **kw)[0] for dl in gd for x in gx)
        res.append((w, b))
    return min(res)


if __name__ == "__main__":
    print(f"lam=kappa=0.75 (proved; must equal paper +0.01168): {best(lam=0.75)[0]:+.5f}")
    for lam in [0.5, 0.3, 0.2, 0.1, 0.05, 0.0]:
        print(f"lam={lam}: alone {best(lam=lam)[0]:+.5f}   mk=0.75: {best(lam=lam, mk=0.75, cap82=False)[0]:+.5f}", flush=True)
    def thr(**kw):
        lo, hi = 0.0, 0.75
        if best(fine=True, lam=0.0, **kw)[0] >= 0: return None
        for _ in range(10):
            mid = (lo + hi) / 2
            if best(fine=True, lam=mid, **kw)[0] < 0: lo = mid
            else: hi = mid
        return lo
    print("threshold lam* alone (fine grid):", thr())
    for mk in [0.6, 0.7, 0.75]:
        print(f"threshold lam* with mk={mk}:", thr(mk=mk, cap82=False), flush=True)
