"""Quantify plain-capacity variants in the weighted model at ell=0.1773, B=41/47.
cP: plain-witness capacity rate, cN: numerator capacity rate (both = 1/(6 kappa) in P78; a better
exchange rate in lem:plain improves BOTH).  LS=True adds the (proved) large-sieve capacity 1/2-2m."""
import numpy as np
import gmodel
from gmodel import E

ell = 0.1773; B = 41 / 47


def count_ls(dl, x, L, mk=0.5, cP=2/9, alpha=5/6, cap82=True, DH=False, ts=np.linspace(1, 1.5, 26), LS=False):
    q = dl * x; best = 9
    for t in ts:
        worst = -9
        for r in np.linspace(max(0, t - 0.5), t, 61):
            m = t - r
            zcap = max(0, (1 - 2 * m) * cP)
            if LS: zcap = max(zcap, 0.5 - 2 * m)
            zp = min(L, zcap); b = 1 - dl * (2 * m) - 2 * q * zp
            if r < 1:
                zm = min(L, max(0, (1 - r) * mk))
                if cap82: zm = min(zm, max(0, (3 - 2 * r) / 8))
                a = 1 - dl * r - 2 * q * zm; v = min(a, b)
            else:
                c = 1 - dl + (alpha - dl) * (r - 1); v = min(c, b)
            worst = max(worst, v)
        best = min(best, worst)
    return max(best, 1 - dl)


def worst(b, fine=False, LS=False, **kw):
    gd = np.linspace(0.05, 0.75, 29 if fine else 15)
    gx = (0, 0.125, 0.25, 0.375, 0.5) if fine else (0, 0.25, 0.5)
    orig = gmodel.count
    if LS:
        gmodel.count = lambda dl, x, L, **k: count_ls(dl, x, L, LS=True, **k)
    try:
        w = -9
        for dl in gd:
            for x in gx:
                e, _ = E(ell, b, B, dl, x, **kw)
                w = max(w, e)
    finally:
        gmodel.count = orig
    return w


def best(fine=False, LS=False, **kw):
    bs = np.linspace(0.0, 0.15, 13 if fine else 7)
    return min((worst(b, fine=fine, LS=LS, **kw), b) for b in bs)


if __name__ == "__main__":
    k35 = 35 / 47; c35 = 1 / (6 * k35)
    print(f"baseline (cP=cN=2/9, kappa=3/4): {best()[0]:+.5f}")
    print(f"kappa=35/47 (beta*=41/47): cP=cN={c35:.4f}: {best(cP=c35, cN=c35, kappa=k35)[0]:+.5f}")
    print(f"large-sieve option for short witnesses (proved): {best(LS=True)[0]:+.5f}")
    # threshold for joint c with fine grid (bisection)
    def f(c, **kw):
        return best(fine=True, cP=c, cN=c, **kw)[0]
    print("fine-grid baseline:", f"{best(fine=True)[0]:+.5f}")
    lo, hi = 0.5, 0.8
    for _ in range(12):
        mid = (lo + hi) / 2
        if f(mid) < 0: hi = mid
        else: lo = mid
    print(f"threshold joint c (cP=cN=c, mk=1/2): c* ~ {hi:.4f}; in exchange-rate terms a* = 1/(kappa c*) = {1/(0.75*hi):.3f}, "
          f"crossing cost lambda* = 1/(6 c*) = {1/(6*hi):.3f}")
    for mk in [0.6, 0.7, 0.75, 0.8]:
        lo, hi = 0.15, 1.2
        if f(hi, mk=mk, cap82=False) >= 0:
            print(f"mk={mk}: no c<=1.2 suffices"); continue
        for _ in range(12):
            mid = (lo + hi) / 2
            if f(mid, mk=mk, cap82=False) < 0: hi = mid
            else: lo = mid
        print(f"mk={mk}: threshold joint c* ~ {hi:.4f} (a*={1/(0.75*hi):.2f}, lambda*={1/(6*hi):.3f})")
