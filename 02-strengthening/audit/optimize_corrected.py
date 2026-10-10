"""Float model of theta under W(w), with switchable corrections:
  A: low-side loss uses Lemma reflected-energy with powerful part O free: (1+3l'-2M')/4 (when 2M'+l'>=1)
  B: inverse slot capacity z <= min((1+w/2-r)/2, 1-r, (3-2r)/8, S), only r<1 (soft W(w) LP result)
  C: amplification NOT relaxed: e(r)=max(1,(1+5r)/6)
"""
import sys, json, numpy as np
from scipy.optimize import minimize
CP = 2/9; ZETA = 1e-4; Z0 = 17/50
def Rcount(d, x, w, S, fixB, fixC, nt=41, nr=121):
    q = x*d; best = 9
    for t in np.linspace(1, 1.5, nt):
        r = np.linspace(max(t-0.5, 0), t, nr); m = np.minimum(np.maximum(t-r, 0), 0.5)
        if fixB:
            zM = np.clip(np.minimum.reduce([(1+w/2-r)/2, 1-r, (3-2*r)/8, np.full_like(r, S)]), 0, None)
            bI = np.where(r < 1, 1-d*r-2*q*zM, 9)
        else:
            zM = np.clip(np.minimum.reduce([(1+w-r)/2, (3-2*r)/8, np.full_like(r, S)]), 0, None)
            bI = np.where(r < 1+w, 1-d*r-2*q*zM, 9)
        ww = 0.0 if fixC else w
        bA = np.maximum(1, (1+5*r/(1+ww))/6) - d*r
        zP = np.minimum(CP*(1-2*m), S)
        bP = np.where(m < 0.5, 1-2*d*m-2*q*zP, 9); bP0 = 1-2*d*m
        best = min(best, np.min(np.vstack([bI, bA, bP, bP0]), axis=0).max())
    return best
def elow(lx, ly, ell, w, fixA):
    M = lx+ly; worst = -9
    for d in np.linspace(0, ell, 41):
        Mp, lp = M-2*d, ell-d
        if fixA:
            e = max(0, (1+3*lp-2*Mp)/4 if 2*Mp+lp >= 1 else 0, Mp+lp-1-w)
        else:
            e = max(0, (1+2*lp-2*Mp)/3, Mp+lp-1-w)
        worst = max(worst, -d + e/2)
    return worst
def parts(p, w, dg, xg, nt, nr, fix, d0=1/50):
    lx, ly, ell = p; h = 1-lx+ell; b = ly-lx; M = lx+ly; S = ell/(h+ZETA)
    low = 5/6+lx/6-ell/6+0.5*max(0, b/6, 2*b-ly)+elow(lx, ly, ell, w, fix[0])
    hi = max((1+d)*(1-ly)/2-ell/2+x*d*ell+(h+ZETA)*(Rcount(d, x, w, S, fix[1], fix[2], nt, nr)+d/2-1/6) for d in dg for x in xg)
    fl = (1+d0)*(1-ly)/2-ell/2+(d0/2)*ell+(h+ZETA)*(1+d0/2-1/6)
    small = h*(Z0-1/6)-ly/2+0.0126
    gram = 0 if b <= 0 else max(0, -(ly-ell-11*b/6))
    return low, hi, fl, small + 0*gram, gram
def theta(p, w, dg, xg, nt, nr, fix):
    if min(p) < 0 or p[2] > p[0]: return 9
    low, hi, fl, small, gram = parts(p, w, dg, xg, nt, nr, fix)
    return max(low, hi, fl) + 10*max(0, small) + 10*gram
if __name__ == "__main__":
    w = float(sys.argv[1]); fix = tuple(c == "1" for c in sys.argv[2])
    dgc = np.linspace(1/50, 5/6, 30); xgc = np.linspace(0, 0.5, 9)
    starts = [(17/48, 23/48, 1/6), (0.36+0.45*w, 0.48-0.25*w, 0.167+0.7*w), (0.4+0.25*w, 0.4+0.25*w, 0.2+0.55*w), (0.38, 0.45, 0.19)]
    best = min((minimize(theta, s, args=(w, dgc, xgc, 17, 51, fix), method="Nelder-Mead",
                         options=dict(xatol=1e-5, fatol=1e-7, maxiter=400)) for s in starts), key=lambda r: r.fun)
    p = best.x
    low, hi, fl, small, gram = parts(p, w, np.linspace(1/50, 5/6, 100), np.linspace(0, .5, 21), 41, 121, fix)
    print(json.dumps(dict(w=w, fixes_ABC=sys.argv[2], theta=max(low, hi, fl), low=low, high=hi, floor=fl, lx=p[0], ly=p[1], ell=p[2])))
