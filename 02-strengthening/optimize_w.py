"""
Conservative-model optimum of theta under W(w).
Usage:  python optimize_w.py <w>     -> prints one JSON line
Model = exactly the one certified by conditional_certificate.py:
  low:  5/6 + lx/6 - ell/6 + max(0,b/6,2b-ly)/2 + e/2,  e = max(0,(1+2ell-2M)/3, M+ell-1-w)
  high: E = (1+d)(1-ly)/2 - s0 - ell/2 + x d ell + (h+zeta)(R + d/2 - 1/6)
  R(d,x) = min_t max_r min{ inverse-with-slots (r<1+w, z<=min((1+w-r)/2,(3-2r)/8,S)),
                            amplified e(r)-d r, plain-with-slots, plain }
  floor: R=1 at delta0=1/50.
"""
import sys, json, numpy as np
from scipy.optimize import minimize
CP = 2/9; ZETA = 1e-4; Z0 = 17/50
def Rcount(d, x, w, S, nt=41, nr=121):
    q = x*d; best = 9
    for t in np.linspace(1, 1.5, nt):
        r = np.linspace(max(t-0.5, 0), t, nr); m = np.minimum(np.maximum(t-r, 0), 0.5)
        zM = np.clip(np.minimum.reduce([(1+w-r)/2, (3-2*r)/8, np.full_like(r, S)]), 0, None)
        bI = np.where(r < 1+w, 1-d*r-2*q*zM, 9)
        bA = np.maximum(1, (1+5*r/(1+w))/6) - d*r
        zP = np.minimum(CP*(1-2*m), S)
        bP = np.where(m < 0.5, 1-2*d*m-2*q*zP, 9); bP0 = 1-2*d*m
        best = min(best, np.min(np.vstack([bI, bA, bP, bP0]), axis=0).max())
    return best
def parts(p, w, dg, xg, nt, nr, d0=1/50):
    lx, ly, ell = p; h = 1-lx+ell; b = ly-lx; M = lx+ly; S = ell/(h+ZETA)
    e = max(0, (1+2*ell-2*M)/3, M+ell-1-w)
    low = 5/6+lx/6-ell/6+0.5*max(0, b/6, 2*b-ly)+0.5*e
    hi = max((1+d)*(1-ly)/2-ell/2+x*d*ell+(h+ZETA)*(Rcount(d, x, w, S, nt, nr)+d/2-1/6) for d in dg for x in xg)
    fl = (1+d0)*(1-ly)/2-ell/2+(d0/2)*ell+(h+ZETA)*(1+d0/2-1/6)
    small = h*(Z0-1/6)-ly/2+0.0126
    return low, hi, fl, small
def theta(p, w, dg, xg, nt, nr, d0=1/50):
    if min(p) < 0: return 9
    low, hi, fl, small = parts(p, w, dg, xg, nt, nr, d0)
    return max(low, hi, fl) + 10*max(0, small)
if __name__ == "__main__":
    w = float(sys.argv[1]); d0 = float(sys.argv[2]) if len(sys.argv) > 2 else 1/50
    dgc = np.linspace(d0, 5/6, 40); xgc = np.linspace(0, 0.5, 11)
    starts = [(0.36+0.45*w, 0.48-0.25*w, 0.167+0.7*w), (0.4+0.25*w, 0.4+0.25*w, 0.2+0.55*w)]
    best = min((minimize(theta, s, args=(w, dgc, xgc, 21, 61, d0), method="Nelder-Mead",
                         options=dict(xatol=1e-6, fatol=1e-8, maxiter=500)) for s in starts),
               key=lambda r: r.fun)
    p = best.x
    dgf = np.linspace(d0, 5/6, 120); xgf = np.linspace(0, 0.5, 26)
    low, hi, fl, small = parts(p, w, dgf, xgf, 41, 121, d0)
    print(json.dumps(dict(w=w, delta0=d0, theta_fine=max(low, hi, fl), low=low, high=hi, floor=fl,
                          lx=p[0], ly=p[1], ell=p[2])))
