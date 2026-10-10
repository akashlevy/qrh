# Independent check of the slot-free high side at (lx,ly,ell)=(989/5000,3141/5000,0), sigma0=1733/2000
from fractions import Fraction as F
import sys
sys.argv = [sys.argv[0], "lx=989/5000", "ly=3141/5000", "ell=0", "sigma0=1733/2000", "theta=1", "dropq=1", "wE=0"]
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import exact_cert_theta as C   # module-level checks run on import (cheap); certify() only under __main__
lx, ly, ell, s0 = F(989,5000), F(3141,5000), F(0), F(1733,2000)
h = 1 - lx + ell; z0 = F(17,50); zeta = F(1,10000); dmin = F(1,100); al = F(5,6)
def R0(d):          # closed form min_t max(L(t), 1-2 d t/3), t* = 5/(5-2d)
    return (15 - 16*d) / (15 - 6*d)
def R0_brute(d, n=3000):
    best = None
    for i in range(n+1):
        t = 1 + F(i, 2*n)
        L = 1 - d + (al - d)*(t - 1)
        v = max(min(L, F(1)), 1 - 2*d*t/3)
        best = v if best is None else min(best, v)
    return best
# 1. closed form vs checker's exact R_t (S=0) at a grid of delta
mx = F(0)
for k in range(1, 84):
    d = F(k, 100)
    if d > al: break
    tstar = 5/(5 - 2*d)
    rc = C.R_t(d, F(0), tstar)
    mx = max(mx, abs(rc - R0(d)))
print("max |R_t(t*) - R0| over delta grid:", float(mx))
print("R0 brute vs closed at 0.3:", float(R0_brute(F(3,10))), float(R0(F(3,10))))
# 2. high exponent E(d) = a(1-ly) - s0 + h(z0-1/6) + d(R + delta/2 - z0), ell=0
def E(d, delta, R):
    a = (1 + delta)/2
    return a*(1-ly) - s0 + h*(z0 - F(1,6)) + d*(R + delta/2 - z0)
worst = (None, None)
N = 4000
for k in range(N+1):
    delta = F(1,50) + (al - F(1,50))*F(k, N)
    for d in (dmin, h + zeta):
        e = E(d, delta, R0(delta))
        if worst[0] is None or e > worst[0]:
            worst = (e, (float(delta), float(d)))
print("max E over nonfloor bins (grid):", float(worst[0]), "at (delta,d) =", worst[1])
# derivative check: E(h+zeta, delta) concave/convex? print a few values
for dl in (F(1,50), F(1,20), F(1,10), F(1,5), F(3,10), F(1,2), F(5,6)):
    print("  delta=%.3f  R0=%.5f  E(h+zeta)=%.6f" % (float(dl), float(R0(dl)), float(E(h+zeta, dl, R0(dl)))))
# 3. floor bin (R=1, delta0=1/50) at d = h+zeta
print("floor E(h+zeta):", float(E(h+zeta, F(1,50), F(1))))
# 4. small rows: h(z0-1/6) - ly/2 + (63/50) dmin ; paper's cruder 2*dmin
print("small rows (63/50 dmin):", float(h*(z0-F(1,6)) - ly/2 + F(63,50)*dmin), " (2 dmin):", float(h*(z0-F(1,6)) - ly/2 + 2*dmin))
# 5. low side and C(s)
print("C(sigma0)=", float(s0 - F(5,6) + lx/3), " low exponent lx/2=", float(lx/2), " margin=", float(s0 - F(5,6) + lx/3 - lx/2))
# 6. principal margins
print("m_w=ly/20=", float(ly/20), " m_z=h/600=", float(h/600))
# 7. local Euler region two at x_r = sigma0 (u=1 and p|u exponents)
xr, wr, zr = s0, F(19,20), F(33,200)
good = max(4-6*xr-6*zr, 1-xr-wr-6*zr, -xr-wr, -xr-6*zr, -6*zr-wr)
ram = max(1-xr-wr, F(3,2)-3*xr, 2-3*xr-wr, 2-4*xr, F(5,2)-4*xr-wr, 3-6*xr)
print("region-two defect exponents at x_r=sigma0: good-prime", float(good), " ramified", float(ram))
# smallest x_r for which good-prime defect < -1 at z_r=33/200, w_r=19/20
print("good-prime needs x_r >", float((5 - 6*zr)/6), "(R term) and", float(2 - wr - 6*zr), "(second-family term)")
# 8. theta-pairing condition (low side) b/6 <= 2 s(M)
M = lx + ly; b = ly - lx
print("theta piece b/6+max(M-1,7M/3-2,-M/6) =", float(b/6 + max(M-1, F(7,3)*M-2, -M/6)))
