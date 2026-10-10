import numpy as np, sys
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
from zomega import *
from scipy.special import j0, gamma
from scipy.integrate import quad
# (a) Poisson normalization with self-dual measure: sum_m f(m) e(-m t/u) = sum_{j = t mod u} fhat(j/u)
R = 60
a, b = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij"); a = a.ravel(); b = b.ravel()
z = cval(a, b)
c = 0.05
f = lambda x: np.exp(-np.pi * c * np.abs(x) ** 2)
fhat = lambda xi: 2 / (np.sqrt(3) * c) * np.exp(-(4 / 3) * np.pi * np.abs(xi) ** 2 / c)
e = lambda x: np.exp(4j * np.pi * np.imag(x) / np.sqrt(3))
u = cval(2, 3); t = cval(1, -1)   # u = 2+3w (norm 7)
lhs = np.sum(f(z) * e(-z * t / u))
# j = t + u*l
rhs = np.sum(fhat((t + u * z) / u))
print("Poisson check:", lhs, rhs)
# (b) ghat(0) formula check
def bump(t, lo, hi):
    t = np.asarray(t, float); out = np.zeros_like(t); ins = (t > lo) & (t < hi); uu = (t[ins] - lo) / (hi - lo); out[ins] = np.exp(-1 / (uu * (1 - uu)) + 4); return out
Om = lambda t: bump(t, 0.5, 1.0)
rho = np.linspace(np.sqrt(0.5), 1, 4001)
def H(t):
    # H(t) = Omega^(eta), |eta|^2 = t : (2/sqrt3) * 2pi int Omega(rho^2) J0(2 pi rho (2/sqrt3) sqrt t) rho drho
    k = 2 * np.pi * (2 / np.sqrt(3)) * np.sqrt(t)
    return (2 / np.sqrt(3)) * 2 * np.pi * np.trapezoid(Om(rho ** 2) * j0(k * rho) * rho, rho)
H0 = H(0.0); print("H(0) =", H0, " expected (2/sqrt3) pi int Om =", (2 / np.sqrt(3)) * np.pi * quad(lambda x: Om(x)[()] if np.ndim(Om(x)) else Om(x), 0.5, 1)[0])
# int_0^inf H(t) t^{-5/6} dt numerically (t up to large; H Schwartz)
ts = np.concatenate([np.linspace(1e-8, 1, 20001), np.linspace(1, 400, 200001)[1:]])
Hs = np.array([H(t) for t in ts[::1]]) if False else None
# substitute t = v^6: int H(t) t^{-5/6} dt = 6 int H(v^6) dv
vs = np.linspace(0, 400 ** (1 / 6), 60001)
Hv = np.array([H(v ** 6) for v in vs])
I_num = 6 * np.trapezoid(Hv, vs)
IOm = quad(lambda x: float(Om(np.array([x]))[0]) * x ** (-1 / 6), 0.5, 1)[0]
I_th = (2 / np.sqrt(3)) ** (2 / 3) * np.pi ** (2 / 3) * gamma(1 / 6) / gamma(5 / 6) * IOm
print("int H t^{-5/6}: numeric %.6f  formula %.6f" % (I_num, I_th))
ghat0 = (2 * np.pi / np.sqrt(3)) / 6 * I_th
print("ghat(0) =", ghat0)
# (c) lattice sum (1/6) sum_{a != 0} H(|a|^12/T) vs (1/6) ghat0 T^{1/6}
from scipy.interpolate import interp1d
Hint = interp1d(vs ** 6, Hv, bounds_error=False, fill_value=0.0)
for T in [1e3, 1e5, 1e7, 1e9]:
    Rr = int(2 * (400 * T) ** (1 / 12)) + 3
    aa, bb = np.meshgrid(np.arange(-Rr, Rr + 1), np.arange(-Rr, Rr + 1), indexing="ij")
    N = norm(aa, bb).astype(float).ravel(); N = N[N > 0]
    s = np.sum(Hint(N ** 6 / T)) / 6
    print("T=%.0e  lattice %.5f  main %.5f  diff %.5f  -H(0)/6 = %.5f" % (T, s, ghat0 * T ** (1 / 6) / 6, s - ghat0 * T ** (1 / 6) / 6, -H0 / 6))
