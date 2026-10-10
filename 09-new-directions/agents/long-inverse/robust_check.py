# Brute-force check (x = 0): adding a hypothetical amplified MIXED count for (M_r S_m),
#   f_mix = max(1,(1+5s)/6) - delta s,  s = r + m  (sixth-power copies applied to the product,
#   assuming a full-family diagonal mixed moment for H >= (DN)^{1+c}),
# does not change R(delta,0) = 1 - 2 delta alpha/(3 alpha - delta), alpha = 5/6.
import numpy as np
def R(delta, mix=False, alpha=5/6):
    best = 9
    for t in np.linspace(1, 1.5, 501):
        r = np.linspace(max(t-0.5,0), t, 1201)[:, None]
        m = np.linspace(0, 0.5, 601)[None, :]
        ok = (r + m >= t - 1e-12)
        fA = np.where(r >= 1, 1 + alpha*(r-1), 1) - delta*r
        fP = 1 - 2*delta*m
        F = np.minimum(fA, fP)
        if mix:
            s = r + m
            F = np.minimum(F, np.where(s >= 1, 1 + alpha*(s-1), 1) - delta*s)
        v = np.where(ok, F, -9).max()
        best = min(best, v)
    return best
for d in [0.3, 0.39, 0.45]:
    print(d, round(R(d), 5), round(R(d, True), 5), 'closed form', round(1 - 2*d*(5/6)/(2.5 - d), 5))
