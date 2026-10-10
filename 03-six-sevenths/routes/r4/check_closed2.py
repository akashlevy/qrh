import numpy as np
from energy_lp3 import emax
def closed4(Mp, lp, w, count):
    a = (1 + 3*lp - 2*Mp)/4 if count == "paper" else (1 + 2*lp - 2*Mp)/3
    t3 = (Mp + 4*lp - 2)/3 - w if 2*Mp - 1 <= lp else -9
    return max(0.0, a, Mp + lp - 1 - w, t3)
rng = np.random.default_rng(3)
for count in ("paper", "refined"):
    over = under = 0; argo = argu = None; nf = 0
    for _ in range(6000):
        Mp = rng.uniform(0.3, 1.3); lp = rng.uniform(0, 0.7); w = float(rng.choice([0, 0.05, 0.1, 1/7, 0.2, 0.3, 0.5]))
        e = emax(Mp, lp, w, count)
        if e < -50: nf += 1; continue
        e = max(e, 0.0); cf = closed4(Mp, lp, w, count)
        if e - cf > over: over = e - cf; argo = (Mp, lp, w, e, cf)
        if cf - e > under: under = cf - e; argu = (Mp, lp, w, e, cf)
    print(count, "infeasible", nf, "LP>closed4 by", over, argo, " closed4>LP by", under, argu)
