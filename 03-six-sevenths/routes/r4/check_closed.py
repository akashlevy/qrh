import numpy as np
from energy_lp3 import emax, closed
rng = np.random.default_rng(2)
for count in ("paper", "refined"):
    over = 0; under = 0; infeas = 0; argo = argu = None
    for _ in range(4000):
        Mp = rng.uniform(0.3, 1.2); lp = rng.uniform(0, 0.6); w = float(rng.choice([0, 0.05, 0.1, 1/7, 0.2, 0.3]))
        e = emax(Mp, lp, w, count)
        if e < -50:
            infeas += 1; continue
        e = max(e, 0.0); cf = closed(Mp, lp, w, count)
        if e - cf > over: over = e - cf; argo = (Mp, lp, w, e, cf)
        if cf - e > under: under = cf - e; argu = (Mp, lp, w, e, cf)
    print(count, "infeasible:", infeas, " LP>closed by", over, argo, "\n   closed>LP by", under, argu)
