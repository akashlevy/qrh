import numpy as np, sys, json, itertools
sys.path.insert(0, '.')
import qmodel as Q
kappa = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
cfg = Q.cfg_quartic("A", kappa=kappa, dmax=kappa)
m = Q.Model(cfg, kappa=kappa, dgrid=61, xgrid=6, Sgrid=np.linspace(0, 0.6, 31), nt=61, refine=10)
best = (9, None); rows = []
for lx in np.linspace(0.0, 1.0, 41):
    for ly in np.linspace(0.05, 1.1, 43):
        for ell in [0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3]:
            th = m.theta_at((lx, ly, ell))
            if th < best[0]:
                best = (th, (lx, ly, ell), m.parts((lx, ly, ell)))
print("kappa", kappa, "grid min theta:", best[0], "at", best[1], "parts(low,high,floor,small)", np.round(best[2], 4))
r = m.optimize(starts=[best[1]])
print(json.dumps(r))
