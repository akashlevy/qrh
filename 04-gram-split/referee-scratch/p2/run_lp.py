import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
sys.path.insert(0, "/private/tmp/claude-501/-Users-akashlevy-Documents-riemann/eab6979c-f095-4107-a475-7e51d88a9fe7/scratchpad/p2")
import model_r4 as M
import tsm_dropq as T
from scipy.optimize import minimize_scalar
cfg = dict(M.DEFAULT); cfg.update(dropq=True, dgrid=321, xgrid=11)
# slot-free: S = 0
sg, g = T.lp_for_S(0.0, cfg, "LSmarked", return_geom=True)
print("S=0 (slot-free), LSmarked, P^2/Y dropped: theta=%.6f" % sg, {k: round(v, 5) for k, v in g.items()})
# with slots allowed, scan S
cache = {}
Ss = np.concatenate([[0.0, 0.005, 0.01, 0.02], np.linspace(0.03, 0.6, 58)])
vals = [T.lp_for_S(S, cfg, "LSmarked", Rcache=cache) for S in Ss]
k = int(np.argmin(vals)); print("best over S grid: S=%.4f theta=%.6f" % (Ss[k], vals[k]))
for S, v in list(zip(Ss, vals))[:8]: print("   S=%.3f theta=%.6f" % (S, v))
# floor delta0 = 1/1000 variant (provable floor 1/2+eta), slot-free
cfg2 = dict(cfg); cfg2.update(delta0=1/1000)
sg2, g2 = T.lp_for_S(0.0, cfg2, "LSmarked", return_geom=True)
print("S=0, delta0=1/1000:", round(sg2, 6), {k: round(v, 5) for k, v in g2.items()})
# sensitivity: paper's t=1-only zero-slot count (R = max(1-delta, 1-2delta/3) = 1-2delta/3) -> tmax=1
cfg3 = dict(cfg); cfg3.update(tmax=1.0)
print("S=0, t=1 only (paper's stated zero-slot count):", round(T.lp_for_S(0.0, cfg3, "LSmarked"), 6))
# sensitivity: trivial floor bound removed (hypothetically no floor bin)
