"""What survives of Idea 1 if the theta pairing gets no saving (finding F1 of V_gaps.md):
theta part bounded by Cauchy-Schwarz (Gram piece b/6 kept), but the P^2/Y Gram piece removed by (LoA6)."""
import sys
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split")
import numpy as np
import model_r4 as M
import theta_split_model as T
for name, over in [("paper (b/6 and P^2/Y)", {}),
                   ("LoA6 only: b/6 kept, P^2/Y dropped", {"gram_q": 0.0}),
                   ("LoA6 only, caseA refined", {"gram_q": 0.0, "caseA": "refined"})]:
    sg, g = M.theta(over) if hasattr(M, "theta") else (None, None)
    print(name, round(sg, 6), {k: round(v, 4) for k, v in g.items()})
# slot-free (ell = 0) with theta by CS and P^2/Y dropped: S -> 0 limit
cfg = dict(M.DEFAULT); cfg.update({"gram_q": 0.0})
for S in [1e-4, 0.01, 0.05]:
    print("S=%g" % S, round(M.lp_for_S(S, cfg), 6))
