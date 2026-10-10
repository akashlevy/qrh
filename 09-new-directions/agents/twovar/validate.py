"""Validation of kvar_lp.py at k = 1 against template_lp.py (sextic) and p78_cubic_lp.py (cubic)."""
import sys
from fractions import Fraction as Fr
sys.path.insert(0, "../new-architectures"); sys.path.insert(0, "../alt-detect")
import kvar_lp as K
import template_lp as T
import p78_cubic_lp as P

def fr(x): return Fr(x).limit_denominator(400)

rows = []
# sextic P78 GRH floor
c = K.default_cfg("sextic1"); c.update(ymax=3.0)
s = K.solve(c); t = T.solve({})
rows.append(("sextic k=1 (P78) GRH floor", s["theta"], t["theta"], "13/15"))
# sextic sensitivity: wall relaxed by 0.1
c2 = dict(c, wE=0.1); s2 = K.solve(c2); t2 = T.solve({"wE": 0.1})
rows.append(("sextic k=1, W(0.1)", s2["theta"], t2["theta"], "0.865"))
# cubic one-variable template, cB = 2,3,4, with/without third piece
for cB in (2, 3, 4):
    for third in (True, False):
        c = K.default_cfg("cubic1"); c.update(cB=cB, third=third, ymax=3.0)
        s = K.solve(c); p = P.solve(dict(P.CUBIC, cB=cB, third=third))
        rows.append((f"cubic k=1 cB={cB} third={third}", s["theta"], p["theta"], ""))
for lab, a, b, tgt in rows:
    ok = "OK" if abs(a - b) < 1e-9 else "MISMATCH"
    print(f"{lab:36s} kvar={a:.9f} ({fr(a)})  reference={b:.9f}  {ok} {tgt}")
