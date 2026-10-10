"""Closed forms (MODEL) vs LP.
(1) Low-side base: signal/low <= Z^{beta-1/2} P^{1/2} prod_i (ZP)^{c_i-1/2}  =>  beta > 1/2 - ell/2 + sum_i (1/2-c_i)(1+ell).
(2) Unmarked, perfect energy, GRH rows, Gram slope = c_i:
        theta_k = 1/2 + k(1/2 - 2c/3)  for equal c_i  (cubic: 1/2 + 5k/18).
    For unequal c_i the partial-copy rows can bind and the LP is higher than
        1/2 + sum_i (1/2 - c_i) + k c_min / 3 (printed as "closed" below)."""
from fractions import Fraction as Fr
import kvar_lp as K

def closed(cs):
    k = len(cs); return 0.5 + sum(0.5 - c for c in cs) + k * min(cs) / 3

cases = [("cubic k=1", [1/3]), ("cubic k=2", [1/3, 1/3]), ("cubic k=3", [1/3]*3),
         ("sextic k=1 (P78 copies)", [1/6]), ("sextic (1,1,4) k=3", [1/6, 1/6, 1/3]),
         ("sextic (1,2)+ray k=2", [1/6, 1/3]), ("mixed (1/6,1/2) k=2", [1/6, 1/2])]
for lab, cs in cases:
    k = len(cs)
    cfg = K.default_cfg({1: "cubic1", 2: "cubic2", 3: "cubic3"}[k])
    cfg.update(c=cs, gb=None, gq=0.0, ls=3, caseA=None, Eperfect=True, third=False, marks=False,
               margin=0.0, margin_y=0.0, ymax=40, lmax=40, Fge1=False)
    s = K.solve(cfg)
    print(f"{lab:28s} LP={s['theta']:.6f}  closed={closed(cs):.6f} ({Fr(closed(cs)).limit_denominator(100)})  "
          f"low-base={0.5+sum(0.5-c for c in cs):.4f}  l={[round(v,3) for v in s['l']]} y={[round(v,3) for v in s['y']]}")
