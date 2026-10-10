"""Discriminating tests of ../certificate.py: the paper geometry/margin passes, a too-large margin and a too-aggressive ell fail.
Expected output: True / [False] / [False].  (Rescued from the Oct-8 session scratchpad `sanity.py`.)"""
import sys, importlib, re, os
from fractions import Fraction as F
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "certificate.py")).read()
def run(ell, sigma0, margin):
    s=re.sub(r"(?m)^ell   = .*", f"ell   = F({ell.numerator},{ell.denominator})", src)
    s=re.sub(r"(?m)^sigma0 = .*", f"sigma0 = F({sigma0.numerator},{sigma0.denominator})", s)
    s=re.sub(r"(?m)^MARGIN = .*", f"MARGIN = F({margin.numerator},{margin.denominator})", s)
    s=s.replace('if __name__ == "__main__":','if False:')
    g={}; exec(s,g)
    return {n:ok for n,ok,_ in g["checks"]}
# 1) paper geometry, paper margin -> must pass
r=run(F(1,6),F(7,8),F(49,440640)); print("paper geometry, margin 49/440640:", all(r.values()))
# 2) paper geometry, too-large margin 3e-4 -> endpoint must fail
r=run(F(1,6),F(7,8),F(3,10000)); print("paper geometry, margin 3e-4 endpoint:", [ok for n,ok in r.items() if n.startswith("endpoint")])
# 3) too-aggressive ell (margin goes negative ~ -1e-5) -> must fail
ell=F(1,6)+F(1,5000); s0=F(11,12)-ell/4
r=run(ell,s0,F(0)); print("ell=1/6+1/5000, sigma0=11/12-ell/4, margin 0 endpoint:", [ok for n,ok in r.items() if n.startswith("endpoint")])
