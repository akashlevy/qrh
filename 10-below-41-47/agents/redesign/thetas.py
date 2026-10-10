import numpy as np
from rmodel import des_paper, des_mollifier, des_Lambda, des_product, des_sparseLS, worst
from frac import des_frac
def theta_reach(design,lo=0.0,hi=0.21,bs=np.linspace(0,0.175,8)):
    for _ in range(13):
        mid=(lo+hi)/2; B=11/12-mid/4
        w=min(worst(design,mid,b,B)[0] for b in bs)
        if w<=0: lo=mid
        else: hi=mid
    return lo,11/12-lo/4
for name,f in [('paper weighted',des_paper),('(a1) mollifier th=0.01',des_mollifier(0.01)),
               ('(a2) Lambda',des_Lambda),('(b) product t=n0',des_product('n0')),
               ('(b) product t=n0 cP=1',des_product('n0',cP=1.0)),
               ('HYP frac diag cP=1',des_frac(1.0)),('HYP sparse LS on C',des_sparseLS)]:
    l,t=theta_reach(f); print(f"{name:26s} ell_max={l:.4f} theta={t:.5f}")
