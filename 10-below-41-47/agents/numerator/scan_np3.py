import numpy as np
from jmodel import worst
ell=0.1773; B=41/47
gd=np.linspace(0.30,0.60,31); gx=(0.3,0.4,0.45,0.5)
for key,par in [('NP3','etaP'),('NM3','etaM')]:
    best=None
    for b in [0.0,0.05,0.1]:
        lo,hi=0.0,0.6
        for _ in range(9):
            mid=(lo+hi)/2; w=worst(ell,b,B,use=('num4','pw',key),grid_d=gd,grid_x=gx,**{par:mid})
            if w[0]<0: lo=mid
            else: hi=mid
        print(f"{key} b={b:.2f} threshold {par}*={lo:.4f}",flush=True)
