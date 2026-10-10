import numpy as np
from jmodel import worst
ell=0.1773; B=41/47
gd=np.linspace(0.30,0.60,31); gx=(0.3,0.4,0.45,0.5)
for b in [0.0,0.025,0.05,0.075]:
    lo,hi=0.0,0.15
    for _ in range(9):
        mid=(lo+hi)/2; w=worst(ell,b,B,use=('num4','pw','N6'),eta6=mid,grid_d=gd,grid_x=gx)
        if w[0]<0: lo=mid
        else: hi=mid
    w0=worst(ell,b,B,use=('num4','pw'),grid_d=gd,grid_x=gx)
    print(f"b={b:.3f}: paper worst={w0[0]:+.5f} {w0[1]};  N6 threshold eta*={lo:.4f}",flush=True)
