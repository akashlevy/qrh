import numpy as np
from jmodel import worst
ell=0.1773; B=41/47
def W(eta):
    return min((worst(ell,b,B,use=('num4','pw','N6'),eta6=eta,grid_d=np.linspace(0.25,0.65,17),grid_x=(0.25,0.4,0.5)),b) for b in np.linspace(0.0,0.12,13))
lo,hi=0.06,0.12
for _ in range(10):
    mid=(lo+hi)/2; w=W(mid)
    if w[0][0]<0: lo=mid
    else: hi=mid
print("eta6 threshold ~",lo,hi, W(lo))
