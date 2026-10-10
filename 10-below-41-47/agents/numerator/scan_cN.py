import numpy as np
from gmodel import worst
ell=0.1773; B=41/47
for cN in [2/9,0.3,0.4,0.5,0.7,1.0,1.5,2.0]:
    res=min((worst(ell,b,B,cN=cN,grid_d=np.linspace(0.05,0.75,15),grid_x=(0,0.25,0.5)),b) for b in np.linspace(0.0,0.15,7))
    print(f"cN={cN:.3f} worst={res[0][0]:+.5f} arg={tuple(round(float(v),3) for v in res[0][1])} b={res[1]:.3f}")
