import numpy as np
from jmodel import worst
ell=0.1773; B=41/47
for eta in [0.0,0.02,0.04,0.06,0.08,0.10]:
    res=min((worst(ell,b,B,use=('num4','pw','N6'),eta6=eta),b) for b in np.linspace(0.0,0.15,7))
    print(f"eta6={eta:.2f} worst={res[0][0]:+.5f} arg={res[0][1]} b={res[1]:.3f}",flush=True)
