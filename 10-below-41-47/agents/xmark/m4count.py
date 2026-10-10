"""Count-based bound for the Moebius 4th moment: S4(r) = max_{d'} [R(d',x=0) + 2 d' r] (paper counts, physical supply,
geometry ell=0.1773, b=0.075), and the resulting Prop-X1 bound sum|MQ|^2 <~ U^{(1+S4)/2}."""
import numpy as np, json
from xm_model import count
ell=0.1773; b=0.075; M=1-ell; lx=(M-b)/2; h=1-lx+ell; L=ell/h
for r in (0.66,0.69,0.75,0.81):
    best=(-9,None)
    for d in np.linspace(0.02,0.75,74):
        v=count(d,0.0,L)+2*d*r
        if v>best[0]: best=(v,round(float(d),3))
    print(json.dumps(dict(r=r,S4=round(best[0],4),at_delta=best[1],propX1_exponent=round((1+best[0])/2,4),
                          P78_rowbranch_at_z022=round(r+0.44,4))))
