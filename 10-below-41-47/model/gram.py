import numpy as np
from gmodel import count
zc=17/50; kappa=0.75
def worst_g(ell,b,g,M=None):
    if M is None: M=1-ell
    lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    B=5/6+lx/6-ell/6+g*b/12
    w=-9
    for dl in np.linspace(0.05,0.75,15):
        for x in (0,0.25,0.5):
            q=dl*x; a=(1+dl)/2; Rv=count(dl,x,L)
            base=a-B+h*(zc-1/6)-ell/2+q*ell
            Eold=base-a*ly+d*(Rv+dl/2-zc)
            n0=1-ly/d; zN=min(L,max(0,(1-2*n0)/(6*kappa)))
            Ewt=base-ly/2+d*(Rv+(1-Rv)/4-q*zN/2-zc)
            w=max(w,min(Eold,Ewt))
    return w,B
# find best theta for each Gram factor g: bisection on ell with best b (gram domination ly-ell-11b/6>0 ignored->report)
for g in [1.0,0.9,0.85,0.8,0.75,0.5]:
    best=None
    for ell in np.linspace(0.15,0.2,26):
        for b in np.linspace(0,0.3,31):
            w,B=worst_g(ell,b,g)
            if w<0 and (best is None or B<best[0]): best=(B,ell,b,w)
    print(f"g={g}: best theta={best[0]:.5f} ell={best[1]:.4f} b={best[2]:.3f} (41/47={41/47:.5f})")
