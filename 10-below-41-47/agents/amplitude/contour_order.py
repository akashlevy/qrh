# Compare max_x min(Eold,Ewt) (contour chosen per amplitude class) vs min(max_x Eold, max_x Ewt) (per delta)
import numpy as np
from gmodel import count
zc=17/50
def Es(ell,b,B,dl,x,kappa=0.75):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    q=dl*x; a=(1+dl)/2; Rv=count(dl,x,L)
    base=a-B+h*(zc-1/6)-ell/2+q*ell
    Eold=base-a*ly+d*(Rv+dl/2-zc)
    n0=1-ly/d; zN=min(L,max(0,(1-2*n0)/(6*kappa)))
    Ewt=base-ly/2+d*(Rv+(1-Rv)/4-q*zN/2-zc)
    return Eold,Ewt
ell=0.1773; B=41/47
xs=np.linspace(0,0.5,11); ds=np.linspace(0.05,0.75,29)
for b in [0.0,0.05,0.075,0.1]:
    A=np.array([[Es(ell,b,B,dl,x) for x in xs] for dl in ds])  # (nd,nx,2)
    perclass=np.max(np.min(A,axis=2),axis=1)   # per delta: max_x min
    perdelta=np.min(np.max(A,axis=1),axis=1)   # per delta: min over contour of max_x
    i1=np.argmax(perclass); i2=np.argmax(perdelta)
    print(f"b={b}: per-class worst {perclass.max():+.5f} at d={ds[i1]:.3f}; per-delta worst {perdelta.max():+.5f} at d={ds[i2]:.3f}")
