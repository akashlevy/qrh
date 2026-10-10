# Idealized high-side scenarios at a given ell; report worst bin over b (b>=0)
import numpy as np
from wmodel import R
zc=17/50
def scen(ell,b,kind,kappa=0.75):
    lx=(1-ell-b)/2; ly=(1-ell+b)/2; h=1-lx+ell; B=11/12-ell/4
    ds=np.linspace(1/50,3/4,300); xs=np.linspace(0,0.5,51)
    D,X=np.meshgrid(ds,xs); q=D*X; a=(1+D)/2; d=h
    Rv=R(D,X)
    if kind.endswith('DH'): Rv=1-D
    base=a-B+h*(zc-1/6)-ell/2+q*ell
    Eold=base-a*ly+d*(Rv+D/2-zc)
    Ewt=base-ly/2+d*((3*Rv+1)/4-zc)-q*(2*ly-d)/(12*kappa)
    Eideal=base-ly/2+d*(Rv-zc)     # numerator O(1) on every bad row
    if kind.startswith('paper'): E=np.minimum(Eold,Ewt)
    elif kind.startswith('ideal'): E=np.minimum(np.minimum(Eold,Ewt),Eideal)
    return E.max()
for kind in ['paper','paperDH','ideal','idealDH']:
    for ell in [0.167,0.1773,0.19,0.2]:
        bs=np.linspace(0,0.25,51); v=[scen(ell,b,kind) for b in bs]; j=int(np.argmin(v))
        print(f"{kind:8s} ell={ell:.4f} best b={bs[j]:.3f} worst E={v[j]:+.5f}")
