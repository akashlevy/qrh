import numpy as np
from wmodel import R
zc=17/50
def worstE(M,ell,b,kappa=0.75,B=None):
    lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell
    if B is None: B=2/3+M/4
    ds=np.linspace(1/50,3/4,300); xs=np.linspace(0,0.5,51)
    D,X=np.meshgrid(ds,xs); q=D*X; a=(1+D)/2; d=h; Rv=R(D,X)
    base=a-B+h*(zc-1/6)-ell/2+q*ell
    Eold=base-a*ly+d*(Rv+D/2-zc)
    zP=max(0,(2*ly-d)/(6*kappa)/d)  # in U units: (1-2n0)/(6k), n0=1-ly/d
    Ewt=base-ly/2+d*((3*Rv+1)/4-zc)-q*(2*ly-d)/(12*kappa)
    E=np.minimum(Eold,Ewt)
    return E.max()
M=116/141
for ell in [0.1773,0.19,0.2,0.22,0.25,0.3]:
    bs=np.linspace(0,0.3,61); v=[worstE(M,ell,b) for b in bs]; j=int(np.argmin(v))
    print(f"M={M:.4f} ell={ell:.4f} wall excess w={M+ell-1:+.4f} best b={bs[j]:.3f} worstE={v[j]:+.5f}")
print('--- vary ell below the wall at fixed M')
for ell in [0.0,0.05,0.1,0.13,0.15,0.17]:
    bs=np.linspace(0,0.3,61); v=[worstE(M,ell,b) for b in bs]; j=int(np.argmin(v))
    print(f"M={M:.4f} ell={ell:.4f} best b={bs[j]:.3f} worstE={v[j]:+.5f}")
