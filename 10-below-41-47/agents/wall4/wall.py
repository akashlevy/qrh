import numpy as np
from gmodel import count
zc=17/50
def worstW(ell,M,b,B,kappa=0.75,**kw):
    lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    w=-9
    for dl in np.linspace(0.05,0.75,15):
        for x in (0,0.25,0.5):
            q=dl*x; a=(1+dl)/2; Rv=count(dl,x,L,**kw)
            base=a-B+h*(zc-1/6)-ell/2+q*ell
            Eold=base-a*ly+d*(Rv+dl/2-zc)
            n0=1-ly/d; zN=min(L,max(0,(1-2*n0)/(6*kappa)))
            Ewt=base-ly/2+d*(Rv+(1-Rv)/4-q*zN/2-zc)
            w=max(w,min(Eold,Ewt))
    return w
B=41/47
for ell in [0.1773,0.185,0.19,0.195,0.2]:
    M=12*(B-5/6)+2*ell   # lossless beyond the wall: B=5/6+M/12-ell/6
    v=min((worstW(ell,M,b,B),b) for b in np.linspace(0,0.15,7))
    print(f"ell={ell:.4f} M={M:.4f} wall excess w={M+ell-1:+.4f} worst={v[0]:+.5f} b={v[1]:.3f}")
