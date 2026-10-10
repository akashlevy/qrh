"""Test: bound the DENOMINATOR (Moebius sum of length Z = U^{1/d}, s-weight Z^s e^{(s+z-1)^2}) by the amplified
inverse moment e(r)=(1+5r)/6 (lem:inverse-amplification) in a joint level-set LP with the numerator, instead of
the pointwise bin bound Z^a.  Exponents in U-units; 'den' = change relative to the paper's pointwise denominator."""
import numpy as np
from jmodel import wit_bounds
zc=17/50
def joint(Rs,dl,q,L,n0,rD,kappa=0.75,use_den=True):
    zN=min(L,max(0,(1-2*n0)/(6*kappa))); c1=1-2*q*zN; lN=dl*n0/2; lD=dl*rD; eD=(1+5*rD)/6
    if not use_den: eD=99
    best=-9
    for c in {Rs,c1-4*lN,eD-lD,0.0}:
        if c>Rs or c<0: continue
        a=min(lN,(c1-c)/4); bD=min(lD,eD-c)
        if a<0 or bD<0: continue
        best=max(best,c+a+bD/2)
    return best-Rs-lD/2
def E(ell,b,B,dl,x,use_den,ts=np.linspace(1,1.5,26)):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h; n0=1-ly/d; rD=1/d
    q=dl*x; a=(1+dl)/2; base=a-B+h*(zc-1/6)-ell/2+q*ell
    V=9; Rold=9
    for t in ts:
        wv=-9; wR=-9
        for r in np.linspace(max(0,t-0.5),t,41):
            Rs=wit_bounds(dl,q,L,t,r); wR=max(wR,Rs)
            wv=max(wv,Rs+joint(Rs,dl,q,L,n0,rD,use_den=use_den))
        V=min(V,wv); Rold=min(Rold,wR)
    # old (pointwise numerator) branch: denominator moment could also be used there; keep paper's
    Eold=base-a*ly+d*(Rold+dl/2-zc)
    return min(Eold,base-ly/2+d*(V-zc))
def worst(ell,b,B,use_den,gd=np.linspace(0.3,0.7,41),gx=(0.25,0.4,0.5)):
    return max((E(ell,b,B,dl,x,use_den),round(dl,3),x) for dl in gd for x in gx)
if __name__=="__main__":
    ell=0.1773; B=41/47
    for b in [0.0,0.05,0.1]:
        print(f"b={b}: paper {worst(ell,b,B,False)}  with denominator moment {worst(ell,b,B,True)}",flush=True)
