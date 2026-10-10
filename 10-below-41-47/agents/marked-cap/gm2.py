"""gmodel variant with marked capacity z_M(r)=min(L,(1+wM-r)/cM) [paper cM=2,wM=0, cap82],
plain capacity z_P(m)=min(L,(1-2m)*cP), long slope alpha; also optional 'xinv' = marks act as full
Mobius length when x=1/2 (z_M(r)=1-r) etc."""
import numpy as np
zc=17/50
def count(dl,x,L,cM=2.0,wM=0.0,cP=2/9,alpha=5/6,cap82=True,ts=np.linspace(1,1.5,26),nr=61,zfun=None):
    q=dl*x; best=9
    for t in ts:
        worst=-9
        for r in np.linspace(max(0,t-0.5),t,nr):
            m=t-r
            zp=min(L,max(0,(1-2*m)*cP)); b=1-dl*(2*m)-2*q*zp
            if r<1:
                if zfun is not None: zm=min(L,max(0,zfun(r)))
                else: zm=min(L,max(0,(1+wM-r)/cM))
                if cap82: zm=min(zm,max(0,(3-2*r)/8))
                a=1-dl*r-2*q*zm; v=min(a,b)
            else:
                c=1-dl+(alpha-dl)*(r-1); v=min(c,b)
            worst=max(worst,v)
        best=min(best,worst)
    return max(best,1-dl)
def E(ell,b,B,dl,x,kappa=0.75,**kw):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    q=dl*x; a=(1+dl)/2; Rv=count(dl,x,L,**kw)
    base=a-B+h*(zc-1/6)-ell/2+q*ell
    Eold=base-a*ly+d*(Rv+dl/2-zc)
    n0=1-ly/d; zN=min(L,max(0,(1-2*n0)/(6*kappa)))
    sig=(1-Rv)/4-2*q*zN/4
    Ewt=base-ly/2+d*(Rv+sig-zc)
    return min(Eold,Ewt),Rv
def worst(ell,b,B,grid_d=np.linspace(0.05,0.75,15),grid_x=(0,0.25,0.5),**kw):
    w=-9;arg=None
    for dl in grid_d:
        for x in grid_x:
            e,Rv=E(ell,b,B,dl,x,**kw)
            if e>w: w=e;arg=(round(dl,3),x,round(Rv,4))
    return w,arg
def best(ell=0.1773,B=41/47,bs=np.linspace(0.0,0.15,7),**kw):
    return min((worst(ell,b,B,**kw),b) for b in bs)
if __name__=="__main__":
    for kw in [dict(),dict(cM=2,wM=0.1,cap82=False),dict(cM=2,wM=0.13,cap82=False),dict(cM=2,wM=0.15,cap82=False),dict(cM=2,wM=0.2,cap82=False),
               dict(cM=1.5,cap82=False),dict(cM=1.4,cap82=False),dict(cM=1.25,cap82=False),dict(cM=1,cap82=False),
               dict(cM=1.25,cP=0.4,cap82=False),dict(cM=1.4,cP=0.4,cap82=False),dict(cM=2,wM=0.1,cP=0.4,cap82=False),
               dict(wM=0.0,cap82=False),dict(cM=2,wM=0.2,cP=0.4,cap82=False)]:
        (w,arg),b=best(**kw)
        print(kw,f"{w:+.5f}",arg,b)
