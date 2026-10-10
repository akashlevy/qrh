"""General exponent model of the weighted-numerator high side.
Knobs: mk (marked prime capacity: z<= (1-r)*mk ... paper mk=1/2), cP (plain capacity rate 1/(6kappa) paper),
alpha (long-witness slope), cN (numerator moment capacity rate), Nmom (numerator moment order: 4 paper),
DH (perfect counts), supply-limited by L=ell/h."""
import numpy as np
zc=17/50
def count(dl,x,L,mk=0.5,cP=2/9,alpha=5/6,cap82=True,DH=False,ts=np.linspace(1,1.5,26)):
    if DH: return 1-dl
    q=dl*x; best=9
    for t in ts:
        worst=-9
        for r in np.linspace(max(0,t-0.5),t,61):
            m=t-r
            # plain witness bound (m may be 0 -> trivial)
            zp=min(L,max(0,(1-2*m)*cP)); b=1-dl*(2*m)-2*q*zp
            if r<1:
                zm=min(L,max(0,(1-r)*mk)); 
                if cap82: zm=min(zm,max(0,(3-2*r)/8))
                a=1-dl*r-2*q*zm
                v=min(a,b)
            else:
                c=1-dl+(alpha-dl)*(r-1)
                v=min(c,b)
            worst=max(worst,v)
        best=min(best,worst)
    return max(best,1-dl)
def E(ell,b,B,dl,x,kappa=0.75,cN=None,Nmom=4,**kw):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    q=dl*x; a=(1+dl)/2; Rv=count(dl,x,L,**kw)
    base=a-B+h*(zc-1/6)-ell/2+q*ell
    Eold=base-a*ly+d*(Rv+dl/2-zc)
    n0=1-ly/d
    if cN is None: cN=1/(6*kappa)
    zN=min(L,max(0,(1-2*n0)*cN)) if Nmom==4 else min(L,max(0,(1-(Nmom//2)*n0)*cN))
    # Holder with Nmom-th moment: sum_C |S| <= #C^{1-1/k} (sum|S|^k|Qsel|^2)^{1/k} * U^{qL - 2 q zN/k}
    k=Nmom
    sig=(1-Rv)/k-2*q*zN/k
    Ewt=base-ly/2+d*(Rv+sig-zc)
    return min(Eold,Ewt),Rv
def worst(ell,b,B,grid_d=np.linspace(0.05,0.75,15),grid_x=(0,0.25,0.5),**kw):
    w=-9;arg=None
    for dl in grid_d:
        for x in grid_x:
            e,Rv=E(ell,b,B,dl,x,**kw)
            if e>w: w=e;arg=(dl,x,Rv)
    return w,arg
if __name__=="__main__":
    import sys
    ell=0.1773; B=41/47
    for name,kw in [('paper',{}),('mk=1',dict(mk=1,cap82=False)),('cP=1',dict(cP=1)),('DIAG',dict(mk=1,cP=1,cap82=False)),
                    ('DIAG+cN=1',dict(mk=1,cP=1,cap82=False,cN=1)),('alpha=0.7',dict(alpha=0.7)),('DH',dict(DH=True)),('N6',dict(Nmom=6))]:
        res=min((worst(ell,b,B,**kw),b) for b in np.linspace(0.0,0.15,7))
        print(f"{name:10s} worst={res[0][0]:+.5f} at delta,x,R={tuple(round(v,3) for v in res[0][1])} b={res[1]:.3f}")
