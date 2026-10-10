"""Vectorized version of gm2 (same formulas as gmodel.py), for fine grids.
count: R = min_t max_{r in [t-1/2,t]} min(witness_inv_or_long(r), plain(t-r)), floored at 1-delta."""
import numpy as np
zc=17/50
def count_v(dl,x,L,cM=2.0,wM=0.0,cP=2/9,alpha=5/6,cap82=True,zfun=None,nt=101,nr=201):
    q=dl*x
    t=np.linspace(1,1.5,nt)[:,None]; s=np.linspace(0,1,nr)[None,:]
    r=np.maximum(0,t-0.5)+(t-np.maximum(0,t-0.5))*s; m=t-r
    zp=np.clip((1-2*m)*cP,0,L); b=1-dl*2*m-2*q*zp
    zm=np.clip(zfun(r) if zfun is not None else (1+wM-r)/cM,0,L)
    if cap82: zm=np.minimum(zm,np.maximum(0,(3-2*r)/8))
    a=1-dl*r-2*q*zm
    c=1-dl+(alpha-dl)*(r-1)
    w=np.where(r<1,a,c)
    v=np.minimum(w,b)
    return max(float(np.min(np.max(v,axis=1))),1-dl)
def E(ell,b,B,dl,x,kappa=0.75,wall=0.0,**kw):
    M=1-ell+wall; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    q=dl*x; a=(1+dl)/2; Rv=count_v(dl,x,L,**kw)
    base=a-B+h*(zc-1/6)-ell/2+q*ell
    Eold=base-a*ly+d*(Rv+dl/2-zc)
    n0=1-ly/d; zN=min(L,max(0,(1-2*n0)/(6*kappa)))
    sig=(1-Rv)/4-2*q*zN/4
    Ewt=base-ly/2+d*(Rv+sig-zc)
    return min(Eold,Ewt),Rv
def worst(ell,b,B,gd,gx,**kw):
    w=-9;arg=None
    for dl in gd:
        for x in gx:
            e,Rv=E(ell,b,B,dl,x,**kw)
            if e>w: w=e;arg=(round(float(dl),3),x,round(Rv,4))
    return w,arg
def best(ell=0.1773,B=41/47,bs=np.linspace(0,0.10,21),gd=np.linspace(0.2,0.75,56),gx=(0,0.2,0.3,0.4,0.45,0.5),**kw):
    return min((worst(ell,b,B,gd,gx,**kw)+(round(float(b),3),) for b in bs),key=lambda z:z[0])
