import numpy as np
def anat(dl,x,L=1.0,mk=0.5,cP=2/9,alpha=5/6,cap82=True,ts=np.linspace(1,1.5,501)):
    q=dl*x; best=(9,None)
    for t in ts:
        worst=(-9,None)
        for r in np.linspace(max(0,t-0.5),t,2001):
            m=t-r
            zp=min(L,max(0,(1-2*m)*cP)); b=1-dl*(2*m)-2*q*zp
            if r<1:
                zm=min(L,max(0,(1-r)*mk))
                if cap82: zm=min(zm,max(0,(3-2*r)/8))
                a=1-dl*r-2*q*zm; v=min(a,b); kind='inv' if a<b else 'plain'
            else:
                c=1-dl+(alpha-dl)*(r-1); v=min(c,b); kind='long' if c<b else 'plain'
            if v>worst[0]: worst=(v,(t,r,m,kind,a if r<1 else c,b,zm if r<1 else 0,zp))
        if worst[0]<best[0]: best=worst
    return best
if __name__=="__main__":
    ell=0.1773;b=0.075;M=1-ell;lx=(M-b)/2;h=1-lx+ell;L=ell/h
    print("L=",L)
    for dl in [0.3,0.4,0.47,0.55,0.6]:
        for kw in [dict(),dict(mk=0.8,cP=0.4,cap82=False),dict(mk=0.75,alpha=0.75,cap82=False),dict(mk=1,cap82=False)]:
            v,info=anat(dl,0.5,L=L,**kw)
            print(dl,kw,round(v,4),tuple(round(z,4) if isinstance(z,float) else z for z in info))
