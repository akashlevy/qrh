"""Anatomy of the count at the worst bins: which t, which adversary split (r,m), which witness binds."""
import numpy as np
from gmodel import E, worst
zc=17/50
def count_detail(dl,x,L,mk=0.5,cP=2/9,alpha=5/6,cap82=True,ts=np.linspace(1,1.5,26)):
    q=dl*x; best=(9,None)
    for t in ts:
        worst=(-9,None)
        for r in np.linspace(max(0,t-0.5),t,61):
            m=t-r
            zp=min(L,max(0,(1-2*m)*cP)); b=1-dl*(2*m)-2*q*zp
            if r<1:
                zm=min(L,max(0,(1-r)*mk))
                if cap82: zm=min(zm,max(0,(3-2*r)/8))
                a=1-dl*r-2*q*zm; v=min(a,b); wit='a' if a<b else 'b'
            else:
                c=1-dl+(alpha-dl)*(r-1); v=min(c,b); wit='c' if c<b else 'b'
            if v>worst[0]: worst=(v,(t,r,m,zp,wit,a if r<1 else c,b))
        if worst[0]<best[0]: best=worst
    return best
ell=0.1773;B=41/47
for b in [0.05,0.075]:
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; L=ell/h
    print('b',b,'L',round(L,4))
    for dl in [0.3,0.4,0.45,0.5,0.55,0.6]:
        for x in [0.25,0.5]:
            e,R=E(ell,b,B,dl,x)
            v,info=count_detail(dl,x,L)
            print(f" dl={dl} x={x} E={e:+.4f} R={v:.4f} t={info[0]:.2f} r={info[1]:.3f} m={info[2]:.3f} zp={info[3]:.3f} bind={info[4]} a/c={info[5]:.4f} b={info[6]:.4f}")
