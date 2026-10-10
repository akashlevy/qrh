import numpy as np
from gmodel import count
ell=0.1773; b=0.075; M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; L=ell/h
print("L=",L)
def detail(dl,x,L,mk=0.5,cP=2/9,alpha=5/6):
    q=dl*x; best=(9,None)
    for t in np.linspace(1,1.5,101):
        worst=(-9,None)
        for r in np.linspace(max(0,t-0.5),t,201):
            m=t-r
            zp=min(L,max(0,(1-2*m)*cP)); bb=1-dl*(2*m)-2*q*zp
            if r<1:
                zm=min(L,max(0,(1-r)*mk)); zm=min(zm,max(0,(3-2*r)/8))
                a=1-dl*r-2*q*zm; v=min(a,bb); which='M' if a<bb else 'P'
            else:
                c=1-dl+(alpha-dl)*(r-1); v=min(c,bb); which='L' if c<bb else 'P'
            if v>worst[0]: worst=(v,(r,m,which,zm if r<1 else 0,zp))
        if worst[0]<best[0]: best=(worst[0],(t,)+worst[1])
    return best
for dl in [0.3,0.4,0.47,0.55,0.6]:
    for x in [0,0.25,0.5]:
        print(dl,x,detail(dl,x,L))
