# independent brute-force check of the joint plain lever c (cP = cN = c) at ell=0.1773, B=41/47, kappa=3/4
import numpy as np, sys
sys.path.insert(0, '../plain-cap')
zc=17/50
def count_bf(dl,x,L,c,mk=0.5,alpha=5/6,nt=2001,nr=2001):
    q=dl*x
    t=np.linspace(1,1.5,nt)[:,None]; s=np.linspace(0,1,nr)[None,:]
    r=(t-0.5)+0.5*s; m=t-r
    zp=np.clip((1-2*m)*c,0,L); b=1-2*dl*m-2*q*zp
    zm=np.clip((1-r)*mk,0,L); zm=np.minimum(zm,np.maximum(0,(3-2*r)/8))
    a=np.where(r<1,1-dl*r-2*q*zm,1-dl+(alpha-dl)*(r-1))
    v=np.minimum(a,b)
    return max(float(np.min(np.max(v,axis=1))),1-dl)
def worst(c,b,ell=(1-12*(41/47-5/6))/3,B=41/47):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h; n0=1-ly/d
    w=-9
    for dl in np.arange(0.40,0.56,0.01):
        x=0.5; q=dl*x; a=(1+dl)/2; R=count_bf(dl,x,L,c)
        zN=min(L,max(0,(1-2*n0)*c))
        base=a-B+h*(zc-1/6)-ell/2+q*ell-ly/2
        E=base+d*(R+min((1-R)/4-q*zN/2,dl*n0/2)-zc)
        w=max(w,E)
    return w
for c in [0.62,0.655,0.70,0.73]:
    print(c, min((worst(c,b),b) for b in np.linspace(0,0.1,11)))
