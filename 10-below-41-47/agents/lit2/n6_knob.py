"""lit2: thresholds. (a) sixth-moment loss lambda: sum_u |S_u(n0)|^6 << U^{1+lambda};
(b) constant zero-density exponent A on sigma>=0.5 combined with Ingham (min).
Uses gmodel.count (paper counts)."""
import numpy as np, gmodel as G
ell=0.1773; B=41/47; kappa=0.75
cache={}
def R(dl,x,L):
    k=(round(dl,4),round(x,4),round(L,5))
    if k not in cache: cache[k]=G.count(dl,x,L)
    return cache[k]
def E(b,dl,x,lam6=None,A=None):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    q=dl*x; a=(1+dl)/2; Rv=R(dl,x,L)
    if A is not None: Rv=min(Rv,max(1-dl,min(3/(2-a),A)*(1-a)))
    base=a-B+h*(G.zc-1/6)-ell/2+q*ell
    Eold=base-a*ly+d*(Rv+dl/2-G.zc)
    n0=1-ly/d; zN=min(L,max(0,(1-2*n0)/(6*kappa))); sig4=(1-Rv)/4-2*q*zN/4
    sig=sig4
    if lam6 is not None:
        z6=min(L,max(0,(1-3*n0)/(6*kappa)))
        sig=min(sig4,(1-Rv+lam6)/6-2*q*z6/6)
    Ewt=base-ly/2+d*(Rv+sig-G.zc)
    return min(Eold,Ewt)
def worst(**kw):
    best=9
    for b in np.linspace(0,0.15,11):
        w=max(E(b,dl,x,**kw) for dl in np.linspace(0.05,0.75,15) for x in np.linspace(0,0.5,6))
        best=min(best,w)
    return best
if __name__=="__main__":
    print("n0 at b=0.06:",1-((1-ell+0.06)/2)/(1-(1-ell-0.06)/2+ell))
    for lam in [0.0,0.05,0.08,0.10,0.12,0.15,0.2,0.32,0.433]:
        print(f"lam6={lam:.3f} worst={worst(lam6=lam):+.5f}")
    for A in [2.14,2.13,2.12,2.11]:
        print(f"A={A} (with Ingham) worst={worst(A=A):+.5f}")
