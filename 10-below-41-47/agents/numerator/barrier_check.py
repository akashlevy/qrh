"""Exponent-level check of the saturation barrier: the flat adversarial configuration
   |S_{chibar_u}(n0)| = U^{sig*} on the witness class C (#C = U^{Rs}), sig* = min((1-Rs)/4 - q zN/2, dl n0/2),
   all other plain polynomials of chibar_u of size U^{o(1)}, witnesses at their lower bounds,
satisfies every available black-box inequality.  Prints the minimum slack over the grid (must be >= 0)."""
import numpy as np
from jmodel import wit_bounds
kap=0.75
def check(ell,b,dls=np.linspace(0.3,0.6,13),xs=np.linspace(0,0.5,6),ts=np.linspace(1,1.5,11)):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; L=ell/h; n0=1-ly/h
    zN=min(L,max(0,(1-2*n0)/(6*kap)))
    worst={}
    def upd(k,v):
        worst[k]=min(worst.get(k,9),v)
    for dl in dls:
        for x in xs:
            q=dl*x
            for t in ts:
                for r in np.linspace(max(0,t-0.5),t,26):
                    m=t-r; Rs=wit_bounds(dl,q,L,t,r)
                    z1=min(L,(1-n0)/(6*kap))
                    # level-set LP over numerator-only inputs: lines Rs, 1-4l-2q zN, 1-2l-2q z1; l<=dl n0/2
                    lm=dl*n0/2; cands=[0,lm,(1-2*q*zN-Rs)/4,(1-2*q*z1-Rs)/2,((1-2*q*zN)-(1-2*q*z1))/2]
                    cands=[c for c in cands if 0<=c<=lm]
                    lam=max(cands,key=lambda l:l+min(Rs,1-4*l-2*q*zN,1-2*l-2*q*z1))
                    cnt=min(Rs,1-4*lam-2*q*zN,1-2*lam-2*q*z1)
                    sig=lam; Rs=cnt   # adversary class: #C=U^cnt, |S|=U^lam
                    # (P1) zero-slot lem:plain, numerator with itself / reflection partner
                    upd('P1 S(n0)S(n0)',1-(Rs+4*sig))
                    upd('P1 S(n0)S(other)',1-(Rs+2*sig))
                    # (P2) positive slots, numerator pair, capacity zN
                    upd('P2 S^2 Q',1-(Rs+4*sig+2*q*zN))
                    # (P2) numerator x short other plain (n'=0) with max capacity
                    upd('P2 S(n0)Q',1-(Rs+2*sig+2*q*z1))
                    # (P1/P2) numerator x plain witness (same presentation allowed), capacity (1-n0-m)/(6k)
                    if m<=0.5:
                        zj=min(L,max(0,(1-n0-m)/(6*kap)))
                        upd('P2 S(n0)S_m Q',1-(Rs+2*sig+dl*m+2*q*zj))
                    # (B) pointwise bin bound
                    upd('B pointwise',dl*n0/2-sig)
                    # (LS) sextic large sieve on numerator x inverse witness and on S^3
                    upd('LS S M_r',max(1,n0+r,2*(1+n0+r)/3)-(Rs+2*sig+dl*r))
                    upd('LS S^3',max(1,3*n0,2*(1+3*n0)/3)-(Rs+6*sig))
    return worst
if __name__=="__main__":
    for b in [0.0,0.05,0.075,0.1]:
        w=check(0.1773,b)
        print(f"b={b}: "+"; ".join(f"{k}: {v:+.4f}" for k,v in w.items()))
