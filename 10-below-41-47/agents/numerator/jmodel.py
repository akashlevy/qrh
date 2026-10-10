"""Joint witness/numerator model.  For each detector parameter t (ours) the adversary picks the
witness split (r, m=t-r).  For that class: count exponent Rs = min(available witness bounds),
numerator exponent sigma(Rs, r, m) = min over available numerator inputs (each via the level-set LP).
High exponent uses max over splits of (Rs + sigma); we minimize over t.
Inputs:
  num4 : paper 4th moment with selected primes (2n0 + 6k z <= 1)                      [PROVED]
  pw   : pointwise bin bound sigma <= delta*n0/2                                        [PROVED]
  NM   : HYPOTHETICAL joint numerator x inverse-witness moment, diag for n0+r+cap*z<=1
  NMamp: HYPOTHETICAL NM extended past diagonal by sixth-power amplification e(A)=(1+5A)/6
  N6   : HYPOTHETICAL sixth moment sum|S|^6 << U^{1+eta6}
"""
import numpy as np
zc=17/50
def wit_bounds(dl,q,L,t,r,mk=0.5,cP=2/9,alpha=5/6,cap82=True):
    m=t-r
    zp=min(L,max(0,(1-2*m)*cP)); b=1-dl*(2*m)-2*q*zp
    if r<1:
        zm=min(L,max(0,(1-r)*mk))
        if cap82: zm=min(zm,max(0,(3-2*r)/8))
        a=1-dl*r-2*q*zm; v=min(a,b)
    else:
        c=1-dl+(alpha-dl)*(r-1); v=min(c,b)
    return max(v,1-dl)
def sigma_class(Rs,dl,q,L,n0,r,m,kappa=0.75,use=('num4','pw'),NMcap=2.0,eta6=0.10,etaP=0.1,etaM=0.1):
    # level-set LP: maximize lam + min(Rs, c_i - k_i lam) over 0<=lam<=lam_max (concave, check breakpoints)
    lines=[(Rs,0.0)]
    zN=min(L,max(0,(1-2*n0)/(6*kappa)))
    if 'num4' in use: lines.append((1-2*q*zN,4.0))
    if 'N6' in use: lines.append((1+eta6,6.0))
    if 'NP3' in use: lines.append((1+etaP-dl*m,4.0))      # sum|S|^4|S_m|^2 << U^{1+etaP}, |S_m|^2>=U^{dl m} on class
    if 'NM3' in use and r<1e9: lines.append((1+etaM-dl*r,4.0))  # sum|S|^4|M_r|^2 << U^{1+etaM}
    if 'NM' in use and r<1:
        A=n0+r
        if A<=1:
            zk=min(L,max(0,(1-A)/NMcap)); lines.append((1-dl*r-2*q*zk,2.0))
        elif 'NMamp' in use:
            lines.append(((1+5*A)/6-dl*r,2.0))
    lmax=dl*n0/2 if 'pw' in use else 0.5
    cands={0.0,lmax}
    for i in range(len(lines)):
        for j in range(i+1,len(lines)):
            (c1,k1),(c2,k2)=lines[i],lines[j]
            if k1!=k2:
                l=(c1-c2)/(k1-k2)
                if 0<=l<=lmax: cands.add(l)
    best=max(l+min(c-k*l for c,k in lines) for l in cands)
    return best-Rs
def E(ell,b,B,dl,x,use=('num4','pw'),ts=np.linspace(1,1.5,26),**kw):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    q=dl*x; a=(1+dl)/2; n0=1-ly/d
    base=a-B+h*(zc-1/6)-ell/2+q*ell
    Rold=9; V=9
    for t in ts:
        wR=-9; wV=-9
        for r in np.linspace(max(0,t-0.5),t,41):
            Rs=wit_bounds(dl,q,L,t,r)
            wR=max(wR,Rs)
            s=sigma_class(Rs,dl,q,L,n0,r,t-r,use=use,**kw)
            wV=max(wV,Rs+s)
        Rold=min(Rold,wR); V=min(V,wV)
    Eold=base-a*ly+d*(Rold+dl/2-zc)
    Ewt=base-ly/2+d*(V-zc)
    return min(Eold,Ewt),Rold,V-Rold
def worst(ell,b,B,grid_d=np.linspace(0.05,0.75,15),grid_x=(0,0.25,0.5),**kw):
    w=-9;arg=None
    for dl in grid_d:
        for x in grid_x:
            e,R,s=E(ell,b,B,dl,x,**kw)
            if e>w: w=e;arg=(round(dl,3),x,round(R,3),round(s,4))
    return w,arg
if __name__=="__main__":
    import sys
    ell=0.1773; B=41/47
    cases=[('paper (num4+pw)',dict(use=('num4','pw'))),
           ('+N6 eta=0.10',dict(use=('num4','pw','N6'),eta6=0.10)),
           ('+NM cap2',dict(use=('num4','pw','NM'),NMcap=2.0)),
           ('+NM cap1',dict(use=('num4','pw','NM'),NMcap=1.0)),
           ('+NM cap2 +amp',dict(use=('num4','pw','NM','NMamp'),NMcap=2.0)),]
    for name,kw in cases:
        res=min((worst(ell,b,B,**kw),b) for b in np.linspace(0.0,0.15,7))
        print(f"{name:18s} worst={res[0][0]:+.5f} at (delta,x,R,sigma)={res[0][1]} b={res[1]:.3f}",flush=True)
