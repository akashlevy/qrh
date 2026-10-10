"""Amplitude-agent extension of gmodel (weighted-numerator high side, count LP).
Adds knobs for: supply (physical/auxiliary), auxiliary amplitude used only in counts,
amplitude caps (slot redesign), 'mean' (moment) treatment of marks on the high side,
zero-dominance hypothesis ZD(x0), prime-only witnesses, per-class contour choice."""
import numpy as np
zc=17/50
TS=np.linspace(1,1.5,201)
def count(dl,x,L,mk=0.5,cP=2/9,alpha=5/6,cap82=True,DH=False,cQ=0.0,ts=TS,nr=161,detail=False):
    """Row-count exponent R(delta,x) as min over t of max over adversarial split (r,m=t-r).
    L: slot supply (U-units) usable at amplitude q=delta*x.
    cQ: capacity of an extra prime-only witness (sum_u |Q_z|^2 << U for z<=cQ); 0 = off."""
    if DH: return (1-dl,None) if detail else 1-dl
    q=dl*x
    extra = 1-2*q*min(L,cQ) if cQ>0 else 9.0
    T=ts[:,None]; s=np.linspace(0,1,nr)[None,:]
    r=np.maximum(0,T-0.5)+(T-np.maximum(0,T-0.5))*s; m=T-r
    zp=np.minimum(L,np.maximum(0,(1-2*m)*cP)); b=1-dl*2*m-2*q*zp
    zm=np.minimum(L,np.maximum(0,(1-r)*mk))
    if cap82: zm=np.minimum(zm,np.maximum(0,(3-2*r)/8))
    a=np.where(r<1,1-dl*r-2*q*zm,1-dl+(alpha-dl)*(r-1))
    v=np.minimum(np.minimum(a,b),extra)
    j=np.argmax(v,axis=1); worst=v[np.arange(len(ts)),j]
    i=int(np.argmin(worst)); best=max(worst[i],1-dl)
    if detail: return best,(float(ts[i]),float(r[i,j[i]]),float(m[i,j[i]]))
    return best

def geom(ell,b,M=None):
    if M is None: M=1-ell
    lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell
    return lx,ly,h

def E(ell,b,B,dl,x,kappa=0.75,xcount=None,Lc=None,mean_marks=False,ZD=None,**kw):
    """High-side exponent min(E_old,E_wt) at d=h.
    xcount: amplitude usable by the counts (>= physical x if auxiliary windows); default x.
    Lc: supply for counts (default physical L=ell/h).
    mean_marks: bound the high side by a moment in the marks instead of pointwise (see NOTES A2).
    ZD: x0 such that rows with x>=x0 have DH counts (hypothetical zero-dominance rigidity)."""
    lx,ly,h=geom(ell,b); d=h; L=ell/h
    q=dl*x; a=(1+dl)/2
    xc = x if xcount is None else max(x,xcount)
    Rv=count(dl,xc,L if Lc is None else Lc,**kw)
    if ZD is not None and x>=ZD: Rv=1-dl
    base=a-B+h*(zc-1/6)-ell/2+q*ell
    Eold=base-ly*a+d*(Rv+dl/2-zc)
    n0=1-ly/d; zN=min(L,max(0,(1-2*n0)/(6*kappa)))
    Ewt=base-ly/2+d*(Rv+(1-Rv)/4-q*zN/2-zc)
    if mean_marks:
        # replace pointwise |Q_all|<=U^{qL} on C by Cauchy-Schwarz with sum_u|Q_z|^2<<U (z<=3/8, marked r=0):
        # sum_C |N||Q_all| <= U^{q(L-z)} (sum_C|N|^2)^{1/2} (sum_u |Q_z|^2)^{1/2}; numerator 2nd moment via 4th.
        z=min(L,3/8)
        sig_mean = q*(L-z) + (Rv/2+1/2)/2 + 1/2 - Rv   # U-exponent of sum_C|N||Q_all| minus R, minus qL handled below
        # E_wt uses d*(R+sigma) with sigma=(1-R)/4-q zN/2 and +q*ell (=d*qL) in base.
        # mean version: total U-exponent = q(L-z) + R/4+1/4 + 1/2 ; subtract d*qL already in base:
        Ewt_mean = base - q*ell - ly/2 + d*(q*(L-z) + Rv/4+1/4+1/2 - zc)
        Ewt=min(Ewt,Ewt_mean)
    return min(Eold,Ewt)

def worst(ell,b,B,ds=np.linspace(0.05,0.75,29),xs=np.linspace(0,0.5,11),**kw):
    w=-9; arg=None
    for dl in ds:
        for x in xs:
            e=E(ell,b,B,dl,x,**kw)
            if e>w: w=e; arg=(round(dl,3),round(x,3))
    return w,arg

def best_b(ell,B,bs=np.linspace(0,0.15,7),**kw):
    return min((worst(ell,b,B,**kw)+(b,)) for b in bs)

if __name__=="__main__":
    ell=0.1773; B=41/47
    print("baseline:",best_b(ell,B))
