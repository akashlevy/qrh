"""(b) upper envelope: diagonal moments of every real order k<=1/n with prime capacity
   sum_u |S~|^{2k}|Q_z|^2 << U  for  k*n0 + z/cP <= 1, z <= L.   sigma = min_k [(1-R)/(2k) - q z/k].
   A product-split numerator S(n1)S(n2), n1+n2=n0, can at best realize these (granting the
   nonexistent higher product moments).  Also: theta reached by each design on the wall."""
import numpy as np
from rmodel import pieces, E_sigma, sigma_wt, des_paper, des_fractional_ideal, des_sparseLS, worst, best_b
def frac_sigma(P,cP):
    best=9
    for k in np.linspace(1.0,1/P['n0'],200):
        z=min(P['L'],max(0.0,(1-k*P['n0'])*cP))
        best=min(best,(1-P['R'])/(2*k)-P['q']*z/k)
    return best
def des_frac(cP):
    return lambda P: min(P['Eold'],E_sigma(P,frac_sigma(P,cP)))
def theta_reach(design,lo=0.15,hi=0.21,bs=np.linspace(0,0.2,11)):
    # largest ell on the wall with worst<=0 ; theta=11/12-ell/4
    for _ in range(14):
        mid=(lo+hi)/2; B=11/12-mid/4
        w=min(worst(design,mid,b,B)[0] for b in bs)
        if w<=0: lo=mid
        else: hi=mid
    return lo,11/12-lo/4
if __name__=="__main__":
    ell=0.1773;B=41/47
    for name,f in [('frac diag moments cP=2/9',des_frac(2/9)),('frac diag moments cP=1',des_frac(1.0))]:
        w,arg,b=best_b(f,ell,B); print(f"{name:28s} worst@41/47={w:+.5f} {arg} b={b:.3f}")
    for name,f in [('paper',des_paper),('frac diag cP=2/9',des_frac(2/9)),('frac diag cP=1',des_frac(1.0))]:
        l,t=theta_reach(f); print(f"{name:20s} ell_max={l:.4f} theta={t:.6f}  (41/47={41/47:.6f})")
