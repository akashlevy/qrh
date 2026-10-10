import numpy as np
from wmodel import R
zc=17/50; kappa=0.75
ell=0.1773; M=1-ell; B=41/47
for b in [0.05,0.075,0.1]:
    lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; n0=1-ly/d; zN=(1-2*n0)/(6*kappa)
    print(f"b={b}: h={h:.4f} n0={n0:.4f}")
    for x in [0.5,0.3]:
      for dl in [0.2,0.3,0.4,0.47,0.55,0.6,0.65]:
        q=dl*x; a=(1+dl)/2
        # weighted: a-B+h(zc-1/6)-ly/2-ell/2+q ell + d((3R+1)/4-zc) - q(2ly-d)/(12k) <=0
        rest=a-B+h*(zc-1/6)-ly/2-ell/2+q*ell+d*(0.25-zc)-q*(2*ly-d)/(12*kappa)
        Rw=-rest/(0.75*d)
        # old: a-B+h(zc-1/6)-a ly-ell/2+q ell+d(R+dl/2-zc)<=0
        rest2=a-B+h*(zc-1/6)-a*ly-ell/2+q*ell+d*(dl/2-zc)
        Ro=-rest2/d
        print(f"  x={x} delta={dl:.2f} R_now={float(R(dl,x)):.4f} need(wt)={Rw:.4f} need(old)={Ro:.4f} DH={1-dl:.3f} Ingham3(1-d)/(3-d)={3*(1-dl)/(3-dl):.4f}")
