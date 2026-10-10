# Allowed average numerator exponent sigma on bad rows: sum_C |S_u| <= U^{R+sigma}
import numpy as np
from wmodel import R
zc=17/50
def need(ell,b,kappa=0.75):
    lx=(1-ell-b)/2; ly=(1-ell+b)/2; h=1-lx+ell; B=11/12-ell/4; d=h
    n0=1-ly/d; zP=(1-2*n0)/(6*kappa)
    print(f"ell={ell} b={b:.3f} lx={lx:.4f} ly={ly:.4f} h={h:.4f} n0={n0:.4f} zP={zP:.4f} L=ell/h={ell/h:.4f}")
    for dl in [0.2,0.3,0.35,0.4,0.45,0.47,0.5,0.55,0.6,0.7,0.75]:
        for x in [0,0.25,0.5]:
            q=dl*x; a=(1+dl)/2; Rv=R(dl,x)
            base=a-B+h*(zc-1/6)-ell/2+q*ell
            # E_sigma = base - ly/2 + d(R+sigma-zc) <= 0
            smax=(-(base-ly/2)/d)-Rv+zc
            # weighted Holder sigma: (3R+1)/4 - q zP/2 - R   (in U units) [selected primes inside]
            sw=(1-Rv)/4-q*zP/2
            # old pointwise numerator: Z^{-a ly}U^{delta/2} vs Z^{-ly/2}: sigma_old = delta/2 - (a-1/2)ly/d
            so=dl/2-(a-0.5)*ly/d
            print(f"  delta={dl:.2f} x={x:.2f} R={Rv:.4f} sigma_allowed={smax:+.4f} sigma_wt={sw:+.4f} sigma_old={so:+.4f}")
need(0.1773,0.075)
