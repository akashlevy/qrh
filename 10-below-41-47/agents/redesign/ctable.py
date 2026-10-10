"""Direction (c): dimension count per failing bin at ell=0.1773, B=41/47."""
import numpy as np
from rmodel import pieces, sigma_wt, sigma_allowed, E_sigma, des_paper
ell=0.1773; B=41/47
for b in [0.06,0.075]:
    P0=pieces(ell,b,B,0.45,0.5)
    print(f"b={b}: lx={P0['lx']:.4f} ly={P0['ly']:.4f} h={P0['h']:.4f} ly/h={P0['ly']/P0['h']:.4f} n0={P0['n0']:.4f} L={P0['L']:.4f}")
    print(" delta  x    R      sig*    sig_wt  E_paper  thD    R+thD  annih? (ly/d-R)/2  #C(Z) Y(Z)")
    for dl in np.arange(0.2,0.76,0.05):
        for x in (0.25,0.5):
            P=pieces(ell,b,B,dl,x)
            e=des_paper(P)
            if e<=0: continue
            sw=sigma_wt(P); ss=sigma_allowed(P); thD=max(0,(sw-ss)/(dl/2))
            ok = P['R']+thD < P['ly']/P['d']
            print(f" {dl:.2f} {x:.2f} {P['R']:.4f} {ss:+.4f} {sw:+.4f} {e:+.4f} {thD:.4f} {P['R']+thD:.4f} {str(ok):5s} {(P['ly']/P['d']-P['R'])/2:+.4f}  {P['R']*P['h']:.3f} {P['ly']:.3f}")
