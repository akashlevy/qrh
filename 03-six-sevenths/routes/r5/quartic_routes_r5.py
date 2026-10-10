# R5 re-run of the prior crude order-n model (strengthening/quartic_routes.py), adding the BGL (AB)^{2/3}
# quartic large sieve (BGL Thm 1.3 holds for every n>3, incl. n=4) and the dFDH Conj 1.3 shape.
import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize
from order_n_engine import Rtable, dgrid, xgrid
R4=Rtable(4)
def theta(energy, n=4, R=R4):
    D=dgrid[:,None]; X=xgrid[None,:]
    def f(p):
        lx,ly,ell=p; h=1-lx+ell; b=ly-lx; M=lx+ly
        if min(lx,ly,ell)<0: return 9
        e=energy(M,ell)
        low=1-1/n+lx/n-ell/n+0.5*max(0,b/n,2*b-ly)+0.5*e
        hi=((1+D)*(1-ly)/2-ell/2+X*D*ell+h*(R+D/2-1/n)).max()
        pen=10*max(0,1e-3-(ly/2-h*(17/50-1/n)))
        return max(low,hi)+pen
    starts=[(0.33,0.5,0.17),(0.25,0.35,0.1),(0.3,0.4,0.15),(0.2,0.3,0.05),(0.3,0.45,0.0),(0.35,0.45,0.02)]
    r=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-10,fatol=1e-12,maxiter=20000)) for s in starts),key=lambda r:r.fun)
    return r.fun,r.x
quad_dual=lambda M,l: max(0,(1+2*l-2*M)/3, M+l-1)
def quartic_conj(M,l):   # dFDH Conj 1.3 / (1.8) n=4: A+B+A^{3/4}B^{1/2}+A^{1/2}B^{3/4}
    T=2*M+l-1
    return max(0, T-M, 0.75*M+0.5*T-M, 0.5*M+0.75*T-M, (1+2*l-2*M)/3)
def quartic_bgl(M,l):    # BGL Thm 1.3, n=4: A+B+(AB)^{2/3}
    T=2*M+l-1
    return max(0, T-M, (2/3)*(M+T)-M, (1+2*l-2*M)/3)
def quartic_34(M,l):     # prior script's "known" (AB)^{3/4}
    T=2*M+l-1
    return max(0, T-M, 0.75*(M+T)-M, (1+2*l-2*M)/3)
for name,en in [("theta route (Patterson quartic conj.)",quad_dual),
                ("Eisenstein + dFDH Conj 1.3",quartic_conj),
                ("Eisenstein + BGL (AB)^{2/3} [proved]",quartic_bgl),
                ("Eisenstein + (AB)^{3/4} [prior script]",quartic_34)]:
    th,p=theta(en); print(f"{name:42s} theta={th:.5f}  (lx,ly,ell)={np.round(p,4)}")
print("6/7 =",6/7)
