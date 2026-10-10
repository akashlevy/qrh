import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize
from order_n_engine import Rtable, dgrid, xgrid
R4=Rtable(4)
def theta(energy):
    D=dgrid[:,None]; X=xgrid[None,:]; n=4
    def f(p):
        lx,ly,ell=p; h=1-lx+ell; b=ly-lx; M=lx+ly
        if min(lx,ly,ell)<0: return 9
        e=energy(M,ell)
        low=1-1/n+lx/n-ell/n+0.5*max(0,b/n,2*b-ly)+0.5*e
        hi=((1+D)*(1-ly)/2-ell/2+X*D*ell+h*(R4+D/2-1/n)).max()
        pen=10*max(0,1e-3-(ly/2-h*(17/50-1/n)))
        return max(low,hi)+pen
    r=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-9,fatol=1e-11,maxiter=8000)) for s in [(0.33,0.5,0.17),(0.25,0.35,0.1),(0.3,0.4,0.15),(0.2,0.3,0.05)]),key=lambda r:r.fun)
    return r.fun,r.x
quad_dual=lambda M,l: max(0,(1+2*l-2*M)/3, M+l-1)                     # theta route (quadratic dual twist), as in paper
def quartic_dual(M,l):   # Eisenstein route: dual twist quartic; conjectured-optimal quartic LS  A+B+A^{3/4}B^{1/2}+A^{1/2}B^{3/4}
    T=2*M+l-1                                                       # marked dual length (all marks active)
    return max(0, T-M, 0.75*M+0.5*T-M, 0.5*M+0.75*T-M, (1+2*l-2*M)/3)
def quartic_dual_known(M,l):  # known quartic LS (Gao-Zhao-type): A+B+(AB)^{3/4}
    T=2*M+l-1
    return max(0, T-M, 0.75*(M+T)-M, (1+2*l-2*M)/3)
for name,en in [("quartic, theta route (needs Patterson quartic conj.)",quad_dual),
                ("quartic, Eisenstein route + conjectural optimal quartic LS",quartic_dual),
                ("quartic, Eisenstein route + known quartic LS (AB)^{3/4}",quartic_dual_known)]:
    th,p=theta(en); print(f"{name:60s} theta={th:.5f}  (lx,ly,ell)={np.round(p,3)}")
