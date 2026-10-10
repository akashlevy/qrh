# Re-run of the R5 heuristic order-n model (wf67/r5/quartic_routes_r5.py) with the CORRECTED Kubota reflection:
#   conductor N(h)^3  => reflected columns T3 = 3M + l - 1 (instead of 2M + l - 1);
#   original columns T0 = 1 - l (both sides quartic twists, so take the better of the two).
import sys, numpy as np, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'/private/tmp/claude-501/-Users-akashlevy-Documents-riemann/eab6979c-f095-4107-a475-7e51d88a9fe7/scratchpad/wf67/r5')
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
    starts=[(0.33,0.5,0.17),(0.25,0.35,0.1),(0.3,0.4,0.15),(0.2,0.3,0.05),(0.3,0.45,0.0),(0.35,0.45,0.02),(0.4,0.6,0.0),(0.1,0.2,0.0)]
    r=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-10,fatol=1e-12,maxiter=20000)) for s in starts),key=lambda r:r.fun)
    return r.fun,r.x
def LS_dF(M,T):  return max(0, T-M, 0.75*M+0.5*T-M, 0.5*M+0.75*T-M)       # de Faveri n=4
def LS_BGL(M,T): return max(0, T-M, (2/3)*(M+T)-M)                      # BGL (AB)^{2/3}
def LS_quad(M,T):return max(0, T-M)                                     # quadratic (perfect)
def make(LS, k, use_direct):
    def en(M,l):
        Tk=k*M+l-1; T0=1-l
        e=LS(M,Tk)
        if use_direct: e=min(e,LS(M,T0))
        return max(e,(1+2*l-2*M)/3)
    return en
rows=[("R5 assumption: conductor N(h)^2, deFaveri, reflected only",make(LS_dF,2,False)),
      ("R5 assumption: conductor N(h)^2, BGL, reflected only",make(LS_BGL,2,False)),
      ("CORRECT: conductor N(h)^3, deFaveri, min(direct,reflected)",make(LS_dF,3,True)),
      ("CORRECT: conductor N(h)^3, BGL, min(direct,reflected)",make(LS_BGL,3,True)),
      ("CORRECT: conductor N(h)^3, deFaveri, reflected only",make(LS_dF,3,False)),
      ("direct quartic LS only (no reflection), deFaveri",lambda M,l: max(LS_dF(M,1-l),(1+2*l-2*M)/3)),
      ("hypothetical: N(h)^3 conductor but quadratic dual twist",make(LS_quad,3,False)),
      ("theta route (paper's cubic analogue: N(h)^2, quadratic)",make(LS_quad,2,False))]
for name,en in rows:
    th,p=theta(en); print(f"{name:62s} theta={th:.5f}  (lx,ly,ell)={np.round(p,4)}  M={p[0]+p[1]:.4f}")
print("6/7 =",6/7, " 7/8 =",7/8)
