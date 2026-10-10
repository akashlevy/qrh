import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize
import order_n_engine as E
def Rtable(n,dgrid,xgrid,cM=0.5,cP=2/9,nt=31,nr=61):
    alpha=(n-1)/n; R=np.zeros((len(dgrid),len(xgrid)))
    for i,d in enumerate(dgrid):
        for j,x in enumerate(xgrid):
            q=x*d; best=9
            for t in np.linspace(1,1.5,nt):
                r=np.linspace(max(t-0.5,0),t,nr); m=np.minimum(np.maximum(t-r,0),0.5)
                b1=np.where(r<1,1-d*r-2*q*cM*(1-r),9); b2=np.maximum(1,1+alpha*(r-1))-d*r
                b3=np.where(m<0.5,1-2*d*m-2*q*cP*(1-2*m),9); b4=1-2*d*m
                best=min(best,np.min(np.vstack([b1,b2,b3,b4]),axis=0).max())
            R[i,j]=best
    return R
def solve(n,energy,dgrid,xgrid,R,fix_ell=None):
    D=dgrid[:,None]; X=xgrid[None,:]
    def parts(p):
        if fix_ell is None: lx,ly,ell=p
        else: lx,ly=p; ell=fix_ell
        h=1-lx+ell; b=ly-lx; M=lx+ly
        e=energy(M,ell)
        low=1-1/n+lx/n-ell/n+0.5*max(0,b/n,2*b-ly)+0.5*e
        hi=((1+D)*(1-ly)/2-ell/2+X*D*ell+h*(R+D/2-1/n)).max()
        pen=10*max(0,1e-3-(ly/2-h*(17/50-1/n)))
        return low,hi,pen
    def f(p):
        if min(p)<0: return 9
        low,hi,pen=parts(p); return max(low,hi)+pen
    starts=[(0.33,0.5,0.17),(0.3,0.4,0.0),(0.31,0.4,0.01),(0.25,0.35,0.1)] if fix_ell is None else [(0.31,0.4),(0.3,0.45),(0.35,0.35)]
    r=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-11,fatol=1e-13,maxiter=40000)) for s in starts),key=lambda r:r.fun)
    return r.fun,r.x,parts(r.x)
bgl=lambda M,l: max(0, M+l-1, (2/3)*(3*M+l-1)-M, (1+2*l-2*M)/3)
quad=lambda M,l: max(0,(1+2*l-2*M)/3, M+l-1)
for (nd,nx) in [(61,21),(151,51)]:
    dg=np.linspace(0,0.75,nd); xg=np.linspace(0,0.5,nx); R4=Rtable(4,dg,xg); R6=Rtable(6,dg,xg)
    print(f"grid {nd}x{nx}")
    print("  n=4 Eis+BGL free ell :", solve(4,bgl,dg,xg,R4))
    print("  n=4 Eis+BGL ell=0    :", solve(4,bgl,dg,xg,R4,fix_ell=0.0))
    print("  n=6 theta   ell=0    :", solve(6,quad,dg,xg,R6,fix_ell=0.0))
    print("  n=6 theta free ell   :", solve(6,quad,dg,xg,R6))
