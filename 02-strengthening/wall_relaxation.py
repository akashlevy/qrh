import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize
dgrid=np.linspace(0.0,0.75,61); xgrid=np.linspace(0,0.5,21)
def Rtable(w,cP=2/9):
    R=np.zeros((len(dgrid),len(xgrid)))
    for i,d in enumerate(dgrid):
        for j,x in enumerate(xgrid):
            q=x*d; best=9
            for t in np.linspace(1,1.5,31):
                r=np.linspace(max(t-0.5,0),t,61); m=np.minimum(np.maximum(t-r,0),0.5)
                zM=np.clip(np.minimum((1+w-r)/2,(3+3*w-2*r)/8),0,None)
                b1=np.where(r<1+w,1-d*r-2*q*zM,9)
                b2=np.maximum(1,(1+5*r/(1+w))/6)-d*r
                b3=np.where(m<0.5,1-2*d*m-2*q*cP*(1-2*m),9); b4=1-2*d*m
                best=min(best,np.min(np.vstack([b1,b2,b3,b4]),axis=0).max())
            R[i,j]=best
    return R
def theta(w):
    R=Rtable(w); D=dgrid[:,None]; X=xgrid[None,:]
    def f(p):
        lx,ly,ell=p; h=1-lx+ell; b=ly-lx; M=lx+ly
        if min(lx,ly,ell)<0: return 9
        e=max(0,(1+2*ell-2*M)/3,M+ell-1-w)
        low=5/6+lx/6-ell/6+0.5*max(0,b/6,2*b-ly)+0.5*e
        hi=((1+D)*(1-ly)/2-ell/2+X*D*ell+h*(R+D/2-1/6)).max()
        pen=10*max(0,1e-3-(ly/2-h*(17/50-1/6)))
        return max(low,hi)+pen
    r=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-9,fatol=1e-11,maxiter=8000)) for s in [(0.355,0.478,0.167),(0.4,0.45,0.25),(0.4,0.5,0.35)]),key=lambda r:r.fun)
    return r.fun,r.x
for w in [0,0.05,0.1,0.2,0.3,0.5]:
    th,p=theta(w); print(f"wall relaxation w={w:.2f}: theta={th:.5f}  (lx,ly,ell)={np.round(p,3)}")
