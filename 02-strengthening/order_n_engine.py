import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize
dgrid=np.linspace(0.0,0.75,61); xgrid=np.linspace(0,0.5,21)
def Rtable(n,cM=0.5,cP=2/9):
    alpha=(n-1)/n; R=np.zeros((len(dgrid),len(xgrid)))
    for i,d in enumerate(dgrid):
        for j,x in enumerate(xgrid):
            q=x*d; best=9
            for t in np.linspace(1,1.5,31):
                r=np.linspace(max(t-0.5,0),t,61); m=np.minimum(np.maximum(t-r,0),0.5)
                b1=np.where(r<1,1-d*r-2*q*cM*(1-r),9); b2=np.maximum(1,1+alpha*(r-1))-d*r
                b3=np.where(m<0.5,1-2*d*m-2*q*cP*(1-2*m),9); b4=1-2*d*m
                best=min(best,np.min(np.vstack([b1,b2,b3,b4]),axis=0).max())
            R[i,j]=best
    return R
def theta(n,R=None,grh=False):
    D=dgrid[:,None]; X=xgrid[None,:]
    def f(p):
        lx,ly,ell=p; h=1-lx+ell; b=ly-lx; M=lx+ly
        if min(lx,ly,ell)<0: return 9
        e=max(0,(1+2*ell-2*M)/3,M+ell-1)
        low=1-1/n+lx/n-ell/n+0.5*max(0,b/n,2*b-ly)+0.5*e
        if grh: hi=(1-ly)/2-ell/2+h*(1-1/n)
        else: hi=((1+D)*(1-ly)/2-ell/2+X*D*ell+h*(R+D/2-1/n)).max()
        pen=10*max(0,1e-3-(ly/2-h*(17/50-1/n)))
        return max(low,hi)+pen
    r=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-9,fatol=1e-11,maxiter=8000)) for s in [(0.35,0.48,0.17),(0.4,0.4,0.2),(0.3,0.5,0.25),(0.25,0.45,0.3)]),key=lambda r:r.fun)
    return r.fun,r.x
if __name__=='__main__':
  for n in [6,4]:
    th,p=theta(n,Rtable(n)); thg,pg=theta(n,grh=True)
    print(f"n={n}: paper-strength lemmas theta={th:.5f} at {np.round(p,3)} | GRH-for-rows ceiling={thg:.5f} at {np.round(pg,3)}")
