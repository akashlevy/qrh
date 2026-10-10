import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize
ALPHA=5/6
dgrid=np.linspace(0.0,0.75,76); xgrid=np.linspace(0,0.5,26)
tg=np.linspace(1,1.5,41)
def Rtable(cM,cP,alpha=ALPHA):
    R=np.zeros((len(dgrid),len(xgrid)))
    for i,d in enumerate(dgrid):
        for j,x in enumerate(xgrid):
            q=x*d; best=9
            for t in tg:
                r=np.linspace(max(t-0.5,0),t,81); m=np.minimum(np.maximum(t-r,0),0.5)
                b1=np.where(r<1,1-d*r-2*q*cM*(1-r),9)
                b2=np.maximum(1,(1+5*r)/6)-d*r
                b3=np.where(m<0.5,1-2*d*m-2*q*cP*(1-2*m),9)
                b4=1-2*d*m
                worst=np.min(np.vstack([b1,b2,b3,b4]),axis=0).max()
                best=min(best,worst)
            R[i,j]=best
    return R
def theta(R,energy_wall=0.0):
    D=dgrid[:,None]; X=xgrid[None,:]
    def f(p):
        lx,ly,ell=p; h=1-lx+ell; b=ly-lx; M=lx+ly
        if min(lx,ly,ell)<0 or ell>ly: return 9
        e=max(0,(1+2*ell-2*M)/3,M+ell-1-energy_wall)
        low=5/6+lx/6-ell/6+0.5*max(0,b/6,2*b-ly)+0.5*e
        hi=((1+D)*(1-ly)/2-ell/2+X*D*ell+h*(R+D/2-1/6)).max()   # includes delta->0 floor (R=1 there)
        pen=10*max(0,1e-3-(ly/2-h*(17/50-1/6)))
        return max(low,hi)+pen
    best=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-9,fatol=1e-11,maxiter=6000)) for s in [(0.355,0.478,0.167),(0.3,0.55,0.15),(0.4,0.4,0.2)]),key=lambda r:r.fun)
    return best.fun,best.x
if __name__=="__main__":
    for cM,cP in [(0.5,2/9),(0.5,1/3),(0.5,1/2),(0.5,1.0),(0.75,2/9),(1.0,2/9),(1.0,1/2),(1.0,1.0)]:
        R=Rtable(cM,cP); th,p=theta(R)
        print(f"c_M={cM:.3f} c_P={cP:.3f}: theta={th:.5f}  (lx,ly,ell)={np.round(p,3)}   R(0.39,1/2)={R[39,25]:.4f}")
