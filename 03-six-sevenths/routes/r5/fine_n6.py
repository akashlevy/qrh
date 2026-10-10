import numpy as np, sys
from scipy.optimize import minimize
def Rtab(n,dg,xg,cM=0.5,cP=2/9,nt=201,nr=401):
    alpha=(n-1)/n; R=np.full((len(dg),len(xg)),9.0)
    ts=np.linspace(1,1.5,nt)
    for i,d in enumerate(dg):
        q=xg[:,None]*d            # (nx,1)
        for t in ts:
            r=np.linspace(max(t-0.5,0),t,nr)[None,:]
            m=np.minimum(np.maximum(t-r,0),0.5)
            b1=np.where(r<1,1-d*r-2*q*cM*(1-r),9); b2=np.maximum(1,1+alpha*(r-1))-d*r
            b3=np.where(m<0.5,1-2*d*m-2*q*cP*(1-2*m),9); b4=1-2*d*m
            v=np.minimum(np.minimum(b1,b2),np.minimum(b3,b4)).max(axis=1)
            R[i]=np.minimum(R[i],v)
    return R
def solve(n,energy,dg,xg,R,starts):
    D=dg[:,None]; X=xg[None,:]
    def f(p):
        lx,ly,ell=p
        if min(p)<0: return 9
        h=1-lx+ell; b=ly-lx; M=lx+ly; e=energy(M,ell)
        low=1-1/n+lx/n-ell/n+0.5*max(0,b/n,2*b-ly)+0.5*e
        hi=((1+D)*(1-ly)/2-ell/2+X*D*ell+h*(R+D/2-1/n)).max()
        pen=10*max(0,1e-3-(ly/2-h*(17/50-1/n)))
        return max(low,hi)+pen
    r=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-11,fatol=1e-13,maxiter=40000)) for s in starts),key=lambda r:r.fun)
    return r.fun,r.x
quad=lambda M,l: max(0,(1+2*l-2*M)/3, M+l-1)
quad_c=lambda M,l: max(0,(1+3*l-2*M)/4,(1+2*l-2*M)/3, M+l-1)
bgl4=lambda M,l: max(0, M+l-1, (2/3)*(3*M+l-1)-M, (1+2*l-2*M)/3)
opt4=lambda M,l: max(0, M+l-1, 0.75*M+0.5*(2*M+l-1)-M, 0.5*M+0.75*(2*M+l-1)-M, (1+2*l-2*M)/3)
dg=np.linspace(0,0.75,151); xg=np.linspace(0,0.5,26)
R6=Rtab(6,dg,xg); R4=Rtab(4,dg,xg)
print("n=6 theta-route  (1+2l-2M)/3      :", solve(6,quad,dg,xg,R6,[(0.356,0.478,0.166),(0.35,0.48,0.17)]))
print("n=6 theta-route  incl (1+3l-2M)/4 :", solve(6,quad_c,dg,xg,R6,[(0.356,0.478,0.166),(0.35,0.48,0.17)]))
print("n=4 Eis + BGL   (free ell)        :", solve(4,bgl4,dg,xg,R4,[(0.3,0.4,0.0),(0.28,0.42,0.01),(0.33,0.5,0.17)]))
print("n=4 Eis + deFaveri (free ell)     :", solve(4,opt4,dg,xg,R4,[(0.3,0.4,0.0),(0.28,0.42,0.01),(0.33,0.5,0.17)]))
print("n=4 theta (Patterson conj)        :", solve(4,quad,dg,xg,R4,[(0.33,0.5,0.17),(0.3,0.45,0.1)]))
