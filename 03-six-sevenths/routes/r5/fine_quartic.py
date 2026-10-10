# High-resolution evaluation of the prior crude order-n model, n=4, ell=0 (x irrelevant: R decreasing in x).
import numpy as np
from scipy.optimize import minimize
n=4; alpha=(n-1)/n
def R0(d, nt=2001, nr=2001):
    best=9.0
    for t in np.linspace(1,1.5,nt):
        r=np.linspace(t-0.5,t,nr)
        f=np.where(r<=1, 1-d*r, 1+alpha*(r-1)-d*r)
        g=1-2*d*(t-r)
        best=min(best, np.minimum(f,g).max())
    return best
dg=np.linspace(0,0.75,751)
Rv=np.array([R0(d,nt=401,nr=801) for d in dg])
def parts(lx,ly,energy):
    h=1-lx; b=ly-lx; M=lx+ly
    e=energy(M)
    low=1-1/n+lx/n+0.5*max(0,b/n,2*b-ly)+0.5*e
    hi=((1+dg)*(1-ly)/2+h*(Rv+dg/2-1/n)).max()
    pen=10*max(0,1e-3-(ly/2-h*(17/50-1/n)))
    return low,hi,pen
bgl=lambda M: max(0, M-1, M-2/3, (1-2*M)/3)
opt=lambda M: max(0, M-1, 0.75*M-0.5, M-0.75, (1-2*M)/3)   # de Faveri n=4 at ell=0: T=2M-1
for name,en in [("BGL",bgl),("deFaveri",opt)]:
    f=lambda p: max(parts(p[0],p[1],en)[:2])+parts(p[0],p[1],en)[2]
    best=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-12,fatol=1e-14,maxiter=40000)) for s in [(0.31,0.40),(0.30,0.40),(0.32,0.39),(0.29,0.41)]),key=lambda r:r.fun)
    lo,hi,_=parts(best.x[0],best.x[1],en)
    print(name, "theta=%.6f"%best.fun, "lx,ly=",np.round(best.x,5), "low=%.6f hi=%.6f"%(lo,hi), " argmax delta=%.4f"%dg[np.argmax((1+dg)*(1-best.x[1])/2+(1-best.x[0])*(Rv+dg/2-1/n))])
print("6/7=%.6f"%(6/7))
