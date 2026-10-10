# Sensitivity of the best achievable sigma0 to each ingredient of the 7/8 argument.
# Needs numpy + scipy.  Each case re-optimizes the geometry (l_x, l_y, ell).
import numpy as np
from scipy.optimize import differential_evolution
Z0=17/50; D0=1/50
ds=np.linspace(D0,5/6,300)[:,None]; xs=np.linspace(0,0.5,101)[None,:]
def model(p, cP=2/9, alpha=5/6, gram=1/6, energy="paper", count="paper", dmax=5/6, floor=True):
    lx,ly,ell=p; h=1-lx+ell; b=ly-lx; M=lx+ly
    if energy=="paper": loss=max(0,(1+2*ell-2*M)/3,M+ell-1)
    else: loss=0.0
    s_low=5/6+lx/6-ell/6+0.5*max(0,gram*b,2*b-ly)+0.5*loss
    d=ds[ds[:,0]<=dmax+1e-12]; x=xs
    if count=="paper":
        D=3-x-4*cP*x; P=(1-x)*(2-4*cP*x); J=(alpha-d)*D+d*P
        R=1-d+(alpha-d)*d*P/(2*J)
        t=1+d*P/(2*J); rs=((2-4*cP*x)*t-(1-2*cP)*x)/D
        cap=np.maximum((1-rs)/2, cP*(1-2*(t-rs))).max()
    else:  # ideal density-hypothesis count
        R=1-d+0*x; cap=0
    s_high=((1+d)*(1-ly)/2-ell/2+x*d*ell+h*(R+d/2-1/6)).max()
    s_floor=(1+D0)*(1-ly)/2-ell/2+D0*ell/2+h*(1+D0/2-1/6) if floor else -9
    dd=d[:,0]; s_mid=((1+dd)*(1-ly)/2-ell/2+dd*ell/2+h*(Z0-1/6)+0.5*(1-2*dd/3+dd/2-Z0)).max() if count=="paper" else -9
    cons=[ell/h-cap, ly/2-h*(Z0-1/6), ly-ell, M-2*ell]
    return max(s_low,s_high,s_floor,s_mid), cons, (s_low,s_high,s_floor,s_mid)
def best(**kw):
    def obj(p):
        s,c,_=model(p,**kw); return s+10*sum(max(0,2e-3-v) for v in c)
    r=differential_evolution(obj,bounds=[(0.0,1.2),(0.05,1.2),(0,0.8)],seed=5,tol=1e-12,maxiter=3000,popsize=40,polish=True)
    s,c,parts=model(r.x,**kw)
    return s,r.x,parts
cases=[("paper lemmas",{}),
 ("c_P=1/(6*0.70)",dict(cP=1/(6*0.70))),
 ("c_P=1/(6*0.5)",dict(cP=1/3)),
 ("alpha=3/4",dict(alpha=3/4)),
 ("alpha=2/3",dict(alpha=2/3)),
 ("Gram exponent 0",dict(gram=0)),
 ("ideal row count R=1-delta",dict(count="ideal")),
 ("GRH for rows (floor only)",dict(count="ideal",dmax=D0)),
 ("ideal energy (no M+l<=1)",dict(energy="ideal")),
 ("ideal energy + ideal count",dict(energy="ideal",count="ideal")),
]
for name,kw in cases:
    s,p,parts=best(**kw)
    print(f"{name:32s} sigma0={s:.5f}  lx={p[0]:.4f} ly={p[1]:.4f} ell={p[2]:.4f}  parts={[round(float(v),4) for v in parts]}")
