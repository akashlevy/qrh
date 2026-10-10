"""Self-consistent bootstrap inside the Part II bookkeeping (float model, corrected per the audit).
Under H(beta): (i) bin ceiling delta <= 2beta-1, (ii) Lemma plain with kappa = 2beta-1 (capacity (1-2m)/(6 kappa)),
allowing kappa < 3/4 HYPOTHETICALLY (the lemma is only proved for kappa in [3/4,1]).
theta(beta) = min over geometry of max(low, high, floor).  A bootstrap closes at beta iff theta(beta) < beta.
Also: idealised row counts  (a) R = 1 - delta (density hypothesis for every row family),
                            (b) no intermediate rows at all (GRH for rows), floor delta0 in {1/50, 0}.
"""
import numpy as np, sys, json
from scipy.optimize import minimize
ZETA=1e-4; Z0=17/50
def Rcount(d,x,S,CP,nt=31,nr=91):
    q=x*d; best=9
    for t in np.linspace(1,1.5,nt):
        r=np.linspace(max(t-0.5,0),t,nr); m=np.minimum(np.maximum(t-r,0),0.5)
        zM=np.clip(np.minimum.reduce([(1-r)/2,1-r,(3-2*r)/8,np.full_like(r,S)]),0,None)
        bI=np.where(r<1,1-d*r-2*q*zM,9)
        bA=np.maximum(1,(1+5*r)/6)-d*r
        zP=np.minimum(CP*(1-2*m),S)
        bP=np.where(m<0.5,1-2*d*m-2*q*zP,9); bP0=1-2*d*m
        best=min(best,np.min(np.vstack([bI,bA,bP,bP0]),axis=0).max())
    return best
def elow(lx,ly,ell):
    M=lx+ly; worst=-9
    for dd in np.linspace(0,ell,21):
        Mp,lp=M-2*dd,ell-dd
        e=max(0,(1+3*lp-2*Mp)/4 if 2*Mp+lp>=1 else 0,Mp+lp-1)
        worst=max(worst,-dd+e/2)
    return worst
def parts(p,mode,beta,dg,xg,d0=1/50):
    lx,ly,ell=p; h=1-lx+ell; b=ly-lx; S=ell/(h+ZETA)
    kap=max(2*beta-1,1e-3); CP=1/(6*kap)
    low=5/6+lx/6-ell/6+0.5*max(0,b/6,2*b-ly)+elow(lx,ly,ell)
    dgs=[d for d in dg if d<=2*beta-1+1e-12]
    if mode=='paper':
        hi=max((1+d)*(1-ly)/2-ell/2+x*d*ell+(h+ZETA)*(Rcount(d,x,S,CP)+d/2-1/6) for d in dgs for x in xg) if dgs else -9
    elif mode=='DH':
        hi=max((1+d)*(1-ly)/2-ell/2+x*d*ell+(h+ZETA)*((1-d)+d/2-1/6) for d in dgs for x in xg) if dgs else -9
    else: hi=-9   # GRH: no intermediate bins
    fl=(1+d0)*(1-ly)/2-ell/2+(d0/2)*ell+(h+ZETA)*(1+d0/2-1/6)
    small=h*(Z0-1/6)-ly/2+0.0126
    gram=0 if b<=0 else max(0,-(ly-ell-11*b/6))
    return low,hi,fl,small,gram
def theta(p,mode,beta,dg,xg,d0):
    if min(p)<0 or p[2]>p[0]: return 9
    low,hi,fl,small,gram=parts(p,mode,beta,dg,xg,d0)
    return max(low,hi,fl)+10*max(0,small)+10*gram
def best(mode,beta,d0=1/50):
    dg=np.linspace(1/50,5/6,28); xg=np.linspace(0,.5,7)
    starts=[(17/48,23/48,1/6),(0.40,0.40,0.20),(0.38,0.45,0.19),(0.36,0.48,0.17)]
    r=min((minimize(theta,s,args=(mode,beta,dg,xg,d0),method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-8,maxiter=500)) for s in starts),key=lambda r:r.fun)
    return r.fun, r.x
if __name__=="__main__":
    out={}
    for mode,betas,d0 in [('paper',[0.875,0.87,0.8667,0.86],1/50),('DH',[0.875,0.87,0.86],1/50),('GRH',[0.87],1/50),('GRH',[0.87],0.0)]:
        for be in betas:
            v,x=best(mode,be,d0); out[f"{mode} beta={be} d0={d0}"]=(round(v,6),[round(t,5) for t in x])
            print(mode,'beta',be,'d0',d0,'theta',round(v,6),np.round(x,5),flush=True)
