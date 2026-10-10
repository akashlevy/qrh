import numpy as np, warnings; warnings.filterwarnings("ignore")
from payoff import theta, dgrid, xgrid, tg
def Rtable(cM,cP,ampl):
    """ampl(r): exponent of sum_u |M_u(U^r)|^2 for r>=1 (paper: max(1,(1+5r)/6))"""
    R=np.zeros((len(dgrid),len(xgrid)))
    for i,d in enumerate(dgrid):
        for j,x in enumerate(xgrid):
            q=x*d; best=9
            for t in np.linspace(1,1.5,41):
                r=np.linspace(max(t-0.5,0),t,81); m=np.minimum(np.maximum(t-r,0),0.5)
                b1=np.where(r<1,1-d*r-2*q*cM*(1-r),9)
                b2=ampl(r)-d*r
                b3=np.where(m<0.5,1-2*d*m-2*q*cP*(1-2*m),9)
                b4=1-2*d*m
                best=min(best,np.min(np.vstack([b1,b2,b3,b4]),axis=0).max())
            R[i,j]=best
    return R
cases=[("paper",0.5,2/9,lambda r:np.maximum(1,(1+5*r)/6)),
("alpha=2/3",0.5,2/9,lambda r:np.maximum(1,(1+2*r)/3)),
("alpha=1/2",0.5,2/9,lambda r:np.maximum(1,(1+r)/2)),
("alpha=0 (ideal)",0.5,2/9,lambda r:np.ones_like(r)),
("alpha=2/3 + perfect caps",1,1,lambda r:np.maximum(1,(1+2*r)/3)),
("alpha=1/2 + perfect caps",1,1,lambda r:np.maximum(1,(1+r)/2)),
("alpha=0 + perfect caps",1,1,lambda r:np.ones_like(r)),
]
for n,cM,cP,a in cases:
    th,p=theta(Rtable(cM,cP,a)); print(f"{n:28s} theta={th:.5f} (lx,ly,ell)={np.round(p,3)}")
