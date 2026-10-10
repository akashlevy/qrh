import numpy as np, warnings; warnings.filterwarnings("ignore")
from functools import lru_cache
ALPHA=5/6
def count_R(delta, x, cM=0.5, cP=2/9, alpha=ALPHA, nt=61, nr=121):
    """R_* = min_t max_{witness (r,m)} min(available count bounds)  (Prop detector-counts, generalized)"""
    q=x*delta; best=9
    for t in np.linspace(1,1.5,nt):
        worst=-9
        for r in np.linspace(max(t-0.5,0),t,nr):
            m=min(max(t-r,0),0.5)        # worst plain length (bounds decrease in m)
            b=[]
            if r<1: b.append(1-delta*r-2*q*cM*(1-r))
            b.append(max(1,(1+5*r)/6)-delta*r)          # amplified inverse, no slots
            if m<0.5: b.append(1-2*delta*m-2*q*cP*(1-2*m))
            b.append(1-2*delta*m)                       # plain, no slots
            worst=max(worst,min(b))
        best=min(best,worst)
    return best
def closed(delta,x,cP=2/9):
    D=3-x-4*cP*x; P=(1-x)*(2-4*cP*x); J=(ALPHA-delta)*D+delta*P
    return 1-delta+(ALPHA-delta)*delta*P/(2*J)
for d,x in [(0.39,0.5),(0.2,0.3),(0.6,0.5),(0.39,0.0)]:
    print(f"delta={d},x={x}: generalized={count_R(d,x):.5f}  paper closed form={closed(d,x):.5f}")
