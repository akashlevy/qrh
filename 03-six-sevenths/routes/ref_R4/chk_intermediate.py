import numpy as np
# W(2/7)+refined geometry and vertex A
for name,(lx,ly,ell) in {"W27":(11489/28328,14840/27149,32673/97961),"A":(4/7,4/7,3/7),"paper":(17/48,23/48,1/6)}.items():
    h=1-lx+ell; sigma0=6/7 if name!="paper" else 7/8
    worst=-9
    for delta in np.linspace(0.02,5/6,400):
        a=(1+delta)/2
        for x in np.linspace(0,0.5,26):
            q=x*delta
            base=a*(1-ly)-sigma0-ell/2+q*ell+h*(0.34-1/6)
            # no-slot count: min over t of max over r of min(fA, 1-delta r (r<1), 1-2delta m)
            best=9
            for t in np.linspace(1,1.5,51):
                rs=np.linspace(t-0.5,t,501); m=t-rs
                fA=np.maximum(1,(1+5*rs)/6)-delta*rs
                fI=np.where(rs<1,1-delta*rs,np.inf)
                fP0=1-2*delta*m
                R=np.max(np.minimum(np.minimum(fA,fI),fP0))
                best=min(best,R)
            for d in (0.01,0.5):
                E=base+d*(best+delta/2-0.34)
                worst=max(worst,E)
    print(name,"max E(d) over d in {0.01,1/2} with NO slots:",worst)
