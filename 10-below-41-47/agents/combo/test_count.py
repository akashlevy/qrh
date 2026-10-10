import numpy as np, time
from combo_model import *
def Rclosed(d,x,kappa=0.75,alpha=5/6):
    c=1/(6*kappa); Dx=3-(1+4*c)*x; Px=(2-4*c*x)*(1-x); J=(alpha-d)*Dx+d*Px
    return 1-d+(alpha-d)*d*Px/(2*J)
dl,x=make_bins(0.01,xs=np.linspace(0,0.5,11))
t0=time.time()
for lev,L in [(PAPER08,0.2104),(PROVED,0.2104),(PAPER08,10.0)]:
    R=count_R(dl,x,L,lev,0.75)
    Rc=Rclosed(dl,x)
    err=R-Rc
    i=np.argmax(np.abs(err))
    print(lev.mark_supply, lev.plainLS, L, 'max|R-Rclosed|=%.2e at d=%.2f x=%.2f R=%.5f Rc=%.5f'%(abs(err[i]),dl[i],x[i],R[i],Rc[i]), 'min err',err.min(),'max err',err.max())
print('time',time.time()-t0)
for k in [0.75, 35/47, 0.72]:
    R=count_R(dl,x,0.21,PROVED,k); Rc=Rclosed(dl,x,k)
    print('kappa',k,'max |R-Rc|',np.abs(R-Rc).max())
