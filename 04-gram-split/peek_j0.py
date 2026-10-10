import numpy as np
from zomega import *
Q=3000; A=np.load("A_Q%d.npy"%Q); ma,mb,mN=np.load("m_Q%d.npy"%Q)
K=A.shape[0]; Ys=np.array([256*2**i for i in range(K)],float)
def bump(t,lo,hi):
    t=np.asarray(t,float); out=np.zeros_like(t); ins=(t>lo)&(t<hi); u=(t[ins]-lo)/(hi-lo); out[ins]=np.exp(-1/(u*(1-u))+4); return out
Om=bump(mN/Q,.5,1)
for i in range(K):
    v=A[i]*Ys[i]**(1/3); T=np.sum(Om*v); 
    print("Y=%6d  sum_m Om v = %+.4f%+.4fi   |.|*Y^{1/6}=%.3f  |.|*Y^{1/2-1/3}... sum_m Om A * Y^{1/2} = %+.4f"%(Ys[i],T.real,T.imag,abs(T)*Ys[i]**(1/6), (np.sum(Om*A[i])*Ys[i]**0.5).real))
