import numpy as np, sys
from referee_A2_calib_lib import *
NU = nus(2048)
nuN = np.array([x[2] for x in NU])
bins = [(7,32),(32,128),(128,512),(512,2048)]
for Q in (6000, 60000, 600000, 2400000):
    a,b,N,r,u = rows_for(Q)
    Om = bump(N/Q, .5, 1)
    res=[]
    for (va,vb,vN) in NU:
        if vN >= 2048: continue
        if Q>=600000 and vN>=512: continue
        J = jacobi((va + vb*r) % N, N).astype(float)
        y = Om*J
        S = np.bincount(u, y, 6)
        res.append((vN, (S**2).sum()/(Om**2*J**2).sum()))
    res=np.array(res)
    line=[]
    for lo,hi in bins:
        s=(res[:,0]>=lo)&(res[:,0]<hi)
        if s.sum(): line.append("[%d,%d): %.3f (n=%d)"%(lo,hi,res[s,1].mean(),s.sum()))
    print("Q=%d rows=%d  Rat('1', f=3, sqfree-norm rows):"%(Q,len(N)), "  ".join(line), flush=True)
