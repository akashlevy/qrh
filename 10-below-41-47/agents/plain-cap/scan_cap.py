"""Scan plain capacity c (cP) and numerator capacity (cN), alone and jointly, and with mk, at ell=0.1773, B=41/47."""
import numpy as np, sys
from gmodel import worst
ell=0.1773; B=41/47
def best(**kw):
    return min((worst(ell,b,B,**kw)[0],b) for b in np.linspace(0.0,0.15,7))
cs=[2/9,0.25,0.3,0.35,0.4,0.5,0.75,1.0]
print("c      cP-only   cN-only   cP=cN=c")
for c in cs:
    r1=best(cP=c)[0]; r2=best(cN=c)[0]; r3=best(cP=c,cN=c)[0]
    print(f"{c:.3f}  {r1:+.5f}  {r2:+.5f}  {r3:+.5f}")
print("with mk (cap82 off for mk>0.5), cP=cN=c")
for mk in [0.5,0.6,0.7,0.75,0.8]:
    print('mk',mk,' '.join(f"c={c:.3f}:{best(mk=mk,cap82=(mk<=0.5),cP=c,cN=c)[0]:+.4f}" for c in [2/9,0.3,0.4,0.5,1.0]))
