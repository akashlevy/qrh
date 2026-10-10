import numpy as np
from wmodel import R
zc=17/50
ell=0.1773; M=1-ell; b=0.075; kappa=0.75; B=41/47
lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h
print(f"lx={lx:.4f} ly={ly:.4f} h={h:.4f} B={B:.6f}")
for x in [0,0.1,0.2,0.3,0.4,0.5]:
    row=[]
    for dl in [0.1,0.2,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75]:
        q=dl*x; a=(1+dl)/2; Rv=R(dl,x)
        base=a-B+h*(zc-1/6)-ell/2+q*ell
        Eold=base-a*ly+d*(Rv+dl/2-zc)
        Ewt=base-ly/2+d*((3*Rv+1)/4-zc)-q*(2*ly-d)/(12*kappa)
        row.append(min(Eold,Ewt))
    print(f"x={x:.1f} "+" ".join(f"{e:+.4f}" for e in row))
print("deltas: 0.1 0.2 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75")
