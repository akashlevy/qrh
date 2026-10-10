import numpy as np, amp_model as A
ell=0.1773; B=41/47
ds=np.linspace(0.05,0.75,15); xs=np.linspace(0,0.5,6)
print("b, L, worst(phys supply), worst(inf supply, aux amp 1/2)")
for b in [0.0,0.1,0.2,0.3,0.35,0.4,0.45]:
    lx,ly,h=A.geom(ell,b)
    w1=A.worst(ell,b,B,ds=ds,xs=xs)[0]; w2=A.worst(ell,b,B,ds=ds,xs=xs,Lc=10.0,xcount=0.5)[0]
    print(f"{b:.2f} L={ell/h:.3f} {w1:+.5f} {w2:+.5f}")
print("threshold for prime-only capacity cQ with infinite aux supply at amp 1/2:")
for cQ in [0.6,0.7,0.8,0.85,0.9,0.95]:
    w,arg,b=A.best_b(ell,B,ds=ds,xs=xs,cQ=cQ,Lc=10.0,xcount=0.5)
    print(f"  cQ={cQ}: worst={w:+.5f} at {arg} b={b:.3f}")
print("amplitude cap x<=xmax (slot redesign that lowers worst amplitude; counts at the capped amplitude):")
for xmax in [0.4,0.3,0.25,0.2,0.15,0.12,0.1]:
    w,arg,b=A.best_b(ell,B,ds=ds,xs=np.linspace(0,xmax,6))
    print(f"  xmax={xmax}: worst={w:+.5f} at {arg} b={b:.3f}")
print("ZD(x0): rows with x>=x0 get DH counts:")
for x0 in [0.5,0.4,0.3,0.2,0.15,0.1]:
    w,arg,b=A.best_b(ell,B,ds=ds,xs=np.linspace(0,0.5,21),ZD=x0)
    print(f"  x0={x0}: worst={w:+.5f} at {arg} b={b:.3f}")
