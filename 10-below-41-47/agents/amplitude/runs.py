import numpy as np, amp_model as A
from amp_model import count, best_b, geom
ell=0.1773; B=41/47
lx,ly,h=geom(ell,0.075); L=ell/h
print(f"L={L:.4f}, max marked cap 7/37={7/37:.4f}")
print("1) supply irrelevance at x=1/2: R(d,1/2) with L phys vs L=inf")
for dl in [0.3,0.4,0.47,0.55,0.6]:
    r1,a1=count(dl,0.5,L,detail=True); r2,a2=count(dl,0.5,10.0,detail=True)
    print(f"  d={dl}: R(L)={r1:.5f} R(inf)={r2:.5f} binding (t,r,m)={tuple(round(v,3) for v in a1)}")
print("2) worst E under scenarios (best b):")
sc=[('baseline',{}),
    ('infinite physical supply (counts)',dict(Lc=10.0)),
    ('aux windows: counts use amp 1/2 for ALL rows, inf supply',dict(xcount=0.5,Lc=10.0)),
    ('mean-treatment of marks on high side (3)',dict(mean_marks=True)),
    ('prime-only witness cap 1/2 (large sieve)',dict(cQ=0.5)),
    ('prime-only witness cap 1/2 + aux supply',dict(cQ=0.5,Lc=10.0,xcount=0.5)),
    ('prime-only witness cap 1 + aux supply (HYP)',dict(cQ=1.0,Lc=10.0,xcount=0.5)),
   ]
for name,kw in sc:
    w,arg,b=best_b(ell,B,**kw)
    print(f"  {name:55s} worst={w:+.5f} at (delta,x)={arg} b={b:.3f}")
