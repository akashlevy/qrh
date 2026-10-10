import numpy as np, amp_model as A
ell=0.1773; B=41/47
ds=np.linspace(0.05,0.75,29); xs=np.linspace(0,0.5,11); bs=np.linspace(0,0.15,7)
def w(mk,**kw): return A.best_b(ell,B,ds=ds,xs=xs,bs=bs,mk=mk,**kw)[0]
for lab,kw in [("unlimited supply (slot powers), cubic branch dropped",dict(Lc=10.0,cap82=False)),
               ("unlimited supply (slot powers), cubic branch 2r+8z<=3 kept",dict(Lc=10.0,cap82=True)),
               ("physical supply L, cubic branch dropped",dict(cap82=False))]:
    lo,hi=0.5,1.0
    if w(hi,**kw)>0: print(lab,": mk=1 not enough, worst",round(w(hi,**kw),5)); continue
    for _ in range(9):
        mid=(lo+hi)/2
        if w(mid,**kw)<0: hi=mid
        else: lo=mid
    print(lab,f": mk* = {hi:.4f}")
