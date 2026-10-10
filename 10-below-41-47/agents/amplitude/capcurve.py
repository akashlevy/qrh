import numpy as np, amp_model as A
ell=0.1773; B=41/47
ds=np.linspace(0.05,0.75,15); xs=np.linspace(0,0.5,6); bs=np.linspace(0,0.15,4)
def w(mk,cP): return A.best_b(ell,B,ds=ds,xs=xs,bs=bs,mk=mk,cP=cP,cap82=(mk<=0.5))[0]
print("marked-capacity factor mk (z<=mk(1-r)) vs minimal plain capacity rate cP (z<=cP(1-2m)) closing 41/47")
for mk in [0.5,0.6,0.7,0.8,0.9,1.0]:
    lo,hi=2/9,1.5
    if w(mk,hi)>0: print(f"  mk={mk}: not closable with cP<=1.5"); continue
    for _ in range(10):
        mid=(lo+hi)/2
        if w(mk,mid)<0: hi=mid
        else: lo=mid
    print(f"  mk={mk}: cP >= {hi:.3f}  (paper 2/9=0.222)")
