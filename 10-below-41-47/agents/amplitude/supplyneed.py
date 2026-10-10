import numpy as np, amp_model as A
ell=0.1773; B=41/47
ds=np.linspace(0.05,0.75,15); xs=np.linspace(0,0.5,6); bs=np.linspace(0,0.15,4)
def w(mk,cP,Lc): return A.best_b(ell,B,ds=ds,xs=xs,bs=bs,mk=mk,cP=cP,cap82=(mk<=0.5),Lc=Lc)[0]
for mk,cP in [(1.0,2/9),(0.9,2/9),(0.8,2/9),(1.0,0.3),(0.75,2/9)]:
    lo,hi=0.2207,1.0
    if w(mk,cP,hi)>0: print(f"mk={mk},cP={cP:.3f}: not closable with L<=1"); continue
    for _ in range(10):
        mid=(lo+hi)/2
        if w(mk,cP,mid)<0: hi=mid
        else: lo=mid
    print(f"mk={mk}, cP={cP:.3f}: counting supply L >= {hi:.4f} (physical 0.2207)")
