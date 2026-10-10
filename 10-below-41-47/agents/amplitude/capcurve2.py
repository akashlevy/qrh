import numpy as np, amp_model as A
ell=0.1773; B=41/47
ds=np.linspace(0.05,0.75,15); xs=np.linspace(0,0.5,6); bs=np.linspace(0,0.15,4)
def w(mk,cP,**kw): return A.best_b(ell,B,ds=ds,xs=xs,bs=bs,mk=mk,cP=cP,cap82=(mk<=0.5),**kw)[0]
for lab,kw in [("infinite supply at physical amplitude",dict(Lc=10.0)),
               ("infinite supply at amplitude 1/2 for all rows (aux windows)",dict(Lc=10.0,xcount=0.5))]:
    print(lab)
    for mk in [0.6,0.7,0.8,0.9,1.0]:
        lo,hi=2/9,1.5
        if w(mk,hi,**kw)>0: print(f"  mk={mk}: not closable with cP<=1.5"); continue
        if w(mk,lo,**kw)<0: print(f"  mk={mk}: closes already at cP=2/9 (worst {w(mk,lo,**kw):+.4f})"); continue
        for _ in range(10):
            mid=(lo+hi)/2
            if w(mk,mid,**kw)<0: hi=mid
            else: lo=mid
        print(f"  mk={mk}: cP >= {hi:.3f}")
    print("  mk=1,cP=2/9:",round(w(1.0,2/9,**kw),5)," mk=1.0,cP=1:",round(w(1.0,1.0,**kw),5))
