import numpy as np, amp_model as A
ell=0.1773; B=41/47
ds=np.linspace(0.05,0.75,29); bs=np.linspace(0,0.15,7)
def bisect(f,lo,hi,it=12):
    # f(lo)<0 (closes), f(hi)>0
    for _ in range(it):
        mid=(lo+hi)/2
        if f(mid)<0: lo=mid
        else: hi=mid
    return lo,hi
cap=lambda xm: A.best_b(ell,B,ds=ds,bs=bs,xs=np.linspace(0,xm,11))[0]
print("amplitude cap threshold x_max:",bisect(cap,0.2,0.3,8))
zd=lambda x0: A.best_b(ell,B,ds=ds,bs=bs,xs=np.unique(np.r_[np.linspace(0,0.5,26),x0-1e-9]),ZD=x0)[0]
print("ZD threshold x0:",bisect(zd,0.2,0.3,8))
cq=lambda c: A.best_b(ell,B,ds=ds,bs=bs,xs=np.linspace(0,0.5,11),cQ=c,Lc=10.0,xcount=0.5)[0]
print("prime-only capacity threshold (with aux supply at amp 1/2):",bisect(lambda c:-cq(c),0.85,0.95,8))
