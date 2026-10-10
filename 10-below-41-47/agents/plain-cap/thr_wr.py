import sys
from model_wr import best
mk=float(sys.argv[1]); lo,hi=float(sys.argv[2]),float(sys.argv[3])
kw={} if mk==0.5 else dict(mk=mk,cap82=False)
for _ in range(6):
    mid=(lo+hi)/2
    if best(lam=mid,**kw)[0]<0: lo=mid
    else: hi=mid
print(f"mk={mk}: lam* in [{lo:.4f},{hi:.4f}]",flush=True)
