"""Threshold absolute cross-mark saving delta* closing 41/47 (cubic branch relaxed by 4*delta, unlimited supply),
and the small-delta slope of the worst exponent."""
import numpy as np, json
from xm_model import worst
gd=np.linspace(0.3,0.65,36); bs=np.linspace(0.0,0.10,11)
def w(**kw): return worst(gd=gd,bs=bs,**kw)
lo,hi=0.15,0.20
for _ in range(7):
    mid=(lo+hi)/2
    if w(delta=mid,cubic='relax',unlimited=True)[0]<0: hi=mid
    else: lo=mid
print(json.dumps(dict(delta_star_relax_unlimited=round(hi,4))),flush=True)
lo,hi=0.15,0.30
for _ in range(7):
    mid=(lo+hi)/2
    if w(delta=mid,cubic='relax',unlimited=False)[0]<0: hi=mid
    else: lo=mid
print(json.dumps(dict(delta_star_relax_physical=round(hi,4),check=w(delta=hi,cubic='relax',unlimited=False))),flush=True)
for d in (0.0,0.005,0.01,0.03):
    print(json.dumps(dict(delta=d,worst=w(delta=d,cubic='keep',unlimited=False))),flush=True)
