import json, numpy as np
from vm import best
B=41/47
for w in [0.0,0.02,0.03,0.035,0.04,0.045]:
    ell=(12*(B-5/6)-1+1+w)  # placeholder
    # on line B: M=12(B-5/6)+2ell ; wall: M+ell=1+w  => ell=(1+w-12(B-5/6))/3
    ell=(1+w-12*(B-5/6))/3
    for lab,kw in [('wall only',dict(wall=w)),('wall+marked',dict(wall=w,wM=w,cap82=False))]:
        res=best(ell=ell,B=B,**kw)
        print(json.dumps(dict(w=w,ell=round(ell,4),cfg=lab,worst=round(res[0],5),arg=res[1],b=res[2])),flush=True)
