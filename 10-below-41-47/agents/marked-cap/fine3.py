import json,sys
from vm import best
cfgs={
 'paper':dict(),
 'mk=1 sat':dict(cM=1,cap82=False),
 'mk .8(1-r),cP .4':dict(zfun=lambda r:0.8*(1-r),cP=0.4,cap82=False),
 'mk .8(1-r),cP .45':dict(zfun=lambda r:0.8*(1-r),cP=0.45,cap82=False),
 'mk sat,cP .4':dict(cM=1,cP=0.4,cap82=False),
 'mk sat,cP .45':dict(cM=1,cP=0.45,cap82=False),
 'mk sat,cP .5':dict(cM=1,cP=0.5,cap82=False),
 'mk .75(1-r),al .75':dict(zfun=lambda r:0.75*(1-r),alpha=0.75,cap82=False),
 'mk sat,al .78':dict(cM=1,alpha=0.78,cap82=False),
 'mk sat,al .76':dict(cM=1,alpha=0.76,cap82=False),
 'mk sat, cP=1':dict(cM=1,cP=1.0,cap82=False),
 'r+2z<=1.13':dict(wM=0.13,cap82=False),
 'r+2z<=1.045':dict(wM=0.045,cap82=False),
}
for k in (sys.argv[1:] or cfgs):
    w,arg,b=best(**cfgs[k]); print(json.dumps(dict(cfg=k,worst=round(w,5),delta=arg[0],x=arg[1],R=arg[2],b=b)),flush=True)
