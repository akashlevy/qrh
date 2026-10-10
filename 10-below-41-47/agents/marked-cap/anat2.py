import numpy as np
def anat(dl,x,L,cM=2.0,wM=0.0,cP=2/9,alpha=5/6,cap82=True,zfun=None,nt=501,nr=1001):
    q=dl*x
    t=np.linspace(1,1.5,nt)[:,None]; s=np.linspace(0,1,nr)[None,:]
    r=np.maximum(0,t-0.5)+(t-np.maximum(0,t-0.5))*s; m=t-r
    zp=np.clip((1-2*m)*cP,0,L); b=1-dl*2*m-2*q*zp
    zm=np.clip(zfun(r) if zfun is not None else (1+wM-r)/cM,0,L)
    if cap82: zm=np.minimum(zm,np.maximum(0,(3-2*r)/8))
    a=1-dl*r-2*q*zm; c=1-dl+(alpha-dl)*(r-1)
    w=np.where(r<1,a,c); v=np.minimum(w,b)
    j=np.argmax(v,axis=1); vals=v[np.arange(nt),j]; i=int(np.argmin(vals))
    rr=r[i,j[i]]; 
    # also report the inverse-marked witness options at that t: best r among r<1 branch
    return dict(R=round(float(vals[i]),4),t=round(float(t[i,0]),4),r_adv=round(float(rr),4),m=round(float(t[i,0]-rr),4),
                zM=round(float(zm[i,j[i]]),4),zP=round(float(zp[i,j[i]]),4),r_plus_2z=round(float(rr+2*zm[i,j[i]]),4))
ell=0.1773;M=1-ell
for b,cfg,kw in [(0.055,'paper',{}),(0.04,'mk sat cP .45',dict(cM=1,cP=0.45,cap82=False)),(0.005,'mk .75 al .75',dict(zfun=lambda r:0.75*(1-r),alpha=0.75,cap82=False)),(0.0,'mk sat',dict(cM=1,cap82=False))]:
    lx=(M-b)/2;h=1-lx+ell;L=ell/h
    for dl in [0.4,0.46,0.5,0.54]:
        print(cfg,'L=%.4f'%L,dl,anat(dl,0.5,L,**kw))
