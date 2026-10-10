"""Parallel scan helper: coarse-to-fine search of min over (ell, b, M<=Mmax) of the worst high exponent."""
import numpy as np, itertools
from multiprocessing import Pool
from wm4 import Mmax, high_worst_fast
B=41/47
COARSE=dict(gd=np.linspace(0.05,0.75,71),nt=41,nr=81)
FINE=dict(gd=np.linspace(0.2,0.75,111),nt=101,nr=201)
def _eval(args):
    ell,b,dM,ekw,hkw,gram=args
    Mm=Mmax(ell,b,B,**ekw)
    if Mm is None: return (9,None)
    M=Mm-dM
    if gram and (M+b)/2-ell-11*b/6<0: return (9,None)
    v,arg=high_worst_fast(ell,M,b,B,**hkw)
    return (v,dict(ell=round(float(ell),4),b=round(float(b),4),M=round(M,4),w=round(M+ell-1,4),arg=arg))
def search(ells,bs,dMs,ekw,hkw=COARSE,gram=False,procs=4):
    jobs=[(e,b,d,ekw,hkw,gram) for e,b,d in itertools.product(ells,bs,dMs)]
    with Pool(procs) as p: res=p.map(_eval,jobs,chunksize=8)
    return min(res,key=lambda z:z[0])
def refine(ekw,gram=False,ells=np.linspace(0.17,0.24,29),bs=np.linspace(0,0.45,46),dMs=(0,0.005,0.01,0.02,0.04)):
    v,a=search(ells,bs,dMs,ekw,COARSE,gram)
    if a is None: return (v,a)
    e0,b0=a['ell'],a['b']
    ells2=np.linspace(e0-0.005,e0+0.005,11); bs2=np.clip(np.linspace(b0-0.02,b0+0.02,9),0,None)
    return search(ells2,bs2,(0,0.0025,0.005,0.01,0.02),ekw,FINE,gram)
