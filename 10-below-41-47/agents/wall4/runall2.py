"""Lighter coarse-to-fine scans for the remaining scenarios (output appended to runall.out)."""
import sys, numpy as np
import scanpar
from scanpar import search, COARSE
FINE2=dict(gd=np.linspace(0.05,0.75,71),nt=101,nr=201)
def refine2(ekw,gram=False,ells=np.linspace(0.17,0.24,29),bs=np.linspace(0,0.45,46),dMs=(0,0.005,0.01,0.02,0.04)):
    v,a=search(ells,bs,dMs,ekw,COARSE,gram)
    if a is None: return (v,a)
    e0,b0=a['ell'],a['b']
    ells2=np.linspace(e0-0.003,e0+0.003,7); bs2=np.clip(np.linspace(b0-0.015,b0+0.015,7),0,None)
    return search(ells2,bs2,(0,0.0025,0.005,0.01),ekw,FINE2,gram)
cases={
 'slope0.5':dict(mode='slope',f=0.5),'slope0.25':dict(mode='slope',f=0.25),'slope0.15':dict(mode='slope',f=0.15),'slope0.1':dict(mode='slope',f=0.1),
 'thetaInf':dict(mode='theta',wT=1.0),
 'thW01':dict(mode='thetaW',wT=1.0,wW=0.01),'thW02':dict(mode='thetaW',wT=1.0,wW=0.02),'thW025':dict(mode='thetaW',wT=1.0,wW=0.025),'thW03':dict(mode='thetaW',wT=1.0,wW=0.03),
 'W043':dict(mode='W',wW=0.043),'paper_gram':dict(mode='paper'),
}
for k in sys.argv[1:]:
    g = k.endswith('_gram'); key=k.replace('_gram','') if g and k.replace('_gram','') in cases else k
    print(k, refine2(cases[key],gram=g), flush=True)
