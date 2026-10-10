import sys, numpy as np
from scanpar import refine
which=sys.argv[1:]
cases={
 'paper':dict(mode='paper'),
 'W040':dict(mode='W',wW=0.040),'W043':dict(mode='W',wW=0.043),'W045':dict(mode='W',wW=0.045),
 'slope0.5':dict(mode='slope',f=0.5),'slope0.15':dict(mode='slope',f=0.15),'slope0.1':dict(mode='slope',f=0.1),'slope0':dict(mode='slope',f=0.0),
 'thetaInf':dict(mode='theta',wT=1.0),'theta0':dict(mode='theta',wT=0.0),'theta02':dict(mode='theta',wT=0.02),'theta03':dict(mode='theta',wT=0.03),
 'thW01':dict(mode='thetaW',wT=1.0,wW=0.01),'thW02':dict(mode='thetaW',wT=1.0,wW=0.02),'thW025':dict(mode='thetaW',wT=1.0,wW=0.025),'thW03':dict(mode='thetaW',wT=1.0,wW=0.03),
}
for k in which:
    g = k.endswith('_gram'); key=k.replace('_gram','')
    print(k, refine(cases[key],gram=g), flush=True)
