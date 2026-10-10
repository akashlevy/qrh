"""Theta-direction reduction (NOTES sec. 5): low excess = max(w - wT, w - b/6, 0), Gram domination dropped
(Corollary D of 04). Also the same with P78's Gram condition l_y - ell - 11b/6 >= 0 kept, for comparison."""
import numpy as np, sys
from wm4 import scanf
B=41/47
H=dict(gd=np.linspace(0.2,0.75,56),nt=61,nr=121)
ells=np.linspace(0.175,0.215,17); bs=np.linspace(0.0,0.40,33); dMs=(0.0,0.005,0.01,0.02)
mode=sys.argv[1]
if mode=='theta':
    for wT in [0.0,0.02,0.03,0.04,0.05,1.0]:
        print('theta wT',wT, scanf(B,ells=ells,bs=bs,dMs=dMs,hkw=H,mode='theta',wT=wT),flush=True)
if mode=='thetagram':
    for wT in [1.0]:
        print('theta+gram wT',wT, scanf(B,ells=ells,bs=bs,dMs=dMs,hkw=H,mode='theta',wT=wT,gram=True),flush=True)
if mode=='Wbig':   # lossless W with b allowed up to 0.4 (no gram) for comparison
    for wW in [0.04,0.045]:
        print('W(bigb)',wW, scanf(B,ells=ells,bs=bs,dMs=dMs,hkw=H,mode='W',wW=wW),flush=True)
