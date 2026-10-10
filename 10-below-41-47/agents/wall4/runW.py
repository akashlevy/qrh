"""Threshold of a lossless wall relaxation W(wW) and of a slope-f excess in the weighted model (target 41/47)."""
import numpy as np, sys
from wm4 import scanf, high_worst_fast
B=41/47
GDm=np.linspace(0.2,0.75,56)
H=dict(gd=GDm,nt=61,nr=121)
ells=np.linspace(0.175,0.205,13); bs=np.linspace(0.0,0.15,13); dMs=(0.0,0.005,0.01)
mode=sys.argv[1]
if mode=='paper':
    print('paper', scanf(B,ells=ells,bs=bs,dMs=dMs,hkw=H,mode='paper'),flush=True)
if mode=='W':
    for wW in [0.035,0.04,0.045,0.05,0.06]:
        print('W',wW, scanf(B,ells=ells,bs=bs,dMs=dMs,hkw=H,mode='W',wW=wW),flush=True)
if mode=='WM':  # shared wall also feeds lem:marked: r+2z<=1+w/2
    for wW in [0.035,0.04,0.045]:
        Hm=dict(H); Hm['wM']=wW/2
        print('W+marked',wW, scanf(B,ells=ells,bs=bs,dMs=dMs,hkw=Hm,mode='W',wW=wW),flush=True)
if mode=='slope':
    for f in [0.5,0.25,0.15,0.1,0.05,0.0]:
        print('slope',f, scanf(B,ells=np.linspace(0.175,0.26,18),bs=bs,dMs=(0.0,0.01,0.02,0.04),hkw=H,mode='slope',f=f),flush=True)
