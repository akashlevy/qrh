import numpy as np
from gmodel import worst
ell=0.1773; B=41/47
def best(**kw):
    return min((worst(ell,b,B,**kw)[0],b) for b in np.linspace(0.0,0.15,7))
for mk in [0.6,0.7,0.8,0.9,1.0]:
    print('mk',mk, ' '.join(f"cP={cp}:{best(mk=mk,cP=cp,cap82=(mk<=0.5))[0]:+.4f}" for cp in [2/9,0.3,0.4,0.5]))
for al in [0.8,0.75,0.7,0.65,0.6]:
    print('alpha',al, f"{best(alpha=al)[0]:+.4f}", 'with mk=.75:', f"{best(alpha=al,mk=0.75,cap82=False)[0]:+.4f}")
