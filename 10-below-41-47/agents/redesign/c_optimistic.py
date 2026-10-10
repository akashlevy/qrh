import numpy as np
from rmodel import pieces, des_paper, sigma_wt, E_sigma, worst
B=41/47; ell=0.1773
def annih_noeuler(P):
    e=des_paper(P)
    if e<=0: return e
    return -1e-9 if P['R']<P['ly']/P['d'] else e
for name,f in [('paper',des_paper),('annihilate, no Euler overhead',annih_noeuler)]:
    res=min((worst(f,ell,b,B)+(b,)) for b in np.linspace(0,0.3,31))
    print(f"{name:32s} worst={res[0]:+.5f} at {res[1]} b={res[2]:.3f}")
# minimal R-reduction needed: the worst bin's R vs ly/h as b varies
for b in [0.0,0.06,0.1,0.15,0.2,0.25,0.3]:
    P=pieces(ell,b,B,0.45,0.5); print(f"b={b:.2f}: R(0.45,1/2)={P['R']:.4f} ly/h={P['ly']/P['h']:.4f} E_paper={des_paper(P):+.4f}")
