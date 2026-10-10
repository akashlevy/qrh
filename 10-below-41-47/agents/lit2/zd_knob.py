"""lit2: plug a family zero-density exponent A(sigma) into the weighted model.
Count exponent R(delta,x) -> min(R_model, A(a)*(1-a)), a=(1+delta)/2 (U-units, DH: A=2 <-> R=1-delta).
Also: (i) required A at each (delta,x) for E<0; (ii) sixth-moment loss threshold lambda6."""
import numpy as np, gmodel as G
ell=0.1773; B=41/47
def Aingham(s): return 3/(2-s)
def Agm(s): return 15/(3+5*s)
def Azeta(s):  # ANTEDB best-known for zeta (Ingham<=0.7, GM 0.7-0.76, then Ivic/Bourgain pieces), crude
    c=[Aingham(s),Agm(s)]
    if s>=0.75: c.append(3/(7*s-4)) if s<=10/13 else c.append(9/(8*s-2))
    if 17/22<=s<=38/49: c.append(2/(9*s-6))
    if s>=41/53: c.append(9/(7*s-1))
    if s>=13/17: c.append(6/(5*s-1))
    return min(c)
def Achen(s): return 7/3                      # Chen 2507.08296 / Gupta-Li, all chars mod q, (qT)-units
def Ahb_quad(s): return 3/(2-s)              # Heath-Brown 1995 quadratic family (Q-aspect) = Ingham form
def AdF6(s):  # de Faveri Cor 1.5 n=6 values quoted in 09 notes: 0.727 at 3/4 -> A=2.908 ; crude constant
    return 2.908
def E_with(ell,b,B,dl,x,Afun=None,**kw):
    M=1-ell; lx=(M-b)/2; ly=(M+b)/2; h=1-lx+ell; d=h; L=ell/h
    q=dl*x; a=(1+dl)/2; Rv=G.count(dl,x,L,**kw)
    if Afun is not None: Rv=min(Rv,max(1-dl,Afun(a)*(1-a)))
    base=a-B+h*(G.zc-1/6)-ell/2+q*ell
    Eold=base-a*ly+d*(Rv+dl/2-G.zc)
    n0=1-ly/d; kappa=0.75; cN=1/(6*kappa)
    zN=min(L,max(0,(1-2*n0)*cN)); sig=(1-Rv)/4-2*q*zN/4
    Ewt=base-ly/2+d*(Rv+sig-G.zc)
    return min(Eold,Ewt),Rv
def worst(Afun,bs=np.linspace(0,0.15,16),dls=np.linspace(0.05,0.75,29),xs=np.linspace(0,0.5,11)):
    best=None
    for b in bs:
        w=-9;arg=None
        for dl in dls:
            for x in xs:
                e,R=E_with(ell,b,B,dl,x,Afun)
                if e>w: w=e;arg=(round(dl,3),round(x,2),round(R,4))
        if best is None or w<best[0]: best=(w,arg,b)
    return best
if __name__=="__main__":
    for name,f in [('none',None),('Ingham',Aingham),('GM',Agm),('zeta-best',Azeta),('Chen 7/3',Achen),('A=2.2',lambda s:2.2),('A=2.15',lambda s:2.15),('A=2.1',lambda s:2.1),('A=2.05',lambda s:2.05),('DH',lambda s:2.0)]:
        w,arg,b=worst(f); print(f"{name:10s} worst={w:+.5f} at (delta,x,R)={arg} b={b:.3f}")
