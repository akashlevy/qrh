# Library: coefficients of Kubota series psi(nu, s, xi; class) for nu = lambda^j * i^e * pi0^k (pi0 primary prime, k=0..3)
import numpy as np, math, cmath, sys
sys.path.insert(0,'.')
from gi import *
dt=np.dtype([('u','<i8'),('v','<i8'),('N','<i8'),('re','<f8'),('im','<f8'),('kl','<i4'),('ki','<i4'),('om','<i4')])
A=np.fromfile('sqf2_2e6.bin',dtype=dt)
U=A['u']; Vv=A['v']; NN=A['N'].astype(np.float64); GT=A['re']+1j*A['im']; KL=A['kl'].astype(np.int64); KI=A['ki'].astype(np.int64)
OM=A['om'].astype(np.int64)
NMAX=2e6
ANG=np.angle(U+1j*Vv)
CLS1=((U-1)%4==0)&(Vv%4==0)     # c = 1 mod 4 (else 1+lambda^3 mod 4)
def qtab(pi0):
    a,b=pi0; q=a*a+b*b
    r=(-a*pow(b%q,q-2,q))%q
    tab=np.full(q,-1,dtype=np.int64)
    for x in range(1,q):
        w=pow(x,(q-1)//4,q); tab[x]={1:0,r:1,q-1:2,(q-r)%q:3}[w]
    return q,r,tab
def sym_prime(x,pi0):
    """(x/pi0)_4 exponent for Gaussian integer x (tuple)"""
    q,r,tab=qtab(pi0); z=(x[0]+x[1]*r)%q
    return int(tab[z])
def Gj(pi0,j):
    a,b=pi0; q,r,tab=qtab(pi0)
    return sum((1j)**((j*int(tab[d]))%4)*cmath.exp(4j*math.pi*d*a/q) for d in range(1,q))
def coeffs(pi0,k,j=0,e=0,omega=0,cls=None):
    """normalized coefficients a(c)=xi(c) g~(nu,c), nu=lambda^j i^e pi0^k, over primary c (optionally class mod 4).
       returns (N array, coef array, class1 bool array)"""
    q,r,tab=qtab(pi0); a,b=pi0
    red=(U%q+(Vv%q)*r)%q; cop=red!=0
    k4=np.where(cop,tab[red],0)
    Cexp=(((q-1)//4)%2)*(((A['N']-1)//4)%2)
    e_pc=(k4+2*Cexp)%4        # (pi0/c)_4
    xi=np.exp(1j*omega*ANG)
    # conj((lambda^j i^e / c)) = i^{-(j*KL + e*KI)}
    lam_ph=(1j)**((-(j*KL+e*KI))%4)
    c1=(1j)**((-k*e_pc)%4)*GT*xi*lam_ph
    Ns=[NN[cop]]; Cs=[c1[cop]]; Ks=[CLS1[cop]]
    # pi0^{k+1} c'
    Gk1=Gj(pi0,k+1); pref=q**k*Gk1/q**((k+1)/2)
    kl0=sym_prime((1,1),pi0); ki0=sym_prime((0,1),pi0)
    lam0=(1j)**((-(j*kl0+e*ki0)*(k+1))%4)
    xi0=(complex(a,b)/abs(complex(a,b)))**(omega*(k+1))
    ph=(1j)**(((k+1)*(k4+e_pc)-k*e_pc)%4)
    c2=pref*xi0*lam0*ph*GT*xi*lam_ph
    P=complex(a,b)**(k+1); z2=(U+1j*Vv)*P
    cls2=((np.round(z2.real).astype(np.int64)-1)%4==0)&(np.round(z2.imag).astype(np.int64)%4==0)
    Ns.append(NN[cop]*q**(k+1)); Cs.append(c2[cop]); Ks.append(cls2[cop])
    Nn=np.concatenate(Ns); Cn=np.concatenate(Cs); Kn=np.concatenate(Ks)
    m=Nn<=NMAX
    Nn,Cn,Kn=Nn[m],Cn[m],Kn[m]
    if cls is not None:
        mm=(Kn==cls); return Nn[mm],Cn[mm],Kn[mm]
    return Nn,Cn,Kn
def coeffs_unit(j=0,e=0,omega=0):
    """nu = lambda^j i^e (no odd part)"""
    xi=np.exp(1j*omega*ANG); lam_ph=(1j)**((-(j*KL+e*KI))%4)
    return NN.copy(), GT*xi*lam_ph, CLS1.copy()
def smooth_sums(Nn,Cn,Xs,sig=0.5,Kn=None):
    lN=np.log(Nn); out=[]
    for X in Xs:
        t=lN-math.log(X); w=np.exp(-t*t/(2*sig*sig))
        out.append(np.sum(Cn*w))
    return np.array(out)
