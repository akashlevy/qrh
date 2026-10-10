# Polar (theta) term test, omega=0: S_{nu,u}(X) = sum_{c = u mod 4} g~(nu,c) V(N(c)/X) ~ tau(nu;u) Vhat(3/4) X^{3/4} + O(X^{1/2+})
# nu in {1, pi, pi^2, pi^3}.  Suzuki: tau(pi^3)=0 ; tau(pi^2;v) = (-1)^{C(pi,v)} g4(1,pi)/N^{3/4} tau(1; pi v)
import numpy as np, math, cmath, sys
sys.path.insert(0,'.')
from conductor_test import *
sig=0.5
Vhat34=math.sqrt(2*math.pi)*sig*math.exp(sig*sig*(0.75**2)/2)
# class mod 4 of primary c: c = 1 mod 4  <=> (u-1)%4==0 and v%4==0
cls_all = ((U-1)%4==0)&(Vv%4==0)
def coeffs_cls(pi0,k):
    # replicate coeffs() but keep track of class of c mod 4
    q,r,tab=qsym_table(pi0); a,b=pi0
    red=(U%q+(Vv%q)*r)%q; cop=red!=0
    k4=np.where(cop,tab[red],0); Cexp=(((q-1)//4)%2)*(((A['N']-1)//4)%2); e_pc=(k4+2*Cexp)%4
    c1=np.where(cop,(1j)**((-k*e_pc)%4)*GT,0)
    Gk1=gauss_local(pi0,k+1); pref=q**k*Gk1/q**((k+1)/2)
    ph=(1j)**(((k+1)*(k4+e_pc)-k*e_pc)%4); c2=pref*ph*GT
    # class of pi0^{k+1} c'
    P=complex(a,b)**(k+1)
    def cls_of(uu,vv):
        z=(uu+1j*vv)
        return ((np.round(z.real).astype(np.int64)-1)%4==0)&(np.round(z.imag).astype(np.int64)%4==0)
    z2=(U+1j*Vv)*P
    cls2=((np.round(z2.real).astype(np.int64)-1)%4==0)&(np.round(z2.imag).astype(np.int64)%4==0)
    N1=NN[cop]; C1=c1[cop]; K1=cls_all[cop]
    N2=NN[cop]*q**(k+1); C2=c2[cop]; K2=cls2[cop]
    Nn=np.concatenate([N1,N2]); Cn=np.concatenate([C1,C2]); Kn=np.concatenate([K1,K2])
    m=Nn<=NMAX
    return Nn[m],Cn[m],Kn[m]
Xs=10**np.arange(2.0,5.01,0.25)
def run(pi0,k):
    Nn,Cn,Kn=coeffs_cls(pi0,k)
    lN=np.log(Nn); out={True:[],False:[]}
    for X in Xs:
        t=lN-math.log(X); w=np.exp(-t*t/(2*sig*sig))
        for cl in (True,False):
            out[cl].append(np.sum(Cn[Kn==cl]*w[Kn==cl])/(Vhat34*X**0.75))
    return {cl:np.array(v) for cl,v in out.items()}
# nu = 1 : use k=0 with any pi0 but that includes pi0 c' terms with k+1=1 ... (k=0: g(1,pi0 c')) -> that is just all c. fine
for pi0 in [primary((1,2)),primary((3,2)),primary((1,4))]:
    q=gnorm(pi0)
    print("pi0",pi0,"q",q, " g~(pi0)=",complex(gauss_local(pi0,1))/math.sqrt(q))
    R={k:run(pi0,k) for k in range(4)}
    for k in range(4):
        for cl in (True,False):
            v=R[k][cl]
            print("  k=%d class %s:"%(k,'1   ' if cl else '1+l3'), " ".join("%7.3f%+7.3fi"%(z.real,z.imag) for z in v[::2]))
    sys.stdout.flush()
print("log10 X:", " ".join("%15.2f"%x for x in np.log10(Xs)[::2]))
