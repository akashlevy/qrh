# Numerical test of the conductor of c -> g~4(pi0^k, c) (Kubota series psi(pi0^k, s, xi)).
# S_flat(X) = sum_{d primary} xi(d)^4 N(d) sum_{c primary} xi(c) g~(pi0^k, c) V(N(c) N(d)^4 / X)
#   (Dirichlet series zeta_lambda(4s-3,xi) psi(pi0^k, s+1/2, xi): entire for omega != 0, clean FE with N(pi0^k)^{1-s})
import numpy as np, math, cmath, sys
sys.path.insert(0,'.')
from gi import *
dt=np.dtype([('u','<i8'),('v','<i8'),('N','<i8'),('re','<f8'),('im','<f8')])
A=np.fromfile('sqf_2e6.bin',dtype=dt)
U=A['u']; Vv=A['v']; NN=A['N'].astype(np.float64); GT=A['re']+1j*A['im']
NMAX=2e6
ang=np.angle(U+1j*Vv)
def qsym_table(pi0):
    a,b=pi0; q=a*a+b*b
    Binv=pow(b%q,q-2,q); r=(-a*Binv)%q
    # chi(x) for x in F_q as exponent
    tab=np.full(q,-1,dtype=np.int64)
    for x in range(1,q):
        w=pow(x,(q-1)//4,q)
        tab[x]={1:0, r:1, q-1:2, (q-r)%q:3}[w]
    return q,r,tab
def gauss_local(pi0,j):
    a,b=pi0; q=a*a+b*b; _,r,tab=qsym_table(pi0)
    # G_j = sum_{d mod pi0} (d/pi0)^j e(d/pi0), d integer reps, e = exp(4 pi i d a/q)
    s=0
    for d in range(1,q):
        s+= (1j)**((j*tab[d])%4) * cmath.exp(4j*math.pi*d*a/q)
    return s
def coeffs(pi0,k,omega):
    """return arrays (N, coef) for the Dirichlet series psi(pi0^k, s+1/2, xi) (normalized), over all primary c<=NMAX"""
    q,r,tab=qsym_table(pi0)
    a,b=pi0
    red=(U%q + (Vv%q)*r)%q
    cop = red!=0
    k4 = np.where(cop, tab[red], 0)       # (c/pi0)_4 exponent
    Cexp = (((q-1)//4)%2) * (((A['N']-1)//4)%2)   # C(pi0,c) mod 2
    # (pi0/c)_4 = (-1)^C (c/pi0)_4  -> exponent
    e_pc = (k4 + 2*Cexp)%4
    xi = np.exp(1j*omega*ang)
    # coprime part: g~(pi0^k,c) = conj((pi0/c))^k g~(1,c)
    c1 = np.where(cop, (1j)**((-k*e_pc)%4) * GT * xi, 0)
    Ns=[NN[cop]]; Cs=[c1[cop]]
    # pi0^{k+1} c' part
    Gk1=gauss_local(pi0,k+1)
    pref = q**k * Gk1 / q**((k+1)/2)
    xi0 = (complex(a,b)/abs(complex(a,b)))**omega
    # [(c'/pi0)(pi0/c')]^{k+1} conj(pi0/c')^k
    ph = (1j)**(((k+1)*(k4+e_pc) - k*e_pc)%4)
    c2 = pref * xi0**(k+1) * ph * GT * xi
    Ns.append(NN[cop]*q**(k+1)); Cs.append(np.where(cop,c2,0)[cop])
    Nall=np.concatenate(Ns); Call=np.concatenate(Cs)
    m=Nall<=NMAX
    return Nall[m], Call[m]
def flat_sum(Nc,Cc,omega,Xs,sig=0.5):
    # d^4 terms: primary d with N(d)^4 <= NMAX
    ds=[]
    R=int((NMAX)**0.125)+2
    for x in range(-R,R+1):
        for y in range(-R,R+1):
            if is_primary((x,y)) and (x*x+y*y)**4<=NMAX:
                ds.append((x*x+y*y, cmath.exp(1j*4*omega*math.atan2(y,x))))
    out=[]
    lN=np.log(Nc)
    for X in Xs:
        tot=0
        for nd,xd in ds:
            t=lN+4*math.log(nd)-math.log(X)
            w=np.exp(-t*t/(2*sig*sig))
            tot+= xd*nd*np.sum(Cc*w)
        out.append(tot)
    return np.array(out)
if __name__=="__main__":
    Xs=10**np.arange(0.0,5.01,0.1)
    pis=[primary((1,2)),primary((3,2)),primary((1,4)),primary((5,2))]
    omega=int(sys.argv[1]) if len(sys.argv)>1 else 1
    res={}
    for pi0 in pis:
        for k in range(4):
            Nc,Cc=coeffs(pi0,k,omega)
            S=flat_sum(Nc,Cc,omega,Xs)
            res[(pi0,k)]=S
            q=gnorm(pi0)
            print("pi0",pi0,"q",q,"k",k, " ".join("%.2f"%v for v in np.log10(np.abs(S)/np.sqrt(Xs)+1e-30)[::5]))
            sys.stdout.flush()
    np.save('cond_res_om%d.npy'%omega, {'Xs':Xs,'res':{str(k):v for k,v in res.items()}}, allow_pickle=True)
