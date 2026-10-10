import numpy as np
from zomega import *
P = primes_upto(3000)
print(len(P), "primes; first:", [(p.a,p.b,p.N,p.kind) for p in P[:6]])
# Gauss sums: |g|^2 = N ; gamma_j for j=1,2,3
def gam(Pr, j):
    if Pr.kind=="split":
        p=Pr.q; x=np.arange(1,p); jj=Pr.symbol(x,np.zeros_like(x))
        return np.sum(Z6[(j*jj)%6]*np.exp(-2j*np.pi*((Pr.b*x)%p)/p))/np.sqrt(Pr.N)
    q=Pr.q; c,d=np.meshgrid(np.arange(q),np.arange(q),indexing="ij"); c=c.ravel(); d=d.ravel()
    nz=(c!=0)|(d!=0); c=c[nz]; d=d[nz]; jj=Pr.symbol(c,d); sgn=-1 if Pr.a<0 else 1
    return np.sum(Z6[(j*jj)%6]*np.exp(-2j*np.pi*sgn*d/q))/np.sqrt(Pr.N)
errs=[]; hd=[]
for Pr in P[:60]:
    g1=gam(Pr,1); g2=gam(Pr,2); g3=gam(Pr,3)
    al=cval(Pr.a,Pr.b)/abs(cval(Pr.a,Pr.b))
    errs.append((abs(abs(g1)-1), abs(g2**3+al)))
    # paper: gamma1 gamma2 = conj(chi_p(4)) gamma3 gamma2^3  -> gamma1 = conj(chi_p(4)) gamma3 gamma2^2
    j4=Pr.symbol(np.array([4]),np.array([0]))[0]
    hd.append(abs(g1 - np.conj(Z6[j4])*g3*g2**2))
print("max ||g1|-1|, |g2^3+alpha|:", np.max(errs,axis=0))
print("max HD residual gamma1 - conj(chi(4)) gamma3 gamma2^2:", max(hd))
# sextic reciprocity: chi_b(a) / chi_a(b) depends only on classes mod 4 (sign)
import itertools
bad=0; vals={}
for P1,P2 in itertools.combinations(P[:40],2):
    j12=P2.symbol(np.array([P1.a]),np.array([P1.b]))[0]   # (P1/P2)_6 = chi_{P2}(P1)
    j21=P1.symbol(np.array([P2.a]),np.array([P2.b]))[0]
    ratio=(j12-j21)%6
    e1=((P1.N-1)//2)%2; e2=((P2.N-1)//2)%2
    vals.setdefault((e1,e2),set()).add(ratio)
print("reciprocity ratio exponents by (eps1,eps2):", vals)
