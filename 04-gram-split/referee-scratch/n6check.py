# sextic model at a split prime p=1 mod 6: chi = sextic Dirichlet char mod p, chi_{p^k}=chi^k
import numpy as np
from sympy import primitive_root
for p in [7,13]:
    r=primitive_root(p); ind={pow(r,i,p):i for i in range(p-1)}
    chi=lambda v,j: np.exp(2j*np.pi*j*ind[v%p]/6)
    for u in [1,2,3]:
        out=[]
        for k in range(0,6):
            if k==0: out.append(1); continue
            M=p**k; s=sum(chi(v,k)*np.exp(2j*np.pi*u*p**3*v/M) for v in range(M) if v%p)
            out.append(round(abs(s),3))
        print(p,u,"|g(u p^3, p^k)| k=0..5:",out, " p^3*sqrt(p)=",round(p**3*np.sqrt(p),3))
