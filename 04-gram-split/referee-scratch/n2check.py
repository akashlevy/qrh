# n=2 rank-one check: g(m,c)=sum_{v mod c} (v/c) e(mv/c); p-part of D(s;mu^3) at p|mu
import numpy as np
from sympy import jacobi_symbol
def g(m,c):
    return sum(jacobi_symbol(v,c)*np.exp(2j*np.pi*m*v/c) for v in range(c) if np.gcd(v,c)==1)
for p in [3,5,7]:
    for u in [1,2]:
        vals=[g(u*p**3,p**k) if k>0 else 1 for k in range(6)]
        print(p,u,[round(complex(x).real,6)+0 for x in vals])
# check P(s)=1+phi(p^2)x-p^3x^2 (x=p^{-2s}) = (1-p x)(1+p^2 x)
