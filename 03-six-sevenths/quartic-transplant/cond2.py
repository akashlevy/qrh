import numpy as np, math, sys
sys.path.insert(0,'.')
from conductor_test import *
# plain (non-flat) smoothed sums S(X) = sum_c xi(c) g~(pi0^k,c) V(N(c)/X), V Gaussian in log, sigma
def plain_sum(Nc,Cc,Xs,sig=0.5):
    lN=np.log(Nc); out=[]
    for X in Xs:
        t=lN-math.log(X); out.append(np.sum(Cc*np.exp(-t*t/(2*sig*sig))))
    return np.array(out)
omega=int(sys.argv[1]); sig=float(sys.argv[2])
Xs=10**np.arange(0.0,5.31,0.1)
# several primes of similar norm for averaging
from gi import primary, gnorm
def primes_near(qmin,qmax,count):
    out=[]
    for q in range(qmin,qmax):
        if q%4!=1 or any(q%p==0 for p in range(2,int(q**0.5)+1)): continue
        for a in range(1,int(q**0.5)+1):
            b2=q-a*a; b=int(round(b2**0.5))
            if b*b==b2: out.append(primary((a,b))); out.append(primary((a,-b))); break
        if len(out)>=count: break
    return out[:count]
for (qmin,qmax) in [(5,6),(13,14),(17,30),(37,60),(97,140)]:
    pis=primes_near(qmin,qmax,6)
    qs=[gnorm(p) for p in pis]
    print("== primes",pis,"norms",qs)
    for k in range(4):
        M=np.zeros(len(Xs))
        for pi0 in pis:
            Nc,Cc=coeffs(pi0,k,omega)
            S=plain_sum(Nc,Cc,Xs,sig)
            M+=np.abs(S)**2/Xs
        M/=len(pis)
        print(" k=%d"%k," ".join("%5.2f"%v for v in np.log10(M)[::3]))
    sys.stdout.flush()
print("log10 X grid:"," ".join("%5.1f"%v for v in np.log10(Xs)[::3]))
