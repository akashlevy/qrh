import sys; sys.path.insert(0,'.')
from gi import *
import numpy as np, math, cmath
P=[tuple(map(int,l.split()[:2])) for l in open('primes_2e4.txt')]
G={}
for l in open('sqf_2e4.txt'):
    u,v,N,re,im=l.split(); G[(int(u),int(v))]=complex(float(re),float(im))
# brute-force check composite squarefree c with N<=300
bad=0; n=0
for c,g in G.items():
    N=gnorm(c)
    if N>400 or N==1: continue
    gb=gauss_brute((1,0),c,P)/math.sqrt(N)
    n+=1
    if abs(gb-g)>1e-8: bad+=1; print("MISMATCH",c,N,g,gb)
print("checked",n,"bad",bad)
# check g^4 = pi^3 conj(pi) for split primes (normalized: gt^4 = pi^3 conj(pi)/N^2)
mx=0
for l in open('primes_2e4.txt'):
    a,b,N,re,im=l.split(); a=int(a);b=int(b);N=int(N)
    if b==0: continue
    gt=complex(float(re),float(im)); pi=complex(a,b)
    mx=max(mx,abs(gt**4 - pi**3*pi.conjugate()/N**2))
print("max |gt^4 - pi^3 pibar/N^2| =",mx)
