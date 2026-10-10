# Check: for primary prime p and primary c' coprime (squarefree), with shift nu coprime to pc',
#   g4(nu, p c') = eps(p,c') * g4(nu, p) * g4(nu * p^k, c')   with eps a sign independent of nu,
# holds for k = 2 (mark -> shift p^2, conductor N(p)^2 on the cofactor), not for k = 1 or 3.
import sys, math, cmath
sys.path.insert(0,'/private/tmp/claude-501/-Users-akashlevy-Documents-riemann/eab6979c-f095-4107-a475-7e51d88a9fe7/scratchpad/wf67q')
from gi import *
def e_im(w): return cmath.exp(2j*math.pi*w.imag)
P=sorted([tuple(map(int,l.split()[:2])) for l in open('/private/tmp/claude-501/-Users-akashlevy-Documents-riemann/eab6979c-f095-4107-a475-7e51d88a9fe7/scratchpad/wf67q/primes_2e4.txt')],key=gnorm)
def residues(c):
    N=gnorm(c); seen=set(); out=[]
    for x in range(N):
        for y in range(N):
            r=gmod((x,y),c)
            if r not in seen: seen.add(r); out.append(r)
    return out
def gauss(c, shift):
    fs,_=factor(c,P); N=gnorm(c); tot=0; cc=complex(*c)
    for r in residues(c):
        k=qsym(r,fs)
        if k is None: continue
        tot+=(1j)**k*e_im(complex(*gmul(shift,r))/cc)
    return tot/math.sqrt(N)
ps=[primary((1,2)),primary((3,2)),primary((1,4))]
cs=[primary((5,2)),(-3,0),primary((1,6)),primary((5,4))]
nus=[(1,0),(2,1),(1,3),(3,4),(7,2)]
for p in ps:
  for c in cs:
    if gnorm(p)==gnorm(c): continue
    pc=primary(gmul(p,c))
    out=[]
    for k in (1,2,3):
        ratios=[]
        for nu in nus:
            if gmod(nu,p)==(0,0) or gmod(nu,c)==(0,0): continue
            pk=(1,0)
            for _ in range(k): pk=gmul(pk,p)
            lhs=gauss(pc,nu); rhs=gauss(p,nu)*gauss(c,gmul(nu,pk))
            ratios.append(lhs/rhs)
        spread=max(abs(r-ratios[0]) for r in ratios)
        out.append("k=%d: ratio=%s spread=%.1e"%(k, complex(round(ratios[0].real,6),round(ratios[0].imag,6)), spread))
    print("p=%s c'=%s :  "%(p,c)+" | ".join(out))
