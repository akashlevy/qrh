# Direct check of which shift appears after Poisson + absorption (quartic, Z[i]).
# B_h column coefficient (from Poisson in u):  mu(z) gamma_{-1}(z) chi_z(h),  chi_z(x) = (x/z)_4,
#   gamma_j(z) = N^{-1/2} sum_x chi_z(x)^j e(x/z),  e(w) = exp(2 pi i Im w)   (as in quartic_absorption_check.py)
# Absorption: mu gamma_{-1} = Phi * conj(alpha) * gamma_1, Phi ray-class.
# Claim: chi_z(h) gamma_1(z) = gt(h^3, z)  (shift h^3),  NOT gt(h, z).
import sys, math, cmath, itertools
sys.path.insert(0,'.')
from gi import *
def e_im(w): return cmath.exp(2j*math.pi*w.imag)
def residues(c):
    N=gnorm(c); seen=set(); out=[]
    for x in range(N):
        for y in range(N):
            r=gmod((x,y),c)
            if r not in seen: seen.add(r); out.append(r)
    return out
P=sorted([tuple(map(int,l.split()[:2])) for l in open('primes_2e4.txt')],key=gnorm)
def gauss(c, j, shift):
    """N^{-1/2} sum_x chi_c(x)^j e(shift*x/c)"""
    fs,_=factor(c,P); N=gnorm(c); tot=0
    cc=complex(*c)
    for r in residues(c):
        k=qsym(r,fs)
        if k is None: continue
        w=complex(*gmul(shift,r))/cc
        tot+=(1j)**((j*k)%4)*e_im(w)
    return tot/math.sqrt(N)
zs=[c for c in [primary((1,2)),primary((3,2)),primary((1,4)),primary((5,2)),(-3,0),primary((1,6))]]
zs+= [primary(gmul(primary((1,2)),primary((3,2))))]
hs=[(2,1),(1,3),(5,0),(3,4),(7,2),(2,5)]
print("z, N(z): ratio  mu*gamma_{-1}*chi_z(h) / (conj(alpha) * G(r)) for r=h^3 (should be h-independent) and r=h")
for z in zs:
    fs,_=factor(z,P); mu=(-1)**len(fs); N=gnorm(z); al=complex(*z)/abs(complex(*z))
    gm1=gauss(z,-1,(1,0))
    r3=[]; r1=[]
    for h in hs:
        k=qsym(h,fs)
        if k is None: continue
        lhs=mu*gm1*(1j)**k
        h3=gmul(h,gmul(h,h))
        G3=gauss(z,1,h3); G1=gauss(z,1,h)
        r3.append(lhs/(al.conjugate()*G3)); r1.append(lhs/(al.conjugate()*G1))
    sp3=max(abs(a-r3[0]) for a in r3); sp1=max(abs(a-r1[0]) for a in r1)
    print(z,N," r=h^3: Phi=%s spread=%.1e |  r=h: spread=%.2f"%(complex(round(r3[0].real,6),round(r3[0].imag,6)),sp3,sp1))
