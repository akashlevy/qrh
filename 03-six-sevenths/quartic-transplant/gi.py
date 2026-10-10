# Gaussian integer helpers (python ints), conventions as in qgauss.c
import math, cmath
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gmod(a,m):
    # remainder of a mod m with rounding
    n=gnorm(m); num=gmul(a,(m[0],-m[1]))
    q=(round(num[0]/n), round(num[1]/n))
    qm=gmul(q,m); return (a[0]-qm[0], a[1]-qm[1])
def is_primary(a): return a[0]%2==1 and a[1]%2==0 and (a[0]+a[1])%4==1
def primary(a):
    for _ in range(4):
        if is_primary(a): return a
        a=(-a[1],a[0])
    raise ValueError(a)
def gpow_mod(a,e,m):
    r=(1,0); a=gmod(a,m)
    while e:
        if e&1: r=gmod(gmul(r,a),m)
        a=gmod(gmul(a,a),m); e>>=1
    return r
UNITS={(1,0):0,(0,1):1,(-1,0):2,(0,-1):3}
def qsym_prime(x,pi):
    """(x/pi)_4 as exponent k (i^k), None if pi|x"""
    if gmod(x,pi)==(0,0): return None
    N=gnorm(pi); w=gpow_mod(x,(N-1)//4,pi)
    for u,k in UNITS.items():
        if gmod((w[0]-u[0],w[1]-u[1]),pi)==(0,0): return k
    raise ValueError
def factor(c, primes):
    """factor primary c into primary primes from list"""
    fs=[]; N=gnorm(c)
    for pi in primes:
        n=gnorm(pi)
        if n>N: break
        while N%n==0 and gmod(c,pi)==(0,0):
            # divide
            num=gmul(c,(pi[0],-pi[1])); c=(num[0]//n, num[1]//n); N//=n; fs.append(pi)
    return fs, c
def qsym(x, c_factors):
    k=0
    for pi in c_factors:
        t=qsym_prime(x,pi)
        if t is None: return None
        k+=t
    return k%4
def gauss_brute(nu, c, primes):
    """g4(nu,c)=sum_{d mod c} (d/c)_4 e(nu d/c), e(z)=exp(4 pi i Re z); nu Gaussian integer or tuple of fractions (num,den)"""
    fs,_=factor(c,primes)
    N=gnorm(c); tot=0
    # complete residue system: x+yi, 0<=x<N, 0<=y<N each class appears N times
    seen={}
    for x in range(N):
        for y in range(N):
            r=gmod((x,y),c)
            if r in seen: continue
            seen[r]=1
    for r in seen:
        k=qsym(r,fs)
        if k is None: continue
        z=gmul(nu,r); # z/c = z*conj(c)/N
        num=gmul(z,(c[0],-c[1])); re=num[0]/N
        tot+= (1j)**k * cmath.exp(4j*math.pi*re)
    return tot
