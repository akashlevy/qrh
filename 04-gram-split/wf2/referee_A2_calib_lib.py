import numpy as np, sys, time
sys.path.insert(0,"/Users/akashlevy/Documents/riemann/gram_split"); sys.path.insert(0,"/Users/akashlevy/Documents/riemann/gram_split/wf2")
from zomega import *
from a2lib import bump, UNITS
def jacobi(a, n):
    a = a % n; n = n.copy(); res = np.ones_like(a); a = a.copy()
    while True:
        act = a != 0
        if not act.any(): break
        # remove factors of 2
        while True:
            ev = act & (a % 2 == 0)
            if not ev.any(): break
            a[ev] //= 2
            r8 = n[ev] % 8
            res[ev] *= np.where((r8 == 3) | (r8 == 5), -1, 1)
        # swap
        sw = act
        aa = a[sw].copy(); nn = n[sw].copy()
        flip = (aa % 4 == 3) & (nn % 4 == 3)
        res[sw] *= np.where(flip, -1, 1)
        a[sw] = nn % aa; n[sw] = aa
    return np.where(n == 1, res, 0)
def sqfree_mask(N):
    Nmax=int(N.max()); s=np.ones(Nmax+1,bool)
    for p in sieve(int(np.sqrt(Nmax))+1):
        p=int(p); s[p*p::p*p]=False
    return s[N]
def rows_for(Q):
    a,b,N = elements_in_annulus(Q/2, Q)
    keep = sqfree_mask(N)          # squarefree norm => m primitive, squarefree, no conj pairs, no inert
    a,b,N = a[keep],b[keep],N[keep]
    r = (-a * np.array([pow(int(x), -1, int(n)) for x,n in zip(b,N)])) % N   # omega = r mod m
    # unit class: m = u R with R primary
    u = np.zeros(len(a),int)
    for i,(ua,ub) in enumerate(UNITS):
        # R = m / u  -> m * conj(u) since |u|=1: conj(ua+ub w) = (ua-ub) - ub w
        x,y = mul(a,b,ua-ub,-ub)
        u[(x%3==1)&(y%3==0)] = i
    return a,b,N,r,u
def nus(NUMAX):
    Pn = primes_upto(NUMAX)
    return [(x[0],x[1],x[2]) for x in squarefree_products(Pn,NUMAX) if x[2]>=7]
if __name__=="__main__":
    # consistency check vs their TW at Q=6000
    sys.path.insert(0,"/Users/akashlevy/Documents/riemann/gram_split/wf2")
    from a2lib import Rows
    R = Rows(6000)
    a,b,N,r,u = rows_for(6000)
    key = {(int(x),int(y)):i for i,(x,y) in enumerate(zip(R.ma,R.mb))}
    ok = np.array([(int(x),int(y)) in key for x,y in zip(a,b)]); a,b,N,r,u = a[ok],b[ok],N[ok],r[ok],u[ok]
    idx = np.array([key[(int(x),int(y))] for x,y in zip(a,b)]); print("rows compared", len(idx))
    NU = nus(128)
    bad=0; tot=0
    for (va,vb,vN) in NU:
        J = jacobi((va + vb*r) % N, N)
        # theirs: prod over p|R of chi_p(nu)^3
        th=[]
        for t in idx:
            e=0; z=False
            for k in R.fac[t]:
                s=R.P[k].symbol(np.array([va]),np.array([vb]))[0]
                if s<0: z=True
                e+=3*s
            th.append(0 if z else (1 if e%6==0 else -1))
        th=np.array(th); bad+=(th!=J).sum(); tot+=len(th)
    print("Jacobi-vs-a2lib quadratic symbol mismatches:", bad, "of", tot)
