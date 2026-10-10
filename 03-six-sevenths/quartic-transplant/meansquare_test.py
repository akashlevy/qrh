# Mean square over squarefree primary rows h, N(h)<=H, of column sums of length X (Gaussian-log weight):
#  (a) actual dual family:  conj(alpha)(c) g~(lambda^-2 h^3, c)   = conj(alpha) g~(lambda^-2,c) chi_c(h)   [cube shift, l=-1]
#  (b) Patterson-bias family: g~(lambda^-2 h, c) = g~(lambda^-2,c) conj(chi_c(h))                         [first-power shift, l=0]
#  (c) spherical cube family: g~(lambda^-2,c) chi_c(h)                                                       [cube shift, l=0]
#  (d) control: random unimodular phases times chi_c(h)
# report R = (1/X) sum_h |B_h|^2  /  (diagonal)  ; orthogonal behaviour ~ 1 + X/H-ish ; bias family should grow like H^{3/4}X^{1/2}/(H) relative
import sys, numpy as np, math
sys.path.insert(0,'.')
from kub import *
from gi import primary, gnorm
rng=np.random.default_rng(1)
P=[tuple(map(int,l.split()[:2])) for l in open('primes_2e6.txt')]
Pn=[(p,gnorm(p)) for p in P]
def sqf_primary_upto(H):
    # squarefree primary h with N(h)<=H, with factorization
    ps=sorted([p for p,n in Pn if n<=H and p[1]!=0],key=gnorm)   # split primes only
    out=[((1,0),[])]
    def rec(cur,fs,start):
        for i in range(start,len(ps)):
            p=ps[i]; n=gnorm(cur)*gnorm(p)
            if n>H: break
            new=gmul(cur,p); out.append((new,fs+[p])); rec(new,fs+[p],i+1)
    rec((1,0),[],0)
    return out
sig=0.5
def run(X,H):
    m=NN<=20*X
    u=U[m]; v=Vv[m]; n=NN[m]; gt=GT[m]; kl=KL[m]; Nc=A['N'][m]
    lam=(1j)**((2*kl)%4)           # g~(lambda^-2,c) = (lambda^2/c)_4 g~(1,c)
    base=gt*lam
    w=np.exp(-(np.log(n/X))**2/(2*sig*sig))
    alb=np.exp(-1j*ANG[m])
    rnd=np.exp(2j*np.pi*rng.random(len(n)))
    hs=[t for t in sqf_primary_upto(H) if t[0]!=(1,0)]
    tabs={}
    acc=np.zeros(4); 
    for h,fs in hs:
        k=np.zeros(len(n),dtype=np.int64); cop=np.ones(len(n),dtype=bool)
        for p in fs:
            if p not in tabs: tabs[p]=qtab(p)
            q,r,tab=tabs[p]
            red=(u%q+(v%q)*r)%q; cop&=(red!=0); k+=np.where(red!=0,tab[red],0)
        # (h/c)_4 = (-1)^{C(h,c)} (c/h)_4
        Ch=((gnorm(h)-1)//4)%2; C=(Ch*(((Nc-1)//4)%2))%2
        e=(k+2*C)%4
        chi=np.where(cop,(1j)**e,0)
        Ba=np.sum(alb*base*chi*w); Bb=np.sum(base*np.conj(chi)*w); Bc=np.sum(base*chi*w); Bd=np.sum(rnd*chi*w)
        acc+=np.array([abs(Ba)**2,abs(Bb)**2,abs(Bc)**2,abs(Bd)**2])
    diag=len(hs)*np.sum(w*w)
    return len(hs), acc/diag
print("ratios (1/X sum_h|B_h|^2)/diag ;  columns: (a) actual cube+angular, (b) first-power shift (biased), (c) cube spherical, (d) random")
for X in (1000,4000):
    for H in (250,1000,4000):
        nh,R=run(X,H)
        print("X=%5d H=%5d #h=%5d :  a=%.3f  b=%.3f  c=%.3f  d=%.3f"%(X,H,nh,*R)); sys.stdout.flush()
