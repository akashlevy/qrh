# Part II bookkeeping, high side idealized to GRH for ALL non-principal rows (floor bin at Re s = 1/2,
# R = 1, no intermediate bins).  Explore which ZERO-BLIND low-side changes would be needed to reach 6/7.
# Variables x = [lx, ly, ell, s, g, e]; minimize s.
from scipy.optimize import linprog
import numpy as np
def solve(gram_b6=True, gram_2b=True, w=0.0, energyA=True, small_rows=False, gram_cond=False, delta0=0.0):
    A=[];b=[]
    def ge(c, r): A.append([-t for t in c]); b.append(-r)
    ge([-1/6,0,1/6,1,-1/2,-1/2], 5/6)                 # low side
    ge([0,0,0,0,1,0],0)
    if gram_b6: ge([1/6,-1/6,0,0,1,0],0)                # g >= b/6
    if gram_2b: ge([2,-1,0,0,1,0],0)                    # g >= 2b - ly = ly - 2lx
    ge([0,0,0,0,0,1],0)
    if energyA: ge([2/4,2/4,-3/4,0,0,1],1/4)            # e >= (1+3ell-2M)/4  (case A, untouched)
    ge([-1,-1,-1,0,0,1],-1-w)                           # e >= M+ell-1-w   (wall, relaxed by w)
    # floor bin (GRH if delta0=0):  s >= (1+d0)(1-ly)/2 - ell/2 + (d0/2) ell + h(1 + d0/2 - 1/6), h = 1-lx+ell
    d0=delta0; c=(1+d0/2-1/6)
    ge([c,(1+d0)/2,1/2-d0/2-c,1,0,0], (1+d0)/2 + c)
    if small_rows:  # h(17/50-1/6) - ly/2 + 63/50*0.01 <= 0
        k=17/50-1/6
        A.append([-k,-1/2,k,0,0,0]); b.append(-k-0.0126)   # k(1-lx+ell) - ly/2 + .0126 <= 0
    if gram_cond:   # ly - ell - 11 b/6 >= 0  (b = ly-lx)
        ge([11/6,1-11/6,-1,0,0,0],0)
    bounds=[(0,None)]*3+[(None,None)]*3
    r=linprog([0,0,0,1,0,0],A_ub=np.array(A),b_ub=np.array(b),bounds=bounds,method='highs')
    return r.fun, r.x
print("6/7 =",6/7, " 13/15 =",13/15)
for kw in [dict(), dict(small_rows=True,gram_cond=True), dict(delta0=1/50),
           dict(gram_b6=False), dict(gram_b6=False,gram_2b=False),
           dict(w=0.2), dict(w=0.2,gram_b6=False), dict(energyA=False,w=1.0),
           dict(energyA=False,w=1.0,gram_b6=False), dict(gram_b6=False,small_rows=True,gram_cond=True)]:
    v,x=solve(**kw); print(kw, round(v,6), np.round(x[:3],4))
