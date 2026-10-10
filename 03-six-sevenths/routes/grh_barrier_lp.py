# LP: Part II bookkeeping with GRH for all non-principal rows (floor bin at 1/2, R=1, q=0).
# Variables: lx, ly, ell, s (=sigma0), g, e.  Minimize s.
# Low side (paper's Prop probe-low with Gram + energy):  s >= 5/6 + lx/6 - ell/6 + g/2 + e/2
#   g >= 0, g >= b/6, g >= 2b - ly  (b = ly - lx)
#   e >= 0, e >= (1+3ell-2M)/4, e >= M+ell-1   (M = lx+ly)   [paper's own energy bound]
# High side, floor bin with GRH (delta=0, R=1): s >= (1-ly)/2 - ell/2 + 5h/6,  h = 1 - lx + ell
from scipy.optimize import linprog
from fractions import Fraction as Fr
import numpy as np, itertools

def solve(energy='paper', extra=False):
    # x = [lx, ly, ell, s, g, e]
    A=[];b=[]
    def ge(coeffs, rhs):  # coeffs.x >= rhs  ->  -coeffs.x <= -rhs
        A.append([-c for c in coeffs]); b.append(-rhs)
    ge([-1/6,0,1/6,1,-1/2,-1/2], 5/6)
    ge([0,0,0,0,1,0],0)
    ge([1/6,-1/6,0,0,1,0],0)      # g >= (ly-lx)/6
    ge([2,-1,0,0,1,0],0)          # g >= 2(ly-lx) - ly = ly - 2lx ... careful
    ge([0,0,0,0,0,1],0)
    if energy=='paper':
        ge([2/4,2/4,-3/4,0,0,1],1/4)   # e >= (1+3ell-2M)/4
    else:
        ge([2/3,2/3,-2/3,0,0,1],1/3)   # e >= (1+2ell-2M)/3
    ge([-1,-1,-1,0,0,1],-1)            # e >= M+ell-1
    # GRH high: s >= (1-ly)/2 - ell/2 + 5(1-lx+ell)/6
    ge([5/6,1/2,1/2-5/6,1,0,0], 1/2+5/6)
    bounds=[(0,None),(0,None),(0,None),(None,None),(None,None),(None,None)]
    if extra:
        ge([-1,1,0,0,0,0],0)  # b>=0
    r=linprog([0,0,0,1,0,0],A_ub=np.array(A),b_ub=np.array(b),bounds=bounds,method='highs')
    return r

for en in ['paper','old']:
    for ex in [False,True]:
        r=solve(en,ex)
        print(en, 'b>=0' if ex else '', r.status, r.fun, np.round(r.x,6))
