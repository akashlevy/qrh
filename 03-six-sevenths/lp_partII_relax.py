# LP lower bounds for Part II with GRH-for-rows high side and various low-side relaxations.
# vars: lx, ly, ell, s, g, e
import numpy as np
from scipy.optimize import linprog
from itertools import product
def solve(use_g=True, use_eA=True, use_eB=True, wB=0.0, high='grh'):
    A=[];b=[]
    def ge(coef, rhs): A.append([-c for c in coef]); b.append(-rhs)
    # low: s >= 5/6 + lx/6 - ell/6 + g/2 + e/2
    ge([-1/6,0,1/6,1,-1/2,-1/2], 5/6)
    ge([0,0,0,0,1,0],0); ge([0,0,0,0,0,1],0)
    if use_g:
        ge([1/6,-1/6,0,0,1,0],0)        # g >= (ly-lx)/6
        ge([2,-1,0,0,1,0],0)            # g >= 2b - ly = ly - 2lx
    if use_eA:
        ge([1/2,1/2,-3/4,0,0,1],1/4)    # e >= (1+3ell-2M)/4  (case A, valid region assumed)
    if use_eB:
        ge([-1,-1,-1,0,0,1],-1-wB)      # e >= M+ell-1-w
    # high, GRH for rows: s >= (1-ly)/2 - ell/2 + 5/6*(1-lx+ell)
    ge([5/6,1/2,1/2-5/6,1,0,0], 1/2+5/6)
    bounds=[(0,1),(0,1),(0,1),(None,None),(0,None),(0,None)]
    r=linprog([0,0,0,1,0,0],A_ub=A,b_ub=b,bounds=bounds,method='highs')
    return r.fun, r.x
for ug,ua,ub in product([True,False],repeat=3):
    v,x=solve(ug,ua,ub)
    print(f"gram={ug!s:5} caseA={ua!s:5} caseB={ub!s:5}: sigma0 >= {v:.6f}  (lx,ly,ell)=({x[0]:.4f},{x[1]:.4f},{x[2]:.4f})  6/7={6/7:.6f}")
print("--- GRH rows, gram+caseA, caseB relaxed by w")
from fractions import Fraction
for w in [0,0.05,0.1,0.2,0.3,0.4,4/7,0.6,0.8,1.0]:
    v,x=solve(True,True,True,wB=w)
    print(f"w={w:.4f}: {v:.6f} (lx,ly,ell)=({x[0]:.4f},{x[1]:.4f},{x[2]:.4f})  frac~{Fraction(v).limit_denominator(2000)}")
