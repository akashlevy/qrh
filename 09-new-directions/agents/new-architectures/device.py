"""Generic 'device' of log-length lam added to the P78 template (GRH for rows):
   low side gains gL*lam (signal/low), rows/signal worsen by gR*lam, wall and case A untouched.
   Implemented by a 4th geometry variable.  Prints theta(lam-optimized) for (gL, gR)."""
import numpy as np
from scipy.optimize import linprog
import template_lp as T

def solve_dev(gL, gR, lam_max=1.0, over=None):
    cfg = dict(T.BASE); cfg.update(over or {})
    c, a, rho, kU, gb, gq, r = (cfg[k] for k in ("c","a","rho","kU","gb","gq","r"))
    A, B = [], []
    def ge(coef, const):  # sigma >= coef.(lx,ly,ell,lam) + const
        A.append(list(coef) + [-1.0]); B.append(-const)
    def le(coef, const):
        A.append(list(coef) + [0.0]); B.append(-const)
    LX=np.array([1,0,0,0.]); LY=np.array([0,1,0,0.]); EL=np.array([0,0,1,0.]); LA=np.array([0,0,0,1.])
    for lam in (0.0,1.0):
        d = lam*EL; Mp = LX+LY-2*d; lp = EL-d
        base = (c*LX - a*EL - d - gL*LA, 1-c)
        gp=[(np.zeros(4),0.0),(gb*(LY-LX),0.0),(gq*(2*(LY-LX)-LY+d),0.0)]
        ep=[(np.zeros(4),0.0),(0.75*lp-0.5*Mp,0.25),((r-1)*(Mp+lp),-1.0),((1/3)*Mp+(4/3)*lp,-2/3)]
        for g in gp:
            for e in ep:
                ge(base[0]+0.5*g[0]+0.5*e[0], base[1]+0.5*g[1]+0.5*e[1])
    ge(-0.5*LY - rho*EL + (1-c)*(-LX+kU*EL) + gR*LA, 0.5+(1-c))
    le(-(LX+LY)+2*EL, 0.05); le(-LY+EL,0.05); le(-LX+EL,0.0)
    res = linprog([0,0,0,0,1], A_ub=np.array(A), b_ub=np.array(B),
                  bounds=[(0,2),(0,3),(0,2),(0,lam_max),(None,None)], method="highs")
    return res.x

if __name__ == "__main__":
    print(" gL     gR   ratio  theta    lam")
    for gL in (1/24, 1/12, 1/6):
        for k in (1, 2, 3, 4, 5, 8, 14):
            x = solve_dev(gL, k*gL)
            print(f"{gL:.4f} {k*gL:.4f} {k:5.1f}  {x[-1]:.5f}  {x[3]:.3f}")
