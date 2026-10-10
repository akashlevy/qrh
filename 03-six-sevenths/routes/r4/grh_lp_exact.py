"""
Exact (rational) LP analysis of the GRH-for-rows regime of Part II with low-side variants.
Variables (lx, ly, ell, s).  Constraints (all must hold):
  LOW  at d in {0, ell}: s >= 5/6 + lx/6 - ell/6 - d + g/2 + e/2  for every pair of pieces
       g pieces: 0, gb*(ly-lx)/6, (2(ly-lx) - ly + d)
       e pieces: 0, caseA (paper (1+3l'-2M')/4 or refined (1+2l'-2M')/3), M'+l'-1-w, (M'+4l'-2)/3 - w
  HIGH (GRH: delta=0, R=1, q=0): s >= (1-ly)/2 - ell/2 + 5h/6,  h = 1 - lx + ell
  side: M - 2 ell >= 0, ly >= ell, lx >= ell
Solved with scipy (HiGHS) then the optimum is verified exactly by an explicit dual certificate
(nonnegative rational multipliers on active constraints, eliminating lx, ly, ell).
"""
import itertools
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog

def rows(gb, caseA, w):
    R = []   # each: (coef vector over (lx,ly,ell,s), rhs, name)   meaning coef.v >= rhs
    for lam in (0, 1):
        # d = lam*ell; M' = lx+ly-2 lam ell; l' = (1-lam) ell
        base = np.array([F(1,6), F(0), -F(1,6) - lam, F(0)])   # terms: lx/6 - ell/6 - d   (moved to rhs side later)
        gp = [("g0", np.array([F(0)]*4)),
              ("gb6", np.array([-F(gb,6), F(gb,6), F(0), F(0)])),
              ("gq", np.array([F(-2), F(1), F(lam), F(0)]))]
        if gb == 0: gp = [gp[0], gp[2]]
        ep = [("e0", np.array([F(0)]*4), F(0))]
        if caseA == "paper":   # (1 + 3 l' - 2 M')/4
            ep.append(("caseA", np.array([-F(1,2), -F(1,2), F(3*(1-lam),4) + lam, F(0)]), F(1,4)))
        elif caseA == "refined":   # (1 + 2 l' - 2M')/3
            ep.append(("caseA", np.array([-F(2,3), -F(2,3), F(2*(1-lam),3) + F(4*lam,3), F(0)]), F(1,3)))
        ep.append(("wall", np.array([F(1), F(1), -2*lam + (1-lam), F(0)]), -1 - w))
        ep.append(("third", np.array([F(1,3), F(1,3), F(-2*lam,3) + F(4*(1-lam),3), F(0)]), -F(2,3) - w))
        for (gn, gv), (en, ev, ec) in itertools.product(gp, ep):
            # s >= 5/6 + base.v + gv.v/2 + (ev.v + ec)/2   ->  s - base.v - gv/2 - ev/2 >= 5/6 + ec/2
            coef = np.array([F(0), F(0), F(0), F(1)]) - base - gv/2 - ev/2
            R.append((coef, F(5,6) + ec/2, f"low d={'ell' if lam else 0} {gn}+{en}"))
    # HIGH GRH: s >= 1/2 - ly/2 - ell/2 + 5/6 (1 - lx + ell)
    R.append((np.array([F(5,6), F(1,2), F(1,2) - F(5,6), F(1)]), F(1,2) + F(5,6), "high GRH"))
    R.append((np.array([F(1), F(1), F(-2), F(0)]), F(0), "M-2ell>=0"))
    R.append((np.array([F(0), F(1), F(-1), F(0)]), F(0), "ly>=ell"))
    R.append((np.array([F(1), F(0), F(-1), F(0)]), F(0), "lx>=ell"))
    return R

def solve(gb, caseA, w):
    R = rows(gb, caseA, w)
    A = -np.array([[float(c) for c in r[0]] for r in R]); b = -np.array([float(r[1]) for r in R])
    res = linprog([0, 0, 0, 1], A_ub=A, b_ub=b, bounds=[(0, 2)]*3 + [(None, None)], method="highs")
    # dual certificate: find y>=0 with sum y_i a_i = e_s, maximize sum y_i rhs_i  (exact via enumeration of active sets)
    x = res.x
    act = [i for i, r in enumerate(R) if abs(float(np.dot([float(c) for c in r[0]], x)) - float(r[1])) < 1e-7]
    best = None
    for S in itertools.combinations(act, 4):
        M = [[R[i][0][k] for i in S] for k in range(4)]
        # solve M y = e_s exactly
        aug = [row[:] + [F(int(k == 3))] for k, row in enumerate(M)]
        ok = True
        for col in range(4):
            piv = next((r for r in range(col, 4) if aug[r][col] != 0), None)
            if piv is None: ok = False; break
            aug[col], aug[piv] = aug[piv], aug[col]
            for r in range(4):
                if r != col and aug[r][col] != 0:
                    f = aug[r][col] / aug[col][col]
                    aug[r] = [a - f*c for a, c in zip(aug[r], aug[col])]
        if not ok: continue
        y = [aug[k][4] / aug[k][k] for k in range(4)]
        if min(y) < 0: continue
        val = sum(yi * R[i][1] for yi, i in zip(y, S))
        if best is None or val > best[0]:
            best = (val, [(R[i][2], yi) for yi, i in zip(y, S) if yi != 0])
    return res.fun, x, best

if __name__ == "__main__":
    for (gb, ca, w) in [(1, "paper", 0), (0, "paper", 0), (0, "refined", 0), (1, "refined", F(1,7)), (1, "paper", F(1,7)),
                        (1, "refined", F(1,10)), (1, "paper", F(4,7)), (0, "refined", F(1,7))]:
        val, x, cert = solve(gb, ca, w)
        print(f"gram_b={gb} caseA={ca:7s} w={str(w):5s}: LP={val:.8f}  geom(lx,ly,ell)={np.round(x[:3],6)}  "
              f"exact={cert[0] if cert else None}")
        if cert:
            for n, y in cert[1]: print(f"        {y} x [{n}]")
