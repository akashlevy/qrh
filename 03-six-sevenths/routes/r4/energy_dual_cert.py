"""
Rigorous upper bound for the refined-count energy LP by weak duality.
For each branch (u-branch x max-branch) the primal is: maximize c.x + c0 s.t. A x <= b(M',l',w), x >= 0,
with b affine in (M', l', w).  Any y >= 0 with A^T y >= c gives  max <= y.b + c0 for ALL parameters.
We collect optimal duals at random parameter points, rationalize them, verify A^T y >= c and y >= 0 EXACTLY,
and record each certified affine bound  y.b(M',l',w) + c0 - M'.  The claim to prove:
   for every branch and every (M',l',w), max E_ref - M' <= max(0, (1+2l'-2M')/3, M'+l'-1-w, (M'+4l'-2)/3-w)
is then reduced to: every certified affine bound we need is one of these four pieces (checked symbolically).
"""
import itertools, numpy as np
from fractions import Fraction as F
from scipy.optimize import linprog
names = "O H A2 A4 Az N0 B0 S4 za v lb".split(); ix = {k: i for i, k in enumerate(names)}; n = len(names)
def vec(**kw):
    r = [F(0)] * n
    for k, c in kw.items(): r[ix[k]] += F(c)
    return r
def add(a, b, s=1): return [x + s * y for x, y in zip(a, b)]
def scal(a, s): return [s * x for x in a]
# rhs entries are affine in (Mp, lp, w): represented as (const, cM, cl, cw)
def branch_data(ub, wb):
    Td = vec(H=2, A2=2, A4=2, Az=2, za=2, N0=-1, B0=-3)          # + Td0, Td0 = -(1+l')
    y = vec(v=1, lb=3)
    c = vec(O=F(1, 6), A2=F(2, 3), A4=F(1, 3))
    c0 = [F(0), F(0), F(0), F(0)]                                  # affine const term of objective
    if wb == "H": c = add(c, vec(H=1))
    else: c = add(c, vec(v=1, lb=1)); c0 = add(c0, [F(0), F(0), F(0), F(-1)])
    c = add(c, vec(S4=-1, Az=-1, B0=-1, za=1, lb=-1))
    uv = {"v": vec(v=1), "za": vec(za=1), "avg": vec(v=F(1, 3), za=F(1, 3))}
    c = add(c, uv[ub], -1)
    c = add(c, add(Td, y, -1), F(-1, 2)); c0 = add(c0, [F(1, 2), F(0), F(1, 2), F(0)])   # -(Td0)/2 = (1+l')/2
    A, b = [], []
    A.append(vec(H=1, O=1)); b.append([F(0), F(1), F(0), F(0)])                        # H + O <= M'
    A.append(vec(A2=2, A4=4, Az=6, O=-1)); b.append([F(0)] * 4)
    A.append(vec(N0=1, B0=1, S4=1, A4=-1)); b.append([F(0)] * 4)
    A.append(vec(A4=1, N0=-1, B0=-1, S4=-1)); b.append([F(0)] * 4)
    A.append(vec(za=1)); b.append([F(0), F(0), F(1), F(0)])                            # za <= l'
    A.append(add(y, Td, -1)); b.append([F(-1), F(0), F(-1), F(0)])                    # y - Td_lin <= -(1+l')
    for o in uv:
        if o != ub: A.append(add(uv[ub], uv[o], -1)); b.append([F(0)] * 4)
    if wb == "H": A.append(add(vec(v=1, lb=1), vec(H=1), -1)); b.append([F(0), F(0), F(0), F(1)])
    else: A.append(add(vec(H=1), vec(v=1, lb=1), -1)); b.append([F(0), F(0), F(0), F(-1)])
    return c, c0, A, b
def evalb(bb, Mp, lp, w): return [r[0] + r[1] * Mp + r[2] * lp + r[3] * w for r in bb]
rng = np.random.default_rng(7)
certs = {}
for ub, wb in itertools.product(("v", "za", "avg"), ("H", "dual")):
    c, c0, A, b = branch_data(ub, wb)
    found = set()
    for _ in range(400):
        Mp, lp, w = rng.uniform(0.3, 1.3), rng.uniform(0, 0.7), rng.choice([0, 0.1, 1/7, 0.3, 0.5, 4/7])
        bb = np.array([float(x) for x in evalb(b, F(Mp), F(lp), F(w))])
        res = linprog(-np.array([float(x) for x in c]), A_ub=np.array([[float(x) for x in r] for r in A]), b_ub=bb,
                      bounds=[(0, None)] * n, method="highs")
        if res.status != 0: continue
        yd = -res.ineqlin.marginals          # duals >= 0
        yq = tuple(F(v).limit_denominator(60) for v in yd)
        found.add(yq)
    good = []
    for yq in found:
        if min(yq) < 0: continue
        # exact dual feasibility: A^T y >= c componentwise
        ok = all(sum(yq[i] * A[i][j] for i in range(len(A))) >= c[j] for j in range(n))
        if not ok: continue
        # affine bound: y.b + c0  -> coefficients (const, M', l', w); subtract M'
        coef = [sum(yq[i] * b[i][k] for i in range(len(A))) + c0[k] for k in range(4)]
        coef[1] -= 1
        good.append(tuple(coef))
    certs[(ub, wb)] = sorted(set(good))
targets = {(F(0), F(0), F(0), F(0)): "0", (F(1, 3), F(-2, 3), F(2, 3), F(0)): "(1+2l'-2M')/3",
           (F(-1), F(1), F(1), F(-1)): "M'+l'-1-w", (F(-2, 3), F(1, 3), F(4, 3), F(-1)): "(M'+4l'-2)/3-w"}
for key, lst in certs.items():
    print(key, "certified affine bounds (const, M', l', w):")
    for co in lst:
        print("    ", [str(x) for x in co], "=", targets.get(co, "OTHER"))
