"""
Independent re-implementation (referee) of the reflected-energy maximization of
Lemma reflected-energy (eq. old-eq:4.6) as used in Lemma probe-row-norm (f=1, rho=1),
with the W(w) relaxation max(H, v+lb) -> max(H, v+lb-w).

Objective is a max of 6 affine pieces (max of 2) + (max of 3), so
  max_x Obj(x) = max_{pieces} max_x piece(x)   (no branch constraints needed).

Variables: O H A2 A4 Az N0 B0 S4 za v lb  (all >= 0)
  A2 : non-slot active primes with j in {1,2,3,5} (valuation >= 2)
  A4 : j=4 primes (valuation >= 4), split N0 + B0 + S4 = A4
  Az : active j=0 primes (valuation >= 6)
Constraints:
  H + O <= M'
  2 A2 + 4 A4 + 6 Az <= O          (true valuations)       [option 'paperA0': 2(A2+A4+Az) <= O]
  za <= l'
  Td = 2H + 2(A2+A4+Az) + 2 za - (1+l') - N0 - 3 B0
  v + 3 lb <= Td
Objective (relative to M'):
  count + max(H, v+lb-w) - (S4+Az) - B0 + za - min(v, za, (v+za)/3) - lb - (Td - v - 3lb)/2 - M'
count: paper O/2 ; refined O/6 + 2A2/3 + A4/3
"""
import numpy as np, itertools
from scipy.optimize import linprog
names = "O H A2 A4 Az N0 B0 S4 za v lb".split(); ix = {k: i for i, k in enumerate(names)}; n = len(names)
def vec(**kw):
    r = np.zeros(n)
    for k, c in kw.items(): r[ix[k]] += c
    return r
def energy(Mp, lp, w=0.0, count="paper", a0="true", ret=False):
    Td = vec(H=2, A2=2, A4=2, Az=2, za=2, N0=-1, B0=-3); Td0 = -(1 + lp)
    A_ub, b_ub, A_eq, b_eq = [], [], [], []
    A_ub.append(vec(H=1, O=1)); b_ub.append(Mp)
    if a0 == "true":
        A_ub.append(vec(A2=2, A4=4, Az=6, O=-1)); b_ub.append(0)
    else:
        A_ub.append(vec(A2=2, A4=2, Az=2, O=-1)); b_ub.append(0)
    A_ub.append(vec(za=1)); b_ub.append(lp)
    A_ub.append(vec(v=1, lb=3) - Td); b_ub.append(Td0)
    A_eq.append(vec(N0=1, B0=1, S4=1, A4=-1)); b_eq.append(0)
    cnt = vec(O=0.5) if count == "paper" else vec(O=1/6, A2=2/3, A4=1/3)
    base = cnt + vec(S4=-1, Az=-1, B0=-1, za=1, lb=-1) - 0.5 * (Td - vec(v=1, lb=3))
    base0 = -0.5 * Td0
    best = (-np.inf, None)
    for mx, mn in itertools.product(("H", "dual"), ("v", "za", "avg")):
        c = base.copy(); c0 = base0
        if mx == "H": c += vec(H=1)
        else: c += vec(v=1, lb=1); c0 -= w
        c -= {"v": vec(v=1), "za": vec(za=1), "avg": vec(v=1/3, za=1/3)}[mn]
        res = linprog(-c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                      bounds=[(0, None)] * n, method="highs")
        if res.status == 0:
            val = -res.fun + c0
            if val > best[0]: best = (val, (mx, mn, dict(zip(names, np.round(res.x, 5)))))
        elif res.status == 3:
            return (np.inf, None) if ret else np.inf
    if best[0] == -np.inf:
        return (None, None) if ret else None   # infeasible: no retained dual dyad
    return (best[0] - Mp, best[1]) if ret else best[0] - Mp

def closed(Mp, lp, w, count):
    a = (1 + 3*lp - 2*Mp)/4 if count == "paper" else (1 + 2*lp - 2*Mp)/3
    return max(0.0, a, Mp + lp - 1 - w, (Mp + 4*lp - 2)/3 - w)

if __name__ == "__main__":
    rng = np.random.default_rng(2024)
    for count in ("paper", "refined"):
        for a0 in ("true", "paperA0"):
            worst = 0; arg = None; ninf = 0; nover = 0; nunder = 0
            for _ in range(4000):
                Mp = rng.uniform(0.2, 1.4); lp = rng.uniform(0, 0.7); w = float(rng.choice([0, 0.05, 0.1, 1/7, 0.2, 2/7, 0.3, 0.5, 4/7]))
                e = energy(Mp, lp, w, count, a0)
                if e is None: ninf += 1; continue
                c = closed(Mp, lp, w, count)
                if e - c > 1e-8: nover += 1
                if c - e > 1e-8: nunder += 1
                if abs(e - c) > worst: worst = abs(e - c); arg = (round(Mp,4), round(lp,4), round(w,4), e, c)
            print(f"count={count:8s} A0-constraint={a0:8s} infeasible={ninf} LP>closed={nover} LP<closed={nunder} max|diff|={worst:.3g} at {arg}")
