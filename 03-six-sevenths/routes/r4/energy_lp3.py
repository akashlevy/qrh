"""
Refined energy LP: separate local types of the frozen powerful row part.
  A2 : active primes with j in {1,2,3,5} (row valuation >= 2; cheapest O-cost 2)   -> energy coefficient 2 (via 2A0)
  A4 : j=4 primes (valuation >= 4, O-cost 4), split N0 (to n), B0 (to b), S4 (small Ramanujan)
  Az : active j=0 primes (valuation >= 6, O-cost 6), counted in S0
count of frozen parts of a given type (O, A2, A4, Az):
  paper  : O/2
  refined: O/6 + 2*A2/3 + A4/3     [sixth-power part O0/6 + radical of the rest; O0 = O - (2A2+4A4)]
"""
import numpy as np
from scipy.optimize import linprog
names = "O H A2 A4 Az N0 B0 S4 za v lb".split()
ix = {k: i for i, k in enumerate(names)}; n = len(names)
def vec(**kw):
    r = np.zeros(n)
    for k, c in kw.items(): r[ix[k]] += c
    return r
def emax(Mp, lp, w=0.0, count="refined", ret_arg=False, Ns=None):
    Ns = 1 + lp if Ns is None else Ns
    Td = vec(H=2, A2=2, A4=2, Az=2, za=2, N0=-1, B0=-3); Td0 = -Ns
    y = vec(v=1, lb=3)
    best = (-99, None)
    for ub in ("v", "za", "avg"):
        for wb in ("H", "dual"):
            c = np.zeros(n); const = 0.0
            c += vec(O=0.5) if count == "paper" else vec(O=1/6, A2=2/3, A4=1/3)
            if wb == "H": c += vec(H=1)
            else: c += vec(v=1, lb=1); const -= w
            c += vec(S4=-1, Az=-1, B0=-1, za=1, lb=-1)
            uv = {"v": vec(v=1), "za": vec(za=1), "avg": vec(v=1/3, za=1/3)}
            c -= uv[ub]; c -= 0.5 * (Td - y); const -= 0.5 * Td0
            A, b = [], []
            def le(row, rhs): A.append(row); b.append(rhs)
            le(vec(H=1, O=1), Mp)
            le(vec(A2=2, A4=4, Az=6, O=-1), 0)
            le(vec(N0=1, B0=1, S4=1, A4=-1), 0); le(vec(A4=1, N0=-1, B0=-1, S4=-1), 0)
            le(vec(za=1), lp)
            le(y - Td, Td0)
            for o in uv:
                if o != ub: le(uv[ub] - uv[o], 0)
            if wb == "H": le(vec(v=1, lb=1) - vec(H=1), w)
            else: le(vec(H=1) - vec(v=1, lb=1), -w)
            res = linprog(-c, A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, None)] * n, method="highs")
            if res.status == 0:
                val = -res.fun + const
                if val > best[0]: best = (val, (ub, wb, {k: round(float(x), 5) for k, x in zip(names, res.x) if abs(x) > 1e-9}))
    return (best[0] - Mp, best[1]) if ret_arg else best[0] - Mp

def closed(Mp, lp, w, count):
    a = (1 + 3*lp - 2*Mp)/4 if count == "paper" else (1 + 2*lp - 2*Mp)/3
    return max(0.0, a, Mp + lp - 1 - w)

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    worst = {("paper"): 0, ("refined"): 0}
    for count in ("paper", "refined"):
        md = 0
        for _ in range(3000):
            Mp = rng.uniform(0.3, 1.2); lp = rng.uniform(0, 0.6); w = rng.choice([0, 0.05, 0.1, 1/7, 0.2, 0.3])
            e = emax(Mp, lp, w, count); cf = closed(Mp, lp, w, count)
            if abs(e - cf) > md:
                md = abs(e - cf); arg = (Mp, lp, w, e, cf)
        print(count, "max |LP - closed form| =", md, "at", arg)
