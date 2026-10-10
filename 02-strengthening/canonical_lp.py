"""
Terminal step of Lemma canonical-moment: is E_ref <= F for every admissible block
when the first invariant is relaxed to F-M-Q-z0 >= -w and case B of
Lemma reflected-energy is relaxed by w (soft: max(H, v+lb-w))?
Variables: F, M, Q, z0, H, O, A0, N0, B0, S0, za, v, lb, el   (V=0, N*=F)
"""
import sys, numpy as np
from scipy.optimize import linprog
names = "F M Q z0 H O A0 N0 B0 S0 za v lb el".split()
ix = {n: i for i, n in enumerate(names)}; n = len(names)
def e(*pairs):
    r = np.zeros(n)
    for k, c in pairs: r[ix[k]] += c
    return r
def run(w, second_margin=0.0, fixF=None, extra=None, verbose=False):
    Td = e(("H", 2), ("A0", 2), ("za", 2), ("F", -1), ("N0", -1), ("B0", -3))
    y = e(("v", 1), ("lb", 3), ("el", 1))
    best = (-9, None)
    for ub in ["v", "za", "avg"]:
        for wb in ["H", "dual"]:
            c = e(("O", .5), ("S0", -1), ("B0", -1), ("za", 1), ("lb", -1), ("el", -2/3)) - (Td - y) / 2 - e(("F", 1))
            const = 0.0
            if wb == "H": c += e(("H", 1))
            else: c += e(("v", 1), ("lb", 1)); const -= w
            uv = {"v": e(("v", 1)), "za": e(("za", 1)), "avg": e(("v", 1/3), ("za", 1/3))}
            c -= uv[ub]
            A, b = [], []
            def le(row, rhs): A.append(row); b.append(rhs)
            le(e(("H", 1), ("O", 1), ("M", -1)), 0)                       # H <= M - O
            le(e(("H", 2), ("A0", 2), ("N0", -1), ("S0", -1), ("B0", -4), ("M", -2), ("O", 1), ("Q", -1)), 0)  # norm ineq (V=0)
            le(e(("N0", 1), ("A0", -1)), 0); le(e(("S0", 1), ("A0", -1)), 0); le(e(("B0", 1), ("A0", -1)), 0)
            le(y - Td, 0)                                                  # retained dyad
            le(e(("za", 1), ("z0", -1)), 0)
            le(e(("M", 1), ("Q", 1), ("z0", 1), ("F", -1)), w)             # relaxed first invariant
            le(e(("M", 3), ("z0", 6), ("F", -4)), -second_margin)          # second invariant
            for o in uv:
                if o != ub: le(uv[ub] - uv[o], 0)
            dual = e(("v", 1), ("lb", 1)); Hh = e(("H", 1))
            if wb == "H": le(dual - Hh, w)
            else: le(Hh - dual, -w)
            if fixF is not None:
                for k, val in fixF.items():
                    le(e((k, 1)), val); le(e((k, -1)), -val)
            else:
                le(e(("F", 1)), 3.0)
            res = linprog(-c, A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, None)] * n, method="highs")
            if res.status == 0:
                val = -res.fun + const
                if val > best[0]: best = (val, (ub, wb, dict(zip(names, np.round(res.x, 4)))))
    return best
if __name__ == "__main__":
    for w in [0.0, 0.1, 0.2, 0.3]:
        val, arg = run(w)
        print(f"w={w}: max (E_ref - F) over all nodes with relaxed invariants = {val:.5f}")
        print("   at", arg)

def run2(w, cons):
    """same LP but with extra linear constraints cons=[(row,rhs)] ; returns max E_ref-F"""
    import numpy as np
    from scipy.optimize import linprog
    Td = e(("H", 2), ("A0", 2), ("za", 2), ("F", -1), ("N0", -1), ("B0", -3))
    y = e(("v", 1), ("lb", 3), ("el", 1))
    best = (-9, None)
    for ub in ["v", "za", "avg"]:
        for wb in ["H", "dual"]:
            c = e(("O", .5), ("S0", -1), ("B0", -1), ("za", 1), ("lb", -1), ("el", -2/3)) - (Td - y) / 2 - e(("F", 1))
            const = 0.0
            if wb == "H": c += e(("H", 1))
            else: c += e(("v", 1), ("lb", 1)); const -= w
            uv = {"v": e(("v", 1)), "za": e(("za", 1)), "avg": e(("v", 1/3), ("za", 1/3))}
            c -= uv[ub]
            A, b = [], []
            def le(row, rhs): A.append(row); b.append(rhs)
            le(e(("H", 1), ("O", 1), ("M", -1)), 0)
            le(e(("H", 2), ("A0", 2), ("N0", -1), ("S0", -1), ("B0", -4), ("M", -2), ("O", 1), ("Q", -1)), 0)
            le(e(("N0", 1), ("A0", -1)), 0); le(e(("S0", 1), ("A0", -1)), 0); le(e(("B0", 1), ("A0", -1)), 0)
            le(y - Td, 0); le(e(("za", 1), ("z0", -1)), 0)
            le(e(("F", 1)), 3.0)
            for row, rhs in cons: le(row, rhs)
            for o in uv:
                if o != ub: le(uv[ub] - uv[o], 0)
            dual = e(("v", 1), ("lb", 1)); Hh = e(("H", 1))
            if wb == "H": le(dual - Hh, w)
            else: le(Hh - dual, -w)
            res = linprog(-c, A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, None)] * n, method="highs")
            if res.status == 0:
                val = -res.fun + const
                if val > best[0]: best = (val, (ub, wb, {k: round(float(x), 4) for k, x in zip(names, res.x) if abs(x) > 1e-9}))
    return best

if __name__ == "__main__":
    print("--- with extra genuine-diagonal constraint M + Q <= F (i.e. r+z<=1 at Q=0) ---")
    for w in [0.1, 0.2, 0.3]:
        cons = [(e(("M", 1), ("Q", 1), ("z0", 1), ("F", -1)), w), (e(("M", 3), ("z0", 6), ("F", -4)), 0.0),
                (e(("M", 1), ("Q", 1), ("F", -1)), 0.0)]
        print(w, run2(w, cons))
    print("--- paper invariants (w=0 relaxation only in case B) ---")
    for w in [0.1, 0.2]:
        cons = [(e(("M", 1), ("Q", 1), ("z0", 1), ("F", -1)), 0.0), (e(("M", 3), ("z0", 6), ("F", -4)), 0.0)]
        print(w, run2(w, cons))

if __name__ == "__main__":
    print("--- candidate invariants: F-M-Q>=0, F-M-Q-z0>=-w*k, 4F-3M-6z0>=0 ---")
    for w in [0.1, 0.2]:
        for k in [0.5, 0.6, 1.0]:
            cons = [(e(("M", 1), ("Q", 1), ("z0", 1), ("F", -1)), w * k), (e(("M", 3), ("z0", 6), ("F", -4)), 0.0),
                    (e(("M", 1), ("Q", 1), ("F", -1)), 0.0)]
            v, a = run2(w, cons)
            print(f"w={w} k={k}: max E-F = {v:.5f}  {a}")
