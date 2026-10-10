"""C.1 add-on: which points (a,a,c) lie in the convex hull of the images g(C0), C0 = {a,b,c > 1/2},
under the (heuristic, k = conductor exponent) reflection group of Z2 (see tits_hull.py)?
Lindelof-on-average on the line Re z' = 1/2 + a for the mu-family (via Landau/positivity) needs
(a, a, 1/2 + eps) in the hull.  We report the smallest a reachable at c = 1/2 + eps."""
import numpy as np, sys
from scipy.optimize import linprog
k = int(sys.argv[1]) if len(sys.argv) > 1 else 3
N = int(sys.argv[2]) if len(sys.argv) > 2 else 12
Sa = np.array([[-1,0,0],[0,1,0],[k,0,1]],float); Sb = np.array([[1,0,0],[0,-1,0],[0,k,1]],float)
Sc = np.array([[1,0,1],[0,1,1],[0,0,-1]],float)
gens = [Sa, Sb, Sc]
elems = {tuple(np.eye(3).ravel()): np.eye(3)}; frontier = [np.eye(3)]
for L in range(N):
    nf = []
    for M in frontier:
        for S in gens:
            P = S @ M; key = tuple(np.round(P, 6).ravel())
            if key not in elems: elems[key] = P; nf.append(P)
    frontier = nf
mats = list(elems.values())
apex = np.array([m @ np.array([0.5, 0.5, 0.5]) for m in mats])         # (n,3)
dirs = np.concatenate([np.stack([m @ e for e in np.eye(3)]) for m in mats])  # (3n,3)
def inhull(x):
    n = len(apex); nd = len(dirs)
    Aeq = np.zeros((4, n + nd)); Aeq[:3, :n] = apex.T; Aeq[:3, n:] = dirs.T; Aeq[3, :n] = 1
    beq = np.concatenate([x, [1.0]])
    r = linprog(np.zeros(n + nd), A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * (n + nd), method="highs")
    return r.status == 0
print("k=%d, %d group elements" % (k, len(mats)))
for c in (0.51, 1.0, 2.0):
    lo, hi = -3.0, 0.5
    if not inhull(np.array([hi, hi, c])): print("c=%.2f: even a=0.5 not in hull?" % c); continue
    for _ in range(40):
        mid = (lo + hi) / 2
        if inhull(np.array([mid, mid, c])): hi = mid
        else: lo = mid
    print("c=%.2f (Re u = %.2f): smallest reachable a = b = %.4f  (Re z' = %.4f)" % (c, c + 0.5, hi, hi + 0.5))
