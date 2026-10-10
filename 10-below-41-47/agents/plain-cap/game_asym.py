"""Thresholds of the move game for asymmetric plain lengths (n1<=n2) at M=1, kappa=3/4.
Conjecture from the (S,D) analysis: capacity = (M - 2*n_min)/(6 kappa) (paper: (M-n1-n2)/(6 kappa))."""
from game import solve
kappa=0.75
def thr(n1,n2,h=0.01):
    m=(n1+n2)/2; D0=n2-n1
    lo=0.0
    for i in range(1,40):
        z=i*h
        if solve(m,z,kappa=kappa,h=h,D0=D0)>1e-9: return (i-1)*h, i*h
    return None
for (n1,n2) in [(0.4,0.4),(0.3,0.4),(0.2,0.4),(0.0,0.4),(0.3,0.5),(0.0,0.6)]:
    print(f"(n1,n2)=({n1},{n2}): feasible up to z in {thr(n1,n2)};  (1-2nmin)/(6k)={(1-2*n1)/(6*kappa):.4f}  paper (1-n1-n2)/(6k)={(1-n1-n2)/(6*kappa):.4f}",flush=True)
