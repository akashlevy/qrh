"""C.1/C.3: Sawin-Whitehead / Heckenberger Weyl-groupoid finiteness for the moment series of Attack C
(uses ../wf2/weyl_groupoid.py; n = 6, node exponents M_ii and edge exponents M_ij as in wf2/A2_theory.md Sec. 4.1).
'finite=False' means the root coefficients exceeded maxcoef (unbounded growth), i.e. the groupoid is infinite
(finite arithmetic root systems of rank <= 5 have small bounded root coefficients)."""
import sys
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/wf2")
from weyl_groupoid import run, objects
cases = [
 ("first moment, reflected (s,m,mu), mu = cubic theta [M=5]", [[4,5,0],[5,4,3],[0,3,5]]),
 ("first moment, unreflected (s,m,A), A = gamma2*abar [M=5]", [[4,5,0],[5,0,1],[0,1,5]]),
 ("first moment, faithful 4-node (t,A,m,s)", [[5,4,0,0],[4,0,1,0],[0,1,0,5],[0,0,5,4]]),
 ("first moment, mu UNWEIGHTED (s,m,mu), mu Dirichlet [M=0]", [[4,5,0],[5,4,3],[0,3,0]]),
 ("first moment, conjugate column orientation (s,m,A) M_AA=1 (control)", [[4,5,0],[5,0,1],[0,1,1]]),
 ("second moment (s1,m1,mu,m2,s2), true pi", [[4,5,0,0,0],[5,4,3,0,0],[0,3,0,3,0],[0,0,3,2,1],[0,0,0,1,2]]),
 ("second moment reduced (m1,mu,m2), true pi (m-nodes gamma1*gamma3 [4] / conj [2])", [[4,3,0],[3,0,3],[0,3,2]]),
 ("second moment, zeta-free model (m-nodes plain [0])", [[0,3,0],[3,0,3],[0,3,0]]),
]
for name, M in cases:
    fin, nb, no, na = run(M, 6, cap=20000, maxcoef=60)
    print("%-80s finite=%-5s bases=%-6d objects=%-3d nonadmissible_seen=%s" % (name, fin, nb, no, na))
