"""Design-space scan of 'modified probes' in the exact 3-node form (s, m, A):
   s : sextic Kubota node whose residue is the row theta rho(m)   [(M_ss, M_sm) = (4,5) paper; (2,1) conjugate]
   m : row variable, no Gauss sum [0]
   A : column variable with arithmetic weight of S-W type M_AA (5 = gamma_2 Patterson theta, the paper; 1 = gamma_4;
       0 = no Gauss sum; 2,4 = sextic Gauss sums)
   m-A edge j : the row twist chi_A(m)^j (paper: j = 1).
Zero detection (the probe's signal): the m-Voronoi (chamber with e_s simple and the m-direction dualised) must turn the
column into a Moebius node (M = 3), i.e. c_A / tau(chi_A^{-j}) = mu(A)*(ray)  (Davenport-Hasse collapse, V_gaps 1.6).
For every variant: groupoid size, whether a Moebius node is met, the residue-adapted representation types, and the best
bound exponent (beyond Q^{5/12}) among representations that do NOT force the signal line, vs CS = M/2."""
import sys
sys.path.insert(0, '.')
import chamber_analysis as CA
import sympy as sp
from groupoid_bfs import bfs
from groupoid_search import Mof

uA, um = CA.uA, CA.um
rows = []
for sor, (Mss, Msm) in [('paper', (4, 5)), ('conj', (2, 1))]:
    for MAA in (5, 1, 0, 2, 4):
        for j in range(1, 6):
            M = [[Mss, Msm, 0], [Msm, 0, j], [0, j, MAA]]
            CA.MODELS['X'] = dict(M=M, names=['s', 'm', 'A'], res={'s': sp.Rational(1, 6)}, free={'m': um, 'A': uA})
            try:
                dist, chi, mij = bfs(M, cap=4000, maxcoef=30)
            except Exception as e:
                continue
            finite = len(dist) < 4000 and max(abs(x) for E in dist for v in E for x in v) < 30
            moeb = any(Mof(E, chi)[k][k] == 3 for E in dist for k in range(3))
            reps = CA.analyse('X', cap=4000, maxcoef=30, verbose=False)
            # signal chamber: P1 forced; detect Moebius at the free A-image
            sigtypes = [d for d in reps.values() if d['P1']]
            sig_moeb = any(3 in d['ftypes'] for d in sigtypes)
            nonforced = [d for d in reps.values() if not d['P1']]
            best = min((min(v for v in (d['triv'], d['ls']) if v is not None) for d in nonforced), default=None)
            short = [d['w'] for d in nonforced if max(d['lnum']) < 0.826 - 1e-9]
            rows.append((sor, MAA, j, len(dist), finite, moeb, sig_moeb, len(reps), best, short))
print("%-6s %4s %2s %6s %6s %8s %10s %6s %10s  %s" % ('s-or', 'M_AA', 'j', '#ch', 'finite', 'Moebius', 'sig=Moeb', '#types', 'best', 'short non-signal images'))
for r in rows:
    print("%-6s %4d %2d %6d %6s %8s %10s %6d %10s  %s" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7],
          'None' if r[8] is None else '%.3f' % r[8], r[9]))
print("CS exponent at M=0.826: 0.413")
