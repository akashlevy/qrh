"""Predicted mixed-diagonal contribution (n1 p = n2 p') to the coherent cross ratio of orig_side.py.
ratio_pred = sum_{p != p'} R^{-1} sum_{m sf,(m,6pp')=1} W(N(p')N(m)/R) W(N(p)N(m)/R) * rho / (L * m2),
rho = density of rows u coprime to pp'm (~1); compared with the measured coherent_ratio."""
import sys, math, json
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/energy-wall")
from eisenstein import eisenstein_primes
import dual_energy as de
from eisenstein import ResidueField
def W(x): return np.where((x > 1) & (x < 2), np.sin(np.pi*(x-1))**2, 0.0)
for U, P, m2 in [(20000,150,0.0999),(40000,150,0.0998),(80000,150,0.0994),(80000,400,0.0994)]:
    R = int(round(U**0.7))
    primes = eisenstein_primes(2*max(R,P)+10)
    fields = {p[1]: ResidueField(*p) for p in primes}
    pnorm = {p[1]: (p[2] if p[0]=="split" else p[2]**2) for p in primes}
    plist = sorted(fields, key=lambda p: pnorm[p])
    sp = [p for p in plist if fields[p].kind=="split" and P < pnorm[p] <= 2*P]
    mA, mB, mN, fac = de.primary_sf_elements(0, 2*R//P+2, fields, plist, pnorm)
    tot = 0.0
    for p in sp:
        for q in sp:
            if p == q: continue
            ok = np.array([all(pnorm[x] != 4 and x not in (p, q) for x in f) for f in fac], bool)
            Nm = mN[ok].astype(float)
            tot += (W(pnorm[q]*Nm/R)*W(pnorm[p]*Nm/R)).sum()/R
    # m = 1 (unit) term: primary_sf_elements starts at N>0 so includes m=1 (N=1)
    print(json.dumps(dict(U=U, P=P, L=len(sp), ratio_pred=round(tot/(len(sp)*m2), 4))))
