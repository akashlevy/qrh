"""check gamma_2 (edge_lib convention) = (lambda/n)_3 * G~(n) (cubic_data convention), lambda = 1 + 2w."""
import sys, math
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
from edge_lib import setup, gauss_vec, Columns
from nlib import load_cd, Prime, W
P, S = load_cd(410000)
fields, pnorm, plist = setup(6000)
prs = {(int(r["a"]), int(r["b"])): Prime(r) for r in P}
bad = 0; worst = 0; tested = 0
for p in plist:
    if pnorm[p] in (3, 4) or pnorm[p] > 6000: continue
    pr = prs[p]
    g2 = gauss_vec(fields[p], 2)
    k = pr.sextic_exp(1, 2)[0]                         # (lambda/pi)_6
    lam3 = np.exp(2j * np.pi * 2 * int(k) / 6)          # (lambda/pi)_3 = (lambda/pi)_6^2
    d = abs(g2 - lam3 * pr.g2t); worst = max(worst, d); bad += d > 1e-9; tested += 1
print(f"primes: {tested} tested (split+inert, N <= 6000), max |gamma_2 - (lambda/pi)_3 G~| = {worst:.2e}, failures {bad}")
# composites: compare beta from edge_lib Columns with conversion of cubic_data sf records
C = Columns(20000, 22000, fields, pnorm, plist, None) if False else None
fields2, pnorm2, plist2 = setup(44000)
C = Columns(20000, 22000, fields2, pnorm2, plist2, None)
key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(S["A"], S["B"]))}
lamk = {}
worst = 0
for i in range(C.n):
    j = key[(int(C.A[i]), int(C.B[i]))]
    g = complex(S["gre"][j], S["gim"][j])
    lam = 1 + 0j
    for p in C.fac[i]:
        k = prs[p].sextic_exp(1, 2)[0]; lam *= np.exp(2j * np.pi * 2 * int(k) / 6)
    nc = C.A[i] + C.B[i] * W
    b2 = lam * g * np.conj(nc) / abs(nc)
    worst = max(worst, abs(b2 - C.beta[i]))
    assert S["mu"][j] == (-1) ** len(C.fac[i])
print(f"composites: {C.n} columns in (20000, 22000]: max |beta_edge - (lambda/n)_3 G~ conj(alpha)| = {worst:.2e}")
