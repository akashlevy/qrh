import sys, numpy as np, cmath, math
sys.path.insert(0, '.')
from theta_slot_test import split_primes, symbols, gauss_sum
n = int(sys.argv[1]); P = int(sys.argv[2]); filt = int(sys.argv[3])
prs = [t for t in split_primes(P, n) if t[0] % filt == 1 % filt][:40]
vals = []
for i,(p1,a1,b1) in enumerate(prs):
    for j,(p2,a2,b2) in enumerate(prs):
        if p1 == p2: continue
        s12 = symbols(np.array([a1]), np.array([b1]), p2, a2, b2, n)[0]   # (pi1/pi2)
        s21 = symbols(np.array([a2]), np.array([b2]), p1, a1, b1, n)[0]   # (pi2/pi1)
        vals.append(np.round(s12*np.conj(s21), 6))
u, c = np.unique(np.array(vals), return_counts=True)
print(list(zip(u, c)))
# classes of pi mod 4 (a mod 4, b mod 4)
print(sorted(set(((a % 4), (b % 4)) for (_, a, b) in prs)))
