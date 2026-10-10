import sys, numpy as np
from sk_read import read
from hls_lib import Family
D = read(sys.argv[1])
F = Family(2000, 10, 7, nP=1, verbose=False)
# map (A,B) of columns -> compare gamma_m at columns: C file doesn't store A,B; compare via sorted multiset per norm
# instead compare prime Gauss sums directly
from hls_lib import gauss_all
bad = 0
for (N, a, b), g in zip(D["primes"][["N", "a", "b"]][:300], D["pg"][:300]):
    if N == 4: continue
    p = (int(a), int(b))
    if p not in F.fields: continue
    ga = gauss_all(F.fields[p])[1:6]
    bad += np.max(np.abs(ga - g)) > 1e-8
print("prime gauss mismatches:", bad)
# composite check: multiset of (N, gamma_1..5) for N in (2000, 4000]
sel = (D["N"] > 2000) & (D["N"] <= 4000)
gC = D["g"][sel]; NC = D["N"][sel]
gP = F.gam[:, 1:6]; NP_ = F.Nn
print("counts", sel.sum(), len(NP_))
key = lambda N, g: tuple([int(N)] + [round(x, 6) for x in np.concatenate([g.real, g.imag])])
s1 = sorted(key(n, g) for n, g in zip(NC, gC)); s2 = sorted(key(n, g) for n, g in zip(NP_, gP))
print("composite multiset equal:", s1 == s2)
