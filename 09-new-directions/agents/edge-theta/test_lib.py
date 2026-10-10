import numpy as np, time
from edge_lib import *
t0=time.time()
fields, pnorm, plist = setup(20000)
rng = np.random.default_rng(0)
bad = 0; tot = 0
for p in plist[:400]:
    F = fields[p]; S = Symb(F)
    A = rng.integers(-500, 500, 50); B = rng.integers(-500, 500, 50)
    e = S(A, B)
    for a, b, k in zip(A, B, e):
        kk = F.sextic_exp((int(a), int(b)))
        kk = -1 if kk is None else kk
        bad += (kk != k); tot += 1
print("symbol LUT mismatches", bad, "/", tot, time.time()-t0)
cols = Columns(2000, 4000, fields, pnorm, plist, None)
print("cols", cols.n, "identity failures at primes:", check_identity(cols, fields, 40))
# check a(n) = conj(gamma_1 beta) at composite n via direct gamma_1? use TM-free check: |a|=1
print("|G| range", np.abs(cols.G).min(), np.abs(cols.G).max())
