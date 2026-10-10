"""Uniform power saving eps (natural U/Z-exponent units) shared by a subset of the 'power-saving' levers:
sM (cross-mark saving in lem:marked), eP (lem:plain edge, plain witness + numerator), eL (amplified inverse
moment beyond the diagonal), sLS (restricted large sieve on bad rows), w (lossless energy wall).
Minimal eps with theta < 41/47.  Also the linear 'theta-rate' of each lever at the proved point."""
import json, itertools, sys, time
from levers import *
bins = BINS_SEARCH
PS = ['sM', 'eP', 'eL', 'sLS', 'w']
out = {}
th0 = theta_star(PROVED, bins)
if 'rates' in sys.argv:
    for k in PS:
        th = theta_star(lev_with({k: 0.005}), bins)
        out['rate_' + k] = (th0 - th) / 0.005
        print(f'theta rate of {k}: {out["rate_"+k]:.4f} per unit (gain at 0.005: {th0-th:.2e})', flush=True)
kmin, kmax = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 and sys.argv[1].isdigit() else (1, 5)
for k in range(kmin, kmax + 1):
    for S in itertools.combinations(PS, k):
        t = time.time()
        f = lambda e: closes41(lev_with({n: e for n in S}), bins)
        e = bisect_lever(f, 0.0, 0.05 if 'w' in S else 0.2, it=14)
        out['+'.join(S)] = e
        print(f'{"+".join(S):22s} eps* = {e}   [{time.time()-t:.0f}s]', flush=True)
json.dump(out, open(f'epsilon_{kmin}_{kmax}.json', 'w'), indent=1)
