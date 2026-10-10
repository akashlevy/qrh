"""(3) theta* with PROVED inputs and with each single small hypothetical improvement; standalone 41/47 thresholds."""
import json, time, sys
from levers import *
bins = BINS_MID if 'mid' in sys.argv else BINS_SEARCH
res = {}
th0 = theta_star(PROVED, bins); res['proved'] = th0
print(f'PROVED theta* = {th0:.7f}', flush=True)
small = {
    'sM': [0.001, 0.005, 0.01, 0.02, 0.05], 'sMnc': [0.0, 0.05], 'mk': [0.55, 0.6], 'eP': [0.001, 0.0025, 0.005, 0.01, 0.02],
    'c': [0.25, 0.3], 'alpha': [5/6 - 0.001, 5/6 - 0.01, 5/6 - 0.02, 5/6 - 0.05], 'eta': [0.2, 0.19, 0.18, 0.15, 0.1],
    'sLS': [0.001, 0.005, 0.01, 0.02], 'w': [0.001, 0.005, 0.01, 0.02], 'g': [0.99, 0.95, 0.9],
}
res['singles'] = {}
for k, vs in small.items():
    for v in vs:
        lev = lev_with({k: v}); t = time.time()
        th = theta_star(lev, bins)
        res['singles'][f'{k}={v:.4f}'] = th
        print(f'{k:6s}={v:.4f}: theta* = {th:.7f}  (gain {th0 - th:.2e})  [{time.time()-t:.0f}s]', flush=True)
json.dump(res, open('singles%s.json' % ('_mid' if 'mid' in sys.argv else ''), 'w'), indent=1)
