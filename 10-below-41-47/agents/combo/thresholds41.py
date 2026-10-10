"""standalone thresholds for theta < 41/47 (kappa = 35/47), all proved options on."""
import json, time
from levers import *
bins = BINS_SEARCH
ranges = {'sM': (0, 0.6), 'sMnc': (0, 0.6), 'mk': (0.5, 1.0), 'eP': (0, 0.2), 'c': (1/(6*35/47), 1.0),
          'alpha': (0.3, 5/6), 'eta': (0.0, 0.34), 'sLS': (0, 0.3), 'w': (0, 0.05), 'g': (0.0, 1.0)}
out = {}
print('proved margin at 41/47:', closes41(PROVED, bins), flush=True)
for k, (lo, hi) in ranges.items():
    t = time.time()
    f = lambda v: closes41(lev_with({k: v}), bins)
    thr = bisect_lever(f, lo, hi, it=16)
    out[k] = thr
    print(f'{k:6s} threshold for 41/47: {thr}   [{time.time()-t:.0f}s]', flush=True)
json.dump(out, open('thresholds41.json', 'w'), indent=1)
