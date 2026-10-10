import json
from levers import *
thr = json.load(open('thresholds41.json'))
P = {'sM': 0, 'eP': 0, 'eL': 0, 'sLS': 0, 'w': 0, 'g': 1.0, 'alpha': 5/6}
out = {}
for S in [('sM','eP','eL','sLS','w'), ('sM','eP','eL','sLS','w','g'), ('eP','eL','sLS','w'), ('sM','eP','eL','sLS','w','g','alpha')]:
    f = lambda x: closes41(lev_with({n: P[n] + x*(thr[n]-P[n]) for n in S}), BINS_SEARCH)
    fs = bisect_lever(f, 0.0, 1.0, it=12)
    out['+'.join(S)] = fs
    print('+'.join(S), 'f* =', round(fs, 4), {n: round(P[n] + fs*(thr[n]-P[n]), 4) for n in S}, flush=True)
json.dump(out, open('subsets_eL.json', 'w'), indent=1)
