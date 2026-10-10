"""(1) Pareto frontier of minimal joint requirements for theta < 41/47 (kappa = 35/47, all proved options on).
Pairs: lever i fixed at a fraction f_i of its standalone threshold, minimal fraction f_j of lever j.
Subsets: minimal uniform fraction f* (every lever in the subset at p + f (t - p))."""
import json, sys, itertools, time
from levers import *
ETA0 = 0.194
bins = BINS_SEARCH
thr = json.load(open('thresholds41.json'))
P = {k: LEVERS[k][1] for k in LEVERS}; P['c'] = 1/(6*35/47)
P['eta'] = ETA0   # anchor eta at its activation point (inactive above ~0.194; best provable 0.34)
SCALE = {k: thr[k] for k in thr}
SCALE['sM'] = thr['sMnc']           # sM with the cubic branch kept has no standalone threshold: use the cap-free one
if SCALE.get('g') is None: SCALE['g'] = 0.0
def val(k, f): return P[k] + f * (SCALE[k] - P[k])
def margin(fr):  # fr: dict lever -> fraction
    return closes41(lev_with({k: val(k, f) for k, f in fr.items()}), bins)
def min_frac(fixed, k, hi=1.0, it=11):
    if margin({**fixed, k: 0.0}) < 0: return 0.0
    if margin({**fixed, k: hi}) >= 0: return None
    lo = 0.0
    for _ in range(it):
        mid = (lo + hi) / 2
        if margin({**fixed, k: mid}) < 0: hi = mid
        else: lo = mid
    return hi
mode = sys.argv[1]
names = sys.argv[2].split(',') if (len(sys.argv) > 2 and mode != 'pairs') else MAIN
out = {}
if mode == 'pairs':
    pairs = [tuple(p.split(':')) for p in sys.argv[2].split(',')]
    for i, j in pairs:
        for fi in (0.25, 0.5, 0.75):
            t = time.time(); fj = min_frac({i: fi}, j)
            out[f'{i}@{fi}|{j}'] = fj
            print(f'{i:5s} at {fi:.2f} ({val(i,fi):.4f}) -> {j:5s} needs f={fj}  ({"%.4f" % val(j,fj) if fj is not None else "-"})  sum={fi+fj if fj is not None else None}  [{time.time()-t:.0f}s]', flush=True)
    json.dump(out, open('pairs_%s.json' % sys.argv[3], 'w'), indent=1)
elif mode == 'subsets':
    kmin, kmax = int(sys.argv[3]), int(sys.argv[4])
    for k in range(kmin, kmax + 1):
        for S in itertools.combinations(names, k):
            t = time.time()
            f = lambda x: margin({n: x for n in S})
            if f(1.0) >= 0: fs = None
            else:
                lo, hi = 0.0, 1.0
                for _ in range(11):
                    mid = (lo + hi) / 2
                    if f(mid) < 0: hi = mid
                    else: lo = mid
                fs = hi
            out['+'.join(S)] = fs
            vals = ', '.join(f'{n}={val(n, fs):.4f}' for n in S) if fs is not None else ''
            print(f'{"+".join(S):34s} f*={fs if fs is None else round(fs,4)}  [{vals}]  [{time.time()-t:.0f}s]', flush=True)
    json.dump(out, open(f'subsets_{kmin}_{kmax}.json', 'w'), indent=1)
