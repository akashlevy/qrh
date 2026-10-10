"""Re-verify key frontier points on the mid/fine bins: margin at 41/47 (kappa = 35/47), exact supply L at the
optimal b, d-range, floor bin and side conditions."""
import json, sys
from levers import *
K = 35/47
def check(name, vals, bins=BINS_MID):
    lev = lev_with(vals)
    m, det = margin_T(lev, TARGET, bins, detail=True)
    b, ell = det['b'], det['ell']
    bh, d2 = B_high(ell, b, lev, K, BINS_FINE, detail=True, exactL=True)
    ok, viol = side_ok(ell, b, lev, kappa=K, report=True)
    dc = dcheck(ell, b, lev, K, BINS_SEARCH, fracs=(1.0, 0.95, 0.8, 0.6))
    dbind = max(dc, key=lambda z: z[1])[0]
    print(f'{name:44s} margin(mid)={m:+.2e}  fine/exact-L margin={bh-TARGET:+.2e}  b={b:.4f} ell={ell:.4f} '
          f'bin=({d2["delta"]:.4f},{d2["x"]:.2f}) floor={d2["floor"]} side_ok={ok} d-binding={dbind}', flush=True)
    return bh - TARGET
cfgs = json.load(open(sys.argv[1])) if len(sys.argv) > 1 else {}
base = {
 'proved (must be > 0)': {},
 'all 7 at f=0.1992': dict(sM=0.0308, eP=0.0120, alpha=0.7835, sLS=0.0055, w=0.0071, g=0.9479),
 '6 (no eta) at f=0.1992': dict(sM=0.0308, eP=0.0120, alpha=0.7835, sLS=0.0055, w=0.0071, g=0.9479),
 'sM+eP+sLS+w at f=0.2568': dict(sM=0.0398, eP=0.0154, sLS=0.0070, w=0.0092),
 'sLS alone 0.02740': dict(sLS=0.0274), 'w alone 0.0359': dict(w=0.0359), 'eP alone 0.0601': dict(eP=0.0601),
 'eL alone 0.0592': dict(eL=0.0592), 'sM alone 0.1549': dict(sM=0.1549), 'eta alone 0.0971': dict(eta=0.0971),
 'alpha alone 0.5832': dict(alpha=0.5832), 'g alone 0.7384': dict(g=0.7384),
}
base.update(cfgs)
if cfgs: base = cfgs
for k, v in base.items():
    check(k, v)
