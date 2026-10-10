"""d-range check: E+B at row dyads d < h (Collet: convex in d, endpoint d = h binds)."""
from levers import *
K = 35/47
cfgs = [('proved @ theta*', PROVED, 0.874886, 0.0871), ('g=0.7385 @41/47', lev_with({'g': 0.7385}), TARGET, None),
        ('w=0.0358 @41/47', lev_with({'w': 0.0358}), TARGET, None), ('sLS=0.0274 @41/47', lev_with({'sLS': 0.0274}), TARGET, None),
        ('eta=0.0972 @41/47', lev_with({'eta': 0.0972}), TARGET, None), ('eP=0.06 @41/47', lev_with({'eP': 0.06}), TARGET, None),
        ('alpha=0.5833 @41/47', lev_with({'alpha': 0.5833}), TARGET, None), ('sMnc=0.155 @41/47', lev_with({'sMnc': 0.155}), TARGET, None)]
for name, lev, T, b in cfgs:
    kap = min(0.75, 2*T - 1)
    if b is None:
        m, det = margin_T(lev, T, BINS_SEARCH, detail=True); b = det['b']
        fl = det['floor']
    else:
        fl = None
    ell = ell_T(b, lev, T)
    dc = dcheck(ell, b, lev, kap, BINS_SEARCH)
    print(f'{name:22s} b={b:.4f} ell={ell:.4f} floor-binding={fl} :', ' '.join(f'{f:.2f}h:{v-T:+.4f}' for f, v in dc), flush=True)
