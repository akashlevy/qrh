"""Is the (relaxed) wall M + ell = 1 + w optimal on the line B_low = T?  Scan M below the wall."""
import numpy as np
from levers import *
T = TARGET; K = 35/47
for name, lev in [('proved', PROVED), ('g=0.74', lev_with({'g': 0.74})), ('w=0.036', lev_with({'w': 0.036})),
                  ('sLS=0.0275', lev_with({'sLS': 0.0275}))]:
    rows = []
    for dm in [0.0, -0.01, -0.02, -0.04]:
        best = (9, None)
        for b in np.linspace(0, 0.2, 41):
            # B_low = 5/6 + M/12 - ell/6 - b(1-g)/12 = T with M + ell = 1 + w + dm
            # => M = 1 + w + dm - ell ;  5/6 + (1+w+dm)/12 - ell/4 - b(1-g)/12 = T
            ell = 4 * (5/6 + (1 + lev.w + dm)/12 - b*(1-lev.g)/12 - T)
            M = 1 + lev.w + dm - ell
            if not side_ok(ell, b, lev, M=M, kappa=K): continue
            v = B_high(ell, b, lev, K, BINS_SEARCH, M=M) - T
            best = min(best, (v, b))
        rows.append(f'dM={dm:+.2f}: {best[0]:+.5f} (b={best[1]})')
    print(f'{name:12s}', ' | '.join(rows), flush=True)
