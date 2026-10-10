"""Validation 1: the 08 certificate (theta = 10499/12000, ell = 0.167, b = 0.086, kappa = 3/4)."""
import numpy as np
from combo_model import *
B08 = 10499/12000
for name, lev in [('08 paper inputs (kappa=3/4, physical supply, no LS)', PAPER08),
                  ('08 + all proved count options (slot powers, sextic LS)', replace(PROVED, kappa=0.75)),
                  ('08 weighted branch only (no E_old)', replace(PAPER08, Eold=False))]:
    for bins, lab in [(BINS_SEARCH, 'delta step .01'), (BINS_FINE, 'delta step .0025, x step .025')]:
        e, det = worstE(0.167, 0.086, B08, lev, 0.75, bins=bins, detail=True, exactL=True)
        print(f'{name:58s} {lab:30s} worst E = {e:+.6f} at delta={det["delta"]:.4f} x={det["x"]:.3f} R={det["R"]:.5f}')
# certificate-style: E_old for delta<=1/3, E_wt for delta>=1/3, margin
