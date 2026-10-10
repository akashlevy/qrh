"""Validation 2: Collet's PR #9 boundary b0 = 11/12 - ell0/4 at her exact geometry and kappa0 = 2 b0 - 1."""
import numpy as np
from scipy.optimize import minimize_scalar
from combo_model import *
ell0 = 0.16712007850681; b0 = 1 - (1/12 + ell0/4); k0 = (5 - 3*ell0)/6
h0 = 0.793947597388136; x0 = 0.373172481118674; y0 = 0.459707440374516
b = y0 - x0
lev = replace(PROVED, Eold=False)
M, lx, ly, h = geom(ell0, b, lev)
print(f'b0={b0:.10f} kappa0={k0:.6f} our geometry lx={lx:.9f} ly={ly:.9f} h={h:.9f} (Collet {x0:.9f} {y0:.9f} {h0:.9f})')
def EB(dl, x=0.5):
    e, _ = high_terms(ell0, b, lev, k0, bins=(np.array([dl]), np.array([x])), R=count_R(np.array([dl]), np.array([x]), ell0/h, lev, k0))
    return float(e[0])
r = minimize_scalar(lambda d: -EB(d), bounds=(0.3, 0.6), method='bounded', options=dict(xatol=1e-10))
print(f'max_delta B_high at x=1/2: {-r.fun:.10f} at delta={r.x:.6f}  (Collet: b0={b0:.10f} at delta0=0.469344858698)')
for name, L2 in [('weighted only', lev), ('with E_old', PROVED)]:
    bh, det = B_high(ell0, b, L2, k0, BINS_FINE, detail=True, exactL=True)
    print(f'{name}: fine-grid B_high={bh:.10f} (delta={det["delta"]}, x={det["x"]}), B_low={B_low(ell0,b,L2):.10f}')
