"""(d) upper bound for any coupling of probes that share rows: per-bin choice of geometry
(b, and ell>=25/141 with B=41/47 fixed), i.e. max_bin min_geometry E.  If this oracle fails,
no linear/joint combination of these probes can succeed in the weighted model."""
import numpy as np
from rmodel import pieces, des_paper
B=41/47
geoms=[(ell,b) for ell in (0.1773,0.18,0.185,0.19,0.2) for b in np.linspace(0,0.2,21)]
worst=-9;arg=None
for dl in np.linspace(0.05,0.75,29):
    for x in (0,0.125,0.25,0.375,0.5):
        e=min((des_paper(pieces(ell,b,B,dl,x)),ell,b) for ell,b in geoms)
        if e[0]>worst: worst=e[0];arg=(round(dl,3),x,e[1],round(e[2],3))
print(f"per-bin oracle over (ell,b): worst={worst:+.5f} at (delta,x,ell*,b*)={arg}")
# same with b only at ell=0.1773
worst=-9
for dl in np.linspace(0.05,0.75,29):
    for x in (0,0.125,0.25,0.375,0.5):
        e=min(des_paper(pieces(0.1773,b,B,dl,x)) for b in np.linspace(0,0.2,41))
        worst=max(worst,e)
print(f"per-bin oracle over b at ell=0.1773: worst={worst:+.5f}")
