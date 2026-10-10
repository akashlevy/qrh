"""Lever definitions on top of the PROVED model (kappa dynamic, slot powers, sextic LS, E_old)."""
from combo_model import *

# name: (builder(value) -> Lev, proved value, description, units)
LEVERS = {
    'sM':    (lambda v: replace(PROVED, sM=v), 0.0,
              'cross-mark saving in lem:marked: excess e=max(0,r+2z-1,(2r+8z-3)/4) lowered by sM', 'U-exp'),
    'sMnc':  (lambda v: replace(PROVED, sM=v, cap82=False), 0.0,
              'as sM, cubic branch 2r+8z<=3 dropped', 'U-exp'),
    'mk':    (lambda v: replace(PROVED, mk=v, cap82=False), 0.5,
              'marked capacity mk(1-r), cubic branch dropped, slot powers', 'rate'),
    'eP':    (lambda v: replace(PROVED, eP=v), 0.0,
              'power saving at lem:plain edge (plain witness AND numerator): capacity (1-2n+6eP)/(6kappa)', 'U-exp'),
    'c':     (lambda v: replace(PROVED, cP=v), None,
              'lem:plain capacity rate c = cP = cN (proved 1/(6kappa))', 'rate'),
    'alpha': (lambda v: replace(PROVED, alpha=v), 5 / 6,
              'long-witness slope (lem:inverse-amplification)', 'slope'),
    'eL':    (lambda v: replace(PROVED, eL=v), 0.0,
              'power saving in the amplified inverse moment beyond the diagonal: U^{max(1,(1+5r)/6 - eL)}', 'U-exp'),
    'eta':   (lambda v: replace(PROVED, eta=v), 0.34,
              'numerator sixth moment excess (best provable ~0.34)', 'U-exp'),
    'sLS':   (lambda v: replace(PROVED, sLS=v), 0.0,
              'restricted large sieve on bad rows: sum_C|S|^2 << U^{(1+R)/2 - q zN - sLS}', 'U-exp'),
    'w':     (lambda v: replace(PROVED, w=v), 0.0,
              'lossless energy wall M+ell <= 1+w (also r+2z <= 1+w/2 in lem:marked)', 'Z-exp'),
    'g':     (lambda v: replace(PROVED, g=v), 1.0,
              'Gram factor P_a^{g/6} (Gram domination enforced)', 'factor'),
}
MAIN = ['sM', 'eP', 'alpha', 'eta', 'sLS', 'w', 'g']


def lev_with(vals):
    """vals: dict name -> value; combine levers on PROVED."""
    lev = PROVED
    for k, v in vals.items():
        if k == 'sM': lev = replace(lev, sM=v)
        elif k == 'sMnc': lev = replace(lev, sM=v, cap82=False)
        elif k == 'mk': lev = replace(lev, mk=v, cap82=False)
        elif k == 'eP': lev = replace(lev, eP=v)
        elif k == 'c': lev = replace(lev, cP=v)
        elif k == 'alpha': lev = replace(lev, alpha=v)
        elif k == 'eL': lev = replace(lev, eL=v)
        elif k == 'eta': lev = replace(lev, eta=(None if v >= 0.34 else v))
        elif k == 'sLS': lev = replace(lev, sLS=v)
        elif k == 'w': lev = replace(lev, w=v)
        elif k == 'g': lev = replace(lev, g=v)
        else: raise KeyError(k)
    return lev
