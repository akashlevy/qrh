"""Validation 4: standalone thresholds of the agents, recomputed (ell = 0.1773 wall, B = 41/47, kappa = 3/4)."""
import numpy as np
from combo_model import *
ELL = (1 - 12*(TARGET - 5/6))/3; K = 0.75; A = PAPER08; bins = BINS_SEARCH
def W(lev, ell=ELL, kappa=K):
    return best_over_b(ell, TARGET, lev, kappa, bins)[0]
def Wwall(w, lev):
    ell = (1 + w - 12*(TARGET - 5/6))/3
    return best_over_b(ell, TARGET, replace(lev, w=w), K, bins)[0]
out = []
def rep(name, v, ref):
    s = f'{name:66s} {("%.4f" % v) if v is not None else "none"}   (agent: {ref})'; print(s, flush=True)
rep('mk* unlimited supply, cap82 off', bisect_lever(lambda m: W(replace(A, mk=m, cap82=False, mark_supply='unlimited')), 0.5, 1.0), 'amplitude 0.744')
rep('joint c* = cP = cN, mk 1/2', bisect_lever(lambda c: W(replace(A, cP=c)), 0.2222, 1.0), 'plain-cap 0.73')
rep('joint c* with mk .75 cap off (physical)', bisect_lever(lambda c: W(replace(A, cP=c, mk=.75, cap82=False)), 0.2222, 1.0), 'plain-cap 0.37')
rep('sixth-moment eta* (min over b)', bisect_lever(lambda e: -W(replace(A, eta=e)), 0.0, 0.3), 'numerator 0.0885 (b=0)')
rep('wall w* alone (no marked feed)', bisect_lever(lambda w: Wwall(w, replace(A, wall_marked=False)), 0.0, 0.08), 'marked-cap 0.041-0.042')
rep('wall w* with lem:marked feed w/2, cap82 off', bisect_lever(lambda w: Wwall(w, replace(A, cap82=False)), 0.0, 0.08), 'marked-cap 0.037-0.038')
rep('wall w* with lem:marked feed w/2, cap82 kept', bisect_lever(lambda w: Wwall(w, A), 0.0, 0.08), '-')
rep('alpha* alone', bisect_lever(lambda a: -W(replace(A, alpha=a)), 0.3, 5/6), '- (0.7 gives +0.0057)')
rep('restricted LS saving sLS* (U-exponent)', bisect_lever(lambda s: -W(replace(A, sLS=s)), 0.0, 0.2), 'redesign: 0.805 -> 0.770 at worst bin')
rep('mk* physical supply cap off', bisect_lever(lambda m: W(replace(A, mk=m, cap82=False)), 0.5, 1.0), 'none (mk=1: +0.0017)')
