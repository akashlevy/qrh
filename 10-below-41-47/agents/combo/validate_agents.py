"""Validation 3: reproduce each agent's thresholds at ell = 0.1773 (wall), B = 41/47, kappa = 3/4 fixed."""
import numpy as np, sys
from combo_model import *
ELL = (1 - 12*(TARGET - 5/6))/3      # 0.17730...
K = 0.75
A = replace(PAPER08)                  # agents' base: physical marked supply, cap82 on, joint lem:plain lever off
NJ = dict(joint=False)
bins = BINS_SEARCH
def W(lev, ell=ELL, B=TARGET, kappa=K, bins=bins):
    return best_over_b(ell, B, lev, kappa, bins)
rows = [
 ('paper (agents: +0.01119 amplitude, +0.01124 marked-cap)', A),
 ('mk=1 saturated, physical supply, cap off (+0.00187)', replace(A, mk=1, cap82=False)),
 ('mk .8, cP .4 [cN fixed] (+0.00032)', replace(A, mk=.8, cap82=False, cP=.4, **NJ)),
 ('mk .8, cP .45 [cN fixed] (-0.00018)', replace(A, mk=.8, cap82=False, cP=.45, **NJ)),
 ('mk sat, cP .5 [cN fixed] (-0.00085)', replace(A, mk=1, cap82=False, cP=.5, **NJ)),
 ('mk .75, alpha .75 (-0.00088)', replace(A, mk=.75, cap82=False, alpha=.75)),
 ('mk sat, alpha .78 (-0.00004)', replace(A, mk=1, cap82=False, alpha=.78)),
 ('mk sat, alpha .76 (-0.00070)', replace(A, mk=1, cap82=False, alpha=.76)),
 ('mk sat, cP 1 [cN fixed] (-0.00828)', replace(A, mk=1, cap82=False, cP=1.0, **NJ)),
 ('r+2z<=1.13, cap off (+0.00231)', replace(A, sM=.13, cap82=False)),
 ('r+2z<=1.045, cap off (+0.00794)', replace(A, sM=.045, cap82=False)),
 ('alpha 0.7 alone (BRIEF +0.0057)', replace(A, alpha=.7)),
 ('DH counts (BRIEF -0.028)', replace(A, DH=True)),
 ('mk=1, UNLIMITED supply, cap82 kept (amplitude +0.0043)', replace(A, mk=1, mark_supply='unlimited')),
 ('kappa=35/47 (plain-cap -0.00002 vs paper)', replace(A)),
 ('sixth moment eta=0 (numerator -0.0116)', replace(A, eta=0.0)),
 ('sLS ideal restricted large sieve (redesign -0.029)', replace(A, sLS=1.0)),
]
for name, lev in rows:
    kap = 35/47 if 'kappa=35/47' in name else K
    e, b = W(lev, kappa=kap)
    print(f'{name:60s} worst={e:+.5f} (b={b:.4f})', flush=True)
