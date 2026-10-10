"""Locate the Symb-LUT vs ResidueField.sextic_exp mismatches reported by edge-theta/test_lib.py."""
import sys
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
import numpy as np
from edge_lib import setup, Symb
fields, pnorm, plist = setup(20000)
rng = np.random.default_rng(0)
bad = []
for p in plist[:400]:
    F = fields[p]; S = Symb(F)
    A = rng.integers(-500, 500, 50); B = rng.integers(-500, 500, 50)
    e = S(A, B)
    for a, b, k in zip(A, B, e):
        kk = F.sextic_exp((int(a), int(b)))
        kk = -1 if kk is None else kk
        if kk != k:
            bad.append((p, F.kind, F.q, int(a), int(b), int(k), kk))
print(len(bad))
for x in bad[:40]: print(x)
