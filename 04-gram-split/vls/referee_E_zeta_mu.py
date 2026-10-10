import sys, numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/gram_split/vls")
from e_lib import *
R = RowsE(6000, ymin_over_Q=4.0, zeta_check=10)
sf = np.where(R.sf)[0]
mu = (-1.0) ** R.nfac[sf]
lhs = R.zetaZ[sf]
rhs = R.psi1[sf] * mu * np.conj(R.alpha[sf]) * np.conj(R.G[2][sf])
print("rows %d: max |zeta - psi1 mu alphabar conj(gamma2)| = %.1e" % (len(sf), np.abs(lhs - rhs).max()))
rhs2 = R.psi1[sf] * np.conj(R.alpha[sf]) * np.conj(R.G[2][sf])
print("without mu: fraction of rows where it differs (|diff|>1e-9): %.3f" % np.mean(np.abs(lhs - rhs2) > 1e-9))
