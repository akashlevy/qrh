import numpy as np, time, sys
from hls_lib import Family
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
from edge_lib import Columns
t0 = time.time()
F = Family(2000, 300, 50, nP=10)
C = Columns(2000, 4000, F.fields, F.pnorm, F.plist, None)
# align columns
key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(C.A, C.B))}
idx = np.array([key[(int(a), int(b))] for a, b in zip(F.A, F.B)])
print("beta match:", np.max(np.abs(F.beta - C.beta[idx])))
print("|gamma_j| ranges:", [(round(np.abs(F.gam[:, j]).min(), 6), round(np.abs(F.gam[:, j]).max(), 6)) for j in range(1, 6)])
# identity gamma_1 gamma_2 = mu alpha G (G from edge_lib)
lhs = F.gam[:, 1] * F.gam[:, 2]; rhs = F.mu * F.alpha * C.G[idx]
print("gamma1 gamma2 = mu alpha G :", np.max(np.abs(lhs - rhs)))
# gamma_5 = conj(gamma_1) * chi(-1)-ish ?
print("gamma5 vs conj gamma1 ratio set:", np.unique(np.round(F.gam[:, 5] / np.conj(F.gam[:, 1]), 6))[:6])
# operator adjointness
rng = np.random.default_rng(0)
a = rng.standard_normal(F.nc) + 1j * rng.standard_normal(F.nc)
Y = rng.standard_normal((F.nk, F.nP)) + 1j * rng.standard_normal((F.nk, F.nP))
print("adjoint:", abs(np.vdot(Y, F.S(a)) - np.vdot(F.Sh(Y), a)))
vals, vecs = F.top(k=3)
print("top eigs", vals, "rows", F.nrows, "cols", F.nc, time.time() - t0)
# brute force
Sfull = np.einsum('kn,pn->kpn', F.Psi.astype(float), F.H).reshape(F.nrows, F.nc)
ev = np.linalg.eigvalsh(Sfull @ Sfull.conj().T)
print("brute top", ev[-3:][::-1])
