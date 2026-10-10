"""Identify the top singular direction of [S_{k,P}] (all marks): correlate the top right singular vector (mark space)
with chi6^j(disc P), chi6^j(P(0)), chi6^j(Res(P,P')) and the top left vector (row space) with chi6^j(disc k), chi6^j(k(0))."""
import sys, math, numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ffdriver import run
from ff_lib import irreducibles, Fq, resultant, deriv
q, kappa, lam, nu = map(int, sys.argv[1:5]); nsub = int(sys.argv[5]) if len(sys.argv) > 5 else 0
F = Fq(q)
marks = irreducibles(lam, q)
rows = None
if nsub:
    rng = np.random.default_rng(3); rows = sorted(rng.choice(q ** kappa, size=nsub, replace=False).tolist())
R = run(q, kappa, nu, lam, marks, threads=3, rows=rows, tag=f"out{q}_{kappa}_{lam}_{nu}")
S = np.where(R["cop"], R["S"], 0); D = np.array([R["d"][R["cop"][:, j], j].sum() for j in range(len(marks))])
U, sv, Vh = np.linalg.svd(S, full_matrices=False)
print(f"q={q} (kappa,lam,nu)=({kappa},{lam},{nu}) rows={S.shape[0]} marks={S.shape[1]} sv^2/meanDiag top6:", np.round(sv[:6]**2/D.mean(), 3),
      " RMT edge", round((1+math.sqrt(S.shape[1]/S.shape[0]))**2, 3))
v = Vh[0].conj(); u = U[:, 0]
def corr(x, y):
    x = x - 0; return abs(np.vdot(x, y)) / (np.linalg.norm(x) * np.linalg.norm(y))
feats = {}
for j in range(6):
    feats[f"chi6^{j}(disc P)"] = np.array([F.chi((-1) ** (lam * (lam - 1) // 2) * resultant(P, deriv(P, q), q), j) for P in marks])
    feats[f"chi6^{j}(P(0))"] = np.array([F.chi(P[0], j) if P[0] else 0 for P in marks])
for name, f in feats.items():
    c = corr(f, v)
    if c > 0.1: print("   mark-vector corr with", name, "=", round(c, 3))
print("   |v| mass profile: max|v|^2*#marks =", round(float(np.max(np.abs(v)**2) * len(v)), 2), "; frac of |v|^2 on top 10% marks =",
      round(float(np.sort(np.abs(v)**2)[::-1][:len(v)//10].sum()), 3))
kpolys = [tuple((int(ki) // q**i) % q for i in range(kappa)) + (1,) for ki in R["kidx"]]
for j in range(6):
    f = np.array([F.chi((-1) ** (kappa * (kappa - 1) // 2) * resultant(k, deriv(k, q), q), j) for k in kpolys])
    g = np.array([F.chi(k[0], j) if k[0] else 0 for k in kpolys])
    for name, x in ((f"chi6^{j}(disc k)", f), (f"chi6^{j}(k(0))", g)):
        c = corr(x, u)
        if c > 0.1: print("   row-vector corr with", name, "=", round(c, 3))
print("   |u| mass: max|u|^2*rows =", round(float(np.max(np.abs(u)**2) * len(u)), 2))
