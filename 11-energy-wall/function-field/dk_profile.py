"""Profile of the twisted series D_{k,P}(u) = sum_nu S_{k,P}(nu) u^nu over squarefree n:
mean_k |S_{k,P}(nu)|^2 / q^nu and the fraction of rows with S = 0, as nu varies (fixed q, kappa, P)."""
import sys, numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ffdriver import run, choose_marks
q = int(sys.argv[1]); kappa = int(sys.argv[2]); lam = int(sys.argv[3]); numax = int(sys.argv[4])
marks = choose_marks(q, lam, 1)
print(f"q={q} kappa={kappa} lam={lam} mark={marks[0]}")
for nu in range(0, numax + 1):
    if q ** (kappa + nu) > 4e9: break
    R = run(q, kappa, nu, lam, marks, tag=f"prof{q}_{kappa}_{lam}_{nu}")
    cop = R["cop"][:, 0]; a2 = R["absS2"][cop, 0].astype(float); d = R["d"][cop, 0]
    print(f" nu={nu:2d}: mean|S|^2/q^nu = {a2.mean()/q**nu:9.5f}   mean d/q^nu = {d.mean()/q**nu:7.4f}   "
          f"rho = {a2.sum()/d.sum()-1:+8.4f}   frac(S=0) = {np.mean(a2==0):6.3f}   max|S|^2/q^nu = {a2.max()/q**nu:8.3f}", flush=True)
