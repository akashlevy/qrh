"""Profile mean_k |S_{k,P}(nu)|^2/q^nu over a random subset of rows (cost per row ~ q^nu).
Prints normalised energy profile and rho(nu) = E/Diag - 1 on the subset; also the theta/control comparison."""
import sys, numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ffdriver import run, choose_marks
q, kappa, lam, numax, nrows = map(int, sys.argv[1:6])
markseed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
budget = float(sys.argv[7]) if len(sys.argv) > 7 else 6e9
marks = choose_marks(q, lam, 1, seed=markseed)
R0 = run(q, kappa, 0, lam, marks, tag=f"pr0_{q}_{kappa}_{lam}")
allrows = [int(r) for r, c in zip(R0["kidx"], R0["cop"][:, 0]) if c]
rng = np.random.default_rng(5)
rows = sorted(rng.choice(allrows, size=min(nrows, len(allrows)), replace=False).tolist())
print(f"q={q} kappa={kappa} lam={lam} P={marks[0]} rows used {len(rows)}/{len(allrows)}  (D-hypotheses: kappa+lam+1={kappa+lam+1}, 2kappa+lam-2={2*kappa+lam-2})")
for nu in range(numax + 1):
    if len(rows) * q ** nu > budget: break
    R = run(q, kappa, nu, lam, marks, tag=f"pr{q}_{kappa}_{lam}_{nu}", rows=rows)
    a2 = R["absS2"][:, 0].astype(float); d = R["d"][:, 0]
    print(f"  nu={nu:2d}: mean|S|^2/q^nu={a2.mean()/q**nu:8.4f}  mean d/q^nu={d.mean()/q**nu:6.4f}  rho={a2.sum()/max(d.sum(),1)-1:+8.4f}  frac0={np.mean(a2==0):5.3f}", flush=True)
