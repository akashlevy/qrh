"""First scan at q = 7: rho = E_P/Diag - 1 for theta (beta) vs random-phase control."""
import sys, json, numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ffdriver import run, choose_marks, stats
q = int(sys.argv[1]) if len(sys.argv) > 1 else 7
grid = [(1,1),(1,2),(1,3),(2,2),(2,3),(2,4),(3,3),(3,4),(3,5),(4,4),(4,5),(4,6)]
for lam in (0, 1, 2):
    for (kappa, nu) in grid:
        if q ** (kappa + nu) > 3e9: continue
        marks = choose_marks(q, lam, 1 if lam <= 1 else 3)
        R = run(q, kappa, nu, lam, marks)
        st = stats(R, coeff_sets=())
        Rc = run(q, kappa, nu, lam, marks, control=11)
        sc = stats(Rc, coeff_sets=())
        print(f"q={q} lam={lam} kappa={kappa} nu={nu}: rows={st['rows']:6d} cols/row={st['cols_per_row']:9.1f} "
              f"rho_theta={np.round(st['rho'],4)} m4={np.round(st['m4'],3)} | rho_ctrl={np.round(sc['rho'],4)} m4c={np.round(sc['m4'],3)}", flush=True)
