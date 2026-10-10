"""Compare, at x = 1/2 (paper slots, S = 0.2055), the JOINT beyond-edge slope lamJ (brute force, s-freedom)
with the equivalent long-Moebius slope kA = 6 lamJ/5 (exact count).  Shows whether joint ~ kA also at x = 1/2."""
import numpy as np, counts_lab as CL, model_sa as MS
CL.install(); cfg = dict(MS.DEFAULT); cfg.update(CL.NEWDEF)
for lam in (0.75,):
    ck = dict(cfg, kA=6 * lam / 5)
    for d in (0.33, 0.37, 0.39, 0.41, 0.45):
        ex = float(MS.R_of(np.array([d]), np.array([0.5]), 0.2055, ck)[0])
        ex0 = float(MS.R_of(np.array([d]), np.array([0.5]), 0.2055, cfg)[0])
        bj = CL.joint_R(d, 0.5, 0.2055, cfg, lamJ=lam, nt=101, nr=1201, nm=401)
        print(f"lam={lam} delta={d}: paper {ex0:.5f} | kA={6*lam/5:.2f} exact {ex:.5f} | joint lam (grid) {bj[0]:.5f} (t={bj[1]:.3f})")
