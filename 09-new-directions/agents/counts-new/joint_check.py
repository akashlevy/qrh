"""Check: a JOINT moment for the product witness M_r S_m with slope lamJ (count 1 + lamJ(s-1) - delta s, s = r+m > 1),
with the detector's freedom r + m >= t, gives R(delta,0) = R0(delta; min(alpha, lamJ)) and the optimal detector is t = 1.
Also: without s-freedom (m = t - r forced) the answer would be different (shows the freedom matters)."""
import numpy as np
import counts_lab as CL, model_sa as MS
cfg = dict(MS.DEFAULT)
for d in (0.27, 0.33, 0.39):
    for lam in (5/6, 0.80, 0.75, 0.70):
        v, t = CL.joint_R(d, 0.0, 0.2, cfg, lamJ=lam, nt=51, nr=601, nm=251)
        print(f"delta={d:.2f} lamJ={lam:.3f}: R={v:.5f} (t*={t:.3f})   R0(min(5/6,lam))={CL.R0(d, min(5/6, lam)):.5f}   R0(5/6)={CL.R0(d):.5f}")
print("no joint (grid check):", [round(CL.joint_R(d, 0.0, 0.2, cfg, None, nt=51, nr=601, nm=251)[0], 5) for d in (0.27, 0.33, 0.39)])
# x = 1/2, paper slots, S = 0.2055
for d in (0.35, 0.39, 0.41):
    base = CL.joint_R(d, 0.5, 0.2055, cfg, None, nt=51, nr=601, nm=251)[0]
    for lam in (5/6, 0.75):
        v, t = CL.joint_R(d, 0.5, 0.2055, cfg, lamJ=lam, nt=51, nr=601, nm=251)
        print(f"x=1/2 delta={d:.2f} lamJ={lam:.3f}: R={v:.5f} (t={t:.3f}) vs paper {base:.5f}")
