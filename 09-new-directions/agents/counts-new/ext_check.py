"""Brute-force check (x = 0, witness freedom r + m >= t) of the closed forms in diag_cf.py:
 joint diagonal extension w (sum_u|M_r S_m|^2 << U for r+m <= 1+w, MVT slope 1 beyond),
 Moebius diagonal up to r1 (then copies slope 5/6), plain sixth moment diagonal up to 3m <= 1+w6."""
import numpy as np, diag_cf as C
def brute(d, kind, p, nt=101, nr=801, nm=401, a=5/6):
    best = np.inf
    for t in np.linspace(1, 1.5, nt):
        r = np.linspace(max(t - .5, 0), t, nr)[:, None]; m = np.linspace(0, .5, nm)[None, :]
        ok = r + m >= t - 1e-12; s = r + m
        if kind == "r1":
            fA = 1 + a * np.clip(r - p, 0, None) - d * r
        else:
            fA = np.maximum(1, 1 + a * (r - 1)) - d * r
        F = np.minimum(fA, 1 - 2 * d * m)
        if kind == "joint":
            F = np.minimum(F, 1 + np.clip(s - 1 - p, 0, None) - d * s)
        if kind == "cube":
            F = np.minimum(F, np.where(3 * m <= 1 + p, 1 - 3 * d * m, np.inf))
        best = min(best, np.where(ok, F, -np.inf).max())
    return best
for d in (0.27, 0.33):
    print(d, "r1=1.0426:", round(brute(d, "r1", 1.0426), 5), round(float(C.R_r1(d, 1.0426)), 5),
          "| joint w=.0657:", round(brute(d, "joint", 0.0657), 5), round(float(C.R_joint(d, 0.0657)), 5),
          "| cube w=.2355:", round(brute(d, "cube", 0.2355), 5), round(float(C.R_cube(d, 0.2355)), 5),
          "| R0", round(C.R0(d), 5))
