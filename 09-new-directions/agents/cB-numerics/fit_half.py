"""Crossing point X_1/2 (window-mean E = 1/2, log-interpolated) per N(h) group; slope c of log X_1/2 vs log N(h)
with bootstrap over rows (h, nu).  usage: fit_half.py run.npz W fams [level]"""
import sys, numpy as np
R = np.load(sys.argv[1]); W = sys.argv[2]; fams = sys.argv[3].split(","); lev = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
Xs = R["Xs"]; Nh = R["Nh"]; v = R[f"valid_{W}"]; lx = np.log(Xs[v])
rng = np.random.default_rng(0)
def cross(E):          # E: rows x X (valid only); smooth over 12-pt window, first downward crossing of lev
    m = np.convolve(E.mean(0), np.ones(12) / 12, mode="same")
    m[:6] = 1; m[-6:] = m[-7]
    k = np.nonzero(m < lev)[0]
    k = [i for i in k if i > 6 and np.all(m[i:i + 12] < lev + 0.15)]
    if not k: return np.nan
    i = k[0]; t = (m[i - 1] - lev) / (m[i - 1] - m[i]); return lx[i - 1] + t * (lx[i] - lx[i - 1])
for f in fams:
    E = (np.abs(R[f"T_{f}_{W}"]) ** 2 / R[f"D_{f}_{W}"])[:, v]
    Ns = [n for n in np.unique(Nh) if n >= float(sys.argv[5] if len(sys.argv) > 5 else 2)]
    xh = {n: cross(E[Nh == n]) for n in np.unique(Nh)}
    def slope(xd):
        pts = [(np.log(n), xd[n]) for n in Ns if np.isfinite(xd[n])]
        if len(pts) < 3: return np.nan
        a = np.array(pts); return np.polyfit(a[:, 0], a[:, 1], 1)[0]
    bs = []
    for _ in range(300):
        xd = {}
        for n in np.unique(Nh):
            rows = E[Nh == n]; xd[n] = cross(rows[rng.integers(0, len(rows), len(rows))])
        bs.append(slope(xd))
    bs = np.array(bs); bs = bs[np.isfinite(bs)]
    print(f"{f:4s} X_half: " + " ".join(f"N{n:.0f}:{np.exp(x):.3g}" for n, x in xh.items()),
          f"| slope c = {slope(xh):.2f}  bootstrap {np.mean(bs):.2f} +- {np.std(bs):.2f} (n={len(bs)})")
