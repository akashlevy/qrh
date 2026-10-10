import sys, numpy as np
sys.path.insert(0, "/home/user/qrh/03-six-sevenths/routes/r4")
import model as M4
cfg = dict(M4.DEFAULT)
S = 0.20547772
def pieces(delta, x, t, r):
    q = x*delta
    kA, w1, w2, cP, mk, w2s = cfg["kA"], cfg["w1"], cfg["w2"], cfg["cP"], cfg["mk"], cfg["w2scale"]
    fA = max(1.0, (1 + 5*kA*r)/6) - delta*r
    zI = max(min(mk*(1+w1-r), w2s*(3+w2-2*r), S), 0)
    fI = 1 - delta*r - 2*q*zI if r < 1+w1 else np.inf
    m = max(t - r, 0)
    zP = min(cP*max(1-2*m, 0), S)
    fP = 1 - 2*delta*m - 2*q*zP if m < 0.5 else np.inf
    fP0 = 1 - 2*delta*m
    return fA, fI, fP, fP0
for delta in [0.3, 0.35, 0.387, 0.39, 0.45, 0.6]:
    for x in [0.0, 0.25, 0.5]:
        best = None
        for t in np.linspace(1, 1.5, 501):
            rs = np.linspace(t-0.5, t, 2001)
            vals = [min(pieces(delta, x, t, r)) for r in rs]
            k = int(np.argmax(vals)); v = vals[k]
            if best is None or v < best[0]:
                best = (v, t, rs[k], pieces(delta, x, t, rs[k]))
        v, t, r, pc = best
        print(f"delta={delta:.3f} x={x:.2f}: R={v:.4f} (DH {1-delta:.4f}, excess {v-(1-delta):+.4f})  t={t:.3f} r*={r:.3f} m*={t-r:.3f}  fA={pc[0]:.4f} fI={pc[1]:.4f} fP={pc[2]:.4f} fP0={pc[3]:.4f}")
