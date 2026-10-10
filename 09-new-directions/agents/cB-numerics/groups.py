"""Group-averaged E = |T|^2/D vs X (mean over rows (h, nu) with the same N(h)), smoothed over
log-X windows of width `win` (in units of the 1/12 grid), with standard errors (rows treated as
independent, window points as correlated: SE from the spread of row-means).
usage: groups.py run.npz W fams [win] [stride]"""
import sys, numpy as np
R = np.load(sys.argv[1])
W = sys.argv[2]; fams = sys.argv[3].split(",")
win = int(sys.argv[4]) if len(sys.argv) > 4 else 12
stride = int(sys.argv[5]) if len(sys.argv) > 5 else 6
Xs = R["Xs"]; Nh = R["Nh"]; v = R[f"valid_{W}"]
kind = R["kind"] if "kind" in R else np.zeros(len(Nh))
idx = np.nonzero(v)[0]
centers = idx[win // 2: len(idx) - win // 2: stride]


def grouped(f):
    E = np.abs(R[f"T_{f}_{W}"]) ** 2 / R[f"D_{f}_{W}"]
    out = {}
    for n in np.unique(Nh):
        sel = Nh == n
        rows = E[sel]
        M = []; S = []
        for c in centers:
            blk = rows[:, c - win // 2: c + win // 2 + 1].mean(1)    # per-row window mean
            M.append(blk.mean()); S.append(blk.std(ddof=1) / np.sqrt(len(blk)) if len(blk) > 1 else np.nan)
        out[n] = (np.array(M), np.array(S), sel.sum())
    return out


if __name__ == "__main__":
    print("X centres: " + " ".join(f"{Xs[c]:8.0f}" for c in centers))
    for f in fams:
        g = grouped(f)
        print(f"== {f} [{W}] mean E (window {win}) +- SE")
        for n, (M, S, k) in g.items():
            print(f"N={n:6.0f} rows={k:3d}: " + " ".join(f"{m:5.2f}({s:4.2f})" for m, s in zip(M, S)))
