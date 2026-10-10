"""Print band-averaged E = |T|^2/D vs X for each family (quick look).  usage: bands.py run.npz [W] [fams]"""
import sys, numpy as np
R = np.load(sys.argv[1])
W = sys.argv[2] if len(sys.argv) > 2 else "lg35"
fams = sys.argv[3].split(",") if len(sys.argv) > 3 else list(R["fams"])
Xs = R["Xs"]; Nh = R["Nh"]; kind = R["kind"]
valid = R[f"valid_{W}"]
prime = (kind // 10) == 1
bands = [(4, 4), (7, 7), (13, 13), (19, 19), (25, 31), (37, 43), (49, 79), (97, 151), (157, 241), (250, 400),
         (400, 650), (650, 1000), (1000, 1600), (1600, 2600), (2600, 4100), (4100, 6600), (6600, 10000)]
cols = [i for i in range(len(Xs)) if valid[i]][::6]
for f in fams:
    E = np.abs(R[f"T_{f}_{W}"]) ** 2 / R[f"D_{f}_{W}"]
    print(f"== {f}  ({W}); rows: N(h) band [n rows]; cols: X")
    print("            " + " ".join(f"{Xs[i]:7.0f}" for i in cols))
    for lo, hi in bands:
        sel = prime & (Nh >= lo) & (Nh <= hi)
        if sel.sum() == 0: continue
        m = E[sel].mean(0)
        print(f"{lo:5d}-{hi:<5d}[{sel.sum():3d}] " + " ".join(f"{m[i]:7.3f}" for i in cols))
