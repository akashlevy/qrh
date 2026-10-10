import sys, numpy as np
from analyze import load, window_stats
d = load(sys.argv[1]); W = window_stats(d)
names = d["names"]
show = ["beta", "rand0", "bP7", "rP7", "bP97", "rP97", "bP409", "X:s100:p78", "XR:s100:p78", "X:s1000:p78", "X:s1000:ones", "X:t31x61:p78x2"]
show = [n for n in show if n in names]
print("F      rows/sec " + " ".join(f"{n[:14]:>15s}" for n in show))
for w in W:
    s = w["sec"][12]
    print(f"{w['F']:<6.3g} {s['rows_per_sector']:8.1f} " + " ".join(f"{s['rms'][n]:7.4f}/{w['rho'][n]:+.4f}" for n in show))
