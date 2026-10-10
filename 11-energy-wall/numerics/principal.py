"""Principal row k = 1 (excluded from all windows): |S_1|^2 / d_1 per design -- the Patterson-bias diagnostic.
wall-spectral s.4.3 predicts for Gauss-sum marks a non-oscillating part of beta X of relative size ~ (L^2/N)^{1/6}/log L,
living on principal rows only, so |S_1|^2/d_1 ~ N (L^2/N)^{1/3}/log^2 L there, vs O(1) (exponential) for random marks."""
import glob, json, numpy as np
from analyze import load
rows = []
for f in sorted(glob.glob("data/energy_N*.npz"), key=lambda f: int(f.split("_N")[1][:-4])):
    d = load(f); N0 = d["meta"]["N0"]
    i1 = int(np.nonzero(d["rowN"] == 1)[0][0])
    r = {n: float(np.abs(d["S"][i1, j]) ** 2 / d["D"][i1, j]) for j, n in enumerate(d["names"]) if n.startswith("X") or n in ("beta",) or n.startswith("rand")}
    rows.append((N0, d["meta"]["cols"], r))
names = [n for n in rows[0][2] if n.startswith("X:")]
print("|S_1|^2/d_1 at the principal row (random-phase expectation ~1)")
print(f"{'design':18s} " + " ".join(f"{'N0='+str(N0):>11s}" for N0, _, _ in rows))
for n in ["beta", "rand0"] + names + [m for m in rows[0][2] if m.startswith("XR:s1000")]:
    print(f"{n:18s} " + " ".join(f"{r[n]:11.2f}" for _, _, r in rows))
json.dump({str(N0): r for N0, _, r in rows}, open("principal.json", "w"), indent=1)
