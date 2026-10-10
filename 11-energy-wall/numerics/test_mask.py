"""check: masked design row sums X:s50:d == sum_{P not | k} c_{d,P} S_{k,P} (per-mark columns g*), and the
masked diagonal == sum_n |beta w X_k|^2 1_{(n,k)=1} computed directly from a Psi row."""
import json, numpy as np
from nlib import Arith, Cols, Rows, mark_chi_m2, gamma_prime, W
from run_energy import slot_primes, psi_rows
Z = np.load("data/energy_test.npz"); names = list(Z["names"]); S, D = Z["S"], Z["D"]; meta = json.loads(str(Z["meta"]))
ar = Arith(20000); cols = Cols(ar, 5000); rows = Rows(ar, 5000)
ps = slot_primes(ar, 50)
assert ps == meta["slots"]["s50"]["primes"]
al = np.array([complex(ar.P[i].a, 0) + ar.P[i].b * W for i in ps]); al /= np.abs(al)
qn = np.array([ar.P[i].N for i in ps], float)
cs = {"p78": np.conj(al) ** 2 * np.array([gamma_prime(ar.P[i], 2) for i in ps]) / np.sqrt(qn),
      "ones": 1 / np.sqrt(qn), "kub": np.conj(np.array([gamma_prime(ar.P[i], 1) for i in ps])) / np.sqrt(qn)}
G = S[:, [names.index(f"g{j}") for j in range(len(ps))]]
mask = np.array([[0.0 if i in set(f) else 1.0 for i in ps] for f in rows.fac])
worst = {}
for d, c in cs.items():
    ref = (G * mask) @ c
    worst[d] = np.abs(ref - S[:, names.index(f"X:s50:{d}")]).max() / np.abs(ref).max()
print("masked row sums vs per-mark recombination (rel. max diff):", worst)
hit = [k for k in range(rows.n) if set(rows.fac[k]) & set(ps)][:40]
M = np.stack([mark_chi_m2(ar.P[i], cols) for i in ps], axis=1)
Psi = psi_rows(rows, cols, hit)
err = 0
for r, k in enumerate(hit):
    X = M @ (cs["p78"] * mask[k])
    d = np.abs(Psi[r]) @ np.abs(cols.beta * X * cols.w) ** 2
    err = max(err, abs(d - D[k, names.index("X:s50:p78")]) / d)
print(f"masked diagonal on {len(hit)} rows: max rel diff {err:.2e}; rows with slot factors: {int((mask.min(1) == 0).sum())} / {rows.n}")
