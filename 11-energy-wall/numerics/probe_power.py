"""positive control for the cross-scale lambda estimator: plant S' = S + eps * c * conj(A_sigma(Y0)) (Y0 = 2^15, an s-range
disjoint from the estimator's Y in [2^16, 2^18], so the plant correlates with them only through the Y-stable rho part)
with c chosen so that the plant has the norm of S at eps = 1, and recompute lambda(rho)."""
import json, math, numpy as np
import probe
from analyze import load
from nlib import Arith
N0, F = 50000, 4
d = load("data/energy_N50000.npz")
keep = d["rowN"] <= N0 / F
rows = dict(A=d["rowA"][keep].astype(np.int64), B=d["rowB"][keep].astype(np.int64), N=d["rowN"][keep])
ar = Arith(410000)
Ys = [2 ** 15, 2 ** 16, 2 ** 17, 2 ** 18]
Aout, Sd, E, pidx = probe.build_A(ar, rows, Ys)
K = N0 / F; x = rows["N"] / K
win = (x > 0.5) & (x <= 1) & (rows["N"] > 1)
V = np.sin(2 * np.pi * (x[win] - 0.5)) ** 2
rc = probe.cls4(rows["A"], rows["B"])[win]
j = d["names"].index("X:s100:p78")
S0 = d["S"][keep][win, j]
def lam(S, classes=range(12)):
    big = Ys[1:]; Tc, Er, fr = {}, {}, {}
    for Y in big:
        a2 = np.abs(Aout[Y][:, win]) ** 2; nz = Aout[(Y, "noise")][:, win]
        fr[Y] = float(np.sum(a2) / np.sum(nz) - 1)
        a2 = np.abs(Aout[Y][list(classes)][:, win]) ** 2; nz = Aout[(Y, "noise")][list(classes)][:, win]
        fr[Y] = float(np.sum(a2) / np.sum(nz) - 1)
        T = np.zeros((12, 12), complex); er = 0.0
        for sg in classes:
            u = Aout[Y][sg][win]; np.add.at(T[sg], rc, V * u * S); er += float(np.sum(V ** 2 * np.abs(u) ** 2 * np.abs(S) ** 2))
        Tc[Y] = T; Er[Y] = er
    out = []
    for a in range(len(big)):
        for b in range(a + 1, len(big)):
            Y1, Y2 = big[a], big[b]
            out.append(float(np.real(np.sum(Tc[Y1] * np.conj(Tc[Y2])))) / math.sqrt(fr[Y1] * fr[Y2] * Er[Y1] * Er[Y2]))
    return float(np.mean(out)), fr
res = {}
for sg0 in range(12):
    plant = np.conj(Aout[Ys[0]][sg0, win])
    if np.linalg.norm(plant) == 0: continue
    plant = plant * np.linalg.norm(S0) / np.linalg.norm(plant)
    r = {}
    for eps in [0, 0.1, 0.3, 1.0]:
        l, fr = lam(S0 + eps * plant, classes=[sg0])
        r[eps] = l
    res[sg0] = r
    print(f"class {sg0:2d}: lambda_sigma0(rho) for eps = 0, 0.1, 0.3, 1: " + "  ".join(f"{v:7.2f}" for v in r.values()) +
          f"   (f at Y=2^16..2^18: {[round(v, 3) for v in fr.values()]})", flush=True)
for eps in [0, 0.1, 0.3, 1.0]:
    v = np.array([r[eps] for r in res.values()])
    print(f"eps={eps:<4}: mean over classes {v.mean():7.2f}  median {np.median(v):7.2f}   (random level 1; need-violation level ~ n_eff = {V.sum()**2/(V**2).sum():.0f})")
json.dump({str(k): {str(e): x for e, x in v.items()} for k, v in res.items()}, open("probe_power.json", "w"), indent=1)
