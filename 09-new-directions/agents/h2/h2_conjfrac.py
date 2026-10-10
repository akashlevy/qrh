"""
h2: explains the uniform excess of t^2 = (T(c)/sd_random)^2 for structured c in h2_coh.py.
For theta, beta(conj n) = conj beta(n) gives O_{conj P, conj Q} = conj O_{P,Q} exactly; all test vectors satisfy
c_{conj P} = conj c_P, so the (P,Q) and (conj P, conj Q) terms of T(c) add coherently: E t^2 = 1 + f,
f = fraction of ordered off-diagonal pairs (P,Q), Q != conj P, whose conjugate pair is also in the mark set.
Random-phase columns break the symmetry (E t^2 = 1).  Reconstructs the mark sets exactly as h2_struct.build does.
"""
import sys, json
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
from edge_lib import setup

def fracs(N0, L0s, nmark, seed=7):
    fields, pnorm, plist = setup(max(2 * N0, N0, 2 * max(L0s)) + 10)
    rng = np.random.default_rng(seed)
    out = {}
    for L0 in L0s:
        ms = [p for p in plist if fields[p].kind == "split" and L0 <= pnorm[p] < 2 * L0]
        if nmark and len(ms) > nmark:
            idx = np.sort(rng.choice(len(ms), nmark, replace=False)); ms = [ms[i] for i in idx]
        S = set(ms)
        conj = lambda p: next(q for q in plist if pnorm[q] == pnorm[p] and q != p and fields[q].kind == "split")
        cj = {p: conj(p) for p in ms}
        tot = good = 0
        for P in ms:
            for Q in ms:
                if Q == P or Q == cj[P]: continue
                tot += 1; good += (cj[P] in S and cj[Q] in S)
        out[L0] = good / tot
    return out

res = {}
for N0, L0s in ((20000, (400, 520, 680, 880, 1150, 1500, 1950, 2550, 3300, 4300)),
                (40000, (450, 600, 780, 1000, 1300, 1700, 2200, 2900, 3700, 4800))):
    f = fracs(N0, L0s, 50)
    r = json.load(open(f"runs/coh_N{N0}.json"))
    blocks = r["blocks"]
    names = [k for k in blocks[0]["theta"] if k != "_sd"]
    fb = np.array([f[b["L0"]] for b in blocks])
    th = np.array([[b["theta"][k] ** 2 for k in names] for b in blocks])
    rn = np.array([[b["random"][k] ** 2 for k in names] for b in blocks])
    # pooled over the 9 structured vectors
    pred = 1 + fb.mean()
    res[N0] = dict(f_by_L0=f, f_mean=float(fb.mean()), pred_t2=float(pred),
                   obs_t2_theta=float(th.mean()), obs_t2_random=float(rn.mean()),
                   obs_by_vec={k: float(th[:, j].mean()) for j, k in enumerate(names)})
    # regression of per-block mean t^2 (over vectors) on f
    y = th.mean(1); A = np.vstack([np.ones_like(fb), fb]).T
    coef = np.linalg.lstsq(A, y, rcond=None)[0]
    res[N0]["fit_t2_vs_f"] = [float(coef[0]), float(coef[1])]
    print(f"N0={N0}: f by L0 {', '.join(f'{k}:{v:.2f}' for k, v in f.items())}")
    print(f"   predicted mean t^2 (theta, structured c) = 1 + <f> = {pred:.3f};  observed theta {th.mean():.3f}, random columns {rn.mean():.3f};"
          f"  per-block fit t^2 = {coef[0]:.2f} + {coef[1]:.2f} f  (prediction 1 + 1 f)")
json.dump(res, open("runs/conjfrac.json", "w"), indent=1)
