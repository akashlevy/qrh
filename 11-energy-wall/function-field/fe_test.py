"""Per-row functional-equation test: for each row k (fixed mark P), v_k(nu) = |S_{k,P}(nu)|^2 / q^nu,
nu = 0..numax.  Report, per row, the best D with v_k(nu) = v_k(D - nu) exactly (all 0<=nu<=D), and the
'polynomial degree' (last nu with v_k(nu) >= 1/q)."""
import sys, numpy as np
from collections import Counter
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ffdriver import run, choose_marks
q, kappa, lam, numax = map(int, sys.argv[1:5])
marks = choose_marks(q, lam, 1)
V = []; S = []
for nu in range(numax + 1):
    R = run(q, kappa, nu, lam, marks, tag=f"fe{q}_{kappa}_{lam}_{nu}")
    cop = R["cop"][:, 0]
    V.append(R["absS2"][:, 0].astype(float) / q ** nu); S.append(R["S"][:, 0])
V = np.array(V).T; S = np.array(S).T; cop = R["cop"][:, 0]; kidx = R["kidx"]
V = V[cop]; S = S[cop]; kidx = kidx[cop]
res = Counter(); exD = Counter()
examples = {}
for i in range(len(V)):
    v = V[i]
    best = None
    for D in range(numax, -1, -1):
        if all(abs(v[n] - v[D - n]) < 1e-9 * max(1, v[n]) for n in range(D + 1)):
            best = D; break
    lastbig = max(n for n in range(numax + 1) if v[n] >= 1.0 / q - 1e-12)
    res[(best, lastbig)] += 1
    if (best, lastbig) not in examples:
        examples[(best, lastbig)] = (int(kidx[i]), np.round(v, 4).tolist())
print(f"q={q} kappa={kappa} lam={lam} P={marks[0]} rows={len(V)}: (exact-FE degree D, last nu with v>=1/q): count")
for key, c in sorted(res.items(), key=lambda t: -t[1]):
    print("  ", key, c, " example row", examples[key][0], examples[key][1])
