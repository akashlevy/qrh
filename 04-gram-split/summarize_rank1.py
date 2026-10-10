import json, numpy as np, sys
d = json.load(open(sys.argv[1])); Ys = np.array(d["Ys"], float); G = np.array(d["G_re"])
K = len(Ys); vals = []; w = []
for i in range(K):
    for j in range(i + 1, K):
        vals.append(G[i, j] * (Ys[i] * Ys[j]) ** (1 / 3)); w.append((i, j))
vals = np.array(vals)
print("off-diag kappa estimates: mean %.4f  sd %.4f  se %.4f  n=%d" % (vals.mean(), vals.std(), vals.std() / np.sqrt(len(vals)), len(vals)))
big = np.array([v for v, (i, j) in zip(vals, w) if i >= K // 2])
print("  pairs with i >= K/2: mean %.4f  sd %.4f  n=%d" % (big.mean(), big.std(), len(big)))
diag = np.array([G[i, i] * Ys[i] for i in range(K)])
# fit diag*Y = c0 + kappa Y^{1/3}
Xd = np.stack([np.ones(K), Ys ** (1 / 3)], 1)
c, *_ = np.linalg.lstsq(Xd, diag, rcond=None)
print("diag fit G_ii Y = %.3f + %.4f Y^{1/3}" % tuple(c))
