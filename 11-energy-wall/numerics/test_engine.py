"""engine check: rowsums() against a direct evaluation with edge_lib sextic symbols (exact psi_k)."""
import sys, time
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
from edge_lib import setup, Symb, Columns
from nlib import Arith, Cols, Rows, rowsums, mark_chi_m2, gamma_prime, W
t0 = time.time()
ar = Arith(410000)
N0 = 3000
cols = Cols(ar, N0); rows = Rows(ar, 4000)
print("cols", cols.n, "rows", rows.n, f"{time.time()-t0:.1f}s")
# compare beta with edge_lib
fields, pnorm, plist = setup(2 * N0 + 10)
C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(cols.A, cols.B))}
idx = np.array([key[(int(a), int(b))] for a, b in zip(C.A, C.B)])
assert C.n == cols.n
print("beta max diff vs edge_lib:", np.abs(cols.beta[idx] - C.beta).max(), " w diff:", np.abs(cols.w[idx] - C.w).max())
rng = np.random.default_rng(0)
coef = np.stack([cols.beta * cols.w, np.exp(2j * np.pi * rng.random(cols.n)) * cols.w], axis=1)
S, D = rowsums(rows, cols, coef, (cols.w ** 2)[:, None], verbose=False)
# direct: psi_k(n) via edge_lib sextic symbols, for 300 random rows (incl. ones with conjugate pairs / inert)
fields2, pnorm2, plist2 = setup(4010)
symb = {}
bad = 0; worst = 0
pick = list(rng.choice(rows.n, 300, replace=False))
special = [i for i in range(rows.n) if len(set(ar.P[j].p for j in rows.fac[i])) < len(rows.fac[i]) or any(ar.P[j].kind == 1 for j in rows.fac[i])]
pick += special[:100]
for i in pick:
    v = np.ones(cols.n)
    for j in rows.fac[i]:
        pr = ar.P[j]; key_p = (pr.a, pr.b)
        if key_p not in symb: symb[key_p] = Symb(fields2[key_p])
        e = symb[key_p](cols.A, cols.B)
        v *= np.where(e < 0, 0, np.where(e % 2 == 1, -1.0, 1.0))
    Sd = v @ coef; Dd = (v * v) @ (cols.w ** 2)
    worst = max(worst, np.abs(Sd - S[i]).max(), abs(Dd - D[i, 0]))
print(f"rowsums vs direct sextic-symbol evaluation on {len(pick)} rows ({len(special[:100])} with inert/conjugate-pair factors): max |diff| = {worst:.2e}")
# marks: chi_P^{-2} vs edge_lib Symb
pr = [p for p in ar.P if p.kind == 0 and p.N > 200][0]
e = Symb(fields2[(pr.a, pr.b)])(cols.A, cols.B)
ref = np.where(e < 0, 0, np.exp(-2j * np.pi * 2 * e.astype(np.int64) / 6))
print("chi_P^-2 max diff:", np.abs(mark_chi_m2(pr, cols) - ref).max())
from edge_lib import gauss_vec
print("gamma_2 / gamma_1 at a prime vs edge_lib:", abs(gamma_prime(pr, 2) - gauss_vec(fields2[(pr.a, pr.b)], 2)), abs(gamma_prime(pr, 1) - gauss_vec(fields2[(pr.a, pr.b)], 1)))
print(f"total {time.time()-t0:.1f}s")
