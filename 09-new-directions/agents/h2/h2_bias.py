"""
h2: mechanism check for the excess variance of T(c) along Kubota-aligned slot vectors (h2_coh.py).
Patterson / Heath-Brown ("Kummer's conjecture for cubic Gauss sums"): sum over primes P of the normalised cubic Gauss
sum gamma_2(P; n) = gamma_2(P) chi_P(n)^{-2} has a conjectured main term ~ c(n) X^{5/6}/log X, with c(n) a cubic-theta
coefficient at index n.  So the slot sum X_c(n) = sum_P c_P chi_P(n)^{-2} with c = gamma_2 (kub) should correlate with
theta-type functions of n, while c = 1, alpha, ... should not.
We measure  r(c, f) = sum_n w X_c(n) conj f(n) / sqrt(sum w |X_c|^2 sum w |f|^2)  for f in a dictionary of unimodular
functions (beta, gamma_2 (= beta alpha), their conjugates, times alpha^j), n over all columns; null |r| ~ 1/sqrt(#n).
"""
import sys, json, math
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/h2")
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
from edge_lib import setup, Symb, Columns, gauss_vec
from eisenstein import cx

Z6 = np.exp(2j * np.pi * np.arange(6) / 6)


def run(N0=20000, L0s=(100, 300, 1000, 3000), out=None):
    fields, pnorm, plist = setup(max(2 * N0, 2 * max(L0s)) + 10)
    C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
    w = C.w
    alpha_n = C.alpha / np.abs(C.alpha)
    beta = C.beta
    gam2 = beta * alpha_n                      # gamma_2(n) up to the twisted-multiplicative convention of edge_lib
    dic = {}
    for j in range(-3, 4):
        dic[f"beta*a^{j}"] = beta * alpha_n ** j
        dic[f"conjbeta*a^{j}"] = np.conj(beta) * alpha_n ** j
    res = dict(N0=N0, ncols=C.n, null=1 / math.sqrt(C.n), rows=[])
    rng = np.random.default_rng(5)
    for L0 in L0s:
        marks = [p for p in plist if fields[p].kind == "split" and L0 <= pnorm[p] < 2 * L0]
        E = np.stack([Symb(fields[p])(C.A, C.B) for p in marks], axis=1)
        chi = np.where(E < 0, 0, Z6[(-2 * E.astype(np.int64)) % 6])
        g2 = np.array([gauss_vec(fields[p], 2) for p in marks])
        al = np.array([cx(p) / abs(cx(p)) for p in marks])
        vecs = {"one": np.ones(len(marks)), "kub=g2": g2, "conj g2": np.conj(g2), "beta_P": g2 * np.conj(al),
                "p78": g2 * np.conj(al) ** 2, "conj p78": np.conj(g2) * al ** 2, "alpha": al,
                "random": np.exp(2j * np.pi * rng.random(len(marks)))}
        for vn, c in vecs.items():
            X = chi @ c
            rs = {k: complex(np.sum(w * X * np.conj(f)) / math.sqrt(np.sum(w * np.abs(X) ** 2) * np.sum(w * np.abs(f) ** 2)))
                  for k, f in dic.items()}
            best = max(rs, key=lambda k: abs(rs[k]))
            rec = dict(L0=L0, marks=len(marks), vec=vn, best=best, best_abs=abs(rs[best]),
                       rms_abs=float(np.sqrt(np.mean([abs(v) ** 2 for v in rs.values()]))),
                       r={k: [v.real, v.imag] for k, v in rs.items()})
            res["rows"].append(rec)
            print(f"L0={L0:<5d} L'={len(marks):<4d} c={vn:9s} max|r| = {abs(rs[best]):.4f} at {best:15s}  rms|r| {rec['rms_abs']:.4f}   (null {res['null']:.4f})", flush=True)
    if out: json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    if "L0s" in args: args["L0s"] = tuple(args["L0s"])
    run(**args)
