"""
h2 T7': is the coherent off-diagonal T(c) = c^H (G - Delta - diag) c special for structured mark vectors c?
Many disjoint mark sets (independent blocks) x 3 dyadic windows (independent after removing the n-diagonal).
For each block and vector c: t(c) = T(c) / sd_random(T), sd over 400 random-phase c on the same block.
Null (no alignment): t ~ N(0,1) (real part), so mean t^2 ~ 1.  Same statistic on random-phase columns as control.
Vectors: one, kub = gamma_2(P) (Kubota-aligned: carries Patterson's bias), conj kub, beta = gamma_2 conj(alpha),
p78 = gamma_2 conj(alpha)^2 (mark-design's coherent choice), conj p78, alpha, alpha^2.
"""
import sys, json, time, math
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/h2")
from h2_struct import build, gram


def run(N0=20000, L0s=(400, 520, 680, 880, 1150, 1500, 1950, 2550, 3300, 4300), nmark=50, nrand=400, out=None):
    t0 = time.time()
    B = build(N0, N0, tuple(L0s), nmark)
    rng = np.random.default_rng(11)
    jw = [j for j, K in enumerate(B["Ks"]) if K >= N0 // 8]
    res = dict(N0=N0, L0s=list(L0s), blocks=[])
    for L0 in L0s:
        sel = np.nonzero(B["mL0"] == L0)[0]
        if len(sel) < 10: continue
        al = B["alpha"][sel]; be = B["betaP"][sel]; g2 = be * al
        vecs = {"one": np.ones(len(sel), complex), "kub": g2, "conjkub": np.conj(g2), "beta": be, "conjbeta": np.conj(be),
                "p78": be * np.conj(al), "conjp78": np.conj(be * np.conj(al)), "alpha": al, "alpha2": al ** 2}
        for j in jw:
            V = B["Vw"][:, j]
            rec = dict(L0=L0, K=B["Ks"][j], marks=len(sel))
            for nm, S, CPm in (("theta", B["SP"][:, sel], B["CP"][:, sel]), ("random", B["SPr"][:, sel], B["CPr"][:, sel])):
                G = gram(S, V)
                Dl = CPm.conj().T @ (B["Dv"][:, j][:, None] * CPm)
                O = G - Dl; O = O - np.diag(np.diag(O))
                Dm = np.mean(np.real(np.diag(G)))
                T = lambda c: float(np.real(c.conj() @ O @ c) / (np.vdot(c, c).real * Dm))
                rr = np.array([T(np.exp(2j * np.pi * rng.random(len(sel)))) for _ in range(nrand)])
                sd = rr.std()
                rec[nm] = {k: T(c) / sd for k, c in vecs.items()}
                rec[nm]["_sd"] = float(sd)
            res["blocks"].append(rec)
        print(f"L0={L0}: " + " ".join(f"K={r['K']}: p78 {r['theta']['p78']:+.2f} kub {r['theta']['kub']:+.2f} one {r['theta']['one']:+.2f}"
                                      for r in res["blocks"][-len(jw):]), flush=True)
    names = list(res["blocks"][0]["theta"].keys()); names.remove("_sd")
    summ = {}
    for nm in ("theta", "random"):
        summ[nm] = {k: dict(mean_t=float(np.mean([b[nm][k] for b in res["blocks"]])),
                            mean_t2=float(np.mean([b[nm][k] ** 2 for b in res["blocks"]])),
                            n=len(res["blocks"])) for k in names}
    res["summary"] = summ
    print("blocks:", len(res["blocks"]), " (null: mean t = 0 +- 1/sqrt(n), mean t^2 = 1 +- sqrt(2/n))")
    for k in names:
        a, b = summ["theta"][k], summ["random"][k]
        print(f"  {k:9s} theta: mean t {a['mean_t']:+.3f}  mean t^2 {a['mean_t2']:.3f}   | random cols: mean t {b['mean_t']:+.3f}  mean t^2 {b['mean_t2']:.3f}")
    print(f"done {time.time()-t0:.1f}s")
    if out: json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    run(**args)
