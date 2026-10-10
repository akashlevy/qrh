"""(H2)/Gram check on the Gram slot (44 marks, N(P) in [300,600)), rows with (P,k)=1 for each entry.
c_PP' = G_PP'/sqrt(G_PP G_P'P'); random model rms ~ n_eff^{-1/2}, n_eff = (sum V)^2/sum V^2.
Rayleigh quotient of the P78-coherent vector c_P = conj(alpha_P)^2 gamma_2(P) N(P)^{-1/2}:
R(c) = c*Gc / (|c|^2 mean diag)  (diagonal-dominance: R = 1 + rho; sup = top eigenvalue / mean diag)."""
import glob, json, math, numpy as np
from analyze import load
from nlib import Arith, gamma_prime, W
out = {}
ar = Arith(410000)
for f in sorted(glob.glob("data/energy_N*.npz"), key=lambda f: int(f.split("_N")[1][:-4])):
    d = load(f); N0 = d["meta"]["N0"]; gram = d["meta"]["gram"]
    J = [d["names"].index(f"g{j}") for j in range(len(gram))]
    Sg = d["S"][:, J]
    Mk = np.array([[g["idx"] not in fs for g in gram] for fs in d["fac"]], float)
    ps = [ar.P[g["idx"]] for g in gram]
    al = np.array([complex(p.a, 0) + p.b * W for p in ps]); al /= np.abs(al)
    c = np.conj(al) ** 2 * np.array([gamma_prime(p, 2) for p in ps]) / np.sqrt([p.N for p in ps])
    rng = np.random.default_rng(3)
    res = []
    for F in [2 ** (j / 2) for j in range(0, 30)]:
        K = N0 / F; x = d["rowN"] / K
        win = (x > 0.5) & (x <= 1) & (d["rowN"] > 1)
        if win.sum() < 200: break
        V = np.sin(2 * np.pi * (x[win] - 0.5)) ** 2
        Sm = Sg[win] * Mk[win]
        G = (V[:, None] * Sm).T @ np.conj(Sm)
        G = G.T          # G[P,P'] = sum V S_P conj S_P'
        dg = np.real(np.diag(G))
        ch = G / np.sqrt(np.outer(dg, dg))
        off = ch[~np.eye(len(dg), dtype=bool)]
        neff = V.sum() ** 2 / (V ** 2).sum()
        Rc = float(np.real(np.conj(c) @ G @ c) / (np.vdot(c, c).real * dg.mean()))
        rr = [np.exp(2j * np.pi * rng.random(len(c))) for _ in range(200)]
        Rr = np.array([float(np.real(np.conj(v) @ G @ v) / (len(v) * dg.mean())) for v in rr])
        lam = float(np.linalg.eigvalsh(G)[-1] / dg.mean())
        mp = (1 + math.sqrt(len(dg) / neff)) ** 2
        res.append(dict(F=F, rows=int(win.sum()), neff=neff, h2_rms=float(np.sqrt(np.mean(np.abs(off) ** 2))),
                        h2_rms_x_sqrt_neff=float(np.sqrt(np.mean(np.abs(off) ** 2)) * math.sqrt(neff)),
                        R_p78=Rc, R_rand_mean=float(Rr.mean()), R_rand_sd=float(Rr.std()), sup=lam, MP_edge=mp))
    out[N0] = res
    print(f"N0={N0}")
    for r in res:
        print(f"  F={r['F']:<6.3g} rows={r['rows']:<6d} h2 rms*sqrt(neff)={r['h2_rms_x_sqrt_neff']:.3f}  R(p78)={r['R_p78']:.3f}  "
              f"R(rand)={r['R_rand_mean']:.3f}±{r['R_rand_sd']:.3f}  sup={r['sup']:.3f} (MP edge {r['MP_edge']:.3f})")
json.dump({str(k): v for k, v in out.items()}, open("h2.json", "w"), indent=1)
