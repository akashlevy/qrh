"""
h2 T8: complete multiplicative-character spectrum, in the row variable k, of the cross-mark products
   g_12(k) = V(N(k)/K) conj D_{P1}(k) D_{P2}(k),
i.e.  F(chi) = sum_k g_12(k) conj chi(k)  for EVERY character chi of (Z[w]/P1P2)^x  (split P1, P2: F_p1^x x F_p2^x).
F(trivial) is the cross term G_12 itself (rows prime to P1P2).  If D_P(k) carried a hidden phase chi_P(k)^a or
psi_k(P) (cocycle / root-number structure), |F| would spike at the corresponding character.  Null: all |F|^2/mean ~ Exp(1).
Reports, aggregated over pairs: the rank of the trivial character, the 36 sextic points, the global max and where it sits.
"""
import sys, json, time, math
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/h2")
from h2_struct import build


def primroot(p):
    fac = [q for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g


def dlog_table(p):
    g = primroot(p); t = np.full(p, -1, np.int64); x = 1
    for e in range(p - 1):
        t[x] = e; x = x * g % p
    return t


def run(N0=20000, L0=40, nmark=14, jwin=(-1, -2), out=None):
    t0 = time.time()
    B = build(N0, N0, (L0,), nmark, keep_psi2=True)
    sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
    from edge_lib import setup
    fields, pnorm, plist = setup(max(2 * N0, 2 * L0) + 10)
    # recover row elements (same construction as build)
    import dual_energy as de
    rA, rB, rN, rfac = de.primary_sf_elements(0, N0, fields, plist, pnorm)
    odd = np.array([all(pnorm[p] != 4 for p in f) for f in rfac], bool)
    keep = np.nonzero(odd)[0]; rA, rB, rN = rA[keep], rB[keep], rN[keep]
    o = np.argsort(rN, kind="stable"); rA, rB, rN = rA[o], rB[o], rN[o]
    assert np.array_equal(rN, B["rN"])
    V = sum(B["Vw"][:, j] for j in jwin)
    marks = B["marks"]; Lp = len(marks)
    xs = []
    for P in marks:
        F = fields[P]; p = F.q
        red = (rA + rB * F.r) % p
        xs.append((dlog_table(p)[red], p))
    res = dict(N0=N0, L0=L0, marks=Lp, pairs=[])
    agg = {"theta": [], "random": []}
    for i in range(Lp):
        for j in range(i + 1, Lp):
            (xi, pi), (xj, pj) = xs[i], xs[j]
            ok = (xi >= 0) & (xj >= 0) & (V > 0)
            rec = dict(i=i, j=j, p1=int(pi), p2=int(pj))
            for nm, S, CPm in (("theta", B["SP"], B["CP"]), ("random", B["SPr"], B["CPr"])):
                dk = B["Psi2"][ok] @ (np.conj(CPm[:, i]) * CPm[:, j])      # exact per-row n-diagonal
                g = V[ok] * (np.conj(S[ok, i]) * S[ok, j] - dk)
                Fm = np.zeros((pi - 1, pj - 1), complex)
                np.add.at(Fm, (xi[ok], xj[ok]), g)
                Fh = np.fft.fft2(Fm)                      # Fh[a,b] = sum g e(-(a x_i/(p_i-1) + b x_j/(p_j-1)))
                A2 = np.abs(Fh) ** 2
                m = A2.mean()
                r0 = A2[0, 0] / m
                rank = float((A2 > A2[0, 0]).mean())          # fraction of characters larger than trivial
                amax = np.unravel_index(A2.argmax(), A2.shape)
                sext = []
                for a in range(6):
                    for b in range(6):
                        if (pi - 1) % 6 == 0 and (pj - 1) % 6 == 0:
                            sext.append(A2[a * (pi - 1) // 6, b * (pj - 1) // 6] / m)
                sext = np.array(sext)
                rec[nm] = dict(triv_over_mean=float(r0), triv_rank=rank, max_over_mean=float(A2.max() / m),
                               argmax=[int(amax[0]), int(amax[1])],
                               argmax_is_sextic=bool(amax[0] % ((pi - 1) // 6) == 0 and amax[1] % ((pj - 1) // 6) == 0),
                               sextic_mean=float(sext.mean()), sextic_max=float(sext.max()),
                               nchar=int((pi - 1) * (pj - 1)), nrows=int(ok.sum()))
                agg[nm].append(rec[nm])
            res["pairs"].append(rec)
    for nm in ("theta", "random"):
        a = agg[nm]
        tr = np.array([x["triv_over_mean"] for x in a]); rk = np.array([x["triv_rank"] for x in a])
        mx = np.array([x["max_over_mean"] for x in a]); nc = np.array([x["nchar"] for x in a])
        sm = np.array([x["sextic_mean"] for x in a]); sx = np.array([x["sextic_max"] for x in a])
        ism = np.array([x["argmax_is_sextic"] for x in a])
        exp_max = np.log(nc) + 0.5772
        summ = dict(pairs=len(a), triv_over_mean_mean=float(tr.mean()), triv_rank_mean=float(rk.mean()),
                    triv_rank_KS_uniform=float(np.max(np.abs(np.sort(rk) - (np.arange(len(rk)) + 0.5) / len(rk)))),
                    max_over_mean_mean=float(mx.mean()), expected_max_iid=float(exp_max.mean()),
                    max_over_expected_max=float((mx / exp_max).mean()), argmax_sextic_frac=float(ism.mean()),
                    sextic_mean=float(sm.mean()), sextic_max_mean=float(sx.mean()))
        res["summary_" + nm] = summ
        print(f"{nm:7s}: pairs {len(a)}  |G12|^2/mean|F|^2 = {tr.mean():.3f} (null 1)  rank of trivial {rk.mean():.3f} (null 0.5, KS {summ['triv_rank_KS_uniform']:.3f})  "
              f"max/mean {mx.mean():.2f} (iid {exp_max.mean():.2f})  argmax sextic {ism.mean():.3f}  sextic pts mean {sm.mean():.3f} max {sx.mean():.2f}", flush=True)
    print(f"done {time.time()-t0:.1f}s")
    if out: json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    if "jwin" in args: args["jwin"] = tuple(args["jwin"])
    run(**args)
