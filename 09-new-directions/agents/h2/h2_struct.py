"""
h2: structure hunt for the cross-mark Gram (hypothesis (H2) of wall-spectral), exact Z[w] arithmetic.

  D_P(k) = sum_{n sf odd primary, N0 < N(n) <= 2N0} beta(n) w(N(n)/N0) psi_k(n) chi_P(n)^{-2}
  G_{P1P2}(K) = sum_k V(N(k)/K) conj D_{P1}(k) D_{P2}(k)      (V = sin^2 bump on (1/2, 1])
  Delta_{P1P2}  = exact n-diagonal  sum_k V sum_n |beta w|^2 conj chi_{P1}^{-2} chi_{P2}^{-2}(n) 1_{(n,k)=1}
  O = G - Delta (off-diagonal).   C = G / sqrt(G11 G22) (normalised).

Tests (theta beta vs random-phase control with identical |w| and identical marks):
  T1  rms|C_12|, mixed-diagonal size |Delta_12|, by window K.
  T2  idea 1: does |O_12|^2 <= O_11 O_22 hold?  fraction of indefinite 2x2 off-diagonal forms.
  T3  phase bias: z = |sum_{i<j} C_ij conj(phi_ij)| / sqrt(sum |C_ij|^2) for structured pair phases phi:
      sextic symbols (P_i/P_j)_6^a (P_j/P_i)_6^b, Groessencharacter angles (alpha_i conj alpha_j)^m, beta(P_i) conj beta(P_j).
  T4  pair invariants: mean |C|^2 by cubic symbol class, conjugate pairs, N(P1 P2) (L0 scan).
  T5  persistence: complex correlation of C_ij across disjoint dyadic K windows, across 12 angular row sectors,
      and across disjoint column halves.  A deterministic (P1,P2)-dependent component would show up here.
  T6  hidden cocycle phases: twist column P at row k by chi_k(P)^a chi_P(k)^b (a,b in Z/6) and recompute rms|C_12|.
  T7  coherent vectors: c^H (G - diag G) c / (|c|^2 mean diag) for structured c vs random c, and the operator norm.
"""
import sys, json, time, math
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
from edge_lib import setup, Symb, Columns, gauss_vec      # noqa
import dual_energy as de                                    # noqa  (path set by edge_lib)
from eisenstein import cx                                   # noqa

Z6 = np.exp(2j * np.pi * np.arange(6) / 6)


def build(N0, KMAX, L0s, nmark, seed=7, block=1500, keep_psi2=False):
    t0 = time.time()
    fields, pnorm, plist = setup(max(2 * N0, KMAX, 2 * max(L0s)) + 10)
    C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
    rA, rB, rN, rfac = de.primary_sf_elements(0, KMAX, fields, plist, pnorm)
    odd = np.array([all(pnorm[p] != 4 for p in f) for f in rfac], bool)
    keep = np.nonzero(odd)[0]
    rA, rB, rN = rA[keep], rB[keep], rN[keep]; rfac = [rfac[i] for i in keep]
    o = np.argsort(rN, kind="stable"); rA, rB, rN = rA[o], rB[o], rN[o]; rfac = [rfac[i] for i in o]
    nr = len(rN)
    rng = np.random.default_rng(seed)
    # marks: split primes, for each L0 the first nmark of norm in [L0, 2L0)
    marks, mL0 = [], []
    for L0 in L0s:
        ms = [p for p in plist if fields[p].kind == "split" and L0 <= pnorm[p] < 2 * L0]
        if nmark and len(ms) > nmark:
            idx = np.sort(rng.choice(len(ms), nmark, replace=False)); ms = [ms[i] for i in idx]
        marks += ms; mL0 += [L0] * len(ms)
    Lp = len(marks)
    print(f"N0={N0}: cols {C.n}, rows {nr}, marks {Lp} {dict(zip(*np.unique(mL0, return_counts=True)))}  {time.time()-t0:.1f}s", flush=True)
    # column marks chi_P(n)^{-2}
    Ecol = np.stack([Symb(fields[p])(C.A, C.B) for p in marks], axis=1)
    chiP2 = np.where(Ecol < 0, 0, Z6[(-2 * Ecol.astype(np.int64)) % 6])
    bw = C.beta * C.w
    u = np.exp(2j * np.pi * rng.random(C.n))                       # random-phase control
    CP = bw[:, None] * chiP2
    CPr = (u * C.w)[:, None] * chiP2
    halfA = C.N <= 1.5 * N0
    # rows: symbols
    rowprimes = sorted({p for f in rfac for p in f}, key=lambda p: pnorm[p])
    par, zer = {}, {}
    for p in rowprimes:
        e = Symb(fields[p])(C.A, C.B); zer[p] = e < 0; par[p] = (e % 2 == 1)
    # chi_P(k) = (k/P)_6 exponents (rows x marks), chi_k(P) = prod_{pi|k} (P/pi)_6 exponents
    yk = np.stack([Symb(fields[p])(rA, rB) for p in marks], axis=1).astype(np.int64)   # -1 if P | k
    MA = np.array([p[0] for p in marks], np.int64); MB = np.array([p[1] for p in marks], np.int64)
    epi = {p: Symb(fields[p])(MA, MB).astype(np.int64) for p in rowprimes}               # (P/pi)_6
    xk = np.zeros((nr, Lp), np.int64); xz = np.zeros((nr, Lp), bool)
    for i, f in enumerate(rfac):
        for p in f:
            e = epi[p]; xz[i] |= e < 0; xk[i] += np.maximum(e, 0)
    xk %= 6
    print(f"row symbols for {len(rowprimes)} primes  {time.time()-t0:.1f}s", flush=True)
    # windows
    Ks = [2 ** j for j in range(3, 40) if 2 ** j <= KMAX]
    Vw = np.zeros((nr, len(Ks)))
    for j, K in enumerate(Ks):
        x = rN / K
        Vw[:, j] = np.where((x > 0.5) & (x <= 1), np.sin(2 * np.pi * (x - 0.5)) ** 2, 0.0)
    Vw[rN == 1] = 0.0
    SP = np.zeros((nr, Lp), complex); SPr = np.zeros((nr, Lp), complex)
    SPA = np.zeros((nr, Lp), complex); SPB = np.zeros((nr, Lp), complex)
    Dv = np.zeros((C.n, len(Ks)))                                     # sum_k V 1_{(n,k)=1}
    Psi2 = np.zeros((nr, C.n), np.float32) if keep_psi2 else None
    for b0 in range(0, nr, block):
        b1 = min(nr, b0 + block)
        Psi = np.empty((b1 - b0, C.n), np.float64)
        for i in range(b0, b1):
            f = rfac[i]
            if not f:
                Psi[i - b0] = 1.0; continue
            pr = np.zeros(C.n, bool); zz = np.zeros(C.n, bool)
            for p in f: pr ^= par[p]; zz |= zer[p]
            Psi[i - b0] = np.where(zz, 0.0, np.where(pr, -1.0, 1.0))
        SPA[b0:b1] = Psi[:, halfA] @ CP[halfA]
        SPB[b0:b1] = Psi[:, ~halfA] @ CP[~halfA]
        SPr[b0:b1] = Psi @ CPr
        Dv += (Psi ** 2).T @ Vw[b0:b1]
        if keep_psi2: Psi2[b0:b1] = Psi ** 2
    SP = SPA + SPB
    print(f"rows done {time.time()-t0:.1f}s", flush=True)
    sect = (np.floor((np.angle(rA + rB * complex(-0.5, math.sqrt(3) / 2)) + np.pi) / (2 * np.pi) * 12).astype(int)) % 12
    # mark data
    g2 = np.array([gauss_vec(fields[p], 2) for p in marks])
    alpha = np.array([cx(p) / abs(cx(p)) for p in marks])
    betaP = g2 * np.conj(alpha)
    e_sym = np.stack([Symb(fields[p])(MA, MB) for p in marks], axis=1).astype(np.int64)   # e_sym[i, j] = (P_i / P_j)_6
    qn = np.array([pnorm[p] for p in marks])
    conjpair = np.zeros((Lp, Lp), bool)
    for i, p in enumerate(marks):
        for j, p2 in enumerate(marks):
            conjpair[i, j] = (i != j) and qn[i] == qn[j]
    return dict(N0=N0, KMAX=KMAX, Ks=Ks, Vw=Vw, SP=SP, SPr=SPr, SPA=SPA, SPB=SPB, CP=CP, CPr=CPr, Dv=Dv,
                sect=sect, rN=rN, xk=xk, xz=xz, yk=yk, mL0=np.array(mL0), marks=marks, betaP=betaP, alpha=alpha,
                e_sym=e_sym, qn=qn, conjpair=conjpair, halfA=halfA, ncols=C.n, nrows=nr, Psi2=Psi2)


def gram(S, V):
    return S.conj().T @ (V[:, None] * S)


def normc(G):
    d = np.sqrt(np.real(np.diag(G)))
    return G / np.outer(d, d)


def zstat(Cu, phi):
    """|sum C conj phi| / sqrt(sum |C|^2) over the given (flattened) unordered pairs; null E z^2 = 1"""
    return float(abs(np.sum(Cu * np.conj(phi))) / math.sqrt(np.sum(np.abs(Cu) ** 2)))


def analyse_block(B, sel, jK, label, out, rng):
    """all tests for the mark subset sel (indices) at window index jK"""
    V = B["Vw"][:, jK]; K = B["Ks"][jK]
    Lp = len(sel)
    iu = np.triu_indices(Lp, 1)
    res = dict(K=K, rows=int((V > 0).sum()), marks=Lp, label=label)
    for nm, S, CPm in (("theta", B["SP"][:, sel], B["CP"][:, sel]), ("random", B["SPr"][:, sel], B["CPr"][:, sel])):
        G = gram(S, V)
        Dl = CPm.conj().T @ (B["Dv"][:, jK][:, None] * CPm)
        O = G - Dl
        C = normc(G)
        d = np.real(np.diag(G))
        On = O / np.sqrt(np.outer(d, d))                 # pure off-diagonal, normalised
        r = {}
        r["rms_C"] = float(np.sqrt(np.mean(np.abs(C[iu]) ** 2)))
        r["rms_O"] = float(np.sqrt(np.mean(np.abs(On[iu]) ** 2)))
        r["rms_Delta_rel"] = float(np.sqrt(np.mean(np.abs((Dl / np.sqrt(np.outer(d, d)))[iu]) ** 2)))
        r["rho_diag_mean"] = float(np.mean(np.real(np.diag(O)) / np.real(np.diag(Dl))))
        # T2: (Off)^{1/2} inequality
        o11 = np.real(np.diag(O))[iu[0]]; o22 = np.real(np.diag(O))[iu[1]]; o12 = O[iu]
        r["frac_ineq_fails"] = float(np.mean(np.abs(o12) ** 2 > np.maximum(o11, 0) * np.maximum(o22, 0)))
        indef = (o11 < 0) | (o22 < 0) | (o11 * o22 < np.abs(o12) ** 2)
        r["frac_2x2_indefinite"] = float(np.mean(indef))
        r["frac_O11_negative"] = float(np.mean(np.real(np.diag(O)) < 0))
        r["median_|O12|/sqrt|O11 O22|"] = float(np.median(np.abs(o12) / np.sqrt(np.abs(o11 * o22) + 1e-300)))
        r["median_|Delta12|/|O12|"] = float(np.median(np.abs(Dl[iu]) / np.abs(o12)))
        es = B["e_sym"][np.ix_(sel, sel)]
        eij = es[iu]; eji = es.T[iu]
        al = B["alpha"][sel]; be = B["betaP"][sel]
        cub = (np.where(eij < 0, -1, (2 * eij) % 6) == 0)
        cp = B["conjpair"][np.ix_(sel, sel)][iu]
        lq = np.log(B["qn"][sel][iu[0]] * B["qn"][sel][iu[1]])
        for tag, Cu in (("G", C[iu]), ("O", On[iu])):
            # T3: phase bias
            zs = {}
            for a in range(6):
                for b in range(6):
                    ph = np.where((eij < 0) | (eji < 0), 0, Z6[(a * eij + b * eji) % 6])
                    zs[f"sym({a},{b})"] = zstat(Cu, ph)
            for m in (1, 2, 3, 4, 5, 6, 12):
                zs[f"alpha^{m}"] = zstat(Cu, (al[iu[0]] * np.conj(al[iu[1]])) ** m)
                zs[f"alpha^-{m}"] = zstat(Cu, (np.conj(al[iu[0]]) * al[iu[1]]) ** m)
            for nmph, ph in (("beta", be[iu[0]] * np.conj(be[iu[1]])), ("conjbeta", np.conj(be[iu[0]]) * be[iu[1]]),
                             ("p78", (be * np.conj(al))[iu[0]] * np.conj((be * np.conj(al))[iu[1]])),
                             ("p78conj", np.conj((be * np.conj(al))[iu[0]]) * (be * np.conj(al))[iu[1]])):
                zs[nmph] = zstat(Cu, ph)
            r["z_" + tag] = zs
            r["z_max_" + tag] = max(zs.values()); r["z_argmax_" + tag] = max(zs, key=zs.get)
            r["z2_mean_" + tag] = float(np.mean(np.array(list(zs.values())) ** 2))
            # T4: pair invariants
            a2 = np.abs(Cu) ** 2; m0 = a2.mean()
            r["cubic1_" + tag] = float(a2[cub].mean() / m0) if cub.any() else None
            r["conjpair_" + tag] = float(a2[cp].mean() / m0) if cp.any() else None
            r["slopeN_" + tag] = float(np.polyfit(lq, np.log(a2 + 1e-300), 1)[0]) if lq.std() > 0 else None
        r["frac_cubic1"] = float(cub.mean())
        # T7: coherent vectors
        Dm = np.mean(d)
        Goff = O - np.diag(np.diag(O))          # pure off-diagonal (k-off-diagonal, P1 != P2)
        def Rq(c): return float(np.real(c.conj() @ Goff @ c) / (np.vdot(c, c).real * Dm))
        cands = {"one": np.ones(Lp, complex), "beta": be, "conjbeta": np.conj(be), "p78": be * np.conj(al),
                 "p78conj": np.conj(be * np.conj(al)), "alpha": al, "conjalpha": np.conj(al)}
        r["T_coh"] = {k: Rq(c) for k, c in cands.items()}
        rr = np.array([Rq(np.exp(2j * np.pi * rng.random(Lp))) for _ in range(200)])
        r["T_rand_mean_sd_max"] = [float(rr.mean()), float(rr.std()), float(np.abs(rr).max())]
        ev = np.linalg.eigvalsh(Goff / Dm)
        r["opnorm_off"] = [float(ev[0]), float(ev[-1])]
        keff = V.sum() ** 2 / (V ** 2).sum()
        r["rmt_off_edge"] = float(2 * math.sqrt(Lp / keff) + Lp / keff)
        res[nm] = r
    out.append(res)
    t, q = res["theta"], res["random"]
    f2 = lambda x: "na" if x is None else f"{x:.2f}"
    print(f"[{label}] K={K:<6d} rows={res['rows']:<5d} L'={Lp:<4d} rms C th {t['rms_C']:.4f} rnd {q['rms_C']:.4f} | rms O th {t['rms_O']:.4f} rnd {q['rms_O']:.4f} | "
          f"Delta/sqrt {t['rms_Delta_rel']:.2e} | ineq fails {t['frac_ineq_fails']:.2f}/{q['frac_ineq_fails']:.2f} "
          f"indef {t['frac_2x2_indefinite']:.2f} | O: zmax th {t['z_max_O']:.2f}({t['z_argmax_O']}) rnd {q['z_max_O']:.2f} "
          f"<z2> {t['z2_mean_O']:.2f}/{q['z2_mean_O']:.2f} | O: cub1 {f2(t['cubic1_O'])}/{f2(q['cubic1_O'])} conj {f2(t['conjpair_O'])}/{f2(q['conjpair_O'])} "
          f"[G conj {f2(t['conjpair_G'])}/{f2(q['conjpair_G'])}] "
          f"| Tcoh p78 {t['T_coh']['p78']:+.3f} rand {t['T_rand_mean_sd_max'][0]:+.3f}±{t['T_rand_mean_sd_max'][1]:.3f} "
          f"| opnorm {t['opnorm_off'][1]:.2f} (rnd {q['opnorm_off'][1]:.2f}, rmt {t['rmt_off_edge']:.2f})", flush=True)


def complex_corr(a, b):
    return complex(np.sum(a * np.conj(b)) / math.sqrt(np.sum(np.abs(a) ** 2) * np.sum(np.abs(b) ** 2)))


def offnorm(S, V, CPm, Dvj):
    """normalised pure off-diagonal (G - Delta)/sqrt(G11 G22) for rows weighted by V, coprimality weights Dvj"""
    G = gram(S, V)
    Dl = CPm.conj().T @ (Dvj[:, None] * CPm)
    d = np.sqrt(np.real(np.diag(G)))
    return (G - Dl) / np.outer(d, d)


def persistence(B, sel, label, out):
    """T5: complex correlations of the off-diagonal O_ij across disjoint windows / sectors / column halves"""
    Lp = len(sel); iu = np.triu_indices(Lp, 1)
    Psi2 = B["Psi2"]
    res = dict(label=label, marks=Lp, windows=[])
    for nm, S, CPm in (("theta", B["SP"][:, sel], B["CP"][:, sel]), ("random", B["SPr"][:, sel], B["CPr"][:, sel])):
        Cs = []
        for j, K in enumerate(B["Ks"]):
            V = B["Vw"][:, j]
            if (V > 0).sum() < 4 * Lp: Cs.append(None); continue
            Cs.append(offnorm(S, V, CPm, B["Dv"][:, j])[iu])
        cw = []
        for j in range(len(Cs) - 1):
            if Cs[j] is None or Cs[j + 1] is None: continue
            c = complex_corr(Cs[j], Cs[j + 1]); cw.append(dict(K=B["Ks"][j], corr_re=c.real, corr_abs=abs(c)))
        res.setdefault("dyadic", {})[nm] = cw
    # sectors and halves at the largest windows
    for j in range(len(B["Ks"]) - 3, len(B["Ks"])):
        V = B["Vw"][:, j]; K = B["Ks"][j]
        w = dict(K=K)
        for nm, S, CPm in (("theta", B["SP"][:, sel], B["CP"][:, sel]), ("random", B["SPr"][:, sel], B["CPr"][:, sel])):
            Csec = []
            for s in range(12):
                Vs = V * (B["sect"] == s)
                Csec.append(offnorm(S, Vs, CPm, Psi2.T @ Vs)[iu])
            Csec = np.array(Csec)                                   # 12 x pairs
            ratio = np.mean(np.abs(Csec.sum(0)) ** 2) / np.mean(np.sum(np.abs(Csec) ** 2, 0))
            # mean pairwise sector correlation
            cc = [complex_corr(Csec[a], Csec[b]).real for a in range(12) for b in range(a + 1, 12)]
            w[nm + "_sector_coherence_ratio"] = float(ratio)      # null ~1, fully deterministic ~12
            w[nm + "_sector_corr_mean"] = float(np.mean(cc))
        hA = B["halfA"]
        for nm, SA, SB in (("theta", B["SPA"][:, sel], B["SPB"][:, sel]),):
            cA = offnorm(SA, V, B["CP"][:, sel] * hA[:, None], B["Dv"][:, j])[iu]
            cB = offnorm(SB, V, B["CP"][:, sel] * (~hA)[:, None], B["Dv"][:, j])[iu]
            c = complex_corr(cA, cB); w["theta_halves_corr"] = [c.real, abs(c)]
        res["windows"].append(w)
        print(f"[{label}] persistence K={K}: sector coherence ratio th {w['theta_sector_coherence_ratio']:.3f} "
              f"rnd {w['random_sector_coherence_ratio']:.3f} (null 1, det 12); mean sector corr th {w['theta_sector_corr_mean']:+.4f} "
              f"rnd {w['random_sector_corr_mean']:+.4f}; column halves corr {w['theta_halves_corr'][0]:+.4f}", flush=True)
    for nm in ("theta", "random"):
        print(f"[{label}] dyadic-window corr {nm}: " + " ".join(f"{d['K']}:{d['corr_re']:+.3f}" for d in res["dyadic"][nm]), flush=True)
    out.append(res)


def twisted(B, sel, jK, label, out):
    """T6: twist column P at row k by chi_k(P)^a chi_P(k)^b and recompute the off-diagonal rms"""
    V = B["Vw"][:, jK]; K = B["Ks"][jK]
    Lp = len(sel); iu = np.triu_indices(Lp, 1)
    xk = B["xk"][:, sel]; yk = B["yk"][:, sel]; bad = B["xz"][:, sel] | (yk < 0)
    res = dict(label=label, K=K)
    for nm, S in (("theta", B["SP"][:, sel]), ("random", B["SPr"][:, sel])):
        tab = np.zeros((6, 6))
        for a in range(6):
            for b in range(6):
                ph = np.where(bad, 0, Z6[(a * xk + b * np.maximum(yk, 0)) % 6])
                C = normc(gram(S * ph, V))
                tab[a, b] = np.sqrt(np.mean(np.abs(C[iu]) ** 2))
        res[nm] = tab.tolist()
        base = tab[0, 0]
        print(f"[{label}] twisted rms/untwisted ({nm}) K={K}: min {tab.min()/base:.3f} max {tab.max()/base:.3f} "
              f"argmax {np.unravel_index(tab.argmax(), tab.shape)}", flush=True)
    out.append(res)


def run(N0=20000, KMAX=None, L0s=(30, 100, 300, 1000, 3000), nmark=60, out=None, seed=7):
    KMAX = KMAX or N0
    rng = np.random.default_rng(seed + 1)
    B = build(N0, KMAX, L0s, nmark, seed, keep_psi2=True)
    res = dict(N0=N0, KMAX=KMAX, L0s=list(L0s), nmark=nmark, ncols=B["ncols"], nrows=B["nrows"], blocks=[], persistence=[], twisted=[])
    jbig = [j for j, K in enumerate(B["Ks"]) if K >= KMAX // 8]
    allsel = np.arange(len(B["marks"]))
    for L0 in L0s:
        sel = np.nonzero(B["mL0"] == L0)[0]
        if len(sel) < 4: continue
        for j in jbig:
            analyse_block(B, sel, j, f"L0={L0}", res["blocks"], rng)
    for j in jbig:
        analyse_block(B, allsel, j, "all", res["blocks"], rng)
    persistence(B, allsel, "all", res["persistence"])
    twisted(B, allsel, len(B["Ks"]) - 1, "all", res["twisted"])
    if out:
        json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    if "L0s" in args: args["L0s"] = tuple(args["L0s"])
    run(**args)
