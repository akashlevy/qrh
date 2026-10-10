"""
mark-design: numerical structure of the cross-mark Gram matrix G beyond the energy wall.

Dual (reflected) energy at the binding dyad, in Gram form (wall-spectral / edge-theta normalisation):
   S_{k,P} = sum_{n sf odd primary, N0 < N(n) <= 2N0} beta(n) w(n) psi_k(n) chi_P(n)^{-2}
   G_{P1 P2}(K) = sum_{N(k) <= K} S_{k,P1} conj S_{k,P2}              (L' x L', PSD)
   E(c) = c^* G c,   D = diag G.
beta = gamma_2 conj(alpha) (P78 theta coefficient), psi_k = chi_k^3 (quadratic rows), marks chi_P^{-2}.

Weight vectors c (dual basis):
   p78   : conj(alpha(P))^2 * gamma_2(P)  -- the dual phase of P78's coherent marks for eta = 1
           (Cor. marked-reflection: zeta_Q contains abar(c)^2 chi_q(sigma_q)^{-2}(-omega_{q,0}),
            omega_{q,0} = -tau^+_{q,2} chi_q(eps_q)^{-2}; single active mark, c_F = 1)
   ones  : 1
   kub   : conj gamma_1(P)       (Kubota-aligned)
   rand  : random unit phases (200 draws)
Questions:
   (1) Is G - D close to low rank?  (eigen-spectrum of Ghat - I, Frobenius share of the top r, stable rank)
   (2) Where does the coherent vector sit in the spectrum?  (Rayleigh quotient vs random / sup; overlap with top
       eigenvectors; gain of the optimal weights c* = G^{-1} loc, i.e. loc^* G^{-1} loc vs Rayleigh)
   (3) Does any restriction of the prime set make G block-structured?
       (a) cubic class of (pi0/P)_3 for a fixed small prime pi0 (= ray class by reciprocity)
       (b) conjugate pairs P, conj(P)
       (c) mutual cubic symbol (P1/P2)_3 = 1 versus != 1, and greedy cubic cliques (size ~ L'^{2/3}?)
Controls: random-phase columns (same |beta w|), and iid complex Gaussian S of the same shape.
"""
import sys, json, time, math
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/energy-wall")
from edge_lib import setup, Symb, Columns, gauss_vec          # noqa: E402
import dual_energy as de                                       # noqa: E402
from eisenstein import cx                                      # noqa: E402


def build(N0, KMAX, L0, seed=3, block=1500):
    t0 = time.time()
    fields, pnorm, plist = setup(max(2 * N0, KMAX, 2 * L0))
    C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
    rA, rB, rN, rfac = de.primary_sf_elements(0, KMAX, fields, plist, pnorm)
    odd = np.array([all(pnorm[p] != 4 for p in f) for f in rfac], bool)
    keep = np.nonzero(odd)[0]
    rA, rB, rN = rA[keep], rB[keep], rN[keep]; rfac = [rfac[i] for i in keep]
    o = np.argsort(rN, kind="stable"); rN = rN[o]; rfac = [rfac[i] for i in o]
    marks = [p for p in plist if fields[p].kind == "split" and L0 <= pnorm[p] < 2 * L0]
    Em = np.stack([Symb(fields[p])(C.A, C.B) for p in marks], axis=1)
    chiP2 = np.where(Em < 0, 0, np.exp(-2j * np.pi * 2 * Em.astype(np.int64) / 6))
    g1 = np.array([gauss_vec(fields[p], 1) for p in marks])
    g2 = np.array([gauss_vec(fields[p], 2) for p in marks])
    al = np.array([cx(p) / abs(cx(p)) for p in marks])
    rng = np.random.default_rng(seed)
    coef = {"theta": C.beta * C.w, "randcol": np.exp(2j * np.pi * rng.random(C.n)) * np.abs(C.beta) * C.w}
    rowprimes = sorted({p for f in rfac for p in f}, key=lambda p: pnorm[p])
    par, zer = {}, {}
    for p in rowprimes:
        e = Symb(fields[p])(C.A, C.B)
        zer[p] = e < 0; par[p] = (e % 2 == 1)
    nr = len(rN)
    S = {k: np.zeros((nr, len(marks)), complex) for k in coef}
    Mcol = {k: (v[:, None] * chiP2) for k, v in coef.items()}
    for b0 in range(0, nr, block):
        b1 = min(nr, b0 + block)
        Psi = np.empty((b1 - b0, C.n))
        for i in range(b0, b1):
            f = rfac[i]
            if not f:
                Psi[i - b0] = 1.0; continue
            pr = np.zeros(C.n, bool); zz = np.zeros(C.n, bool)
            for p in f: pr ^= par[p]; zz |= zer[p]
            Psi[i - b0] = np.where(zz, 0.0, np.where(pr, -1.0, 1.0))
        for k in coef:
            S[k][b0:b1] = Psi @ Mcol[k]
    print(f"N0={N0} cols={C.n} rows={nr} marks={len(marks)} ({time.time()-t0:.1f}s)", flush=True)
    return dict(S=S, rN=rN, marks=marks, pnorm=pnorm, fields=fields, g1=g1, g2=g2, al=al, ncols=C.n)


def cubic_exp(F, x):
    """exponent of (x/pi)_3 in Z/3, from the sextic exponent (x/pi)_3 = (x/pi)_6^2."""
    k = F.sextic_exp(x)
    return None if k is None else (2 * k) % 3


def spectral(G, c_named, rng, nrand=200, S=None):
    """statistics of PSD Gram G against weight vectors (S = the K' x L' matrix, for single-row diagnostics)."""
    D = np.real(np.diag(G)); Db = D.mean()
    Gh = G / np.sqrt(np.outer(D, D))
    Lp = G.shape[0]
    off = Gh[~np.eye(Lp, dtype=bool)]
    ev, U = np.linalg.eigh(Gh)
    mu = np.linalg.eigvalsh(Gh - np.eye(Lp))
    mu2 = np.sort(mu ** 2)[::-1]
    fro = mu2.sum()
    out = dict(D_cv=float(D.std() / Db), off_rms=float(np.sqrt(np.mean(np.abs(off) ** 2))),
               off_max=float(np.abs(off).max()), lam_max=float(ev[-1]), lam_min=float(ev[0]),
               top1_share=float(mu2[0] / fro), top3_share=float(mu2[:3].sum() / fro),
               top10pct_share=float(mu2[:max(1, Lp // 10)].sum() / fro),
               stable_rank=float(fro / mu2[0]), Lp=Lp)
    u1 = U[:, -1]
    out["u1_ipr"] = float(np.sum(np.abs(u1) ** 4))          # 1/Lp if delocalised, ~1 if localised
    R = {}
    Kp_rows = S.shape[0] if S is not None else None
    Ginv = np.linalg.inv(G) if (Kp_rows is None or Kp_rows > 2 * Lp) else None   # G singular/ill-posed if K' <~ L'
    if S is not None:
        rn = np.sum(np.abs(S) ** 2, axis=1)
        out["row_norm_max_over_mean"] = float(rn.max() / rn.mean())
        Sh = S / np.sqrt(rn)[:, None]
        # does one quadratic character (one row s_k) carry the top eigen-direction?  (|<s_k/|s_k|, u1>|^2, D-scaled)
        u1D = U[:, -1] / np.sqrt(D); u1D /= np.linalg.norm(u1D)
        out["u1_best_single_row_overlap"] = float(np.max(np.abs(np.conj(Sh) @ u1D) ** 2))
    for name, c in c_named.items():
        c = c / np.linalg.norm(c)
        Rc = float(np.real(np.conj(c) @ G @ c) / Db)
        ovl = float(np.abs(np.conj(U[:, -1]) @ (c * np.sqrt(D) / np.linalg.norm(c * np.sqrt(D)))) ** 2)
        # optimal weights c* = G^{-1} loc: quality loc^* G^{-1} loc; gain over the coherent choice c = loc
        if Ginv is not None:
            q_opt = float(np.real(np.conj(c) @ Ginv @ c)) * Db
            cs = Ginv @ c
            R[name] = dict(R=Rc, overlap_u1=ovl, gain_opt=q_opt * Rc,
                           cstar_flat=float(np.abs(cs).max() / np.sqrt(np.mean(np.abs(cs) ** 2))))
        else:
            R[name] = dict(R=Rc, overlap_u1=ovl, gain_opt=float("nan"), cstar_flat=float("nan"))
    rr = []
    for _ in range(nrand):
        c = np.exp(2j * np.pi * rng.random(Lp)) / np.sqrt(Lp)
        rr.append(float(np.real(np.conj(c) @ G @ c) / Db))
    rr = np.array(rr)
    R["rand"] = dict(mean=float(rr.mean()), std=float(rr.std()), max=float(rr.max()))
    # best flat-modulus vector found from the top eigenvector's phases (sup over the torus, heuristic)
    cflat = np.exp(1j * np.angle(U[:, -1] / np.sqrt(D))) / np.sqrt(Lp)
    R["flat_top"] = dict(R=float(np.real(np.conj(cflat) @ G @ cflat) / Db))
    R["sup"] = dict(R=float(np.linalg.eigvalsh(G)[-1] / Db))
    out["R"] = R
    return out


def gaussian_control(Kp, Lp, rng):
    return (rng.standard_normal((Kp, Lp)) + 1j * rng.standard_normal((Kp, Lp))) / math.sqrt(2)


def design_blocks(G, labels):
    """rms of normalised off-diagonal entries within equal labels vs across labels."""
    D = np.real(np.diag(G)); Gh = np.abs(G / np.sqrt(np.outer(D, D))) ** 2
    lab = np.asarray(labels)
    same = (lab[:, None] == lab[None, :]) & ~np.eye(len(lab), dtype=bool)
    diff = lab[:, None] != lab[None, :]
    return dict(within=float(np.sqrt(Gh[same].mean())) if same.any() else None,
                across=float(np.sqrt(Gh[diff].mean())) if diff.any() else None,
                n_within=int(same.sum()), n_across=int(diff.sum()))


def run(N0=20000, KMAX=None, L0=800, fracs=(1 / 16, 1 / 4, 1), out=None):
    KMAX = KMAX or N0
    B = build(N0, KMAX, L0)
    rng = np.random.default_rng(11)
    marks, fields = B["marks"], B["fields"]
    Lp = len(marks)
    cvec = {"p78": np.conj(B["al"]) ** 2 * B["g2"], "ones": np.ones(Lp, complex), "kub": np.conj(B["g1"])}
    res = dict(N0=N0, KMAX=KMAX, L0=L0, cols=B["ncols"], marks=Lp, windows=[])
    for f in fracs:
        i = int(np.searchsorted(B["rN"], f * N0, side="right"))
        Kp = i
        win = dict(K_over_N0=f, rows=Kp)
        for key in ("theta", "randcol"):
            S = B["S"][key][:i]
            G = S.conj().T @ S
            win[key] = spectral(G, cvec, rng, S=S)
        Sg = gaussian_control(Kp, Lp, rng)
        win["gauss"] = spectral(Sg.conj().T @ Sg, cvec, rng, S=Sg)
        g = (1 + math.sqrt(Lp / Kp)) ** 2
        win["MP_edges"] = [(1 - math.sqrt(Lp / Kp)) ** 2, g]
        res["windows"].append(win)
        t, r_, gg = win["theta"], win["randcol"], win["gauss"]
        print(f"K/N0={f:<7g} rows={Kp:6d} L'={Lp}  MP=[{win['MP_edges'][0]:.3f},{g:.3f}]")
        for nm, s in (("theta", t), ("randcol", r_), ("gauss", gg)):
            print(f"   {nm:8s} lam=[{s['lam_min']:.3f},{s['lam_max']:.3f}] off_rms={s['off_rms']:.4f} "
                  f"top1={s['top1_share']:.3f} top3={s['top3_share']:.3f} top10%={s['top10pct_share']:.3f} "
                  f"srank={s['stable_rank']:.1f} ipr(u1)={s['u1_ipr']:.3f} Dcv={s['D_cv']:.3f} "
                  f"rowmax/mean={s.get('row_norm_max_over_mean', float('nan')):.2f} u1~row={s.get('u1_best_single_row_overlap', float('nan')):.3f}")
            R = s["R"]
            print(f"      R: p78={R['p78']['R']:.3f} (ovl_u1 {R['p78']['overlap_u1']:.3f}, opt-gain {R['p78']['gain_opt']:.3f},"
                  f" c*flat {R['p78']['cstar_flat']:.2f})  ones={R['ones']['R']:.3f}  kub={R['kub']['R']:.3f}  "
                  f"rand={R['rand']['mean']:.3f}+-{R['rand']['std']:.3f} (max {R['rand']['max']:.3f})  "
                  f"flat_top={R['flat_top']['R']:.3f}  sup={R['sup']['R']:.3f}", flush=True)
    # ---- designs on the largest window
    i = int(np.searchsorted(B["rN"], fracs[-1] * N0, side="right"))
    G = B["S"]["theta"][:i].conj().T @ B["S"]["theta"][:i]
    Gr = B["S"]["randcol"][:i].conj().T @ B["S"]["randcol"][:i]
    des = {}
    # (a) cubic class relative to fixed small primes pi0 (ray class by cubic reciprocity)
    for pi0n in (7, 13):
        pi0 = [p for p in fields if B["pnorm"][p] == pi0n][0]
        lab = [cubic_exp(fields[P], pi0) for P in marks]
        des[f"class_(pi0/P)_3,N(pi0)={pi0n}"] = dict(theta=design_blocks(G, lab), randcol=design_blocks(Gr, lab),
                                                     sizes=[int(np.sum(np.array(lab) == j)) for j in range(3)])
    # (b) conjugate pairs
    idx = {P: j for j, P in enumerate(marks)}
    conj_pairs = []
    for P, j in idx.items():
        a, b = P; Pc = None
        for Q in marks:
            if B["pnorm"][Q] == B["pnorm"][P] and Q != P: Pc = Q
        if Pc is not None and idx[Pc] > j: conj_pairs.append((j, idx[Pc]))
    D = np.real(np.diag(G)); Gh2 = np.abs(G / np.sqrt(np.outer(D, D))) ** 2
    Dr = np.real(np.diag(Gr)); Grh2 = np.abs(Gr / np.sqrt(np.outer(Dr, Dr))) ** 2
    offmask = ~np.eye(Lp, dtype=bool)
    des["conjugate_pairs"] = dict(n=len(conj_pairs),
                                  theta_pair_rms=float(np.sqrt(np.mean([Gh2[a, b] for a, b in conj_pairs]))),
                                  theta_all_rms=float(np.sqrt(Gh2[offmask].mean())),
                                  rand_pair_rms=float(np.sqrt(np.mean([Grh2[a, b] for a, b in conj_pairs]))),
                                  rand_all_rms=float(np.sqrt(Grh2[offmask].mean())))
    # (c) mutual cubic symbol (P1/P2)_3
    sym = np.full((Lp, Lp), -1)
    for a in range(Lp):
        for b in range(Lp):
            if a != b:
                sym[a, b] = cubic_exp(fields[marks[b]], marks[a])
    one = (sym == 0); notone = (sym > 0)
    des["mutual_symbol"] = dict(theta_sym1=float(np.sqrt(Gh2[one].mean())), theta_symne1=float(np.sqrt(Gh2[notone].mean())),
                                rand_sym1=float(np.sqrt(Grh2[one].mean())), rand_symne1=float(np.sqrt(Grh2[notone].mean())),
                                frac_sym1=float(one.sum() / (one.sum() + notone.sum())))
    # greedy cubic cliques (pairwise (P1/P2)_3 = 1), many random orders
    best = []
    for t in range(300):
        order = rng.permutation(Lp); cl = []
        for a in order:
            if all(sym[a, b] == 0 for b in cl): cl.append(a)
        if len(cl) > len(best): best = cl
    cl = best
    if len(cl) >= 2:
        sub = np.ix_(cl, cl); Gc = G[sub]; Dc = np.real(np.diag(Gc))
        cc = (np.conj(B["al"]) ** 2 * B["g2"])[cl]; cc = cc / np.linalg.norm(cc)
        des["greedy_cubic_clique"] = dict(size=len(cl), Lp=Lp, ratio_log=float(math.log(len(cl)) / math.log(Lp)),
                                          R_p78=float(np.real(np.conj(cc) @ Gc @ cc) / Dc.mean()))
    res["designs"] = des
    print(json.dumps(des, indent=1), flush=True)
    if out:
        json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    run(**args)
