"""
Per-pair test of the proposed sufficient statement (Lemma M(w), EDGE_THETA.md s.5), Moebius side:

  T(P1,P2) = sum_{h ~ H} W (h/P1)_6 conj(h/P2)_6 |M_h|^2,   M_h = sum_n a(n) w (h/n)_6,   P1 != P2 slot primes,
  minus the explicit coincidence term C(P1,P2) = sum_h W sum_m (a w)(m P2) conj((a w)(m P1)) 1_{(h, m P1 P2)=1}.
Lemma M(w) asks |T - C| << Z^{2+eps} (count units: #cols^2), i.e. relative to the P1 = P2 diagonal
T(P,P) ~ #h * #cols a saving #cols/#h = 1/(F_m L).  Random model: |T - C| ~ sqrt(#h) * #cols.
We report r = |T - C| / (#cols * sqrt(#h)) (random-model units) and s = |T - C| / #cols^2 (need s << Z^eps),
for the theta partner a = mu conj G and random-phase controls.
"""
import sys, json, math, time
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
import dual_mobius as dm
from edge_lib import setup, Columns
from eisenstein import mul

SCR = "/tmp/claude-0/-home-user-qrh/83fc0a9f-e093-5c20-9a33-d92b5adf4732/scratchpad"


def main(N0=1000, L0=31, Fmax=32, seeds=4, out=None):
    dump = f"{SCR}/pairs_N{N0}_L{L0}.npz"
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        dm.run(N0=N0, L0=L0, Fmax=Fmax, seeds=seeds, dump=dump)
    Z = np.load(dump)
    hN, M2, Xm, names = Z["hN"], Z["M2"], Z["Xm"], list(Z["names"])
    fields, pnorm, plist = setup(max(2 * N0, 2 * L0, 50))
    C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
    marks = [p for p in plist if fields[p].kind == "split" and L0 <= pnorm[p] < 2 * L0]
    Lp = len(marks)
    rng = np.random.default_rng(12345)
    # rebuild coefficient vectors exactly as dual_mobius (same rng order)
    coefs = [C.amob, C.amob * np.conj(C.alpha), C.mu.astype(complex)] + [np.exp(2j * np.pi * rng.random(C.n)) for _ in range(seeds)]
    aw = np.stack(coefs, axis=1) * C.w[:, None]
    # coincidences n_i = m p2, n_i2 = m p1 :  contributes to T(p1, p2) the term aw_i conj(aw_i2) 1_{(h, n_i p1) = 1}
    key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(C.A, C.B))}
    mset = {p: j for j, p in enumerate(marks)}
    coinc = []
    for i, f in enumerate(C.fac):
        for p2 in f:
            if p2 not in mset: continue
            rest = [p for p in f if p != p2]
            m = (1, 0)
            for p in rest: m = mul(m, p)
            for p1 in marks:
                if p1 == p2 or p1 in rest: continue
                i2 = key.get(mul(m, p1))
                if i2 is not None:
                    coinc.append((i, i2, mset[p1], mset[p2]))
    # coprimality of h with n_i p1: need column zero patterns per h -> recompute cheaply from symbols
    from edge_lib import Symb
    symb = {p: Symb(fields[p]) for p in C.primes}
    hA, hB = Z["hA"], Z["hB"]
    res = []
    Fs = [f for f in (1, 2, 4, 8, 16, 32, 64, 128) if f <= Fmax]
    for F in Fs:
        H = F * N0 * L0
        x = hN / H
        Wt = np.where((x > 0.5) & (x <= 1), np.sin(2 * np.pi * (x - 0.5)) ** 2, 0.0)
        sel = Wt > 0
        Wt_s = Wt[sel]; Xs = Xm[sel].astype(np.complex128); M2s = M2[sel]
        nfam = Wt_s.sum()
        # T[j][P1,P2] = sum_h W chi_h(P1) conj chi_h(P2) |M_h|^2
        T = np.einsum("h,ha,hb,hj->jab", Wt_s, Xs, np.conj(Xs), M2s)
        # coincidence terms
        Cc = np.zeros_like(T)
        if coinc:
            A_, B_ = hA[sel], hB[sel]
            for (i, i2, j1, j2) in coinc:
                z = np.ones(sel.sum(), bool)
                for p in C.fac[i]:
                    z &= symb[p](A_, B_) >= 0
                z &= symb[marks[j1]](A_, B_) >= 0 if marks[j1] in symb else (Xs[:, j1] != 0)
                Cc[:, j1, j2] += aw[i] * np.conj(aw[i2]) * np.sum(Wt_s * z)
        off = ~np.eye(Lp, dtype=bool)
        R = np.abs(T - Cc)[:, off]                       # coefficient x pairs
        diagT = np.real(np.einsum("jaa->ja", T)).mean(axis=1)
        unit_rand = C.n * math.sqrt(nfam)                 # random-model size of |T - C|
        r = np.sqrt(np.mean(R ** 2, axis=1)) / unit_rand
        s = np.sqrt(np.mean(R ** 2, axis=1)) / C.n ** 2
        coinc_size = np.sqrt(np.mean(np.abs(Cc[:, off]) ** 2, axis=1)) / C.n ** 2
        row = dict(F=F, H=H, nfam=float(nfam), Fm_counts=float(nfam / (C.n * Lp)),
                   r={names[j]: float(r[j]) for j in range(len(names))},
                   s={names[j]: float(s[j]) for j in range(len(names))},
                   rel_to_diag={names[j]: float(np.sqrt(np.mean(R[j] ** 2)) / diagT[j]) for j in range(len(names))},
                   coinc_over_cols2={names[j]: float(coinc_size[j]) for j in range(len(names))})
        res.append(row)
        rr = np.sqrt(np.mean([row['r'][n] ** 2 for n in names if n.startswith('random')]))
        print(f"F={F:<4} #h~{nfam:<8.0f} |T-C| in random units: theta {row['r']['thetaP']:.3f}  g2 {row['r']['g2P']:.3f}  mu {row['r']['mu']:.3f}  random {rr:.3f} | "
              f"|T-C|/cols^2: theta {row['s']['thetaP']:.3f} (coincidence term/cols^2 {row['coinc_over_cols2']['thetaP']:.2f}) | rel. to diag: theta {row['rel_to_diag']['thetaP']:.4f} need {C.n/nfam:.4f}", flush=True)
    if out:
        json.dump(dict(N0=N0, L0=L0, marks=Lp, cols=C.n, coincidences=len(coinc), res=res), open(out, "w"), indent=1)


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    main(**args)
