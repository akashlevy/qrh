"""
Truth check, dual side: the Moebius-twisted sextic family in the diagonal regime (NOTES s.6, s.20).

  D(H) = sum_{h primary} W(N(h)/H) |Pi_h|^2 |M_h|^2,
  M_h  = sum_{n sf, odd, N0 < N(n) <= 2N0} a(n) w(N(n)/N0) (h/n)_6,
  Pi_h = sum_{P split prime, L0 <= N(P) < 2L0} c_P (h/P)_6      (marks; Pi = 1 if unmarked).

W = sin^2 bump on (1/2, 1] (family h ~ H).  Exact diagonal:
  Ddiag(H) = sum_h W * [ sum_n |a w|^2 1_{(h,n)=1} * sum_P |c_P|^2 1_{(h,P)=1} + cross(h) ],
  cross(h) = sum over n1 P1 = n2 P2, (n1,P1) != (n2,P2)  (exact coincidences of the product index).
Relative off-diagonal delta = D/Ddiag - 1.  Family/length ratio F = H / (N0 * L0) (L0 = 1 unmarked).

By Poisson in h (NOTES s.6), W(w) <=> |delta| << Z^eps / F on the marked family beyond the wall;
the unmarked family is provably within (theta reflection, M <= 1).  We record delta for
  a = mu conj G   (exact Poisson partner of beta = gamma_2 conj alpha: a gamma_1 = conj beta),
  a = mu conj(G alpha) (partner of gamma_2, which HAS the Patterson bias),
  a = mu, and random unit phases (control, several seeds).
"""
import sys, json, time, math
import numpy as np
from edge_lib import setup, Symb, Columns, primary_all, gauss_vec


def run(N0=1000, Fmax=256, L0=0, seeds=6, block=20000, out=None, dump=None):
    t0 = time.time()
    Lfac = L0 if L0 else 1
    Hmax = int(Fmax * N0 * Lfac)
    fields, pnorm, plist = setup(max(2 * N0, 2 * L0, 50))
    C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
    hA, hB, hN = primary_all(Hmax)
    nh = len(hA)
    print(f"N0={N0} L0={L0} cols={C.n} family h<= {Hmax}: {nh}  {time.time()-t0:.1f}s", flush=True)
    rng = np.random.default_rng(12345)
    coefs = {"thetaP": C.amob,                      # mu conj G: Poisson partner of beta (no bias)
             "g2P": C.amob * np.conj(C.alpha),       # mu conj(G alpha): partner of gamma_2 (Patterson-biased)
             "mu": C.mu.astype(complex)}
    for s in range(seeds):
        coefs[f"random{s}"] = np.exp(2j * np.pi * rng.random(C.n))
    names = list(coefs)
    Cm = np.stack([coefs[k] * C.w for k in names], axis=1).astype(np.complex128)     # n x m
    aw2 = (C.w ** 2).astype(np.float64)                                                # |a w|^2 (|a| = 1)
    # marks
    if L0:
        marks = [p for p in plist if fields[p].kind == "split" and L0 <= pnorm[p] < 2 * L0]
        g1 = np.array([gauss_vec(fields[p], 1) for p in marks])
        cP = {"coherent": np.ones(len(marks), complex),
              "random": np.exp(2j * np.pi * rng.random(len(marks))),
              "conj gamma1 (Kubota-aligned)": np.conj(g1)}
        msymb = [Symb(fields[p]) for p in marks]
        print(f"marks {len(marks)} in [{L0},{2*L0})", flush=True)
    else:
        marks = []; cP = {"none": np.ones(1, complex)}; msymb = []
    # symbols for primes dividing columns
    symb = {p: Symb(fields[p]) for p in C.primes}
    pidx = {p: i for i, p in enumerate(C.primes)}
    colfac = [[pidx[p] for p in f] for f in C.fac]
    # cross terms n1 P1 = n2 P2 (only marked): list of (i1, i2, j1, j2) with n_{i1} P_{j1} = n_{i2} P_{j2}
    cross = []
    if L0:
        key = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(C.A, C.B))}
        mset = {p: j for j, p in enumerate(marks)}
        from eisenstein import mul
        for i, f in enumerate(C.fac):
            for p2 in f:
                if p2 not in mset: continue
                # m = n / p2 ;  n_{i2} = m * p1 for p1 != p2 in marks
                rest = [p for p in f if p != p2]
                m = (1, 0)
                for p in rest: m = mul(m, p)
                for p1 in marks:
                    if p1 == p2 or p1 in rest: continue
                    e = mul(m, p1)
                    # make primary-consistent: product of primary elements is primary
                    i2 = key.get(e)
                    if i2 is not None:
                        cross.append((i, i2, mset[p2], mset[p1]))   # n_i P_{p1} = n_{i2} P_{p2}
        print(f"cross coincidences {len(cross)}", flush=True)
    M2 = np.zeros((nh, len(names)))
    XmAll = np.zeros((nh, len(marks)), np.complex64) if (L0 and dump) else None
    d1 = np.zeros(nh)
    PI2 = {k: np.zeros(nh) for k in cP}
    dP = np.zeros(nh)
    crossv = {k: np.zeros((nh, len(names)), complex) for k in cP} if L0 else None
    for s0 in range(0, nh, block):
        s1 = min(nh, s0 + block)
        A, B = hA[s0:s1], hB[s0:s1]
        E = np.zeros((s1 - s0, len(C.primes)), np.int8)
        for p, i in pidx.items():
            E[:, i] = symb[p](A, B)
        X = np.empty((s1 - s0, C.n), np.complex128)
        Z = np.empty((s1 - s0, C.n), bool)
        for j, f in enumerate(colfac):
            e = E[:, f]
            z = (e < 0).any(axis=1)
            k = e.astype(np.int64).sum(axis=1) % 6
            X[:, j] = np.where(z, 0, np.exp(2j * np.pi * k / 6))
            Z[:, j] = ~z
        S = X @ Cm
        M2[s0:s1] = np.abs(S) ** 2
        d1[s0:s1] = Z.astype(np.float64) @ aw2
        if L0:
            Em = np.stack([sy(A, B) for sy in msymb], axis=1)
            Xm = np.where(Em < 0, 0, np.exp(2j * np.pi * Em.astype(np.int64) / 6))
            dP[s0:s1] = (Em >= 0).sum(axis=1)
            if XmAll is not None: XmAll[s0:s1] = Xm
            for k, c in cP.items():
                PI2[k][s0:s1] = np.abs(Xm @ c) ** 2
            # cross(h) = sum over coincidences of b-products: a_{i}w_i c_{p1}... with index identity
            # n_i P_{p1} = n_{i2} P_{p2}:  term = (a w)_i c_{p1} conj((a w)_{i2} c_{p2}) * 1_{(h, n_i p1)=1}
            if cross:
                ci = np.array(cross)
                zz = Z[:, ci[:, 0]] & (Em[:, ci[:, 3]] >= 0)            # (h, n_i P_{p1}) = 1
                for k, c in cP.items():
                    wts = (Cm[ci[:, 0]] * c[ci[:, 3]][:, None]) * np.conj(Cm[ci[:, 1]] * c[ci[:, 2]][:, None])
                    crossv[k][s0:s1] = zz.astype(np.float64) @ wts
        print(f"  block {s0}-{s1} {time.time()-t0:.1f}s", flush=True)
    # windows
    sect = (np.floor((np.angle(hA + hB * complex(-0.5, math.sqrt(3) / 2)) + np.pi) / (2 * np.pi) * 12).astype(int)) % 12
    res = {}
    Fs = [f for f in (1 / 4, 1 / 2, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024) if f <= Fmax]
    for k in cP:
        weight = PI2[k] if L0 else np.ones(nh)
        rows = {}
        for F in Fs:
            H = F * N0 * Lfac
            x = hN / H
            Wt = np.where((x > 0.5) & (x <= 1), np.sin(2 * np.pi * (x - 0.5)) ** 2, 0.0)
            if Wt.sum() == 0: continue
            nfam = float(Wt.sum())
            if L0:
                # diag = sum_n |aw|^2 1_{(h,n)=1} * sum_P |c_P|^2 1_{(h,P)=1} + exact coincidences (|c_P| = 1)
                dg = d1 * dP
                Dd = np.array([np.sum(Wt * (dg + crossv[k][:, j].real)) for j in range(len(names))])
            else:
                Dd = np.full(len(names), np.sum(Wt * d1))
            Dt = (Wt[:, None] * weight[:, None] * M2).sum(axis=0)
            # sector samples: 12 angular sectors of h (independent sub-families)
            dsec = []
            for s_ in range(12):
                Ws = Wt * (sect == s_)
                if L0:
                    Dds = np.array([np.sum(Ws * (dg + crossv[k][:, j].real)) for j in range(len(names))])
                else:
                    Dds = np.full(len(names), np.sum(Ws * d1))
                dsec.append((Ws[:, None] * weight[:, None] * M2).sum(axis=0) / Dds - 1)
            dsec = np.array(dsec)
            ncolsP = C.n * (len(marks) if L0 else 1)
            rows[F] = dict(H=H, nfam=nfam, Fc=nfam / ncolsP,
                           delta={nm: float(Dt[j] / Dd[j] - 1) for j, nm in enumerate(names)},
                           sec_rms={nm: float(np.sqrt(np.mean(dsec[:, j] ** 2))) for j, nm in enumerate(names)},
                           sec_mean={nm: float(np.mean(dsec[:, j])) for j, nm in enumerate(names)})
        res[k] = rows
        print(f"-- marks: {k}", flush=True)
        for F, r in rows.items():
            dl = r["delta"]
            rnd = np.array([dl[n] for n in names if n.startswith("random")])
            print(f"F={F:<7g} H={r['H']:<9g} #h~{r['nfam']:<9.0f} " +
                  "  ".join(f"{n}:{dl[n]:+.4f}" for n in names if not n.startswith("random")) +
                  f"  random rms {np.sqrt(np.mean(rnd**2)):.4f}  1/Fc={1/r['Fc']:.4f} | sector rms theta {r['sec_rms']['thetaP']:.4f} rand {np.sqrt(np.mean([r['sec_rms'][n]**2 for n in names if n.startswith('random')])):.4f}", flush=True)
    if dump:
        np.savez(dump, hA=hA, hB=hB, hN=hN, M2=M2, d1=d1, dP=dP, Xm=XmAll if XmAll is not None else np.zeros(0), names=np.array(names), ncols=C.n)
    if out:
        json.dump(dict(N0=N0, L0=L0, cols=C.n, names=names, res={k: {str(F): v for F, v in r.items()} for k, r in res.items()}),
                  open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    run(**args)
