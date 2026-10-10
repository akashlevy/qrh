"""
Truth check, primal side: reflected (quadratic-row) marked energy beyond the edge, exact diagonal.

  S_k(c) = sum_{n sf odd primary, N0 < N(n) <= 2N0} c(n) w(N(n)/N0) psi_k(n),  psi_k(n) = prod_{pi|k} (n/pi)_6^3,
  rows k sf odd primary, row window k ~ K: weight V(N(k)/K) = sin^2 bump on (1/2, 1].
  rho(K) = sum_k V (|S_k|^2 - d_k) / sum_k V d_k,   d_k = sum_n |c w|^2 1_{(n,k)=1}  (exact n-diagonal).

Coefficients: beta (theta), beta*X (X = sum_{P in slot} chi_P^{-2}, coherent / random slot coefficients),
random unit phases (controls).  Samples: rows split into 12 angular sectors of k (independent sub-families).
Also the cross-mark Gram (hypothesis (H2) of wall-spectral):
  G_{P1P2}(K) = sum_k V S_{k,P1} conj S_{k,P2},  S_{k,P} = S_k(beta chi_P^{-2}),
  normalised c_{P1P2} = G_{P1P2}/sqrt(G_{P1P1}G_{P2P2}); its RMS over P1 != P2 vs #rows.
W(w) needs rho << Z^eps on F = N/K <= Z^w; (H2) needs |c_{P1P2}| <~ L^{-0.93}.
"""
import sys, json, time, math
import numpy as np
from edge_lib import setup, Symb, Columns, gauss_vec


def run(N0=10000, KMAX=None, L0=200, nmark=None, seeds=6, sectors=12, block=2000, out=None):
    t0 = time.time()
    KMAX = KMAX or N0
    fields, pnorm, plist = setup(max(2 * N0, KMAX, 2 * L0))
    C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
    import dual_energy as de
    rA, rB, rN, rfac = de.primary_sf_elements(0, KMAX, fields, plist, pnorm)
    odd = np.array([all(pnorm[p] != 4 for p in f) for f in rfac], bool)
    keep = np.nonzero(odd)[0]
    rA, rB, rN = rA[keep], rB[keep], rN[keep]; rfac = [rfac[i] for i in keep]
    o = np.argsort(rN, kind="stable"); rA, rB, rN = rA[o], rB[o], rN[o]; rfac = [rfac[i] for i in o]
    nr = len(rN)
    print(f"N0={N0}: cols {C.n}, rows {nr} (N(k) <= {KMAX})  {time.time()-t0:.1f}s", flush=True)
    # mark primes
    marks = [p for p in plist if fields[p].kind == "split" and L0 <= pnorm[p] < 2 * L0]
    if nmark: marks = marks[:nmark]
    Em = np.stack([Symb(fields[p])(C.A, C.B) for p in marks], axis=1)                  # cols x L'
    chiP2 = np.where(Em < 0, 0, np.exp(-2j * np.pi * 2 * Em.astype(np.int64) / 6))     # chi_P(n)^{-2}
    rng = np.random.default_rng(7)
    X_coh = chiP2.sum(axis=1) / math.sqrt(len(marks))
    cr = np.exp(2j * np.pi * rng.random(len(marks)))
    X_rnd = (chiP2 * cr).sum(axis=1) / math.sqrt(len(marks))
    coefs = {"beta": C.beta, "beta*X(coh)": C.beta * X_coh, "beta*X(rnd)": C.beta * X_rnd}
    for s in range(seeds):
        coefs[f"rand{s}"] = np.exp(2j * np.pi * rng.random(C.n))
    for s in range(2):
        coefs[f"rand*|X(coh)|{s}"] = np.exp(2j * np.pi * rng.random(C.n)) * np.abs(X_coh)
    names = list(coefs)
    Cm = np.stack([coefs[k] * C.w for k in names], axis=1)                                   # cols x m
    Cabs = np.abs(Cm) ** 2
    # per-mark columns for the (H2) Gram: beta chi_P^{-2} and a random control
    CP = (C.beta * C.w)[:, None] * chiP2
    CPr = (np.exp(2j * np.pi * rng.random(C.n)) * C.w)[:, None] * chiP2
    # row symbols: exponent parity of (n/pi)_6 for each prime pi dividing some row
    rowprimes = sorted({p for f in rfac for p in f}, key=lambda p: pnorm[p])
    par = {}; zer = {}
    for p in rowprimes:
        e = Symb(fields[p])(C.A, C.B)
        zer[p] = e < 0; par[p] = (e % 2 == 1)
    print(f"row symbols for {len(rowprimes)} primes  {time.time()-t0:.1f}s", flush=True)
    S = np.zeros((nr, len(names)), complex); D = np.zeros((nr, len(names)))
    SP = np.zeros((nr, len(marks)), complex); SPr = np.zeros((nr, len(marks)), complex)
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
        S[b0:b1] = Psi @ Cm
        D[b0:b1] = (Psi ** 2) @ Cabs
        SP[b0:b1] = Psi @ CP
        SPr[b0:b1] = Psi @ CPr
    print(f"rows done {time.time()-t0:.1f}s", flush=True)
    sect = (np.floor((np.angle(rA + rB * complex(-0.5, math.sqrt(3) / 2)) + np.pi) / (2 * np.pi) * sectors).astype(int)) % sectors
    res = dict(N0=N0, cols=C.n, rows=nr, marks=len(marks), L0=L0, names=names, windows=[])
    Ks = [2 ** j for j in range(3, 40) if 2 ** j <= KMAX]
    for K in Ks:
        x = rN / K
        V = np.where((x > 0.5) & (x <= 1), np.sin(2 * np.pi * (x - 0.5)) ** 2, 0.0)
        V[rN == 1] = 0.0                                    # principal row recorded separately
        if V.sum() < 3: continue
        nrow = int(((x > 0.5) & (x <= 1)).sum())
        rho = ((V[:, None] * (np.abs(S) ** 2 - D)).sum(0) / (V[:, None] * D).sum(0))
        rho_sec = []
        for s in range(sectors):
            Vs = V * (sect == s)
            if Vs.sum() == 0: continue
            rho_sec.append(((Vs[:, None] * (np.abs(S) ** 2 - D)).sum(0) / (Vs[:, None] * D).sum(0)))
        rho_sec = np.array(rho_sec)
        # (H2) normalised cross-mark Gram
        def corr(SPm):
            G = (V[:, None, None] * SPm[:, :, None] * np.conj(SPm[:, None, :])).sum(0)
            d = np.sqrt(np.real(np.diag(G)))
            c = G / np.outer(d, d)
            off = c[~np.eye(len(d), dtype=bool)]
            return float(np.sqrt(np.mean(np.abs(off) ** 2))), float(np.max(np.abs(off)))
        h2, h2max = corr(SP); h2r, h2rmax = corr(SPr)
        win = dict(K=K, F=N0 / K, nrow=nrow, rho={n: float(rho[j]) for j, n in enumerate(names)},
                   rho_sector_rms={n: float(np.sqrt(np.mean(rho_sec[:, j] ** 2))) for j, n in enumerate(names)},
                   rho_sector_mean={n: float(np.mean(rho_sec[:, j])) for j, n in enumerate(names)},
                   H2_rms=h2, H2_max=h2max, H2_rand_rms=h2r, H2_rand_max=h2rmax)
        res["windows"].append(win)
        rr = np.array([win["rho_sector_rms"][n] for n in names if n.startswith("rand") and "|" not in n])
        print(f"K={K:<7d} F={N0/K:<9.4g} rows~{nrow:<6d} rho: " +
              "  ".join(f"{n}:{win['rho'][n]:+.4f}(sec rms {win['rho_sector_rms'][n]:.3f})" for n in names[:3]) +
              f" | rand sec rms {np.sqrt(np.mean(rr**2)):.3f} | H2 rms {h2:.4f} (rand {h2r:.4f}) max {h2max:.3f}", flush=True)
    S1 = np.abs(S[rN == 1]) ** 2 / D[rN == 1] if (rN == 1).any() else None
    if S1 is not None:
        res["principal_row_ratio"] = {n: float(S1[0, j]) for j, n in enumerate(names)}
        print("principal row |S_1|^2/d_1:", {n: round(float(S1[0, j]), 2) for j, n in enumerate(names[:3])}, flush=True)
    if out:
        json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    run(**args)
