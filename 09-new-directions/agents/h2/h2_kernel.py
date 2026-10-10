"""
h2, idea 1 (a): can a mixed second moment be bounded by the off-diagonals of the two factors (instead of their full
diagonals, as Cauchy-Schwarz does)?  Exact Z[w] test of the obstruction.

  Gamma = Psi^T diag(V) Psi   (N x N, n-space Gram of the K' quadratic rows psi_k in the window)
  Q(x)  = <x, Gamma x> - <x, Dv x>   (off-diagonal form; Dv(n) = sum_k V 1_{(n,k)=1} is the exact diagonal)
Beyond the wall (K' < N) Gamma has a kernel of dimension >= N - K', on which Q = -Dv ~ -K: Q is indefinite with
"negative mass" comparable to the full diagonal, so no Cauchy-Schwarz for Q.

Certificate: starting from an actual mark column x1 = beta w chi_P^{-2}, alternating projections produce a FLAT x2
(|x2(n)| = |x1(n)| for all n) with
   Psi_V x2 = Psi_V x1   (identical row sums: D_2(k) = D_1(k) for every row k)   and   sum_n Dv conj(x1) x2 = 0.
Then O_22 = O_11, mixed diagonal Delta_12 = 0, but O_12 = G_11 (the full energy).  So no inequality
|O_12| <= F(|x1|, |x2|, O_11, O_22, Delta_12) with F(.., 0) << diagonal can exist; arithmetic of chi_P must enter.
"""
import sys, json, time, math
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/h2")
from h2_struct import build


def run(N0=10000, L0=1000, jK=-2, iters=3000, out=None, seed=3):
    t0 = time.time()
    B = build(N0, N0, (L0,), 8, seed, keep_psi2=True)
    Ks = B["Ks"]; V = B["Vw"][:, jK]; K = Ks[jK]
    rows = np.nonzero(V > 0)[0]
    # rebuild the row characters on the window (Psi^2 kept; signs from SP structure are not stored) -> recompute Psi
    sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
    from edge_lib import setup, Symb, Columns
    import dual_energy as de
    fields, pnorm, plist = setup(max(2 * N0, 2 * L0) + 10)
    C = Columns(N0, 2 * N0, fields, pnorm, plist, None)
    rA, rB, rN, rfac = de.primary_sf_elements(0, N0, fields, plist, pnorm)
    odd = np.array([all(pnorm[p] != 4 for p in f) for f in rfac], bool)
    keep = np.nonzero(odd)[0]; rN = rN[keep]; rfac = [rfac[i] for i in keep]
    o = np.argsort(rN, kind="stable"); rN = rN[o]; rfac = [rfac[i] for i in o]
    assert np.array_equal(rN, B["rN"])
    par, zer = {}, {}
    for p in sorted({p for i in rows for p in rfac[i]}, key=lambda p: pnorm[p]):
        e = Symb(fields[p])(C.A, C.B); zer[p] = e < 0; par[p] = (e % 2 == 1)
    Psi = np.empty((len(rows), C.n))
    for t, i in enumerate(rows):
        pr = np.zeros(C.n, bool); zz = np.zeros(C.n, bool)
        for p in rfac[i]: pr ^= par[p]; zz |= zer[p]
        Psi[t] = np.where(zz, 0.0, np.where(pr, -1.0, 1.0))
    sv = np.sqrt(V[rows])
    A = sv[:, None] * Psi                                    # K' x N, E(x) = |A x|^2
    Dv = (Psi ** 2 * V[rows][:, None]).sum(0)                # exact diagonal weights
    Nn = C.n; Kp = len(rows)
    print(f"N0={N0} cols N={Nn}, window K={K}: K'={Kp} rows, mean diag weight {Dv.mean():.2f}  {time.time()-t0:.1f}s", flush=True)
    # spectrum of Gamma and of Q = Gamma - diag(Dv)
    Gam = A.T @ A
    ev = np.linalg.eigvalsh(Gam)
    tol = 1e-8 * ev[-1]
    nker = int((ev < tol).sum())
    evQ = np.linalg.eigvalsh(Gam - np.diag(Dv))
    x1 = B["CP"][:, 0]
    w = np.abs(x1)
    # mass of mark vectors in ker Gamma
    U, s, Vt = np.linalg.svd(A, full_matrices=False)          # row space of A = range(Gamma)
    def kerfrac(x):
        pr = Vt @ x
        return float(1 - np.vdot(pr, pr).real / np.vdot(x, x).real)
    kf = [kerfrac(B["CP"][:, j]) for j in range(B["CP"].shape[1])]
    res = dict(N0=N0, N=Nn, K=K, Kp=Kp, mean_diag=float(Dv.mean()), nker=nker, N_minus_Kp=Nn - Kp,
               Q_min=float(evQ[0]), Q_max=float(evQ[-1]), Q_frac_negative=float((evQ < 0).mean()),
               Q_median_negative=float(np.median(evQ[evQ < 0])), Gamma_top=float(ev[-1]),
               mark_mass_in_kernel=kf, predicted_kernel_mass=1 - Kp / Nn)
    print(f"Gamma: rank-deficiency {nker} (N-K' = {Nn-Kp}); Q eigenvalues in [{evQ[0]:.1f}, {evQ[-1]:.1f}], "
          f"{100*(evQ<0).mean():.1f}% negative, median negative {np.median(evQ[evQ<0]):.1f} (mean diag {Dv.mean():.1f}); "
          f"mark mass in ker: {np.mean(kf):.3f} (1-K'/N = {1-Kp/Nn:.3f})", flush=True)
    # flat counterexample by alternating projections
    M = np.vstack([A.astype(complex), (Dv * np.conj(x1))[None, :]])
    b = np.concatenate([A @ x1, [0.0]])
    MMH = M @ M.conj().T
    L = np.linalg.cholesky(MMH)
    def projB(y):
        r = M @ y - b
        z = np.linalg.solve(L.conj().T, np.linalg.solve(L, r))
        return y - M.conj().T @ z
    rng = np.random.default_rng(seed)
    y = w * np.exp(2j * np.pi * rng.random(Nn))
    hist = []
    for it in range(iters):
        y = projB(y)
        dev = np.sqrt(np.mean((np.abs(y) - w) ** 2) / np.mean(w ** 2))
        y = w * np.exp(1j * np.angle(y))
        if it % 500 == 0 or it == iters - 1:
            cres = np.linalg.norm(M @ y - b) / np.linalg.norm(b)
            hist.append((it, float(dev), float(cres)))
            print(f"  iter {it}: modulus deviation (before re-flattening) {dev:.2e}, constraint residual (flat y) {cres:.2e}", flush=True)
    x2 = projB(y)                       # satisfies constraints exactly; nearly flat
    moddev = float(np.max(np.abs(np.abs(x2) - w) / np.maximum(w, 1e-12)))
    rmsdev = float(np.sqrt(np.mean((np.abs(x2) - w) ** 2) / np.mean(w ** 2)))
    G11 = float(np.linalg.norm(A @ x1) ** 2); G22 = float(np.linalg.norm(A @ x2) ** 2)
    G12 = complex(np.vdot(A @ x1, A @ x2))
    D11 = float(np.sum(Dv * np.abs(x1) ** 2)); D22 = float(np.sum(Dv * np.abs(x2) ** 2))
    D12 = complex(np.sum(Dv * np.conj(x1) * x2))
    O11, O22, O12 = G11 - D11, G22 - D22, G12 - D12
    cert = dict(max_rel_modulus_dev=moddev, rms_rel_modulus_dev=rmsdev, G11=G11, G22=G22, G12=[G12.real, G12.imag],
                Delta11=D11, Delta22=D22, Delta12_abs=abs(D12), O11=O11, O22=O22, O12_abs=abs(O12),
                ratio_O12_over_sqrt_O11O22=abs(O12) / math.sqrt(abs(O11 * O22)), ratio_O12_over_diag=abs(O12) / D11, hist=hist)
    res["flat_counterexample"] = cert
    print(f"certificate: |x2|/|x1| deviation max {moddev:.2e} rms {rmsdev:.2e}; G11 {G11:.1f} G22 {G22:.1f} |G12| {abs(G12):.1f}; "
          f"Delta11 {D11:.1f} |Delta12| {abs(D12):.2e}; O11 {O11:.1f} O22 {O22:.1f} |O12| {abs(O12):.1f}; "
          f"|O12|/sqrt|O11 O22| = {cert['ratio_O12_over_sqrt_O11O22']:.1f}, |O12|/Delta11 = {cert['ratio_O12_over_diag']:.3f}", flush=True)
    if out: json.dump(res, open(out, "w"), indent=1)
    return res


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    run(**args)
