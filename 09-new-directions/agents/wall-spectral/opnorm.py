"""
Route 4 numerics: W(w) must hold uniformly in the slot coefficients a(P) (P78's Lemma reflected-energy
is stated for arbitrary bounded product coefficients).  For marks X(n) = sum_P c_P chi_P(n)^{-2},
   E(c) = sum_k |sum_P c_P S_{k,P}|^2,   S_{k,P} = sum_n theta(n) w(n) psi_k(n) chi_P(n)^{-2},
so  sup_c E(c)/||c||^2 = sigma_max(S)^2.  W for all marks  <=>  sigma_max^2 <~ K' * D * Z^eps,
D = mean |S_{k,P}|^2 (diagonal size).  Random-matrix benchmark: sigma_max^2 ~ D (sqrt K' + sqrt L')^2.
We also record max_P of the per-P energy (single-mark wall, NOTES s.9) and max_k row norms.
"""
import sys, math, time, json
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/energy-wall")
import dual_energy as de
from eisenstein import eisenstein_primes, ResidueField, mul, cx, Z6

def build(N0, KMAX, L0, seed=1):
    t0 = time.time()
    hi = max(2 * N0, KMAX, 2 * L0)
    primes = eisenstein_primes(hi + 10)
    fields = {p[1]: ResidueField(*p) for p in primes}
    pnorm = {p[1]: (p[2] if p[0] == "split" else p[2] ** 2) for p in primes}
    plist = sorted(fields, key=lambda p: pnorm[p])
    cA, cB, cN, cfac = de.primary_sf_elements(N0, 2 * N0, fields, plist, pnorm)
    rA, rB, rN, rfac = de.primary_sf_elements(0, KMAX, fields, plist, pnorm)
    o = np.argsort(rN, kind="stable"); rN = rN[o]; rfac = [rfac[i] for i in o]
    need = sorted({p for f in cfac for p in f}, key=lambda p: pnorm[p])
    g2 = {p: de.gauss2_vec(fields[p]) for p in need}
    theta = np.empty(len(cA), complex)
    for i, f in enumerate(cfac):
        g = 1 + 0j; acc = (1, 0)
        for p in f:
            k = fields[p].sextic_exp(acc) if acc != (1, 0) else 0
            g *= g2[p] * Z6[(4 * k) % 6]; acc = mul(acc, p)
        theta[i] = g * np.conj(complex(cA[i], 0) + cB[i] * complex(-0.5, math.sqrt(3) / 2)) / math.sqrt(cN[i])
    w = np.sin(np.pi * (cN / N0 - 1)) ** 2
    rowprimes = sorted({p for f in rfac for p in f}, key=lambda p: pnorm[p])
    par = {}; zer = {}
    for p in rowprimes:
        e = de.sextic_exp_vec(fields[p], cA, cB)
        zer[p] = e < 0; par[p] = (e % 2 == 1) & (e >= 0)
    Psi = np.empty((len(rfac), len(cA)), np.float32)
    for i, f in enumerate(rfac):
        pr_ = np.zeros(len(cA), bool); zz = np.zeros(len(cA), bool)
        for p in f: pr_ ^= par[p]; zz |= zer[p]
        Psi[i] = np.where(zz, 0.0, np.where(pr_, -1.0, 1.0))
    marks = [p for p in plist if fields[p].kind == "split" and L0 <= pnorm[p] < 2 * L0]
    Cm = np.empty((len(cA), len(marks)), np.complex64)
    for j, p in enumerate(marks):
        e = de.sextic_exp_vec(fields[p], cA, cB).astype(np.int64)
        Cm[:, j] = np.where(e >= 0, np.exp(-2j * np.pi * 2 * e / 6), 0)
    print(f"cols {len(cA)} rows {len(rfac)} marks {len(marks)}  {time.time()-t0:.1f}s", flush=True)
    return dict(theta=theta, w=w, Psi=Psi, Cm=Cm, rN=rN, cN=cN, rfac=rfac, marks=marks)

def analyse(B, coef, fracs, N0, label):
    Psi, Cm, rN = B["Psi"], B["Cm"], B["rN"]
    out = []
    M = (coef * B["w"])[:, None] * Cm                       # ncols x L'
    Sfull = Psi.astype(np.complex64) @ M.astype(np.complex64)   # rows x L'
    for f in fracs:
        i = np.searchsorted(rN, f * N0, side="right")
        S = Sfull[:i]
        D = np.mean(np.abs(S) ** 2)
        sv = np.linalg.svd(S, compute_uv=False)
        Kp, Lp = S.shape
        rm = (math.sqrt(Kp) + math.sqrt(Lp)) ** 2
        perP = np.sum(np.abs(S) ** 2, axis=0) / (Kp * D)
        out.append(dict(K_over_N0=f, rows=Kp, marks=Lp, smax2_over_KD=float(sv[0] ** 2 / (Kp * D)),
                        rmt_over_K=float(rm / Kp), s2_over_KD=float(sv[1] ** 2 / (Kp * D)),
                        max_perP=float(perP.max()), max_row=float(np.max(np.sum(np.abs(S) ** 2, axis=1)) / (Lp * D))))
        r = out[-1]
        print(f"{label:10s} K/N0={f:<7g} rows={Kp:6d} L'={Lp:4d}  s_max^2/(K'D)={r['smax2_over_KD']:.3f} "
              f"(RMT {r['rmt_over_K']:.3f})  s_2^2/(K'D)={r['s2_over_KD']:.3f}  max per-P={r['max_perP']:.3f}  "
              f"max row/(L'D)={r['max_row']:.2f}", flush=True)
    return out

if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    N0 = args.get("N0", 20000); KMAX = args.get("KMAX", N0); L0 = args.get("L0", 500)
    B = build(N0, KMAX, L0)
    fracs = [1 / 64, 1 / 16, 1 / 4, 1 / 2, 1][:: 1]
    fracs = [f for f in fracs if f * N0 <= KMAX]
    rng = np.random.default_rng(2)
    res = {}
    res["theta"] = analyse(B, B["theta"], fracs, N0, "theta")
    res["random"] = analyse(B, np.abs(B["theta"]) * np.exp(2j * np.pi * rng.random(len(B["theta"]))), fracs, N0, "random")
    # adversary: theta coefficients replaced by |theta| * psi_{k0}(n) * chi_{P0}(n)^2  (aligned with one (k,P) cell)
    i0 = np.searchsorted(B["rN"], KMAX // 8)
    adv = np.abs(B["theta"]) * B["Psi"][i0] * np.conj(B["Cm"][:, 0])
    res["adv"] = analyse(B, adv, fracs, N0, "adversary")
    json.dump(dict(N0=N0, KMAX=KMAX, L0=L0, res=res), open(args.get("out", f"opnorm_N{N0}_L{L0}.json"), "w"), indent=1)
