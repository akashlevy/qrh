"""Exponent fits for the truth-check numerics (primal_N*.json, dual_*.json) -> analysis.txt / fits.json."""
import json, glob, math
import numpy as np

out = []; fits = {}
def P(s=""): out.append(s); print(s)

def lsq(X, y):
    A = np.column_stack([np.ones(len(y))] + X)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c
    se = np.sqrt(np.diag(np.linalg.inv(A.T @ A)) * (r @ r) / max(1, len(y) - A.shape[1]))
    return c, se

# ---------------- primal (reflected quadratic-row energy) ----------------
rows = []
for f in sorted(glob.glob("primal_N*.json")):
    d = json.load(open(f))
    for w in d["windows"]:
        if w["nrow"] < 30: continue
        rnd = [w["rho_sector_rms"][n] for n in d["names"] if n.startswith("rand") and "|" not in n]
        rows.append(dict(N0=d["N0"], K=w["K"], F=w["F"], n=w["nrow"] / 12,
                         th=w["rho_sector_rms"]["beta"], thX=w["rho_sector_rms"]["beta*X(coh)"],
                         thXr=w["rho_sector_rms"]["beta*X(rnd)"], rnd=float(np.sqrt(np.mean(np.square(rnd)))),
                         h2=w["H2_rms"], h2r=w["H2_rand_rms"], nrow=w["nrow"], L0=d["L0"], marks=d["marks"]))
P("PRIMAL reflected energy E(K) = sum_k |sum_n c(n) w psi_k(n)|^2, relative off-diagonal rho per sector (12 sectors)")
P(f"{'N0':>6} {'K':>6} {'F=N0/K':>8} {'rows/sec':>8} {'beta':>7} {'betaX':>7} {'betaXr':>7} {'random':>7} {'H2':>7} {'H2rand':>7}")
for r in rows:
    P(f"{r['N0']:>6} {r['K']:>6} {r['F']:>8.3g} {r['n']:>8.0f} {r['th']:>7.3f} {r['thX']:>7.3f} {r['thXr']:>7.3f} {r['rnd']:>7.3f} {r['h2']:>7.4f} {r['h2r']:>7.4f}")
n = np.array([r["n"] for r in rows]); F = np.array([r["F"] for r in rows])
for key, lab in (("th", "beta"), ("thX", "beta*X coherent"), ("thXr", "beta*X random slots"), ("rnd", "random control")):
    y = np.log([r[key] for r in rows])
    c, se = lsq([np.log(n), np.log(F)], y)
    fits[f"primal_{key}"] = dict(const=c[0], exp_rows=c[1], se_rows=se[1], exp_F=c[2], se_F=se[2])
    P(f"fit log rho_rms = a + b log(rows) + c log F   [{lab:22s}]  b = {c[1]:+.3f} +- {se[1]:.3f}   c = {c[2]:+.3f} +- {se[2]:.3f}")
ratio = np.array([r["th"] / r["rnd"] for r in rows] + [r["thX"] / r["rnd"] for r in rows])
P(f"theta-type / random sector rms: mean {ratio.mean():.3f}, median {np.median(ratio):.3f}, max {ratio.max():.3f} over {len(ratio)} windows")
nn = np.array([r["nrow"] for r in rows])
for key, lab in (("h2", "(H2) beta chi_P^-2 cross-mark corr"), ("h2r", "(H2) random control")):
    y = np.log([r[key] for r in rows])
    c, se = lsq([np.log(nn), np.log(F)], y)
    fits[f"H2_{key}"] = dict(const=c[0], exp_rows=c[1], se_rows=se[1], exp_F=c[2], se_F=se[2])
    P(f"fit log H2_rms = a + b log(rows) + c log F   [{lab:36s}]  a = {c[0]:+.3f}  b = {c[1]:+.3f} +- {se[1]:.3f}   c = {c[2]:+.3f} +- {se[2]:.3f}")
P()

# ---------------- dual (Moebius-sextic family, diagonal regime) ----------------
P("DUAL Moebius-sextic family D(H) = sum_h |Pi_h|^2 |sum_n a(n) (h/n)_6|^2, relative off-diagonal delta per sector (12 sectors)")
P("need (= W(w) after Poisson): |delta| <~ Z^eps / Fc,  Fc = #family / #(product columns)")
drows = []
for f in sorted(glob.glob("dual_*.json")):
    d = json.load(open(f))
    for mk, win in d["res"].items():
        for Fs, r in win.items():
            if r["nfam"] < 100: continue
            rnd = [r["sec_rms"][nm] for nm in d["names"] if nm.startswith("random")]
            drows.append(dict(run=f.replace(".json", ""), marks=mk, F=float(Fs), nfam=r["nfam"], Fc=r["Fc"],
                              th=r["sec_rms"]["thetaP"], g2=r["sec_rms"]["g2P"], mu=r["sec_rms"]["mu"],
                              rnd=float(np.sqrt(np.mean(np.square(rnd)))), th_full=r["delta"]["thetaP"]))
P(f"{'run':>22} {'marks':>10} {'F':>6} {'#h':>8} {'1/Fc':>7} {'theta':>7} {'g2':>7} {'mu':>7} {'random':>7} {'th*Fc/12':>8}")
for r in drows:
    P(f"{r['run']:>22} {r['marks'][:10]:>10} {r['F']:>6g} {r['nfam']:>8.0f} {1/r['Fc']:>7.4f} {r['th']:>7.4f} {r['g2']:>7.4f} {r['mu']:>7.4f} {r['rnd']:>7.4f} {r['th']*r['Fc']/12:>8.3f}")
for sel, lab in ((lambda r: r["marks"] == "none", "unmarked"), (lambda r: r["marks"] != "none", "marked (all slot types)")):
    R = [r for r in drows if sel(r)]
    nh = np.array([r["nfam"] / 12 for r in R]); Fc = np.array([r["Fc"] for r in R])
    for key in ("th", "g2", "mu", "rnd"):
        y = np.log([r[key] for r in R])
        c, se = lsq([np.log(nh)], y)
        c2, se2 = lsq([np.log(nh), np.log(Fc)], y)
        fits[f"dual_{lab.split()[0]}_{key}"] = dict(exp_nh=c[1], se=se[1], exp_nh_2=c2[1], exp_Fc_2=c2[2], se_Fc_2=se2[2])
        P(f"[{lab:24s} {key:>4}] log delta_rms = a + b log(#h/sector): b = {c[1]:+.3f} +- {se[1]:.3f};  with log Fc: b = {c2[1]:+.3f}, c(Fc) = {c2[2]:+.3f} +- {se2[2]:.3f}")
    rat = np.array([r["th"] / r["rnd"] for r in R])
    marg = np.array([r["th"] * r["Fc"] / 12 for r in R if r["F"] >= 1])   # per-sector: allowed = 12/Fc
    P(f"[{lab}] theta/random rms ratio: mean {rat.mean():.3f} median {np.median(rat):.3f};  sector delta_theta / allowed(=12/Fc) for F>=1: median {np.median(marg):.3f} max {marg.max():.3f}")
json.dump(fits, open("fits.json", "w"), indent=1)
open("analysis.txt", "w").write("\n".join(out) + "\n")
