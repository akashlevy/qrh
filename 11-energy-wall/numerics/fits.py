"""
numerics: fits of rho_rms against rows/sector, F = N/K and mark length, controls, bias tests, extrapolation.
Reads analysis_windows.json (from analyze.main) and writes fits.json + fits.out (tables for NOTES.md).
"""
import json, math, os, sys
import numpy as np
from analyze import ols

HERE = os.path.dirname(os.path.abspath(__file__))
GEOMS = {  # (M, ell, w) weighted-geometry points (wall4 s.8 / BRIEF)
    "W(0.043) optimum (wall4 s.8)": (0.8511, 0.1915, 0.0426),
    "W(0.040) point": (0.8492, 0.191, 0.0402),
    "BRIEF range low (M=0.82, l=0.18)": (0.82, 0.18, 0.0),
    "BRIEF range high (M=0.86, l=0.20)": (0.86, 0.20, 0.06),
}


def pts_for(stats, sel, ns=12, Fmin=1.0, Fmax=1e9, minrps=8.0, Lof=None):
    P = []
    for N0, W in stats.items():
        for w in W:
            if not (Fmin - 1e-9 <= w["F"] <= Fmax): continue
            s = w["sec"][str(ns)] if str(ns) in w["sec"] else w["sec"][ns]
            if s["rows_per_sector"] < minrps: continue
            for n, v in s["rms"].items():
                if sel(n) and v == v:
                    P.append(dict(N0=int(N0), F=w["F"], rps=s["rows_per_sector"], rms=v, name=n,
                                  L=(Lof(n) if Lof else None), rho=w["rho"][n], mean=s["mean"][n]))
    return P


def fit(P, useL=False):
    y = np.log([p["rms"] for p in P])
    X = [np.log([p["rps"] for p in P]), np.log([p["F"] for p in P])]
    if useL: X.append(np.log([p["L"] for p in P]))
    b, se, cov, s = ols(X, y)
    return dict(coef=b.tolist(), se=se.tolist(), cov=cov.tolist(), resid_sd=s, n=len(y))


def fmt(fr, useL=False):
    b, se = fr["coef"], fr["se"]
    t = f"a={b[0]:+.2f}  b={b[1]:+.3f}±{se[1]:.3f}  c={b[2]:+.3f}±{se[2]:.3f}"
    if useL: t += f"  d={b[3]:+.3f}±{se[3]:.3f}"
    return t + f"  (n={fr['n']}, resid sd {fr['resid_sd']:.2f})"


def margin(fr, geom, useL=False):
    M, l, w = geom
    v = np.array([0, M, w] + ([l] if useL else []))
    b = np.array(fr["coef"]); cov = np.array(fr["cov"])
    e = float(v @ b); s = float(np.sqrt(v @ cov @ v))
    return -e, s      # rho ~ Z^{e}; margin = -e (W needs rho << Z^eps)


def ratio(stats, a, b, ns=12, Fmin=1.0, minrps=8.0):
    """geometric mean over windows of rms(a)/rms(b) (names or name lists pooled in quadrature)."""
    r = []
    for N0, W in stats.items():
        for w in W:
            if w["F"] < Fmin - 1e-9: continue
            s = w["sec"][str(ns)]
            if s["rows_per_sector"] < minrps: continue
            va = [s["rms"][n] for n in a if s["rms"][n] == s["rms"][n]]
            vb = [s["rms"][n] for n in b if s["rms"][n] == s["rms"][n]]
            if va and vb:
                r.append(math.log(np.sqrt(np.mean(np.square(va))) / np.sqrt(np.mean(np.square(vb)))))
    r = np.array(r)
    return float(np.exp(r.mean())), float(r.std() / math.sqrt(max(1, len(r) - 1))), len(r)


def bias(stats, sel, Fmin=1.0, minrps=8.0, sym=frozenset()):
    """z = rho_full / (rms_12/sqrt 12) over windows (sqrt 6 for conjugation-symmetric families): mean, rms (random: 0, 1)."""
    z = []
    for N0, W in stats.items():
        for w in W:
            if w["F"] < Fmin - 1e-9: continue
            s = w["sec"]["12"]
            if s["rows_per_sector"] < minrps: continue
            for n, v in s["rms"].items():
                if sel(n) and v == v and v > 0:
                    z.append(w["rho"][n] / (v / math.sqrt(6 if n in sym else 12)))
    z = np.array(z)
    return float(z.mean()), float(np.sqrt(np.mean(z ** 2))), float(z.mean() / (z.std() / math.sqrt(len(z)))) if len(z) > 2 else 0.0, len(z)


def symmetric_names():
    """names whose row sums satisfy |S_k| = |S_kbar| exactly (conjugation symmetry; checked on the N0=50000 run)."""
    from analyze import load
    d = load(os.path.join(HERE, "data", "energy_N50000.npz"))
    A, B = d["rowA"].astype(int), d["rowB"].astype(int)
    key = {(a, b): i for i, (a, b) in enumerate(zip(A, B))}
    pr = [(i, key[(a - b, -b)]) for i, (a, b) in enumerate(zip(A, B)) if (a - b, -b) in key and key[(a - b, -b)] > i]
    I = np.array([p[0] for p in pr]); Jj = np.array([p[1] for p in pr])
    a = np.abs(d["S"][I]) ** 2; b = np.abs(d["S"][Jj]) ** 2
    rel = np.max(np.abs(a - b) / (a + b + 1e-12), axis=0)
    return {n for n, r in zip(d["names"], rel) if r < 1e-8}


def systematic(stats, sel, Frange, sym):
    """inverse-variance weighted mean of rho_full over windows in Frange (sd per window: rms_12/sqrt(12 or 6))."""
    num = den = 0.0; k = 0
    for N0, W in stats.items():
        for w in W:
            if not (Frange[0] - 1e-9 <= w["F"] <= Frange[1] + 1e-9): continue
            s = w["sec"]["12"]
            if s["rows_per_sector"] < 8: continue
            for n, v in s["rms"].items():
                if sel(n) and v == v and v > 0:
                    sd = v / math.sqrt(6 if n in sym else 12)
                    num += w["rho"][n] / sd ** 2; den += 1 / sd ** 2; k += 1
    return num / den, 1 / math.sqrt(den), k


def main():
    J = json.load(open(os.path.join(HERE, "analysis_windows.json")))
    stats = J["stats"]; metas = J["metas"]
    lines = []
    def pr(s=""): print(s); lines.append(s)
    sym = symmetric_names()
    any_meta = next(iter(metas.values()))
    Ln = any_meta["Ln"]
    slotL = {k: float(np.exp(np.mean(np.log(v)))) for k, v in Ln.items()}
    # 2-slot: geometric mean of product norms
    for k, v in any_meta["slots"].items():
        if "p1" in v:
            slotL[k] = float(np.exp(np.mean(np.log(Ln[k][:len(v["p1"])])) + np.mean(np.log(Ln[k][len(v["p1"]):]))))
    pr(f"runs: N0 = {sorted(int(x) for x in stats)};  windows F = N0/K = 2^(j/2) >= 1, 12 angular row sectors, rows/sector >= 8")
    pr(f"slot mark lengths (geometric mean norm): " + ", ".join(f"{k}: L={v:.0f} ({len(Ln[k])} primes)" for k, v in slotL.items()))
    res = {}
    # ---------------- unmarked + random
    pr("\n== unmarked theta and random-phase control (|c| = |beta w|)")
    for lab, sel in [("beta", lambda n: n == "beta"), ("rand(4 seeds)", lambda n: n.startswith("rand"))]:
        f = fit(pts_for(stats, sel)); res[lab] = f; pr(f"  {lab:18s} {fmt(f)}")
    # ---------------- single marks
    pr("\n== single-mark energies E_P (rows with (P,k)=1), theta beta chi_P^-2 vs control e^{ith} 1_{(n,P)=1}")
    singles = sorted({int(n[2:]) for n in next(iter(stats.values()))[0]["rho"] if n.startswith("bP")})
    pr(f"  {'N(P)':>6s} {'fit theta':60s} | {'control':40s} | theta/control rms")
    for q in singles:
        ft = fit(pts_for(stats, lambda n, q=q: n == f"bP{q}")); fc = fit(pts_for(stats, lambda n, q=q: n == f"rP{q}"))
        rt = ratio(stats, [f"bP{q}"], [f"rP{q}"])
        res[f"bP{q}"] = ft; res[f"rP{q}"] = fc
        pr(f"  {q:6d} b={ft['coef'][1]:+.3f}±{ft['se'][1]:.3f} c={ft['coef'][2]:+.3f}±{ft['se'][2]:.3f} (n={ft['n']}) | "
           f"b={fc['coef'][1]:+.3f}±{fc['se'][1]:.3f} c={fc['coef'][2]:+.3f}±{fc['se'][2]:.3f} | {rt[0]:.3f}±{rt[1]:.3f}")
    Pall = pts_for(stats, lambda n: n.startswith("bP"), Lof=lambda n: int(n[2:]))
    fP = fit(Pall, useL=True); res["single pooled"] = fP
    pr(f"  pooled theta (d = exponent of N(P)):   {fmt(fP, True)}")
    Pc = pts_for(stats, lambda n: n.startswith("rP"), Lof=lambda n: int(n[2:]))
    fPc = fit(Pc, useL=True); res["single pooled control"] = fPc
    pr(f"  pooled control:                        {fmt(fPc, True)}")
    # inside the open window 1 < F < N(P) only
    Pin = [p for p in Pall if p["F"] < p["L"]]
    fPin = fit(Pin, useL=True); res["single pooled, F < N(P)"] = fPin
    pr(f"  pooled theta, only F < N(P) (unproved window): {fmt(fPin, True)}")
    rt = ratio(stats, [f"bP{q}" for q in singles], [f"rP{q}" for q in singles])
    pr(f"  theta/control rms ratio, all P pooled: {rt[0]:.3f} ± {rt[1]:.3f} ({rt[2]} windows)")
    bm = bias(stats, lambda n: n.startswith("bP"), sym=sym); bc = bias(stats, lambda n: n.startswith("rP"), sym=sym)
    pr(f"  bias test z = rho_full/(rms/sqrt12): theta mean {bm[0]:+.3f} rms {bm[1]:.2f} (t={bm[2]:+.2f}, n={bm[3]}); control mean {bc[0]:+.3f} rms {bc[1]:.2f} (t={bc[2]:+.2f})")
    # ---------------- multi-mark designs
    pr("\n== coherent multi-mark core sum N (exact mask 1_{(P,k)=1}); control e^{ith}|beta X w|")
    designs = sorted({n.split(":")[2] for n in next(iter(stats.values()))[0]["rho"] if n.startswith("X:")})
    slots = list(slotL)
    pr(f"  {'slot':8s} {'design':6s} {'theta fit':52s} | {'control fit':36s} | theta/ctrl")
    for sl in slots:
        for d in designs:
            nm = f"X:{sl}:{d}"
            if nm not in next(iter(stats.values()))[0]["rho"]: continue
            ft = fit(pts_for(stats, lambda n, nm=nm: n == nm)); fc = fit(pts_for(stats, lambda n, nm=nm: n == "XR" + nm[1:]))
            rt = ratio(stats, [nm], ["XR" + nm[1:]])
            res[nm] = ft; res["XR" + nm[1:]] = fc
            pr(f"  {sl:8s} {d:6s} b={ft['coef'][1]:+.3f}±{ft['se'][1]:.3f} c={ft['coef'][2]:+.3f}±{ft['se'][2]:.3f} (n={ft['n']}) | "
               f"b={fc['coef'][1]:+.3f}±{fc['se'][1]:.3f} c={fc['coef'][2]:+.3f}±{fc['se'][2]:.3f} | {rt[0]:.3f}±{rt[1]:.3f}")
    for d in ["p78", "gauss", "ones", "kub", "rnd"]:
        Pd = pts_for(stats, lambda n, d=d: n.startswith("X:s") and n.endswith(":" + d), Lof=lambda n: slotL[n.split(":")[1]])
        if Pd:
            f = fit(Pd, useL=True); res[f"design {d} pooled"] = f
            pr(f"  pooled {d:6s} over single slots (d = exponent of L): {fmt(f, True)}")
    Pd = pts_for(stats, lambda n: n.startswith("X:") and (n.endswith(":p78") or n.endswith(":p78x2")), Lof=lambda n: slotL[n.split(":")[1]])
    fD = fit(Pd, useL=True); res["coherent p78 all slots"] = fD
    pr(f"  pooled P78-coherent (p78 + p78x2, all slots):    {fmt(fD, True)}")
    Pdc = pts_for(stats, lambda n: n.startswith("XR:") and (n.endswith(":p78") or n.endswith(":p78x2")), Lof=lambda n: slotL[n.split(":")[1]])
    fDc = fit(Pdc, useL=True); res["coherent p78 control"] = fDc
    pr(f"  its control:                                      {fmt(fDc, True)}")
    Pin = [p for p in Pd if p["F"] < p["L"]]
    pr(f"  (all {len(Pin)}/{len(Pd)} P78-coherent points lie in the unproved window 1 <= F < L)")
    rt = ratio(stats, [n for n in next(iter(stats.values()))[0]["rho"] if n.startswith("X:")],
               [n for n in next(iter(stats.values()))[0]["rho"] if n.startswith("XR:")])
    pr(f"  theta/control rms ratio, all designs pooled: {rt[0]:.3f} ± {rt[1]:.3f} ({rt[2]} windows)")
    for lab, sel in [("p78 (all slots)", lambda n: n.startswith("X:") and "p78" in n), ("ones", lambda n: n.startswith("X:") and n.endswith(":ones")),
                     ("all designs", lambda n: n.startswith("X:")), ("controls", lambda n: n.startswith("XR:"))]:
        bm = bias(stats, sel, sym=sym)
        pr(f"  bias test {lab:16s}: z mean {bm[0]:+.3f}, z rms {bm[1]:.2f}, t = {bm[2]:+.2f} (n={bm[3]})  [z uses sqrt6 for |S_k|=|S_kbar| families]")
    pr(f"\n== conjugation-symmetric families (|S_k| = |S_kbar| exactly): {sorted(n for n in sym if not n.startswith('g'))}")
    pr("== absolute bound on a deterministic off-diagonal component: weighted mean of rho_full (+- 1 s.e.)")
    for lab, sel in [("unmarked beta", lambda n: n == "beta"), ("single marks", lambda n: n.startswith("bP")),
                     ("single controls", lambda n: n.startswith("rP")), ("P78 designs", lambda n: n.startswith("X:") and "p78" in n),
                     ("all designs", lambda n: n.startswith("X:")), ("design controls", lambda n: n.startswith("XR:"))]:
        a = systematic(stats, sel, (1, 4), sym); b = systematic(stats, sel, (5, 64), sym); c = systematic(stats, sel, (90, 1e9), sym)
        res[f"systematic {lab}"] = dict(F1_4=a, F5_64=b, F90=c)
        pr(f"  {lab:16s} F in [1,4]: {a[0]:+.4f} ± {a[1]:.4f} (n={a[2]})   F in [5.7,64]: {b[0]:+.4f} ± {b[1]:.4f}   F >= 90: {c[0]:+.4f} ± {c[1]:.4f}")
    pr("== deep weighted-geometry regime L/F >= 10 only")
    Pd10 = [p for p in Pd if p["L"] / p["F"] >= 10]
    f10 = fit(Pd10, useL=True); res["coherent p78, L/F>=10"] = f10
    pr(f"  P78-coherent designs, L/F >= 10: {fmt(f10, True)}")
    Ps10 = [p for p in Pall if p["L"] / p["F"] >= 10]
    f10s = fit(Ps10, useL=True); res["single, N(P)/F>=10"] = f10s
    pr(f"  single marks, N(P)/F >= 10:      {fmt(f10s, True)}")
    Pw = [p for p in Pd if p["L"] / p["F"] >= 10 and p["L"] ** 2 <= 1.5 * p["N0"]]
    fw = fit(Pw, useL=True); res["coherent p78, weighted-like"] = fw
    pr(f"  P78-coherent, L/F >= 10 and L^2 <= N (weighted-like, {sorted({(p['name'].split(':')[1], p['N0']) for p in Pw})[:3]}...): {fmt(fw, True)}")
    Psw = [p for p in Pall if p["L"] / p["F"] >= 10 and p["L"] ** 2 <= 1.5 * p["N0"]]
    fsw = fit(Psw, useL=True); res["single, weighted-like"] = fsw
    pr(f"  single marks, N(P)/F >= 10 and N(P)^2 <= N: {fmt(fsw, True)}")
    # ---------------- Gram per-mark energies
    gn = [n for n in next(iter(stats.values()))[0]["rho"] if n.startswith("g")]
    if gn:
        f = fit(pts_for(stats, lambda n: n.startswith("g") and n[1:].isdigit())); res["gram per-mark"] = f
        pr(f"\n== per-mark energies E_P for all {len(gn)} slot primes N(P) in [300,600) pooled: {fmt(f)}")
    # ---------------- sector-count method for b (fixed window effects)
    pr("\n== sector-count method: rms over S in {6,12,24,48} sectors at fixed window (window fixed effects)")
    for lab, sel in [("beta", lambda n: n == "beta"), ("rand", lambda n: n.startswith("rand")), ("single marks", lambda n: n.startswith("bP")),
                     ("P78 designs", lambda n: n.startswith("X:") and "p78" in n), ("design controls", lambda n: n.startswith("XR:"))]:
        ys, xs, gid = [], [], []
        g = 0
        for N0, W in stats.items():
            for w in W:
                if w["F"] < 1 - 1e-9: continue
                names = [n for n in w["rho"] if sel(n)]
                for n in names:
                    for S in ("6", "12", "24", "48"):
                        s = w["sec"][S]
                        if s["rows_per_sector"] < 8 or not s["rms"][n] > 0: continue
                        ys.append(math.log(s["rms"][n])); xs.append(math.log(s["rows_per_sector"])); gid.append(g)
                    g += 1
        ys, xs, gid = np.array(ys), np.array(xs), np.array(gid)
        # demean within groups
        yd = ys.copy(); xd = xs.copy()
        for k in np.unique(gid):
            m = gid == k; yd[m] -= ys[m].mean(); xd[m] -= xs[m].mean()
        bb = float(xd @ yd / (xd @ xd)); r = yd - bb * xd
        se = math.sqrt(float(r @ r) / max(1, len(r) - len(np.unique(gid)) - 1) / float(xd @ xd))
        res[f"sector-count b {lab}"] = (bb, se)
        pr(f"  {lab:16s} b = {bb:+.3f} ± {se:.3f}  ({len(ys)} points, {len(np.unique(gid))} windows x names)")
    # ---------------- extrapolation
    pr("\n== extrapolation: rho ~ Z^{bM + cw + d l}; W needs rho << Z^eps; margin = -(bM + cw + dl) (Z-exponent)")
    for gname, geom in GEOMS.items():
        out = []
        for lab, useL in [("beta", False), ("single pooled", True), ("single pooled, F < N(P)", True), ("single, N(P)/F>=10", True),
                          ("single, weighted-like", True), ("coherent p78 all slots", True), ("coherent p78, L/F>=10", True),
                          ("coherent p78, weighted-like", True), ("design p78 pooled", True),
                          ("rand(4 seeds)", False), ("coherent p78 control", True)]:
            m, s = margin(res[lab], geom, useL)
            out.append(f"{lab}: {m:.3f}±{s:.3f}")
        pr(f"  {gname} (M,l,w)={geom}:\n     " + "\n     ".join(out))
    json.dump(res, open(os.path.join(HERE, "fits.json"), "w"), indent=1)
    open(os.path.join(HERE, "fits.out"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
