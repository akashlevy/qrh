"""Tables and fits for the function-field energy campaigns (writes analysis.txt)."""
import json, math, sys, os
import numpy as np
HERE = "/home/user/qrh/11-energy-wall/function-field"
out = []
def P(s=""):
    out.append(s); print(s)

def load(fn):
    p = os.path.join(HERE, fn)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []

# ---------------- large q at fixed degrees ----------------
L = load("campaign_largeq.jsonl")
P("== Large q, fixed degrees (kappa, lam, nu); D = 2kappa+lam-2 (kappa>=3; for kappa=2 the measured degree is lam+1) ==")
P("   rho = E_P/Diag - 1 (theta);  rho*q;  m4 = mean|S|^4/mean(|S|^2)^2;  cross_rms*sqrt(rows) (RMT ~ 1); 's' = row subset")
groups = {}
for r in L:
    groups.setdefault((r["kappa"], r["lam"], r["nu"]), []).append(r)
for key in sorted(groups):
    k, l, n = key
    P(f" (kappa,lam,nu) = {key}:")
    xs, ys, ys1 = [], [], []
    for r in sorted(groups[key], key=lambda t: t["q"]):
        th, co = r["theta"], r["control"]
        rho = th["rho_mean"]; q = r["q"]
        P(f"   q={q:2d} rows={th['rows']:5d}{'s' if r['subset'] else ' '} rho={rho:+.4f} rho*q={rho*q:+7.3f} (1+rho)*q={(1+rho)*q:7.3f} "
          f"m4={np.nanmean(th['m4']):.3f} cross*sqrt(rows)={(th['cross_rms'] or 0)*math.sqrt(th['rows']):.2f} | ctrl rho={np.mean(co['rho']):+.4f}")
        if not r["subset"]:
            xs.append(math.log(q)); ys.append(math.log(abs(rho) + 1e-12)); ys1.append(math.log(1 + rho))
    if len(xs) >= 3:
        b = np.polyfit(xs, ys, 1); b1 = np.polyfit(xs, ys1, 1)
        P(f"   fit (full-row runs): |rho| ~ q^{b[0]:+.2f};  1+rho ~ q^{b1[0]:+.2f}")

# ---------------- fixed q, growing degree ----------------
for fn in ("campaign_q7.jsonl", "campaign_q13.jsonl"):
    R = load(fn)
    if not R:
        continue
    P(""); P(f"== {fn}: fixed q, growing degrees; window = kappa < nu < kappa+lam-2 (genuinely beyond the wall on both sides of the FE) ==")
    th_pts, co_pts = [], []
    for r in R:
        th, co = r["theta"], r["control"]
        rows = th["rows"]
        rr = np.array(th["rho"]); cc = np.array(co["rho"])
        P(f" q={r['q']} (k,l,nu)=({r['kappa']},{r['lam']},{r['nu']}) D={r['D']} win={int(r['window'])} rows={rows:5d}{'s' if r['subset'] else ' '} "
          f"rho per mark {np.round(rr,4).tolist()}  |rho|*sqrt(rows): th {np.sqrt(np.mean(rr**2))*math.sqrt(rows):.2f} ct {np.sqrt(np.mean(cc**2))*math.sqrt(rows):.2f} "
          f"| m4 {np.nanmean(th['m4']):.2f} | cross*sqrt(rows) th {th['cross_rms']*math.sqrt(rows):.2f} ct {co['cross_rms']*math.sqrt(rows):.2f} "
          f"| Ncoh {th['N_coherent']:+.4f} Nrnd {th['N_random']:+.4f} | op/edge {th['opnorm2_over_diag']/th['rmt_edge']:.3f}")
        if r["nu"] <= r["D"]:
            for v in rr: th_pts.append((math.log(rows), math.log(abs(v) + 1e-9)))
            for v in cc: co_pts.append((math.log(rows), math.log(abs(v) + 1e-9)))
    for name, pts in (("theta", th_pts), ("control", co_pts)):
        if len(pts) > 5:
            X = np.array(pts)
            b = np.polyfit(X[:, 0], X[:, 1], 1)
            P(f" pooled fit {name}: |rho_P| ~ rows^{b[0]:+.3f}  (RMT: -0.5);  median |rho|*sqrt(rows) = {np.median(np.exp(X[:,1]+0.5*X[:,0])):.2f}")

A = load("allmarks.jsonl")
if A:
    P(""); P("== All marks of degree lam (true mark family), q = 7 ==")
    for r in A:
        P(" " + json.dumps(r))
open(os.path.join(HERE, "analysis.txt"), "w").write("\n".join(out) + "\n")
