"""aggregate probe_N*.json: R statistics, cross-scale lambda for rho and pi = rho*zeta, need scaling, extrapolation."""
import json, glob, math, numpy as np
from analyze import ols
theta = ["X:s30:p78", "X:s100:p78", "X:s300:p78", "X:s1000:p78", "X:t31x61:p78x2", "bP97", "bP409", "beta"]
ctrl = ["XR:s100:p78", "XR:s1000:p78"]
files = sorted(glob.glob("probe_N*.json"), key=lambda f: int(f[7:-5]))
agg = {}
needpts = []
lines = []
def pr(s=""): print(s); lines.append(s)
pr("rank-one pairing T(u) = sum_k V u(k) S_k over 12 row ray classes mod 4; R = |T|^2 / random-model size (mean 1 if u is unrelated to S)")
pr("lambda = cross-scale estimate (pairs of Y in [2^16, 2^18], disjoint s) of |T(rho)|^2 / random size; need: R <~ R_CS = n_eff/C (Cauchy-Schwarz with W)")
for f in files:
    J = json.load(open(f)); N0 = J["N0"]
    pr(f"\nN0={N0}: windows F = {[w['F'] for w in J['windows']]}, direct kernel check {J['check']:.1e}")
    for w in J["windows"]:
        sf = w["signal_fraction"]["big"]["f"]
        row = {}
        for grp, names in (("theta", theta), ("ctrl", ctrl)):
            for key in ["A(Y=big)", "A*zeta(Y=big)", "cross-scale lambda(rho)", "cross-scale lambda(pi=rho*zeta)", "zeta", "random", "adversary conj(S)"]:
                v = [w["tests"][n][key]["R_mean"] for n in names if n in w["tests"]]
                row[(grp, key)] = v
                agg.setdefault((grp, key), []).extend(v)
        adv = np.mean(row[("theta", "adversary conj(S)")])
        needpts.append((w["rows"], adv))
        pr(f"  F={w['F']:<4g} rows={w['rows']:<6d} f_sig(Y=2^18)={sf:.3f}  theta: R[A]={np.mean(row[('theta','A(Y=big)')]):.2f} "
           f"R[A zeta]={np.mean(row[('theta','A*zeta(Y=big)')]):.2f} lambda(rho)={np.mean(row[('theta','cross-scale lambda(rho)')]):.2f} "
           f"lambda(pi)={np.mean(row[('theta','cross-scale lambda(pi=rho*zeta)')]):.2f} R[zeta]={np.mean(row[('theta','zeta')]):.2f} | "
           f"ctrl lambda(rho)={np.mean(row[('ctrl','cross-scale lambda(rho)')]):.2f} | R_CS={adv:.1f}")
pr("\npooled over all windows and N0 (mean ± s.e.; s.e. from the scatter, windows/designs treated as independent):")
for (grp, key), v in sorted(agg.items()):
    v = np.array(v)
    pr(f"  {grp:5s} {key:34s} {v.mean():7.3f} ± {v.std()/math.sqrt(len(v)):.3f}   (max {v.max():.2f}, n={len(v)})")
r = np.array(needpts)
b, se, cov, s = ols([np.log(r[:, 0])], np.log(r[:, 1]))
pr(f"\nCauchy-Schwarz level R_CS ~ rows^{b[1]:.3f}±{se[1]:.3f} (expected 1: need/truth for |T|^2 grows like K)")
lt = np.array(agg[("theta", "cross-scale lambda(pi=rho*zeta)")]); lr = np.array(agg[("theta", "cross-scale lambda(rho)")])
lam = max(lt.mean() + 2 * lt.std() / math.sqrt(len(lt)), lr.mean() + 2 * lr.std() / math.sqrt(len(lr)))
pr(f"extrapolation: |T|/need ~ sqrt(lambda / R_CS) ~ sqrt(lambda C / n_eff) ~ K^(-{b[1]/2:.3f}); at Q = K = Z^0.851: Z^(-{0.851*b[1]/2:.3f})"
   f"  (lambda <= {lam:.2f} at 2 s.e.)")
open("probe_summary.out", "w").write("\n".join(lines) + "\n")
