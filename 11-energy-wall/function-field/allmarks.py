"""
Full marked energy with ALL monic irreducible marks of degree lam (the true mark family), q = 7:
  N(a) = sum_k | sum_P a(P) q^{-lam/2} 1_{(P,k)=1} S_{k,P} |^2   vs  Diag = sum_P |a(P)|^2 q^{-lam} Diag_P
for coherent a = 1, random phases, a(P) = 'Gauss-aligned' (conj of the normalised cubic Gauss sum phase gamma_2(P)),
and the worst case: ||S||_op^2 / mean Diag_P against the random-matrix edge (1 + sqrt(#marks/#rows))^2.
Also: distribution of rho_P over all marks (SM for every mark), and Aff-orbit structure of the marks.
usage: python3 allmarks.py q kappa lam nu [rowsubset]
"""
import sys, json, math, time
import numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ffdriver import run
from ff_lib import irreducibles, random_irreducibles, Fq, resultant, deriv

q, kappa, lam, nu = map(int, sys.argv[1:5])
nsub = int(sys.argv[5]) if len(sys.argv) > 5 else 0
ctrl = int(sys.argv[6]) if len(sys.argv) > 6 else 0
F = Fq(q)
marks = irreducibles(lam, q) if q ** lam <= 5000 else None
if marks is None:
    raise SystemExit("too many marks")
rows = None
if nsub:
    rng = np.random.default_rng(3)
    rows = sorted(rng.choice(q ** kappa, size=nsub, replace=False).tolist())
t0 = time.time()
R = run(q, kappa, nu, lam, marks, threads=3, rows=rows, control=ctrl, tag=f"all{q}_{kappa}_{lam}_{nu}")
secs = time.time() - t0
S, d, cop, a2 = R["S"], R["d"], R["cop"], R["absS2"]
nk, nm = S.shape
Sc = np.where(cop, S, 0)
E = np.array([a2[cop[:, j], j].sum() for j in range(nm)], float)
D = np.array([d[cop[:, j], j].sum() for j in range(nm)], float)
rho = E / D - 1
# normalised cubic Gauss-sum phase of P: gamma_2(P) = -(-tau_2/sqrt q)^lam * chi6(Res(P,P'))^2  (verify_gauss.py)
tau2 = F.tau(2)
g2 = np.array([-((-tau2 / math.sqrt(q)) ** lam) * F.chi(resultant(P, deriv(P, q), q), 2) for P in marks])
rng = np.random.default_rng(11)
coef = {"coherent": np.ones(nm), "random": np.exp(2j * np.pi * rng.random(nm)), "gauss_aligned": np.conj(g2)}
out = dict(q=q, kappa=kappa, lam=lam, nu=nu, control=bool(ctrl), rows=int(cop[:, 0].sum()) if rows is None else len(R["kidx"]),
           nmarks=nm, subset=rows is not None, secs=secs,
           rho_mean=float(rho.mean()), rho_rms=float(np.sqrt(np.mean(rho ** 2))), rho_min=float(rho.min()), rho_max=float(rho.max()),
           n_distinct_rho=int(len(np.unique(np.round(rho, 10)))))
for name, a in coef.items():
    c = a * q ** (-lam / 2)
    Nval = np.sum(np.abs(Sc @ c) ** 2)
    Ndiag = np.sum(np.abs(c) ** 2 * D)
    out[f"N_{name}"] = float(Nval / Ndiag - 1)
sv = np.linalg.svd(Sc, compute_uv=False)
out["opnorm2_over_meandiag"] = float(sv[0] ** 2 / np.mean(D))
out["sv2_top5_over_meandiag"] = (sv[:5] ** 2 / np.mean(D)).round(4).tolist()
out["rmt_edge"] = float((1 + math.sqrt(nm / nk)) ** 2)
G = Sc.conj().T @ Sc
dg = np.sqrt(np.real(np.diag(G)))
C = np.abs(G / np.outer(dg, dg))[~np.eye(nm, dtype=bool)]
out["cross_rms"] = float(np.sqrt(np.mean(C ** 2))); out["cross_max"] = float(C.max())
out["one_over_sqrt_rows"] = float(1 / math.sqrt(nk))
print(json.dumps(out, indent=1))
with open("/home/user/qrh/11-energy-wall/function-field/allmarks.jsonl", "a") as f:
    f.write(json.dumps(out) + "\n")
