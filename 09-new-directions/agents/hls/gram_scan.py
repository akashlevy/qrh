"""
hls: is (HLS) true?  Top eigenvalue of the hybrid Gram matrix vs (#rows + #cols), with controls,
and Rayleigh quotients of theta-type (Gauss-sum) vectors.

For each config (N0, K, L0, nP):
  lam_top(S S^*) for the arithmetic family xi_{k,P} = psi_k eta_P;
  controls: 'prod' (random signs for psi, random cube roots for eta, same zero pattern),
            'quad' (Heath-Brown family alone: psi_k for #k' ~ #rows quadratic characters, eta = 1),
            'seqbound' = #P * lam_top(Psi Psi^T) (P78 sequential bound, quadratic first);
  test vectors a: random phases, gamma_j (j=1..5), mu*gamma_j, beta = gamma_2 conj(alpha), mu*beta,
            the 'column' extremal conj(xi_{k0,P0}), the top eigenvector;
  rho(a) = sum_{k,P}|S a|^2 / sum_{k,P} sum_n |a_n|^2 |xi_{kP}(n)|^2  (1 = diagonal; theta excess shows as rho >> 1);
  the k = 1 rows (pure cubic sieve) are reported separately.
Usage: nice python3 gram_scan.py '[{"N0":20000,"K":800,"L0":100,"nP":20}, ...]' out.json
"""
import sys, json, time, math
import numpy as np
from hls_lib import Family, rand_control


def diag_norm(F, a):
    ck = np.sum(F.Psi.astype(np.float64) ** 2, axis=0)          # #k coprime to n
    cP = np.sum(np.abs(F.H) ** 2, axis=0)                        # #P coprime to n
    return float(np.sum(np.abs(a) ** 2 * ck * cP))


def rho(F, a):
    Y = F.S(a)
    tot = float(np.sum(np.abs(Y) ** 2))
    k1 = float(np.sum(np.abs(Y[0]) ** 2)) if F.kN[0] == 1 else 0.0
    d = diag_norm(F, a)
    na = float(np.sum(np.abs(a) ** 2))
    return dict(R=tot / na, rho=tot / d, rho_wo_k1=(tot - k1) / d, R_k1=k1 / na)


def run(cfg, seed=0):
    t0 = time.time()
    F = Family(cfg["N0"], cfg["K"], cfg["L0"], nP=cfg.get("nP"), mod4=cfg.get("mod4", False))
    res = dict(cfg=cfg, nc=F.nc, nk=F.nk, nP=F.nP, nrows=F.nrows)
    lam, vec = F.top(k=3)
    res["lam"] = lam.tolist()
    res["lam_over_sum"] = float(lam[0] / (F.nrows + F.nc))
    res["lam_over_max"] = float(lam[0] / max(F.nrows, F.nc))
    # where does the top eigenvector live?  its k=1 share and its overlap with gauss vectors
    a = vec[:, 0]
    Y = F.S(a)
    res["top_k1_share"] = float(np.sum(np.abs(Y[0]) ** 2) / np.sum(np.abs(Y) ** 2))
    res["top_maxrow_share"] = float(np.max(np.abs(Y) ** 2) / np.sum(np.abs(Y) ** 2))
    res["top_overlap"] = {f"g{j}": float(abs(np.vdot(F.gam[:, j], a)) / math.sqrt(F.nc)) for j in range(1, 6)}
    # random-matrix reference (sqrt rows + sqrt cols)^2
    res["rmt"] = float((math.sqrt(F.nrows) + math.sqrt(F.nc)) ** 2)
    # controls
    for kind in ("prod",):
        G = rand_control(F, kind, seed=seed + 11)
        lg, _ = G.top(k=1)
        res[f"lam_{kind}"] = float(lg[0])
    # sequential bound: #P * lam(Psi Psi^T) (quadratic sieve per P), and cubic-first: #k * lam(H H^*)
    from scipy.sparse.linalg import eigsh, LinearOperator
    Pq = F.Psi.astype(np.float64)
    if F.nk <= F.nc:
        lq = np.linalg.eigvalsh(Pq @ Pq.T)[-1]
    else:
        lq = eigsh(LinearOperator((F.nc, F.nc), matvec=lambda x: Pq.T @ (Pq @ x)), k=1, which="LA")[0][0]
    lc = np.linalg.eigvalsh(F.H @ F.H.conj().T)[-1]
    res["lam_quad_perP"] = float(lq)
    res["lam_cubic_perk"] = float(lc)
    res["seq_bound"] = float(min(F.nP * lq, F.nk * lc))
    # Heath-Brown family alone with ~#rows characters
    if cfg.get("quadref", True):
        Kq = cfg["K"] * F.nP
        Q = Family(cfg["N0"], Kq, cfg["L0"], nP=1, gauss=False, verbose=False)
        Q.H = np.where(Q.H == 0, 1.0 + 0j, 1.0 + 0j)
        lqq, _ = Q.top(k=1)
        res["quadref"] = dict(nk=Q.nk, lam=float(lqq[0]), lam_over_sum=float(lqq[0] / (Q.nk + Q.nc)))
    # test vectors
    rng = np.random.default_rng(seed)
    vecs = {}
    for s in range(3):
        vecs[f"rand{s}"] = np.exp(2j * np.pi * rng.random(F.nc))
    for j in range(1, 6):
        vecs[f"g{j}"] = F.gam[:, j]
        vecs[f"mu*g{j}"] = F.mu * F.gam[:, j]
    vecs["beta"] = F.beta
    vecs["mu*beta"] = F.mu * F.beta
    k0 = F.nk // 2; P0 = 0
    vecs["column"] = np.conj(F.Psi[k0].astype(np.float64) * F.H[P0])
    vecs["top"] = a
    res["vectors"] = {k: rho(F, v) for k, v in vecs.items()}
    res["time"] = time.time() - t0
    return res


if __name__ == "__main__":
    cfgs = json.loads(sys.argv[1])
    out = sys.argv[2] if len(sys.argv) > 2 else None
    allres = []
    for cfg in cfgs:
        r = run(cfg)
        allres.append(r)
        v = r["vectors"]
        print(f"N0={cfg['N0']} K={cfg['K']} L0={cfg['L0']}: cols={r['nc']} k={r['nk']} P={r['nP']} rows={r['nrows']}"
              f" | lam={r['lam'][0]:.0f}  lam/(rows+cols)={r['lam_over_sum']:.3f}  lam/max={r['lam_over_max']:.3f}"
              f"  rmt={r['rmt']:.0f}  prod-ctrl={r['lam_prod']:.0f}  seq={r['seq_bound']:.0f}"
              + (f"  HBquad(k'={r['quadref']['nk']}): lam/(k'+cols)={r['quadref']['lam_over_sum']:.3f}" if 'quadref' in r else "")
              + f"  top: k1share={r['top_k1_share']:.3f} maxrow={r['top_maxrow_share']:.3f}", flush=True)
        print("    rho: " + "  ".join(f"{k}={v[k]['rho']:.3f}" + (f"(w/o k1 {v[k]['rho_wo_k1']:.3f})" if k in ('g2','g4','beta','mu*beta','g1','g5') else "")
                                   for k in v), flush=True)
        if out:
            json.dump(allres, open(out, "w"), indent=1)
