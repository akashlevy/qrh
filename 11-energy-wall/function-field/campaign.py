"""
Campaign driver: for each (q, kappa, lam, nu) measure, on the same rows,
  theta  : rho_P = E_P/Diag_P - 1 per mark, 4th-moment ratio m4, cross-mark Gram, N(c) for coherent/random a(P),
           operator norm of [S_{k,P}] vs the random-matrix edge (1 + sqrt(#marks/#rows))^2
  control: the same with i.i.d. uniform sextic phases on the same support (FFRANDOM)
Rows: all squarefree k of degree kappa, or a uniformly random subset when q^(kappa+nu)*#marks > budget.
Output: one JSON line per configuration appended to campaign.jsonl (in this directory).
usage: python3 campaign.py <setname>
"""
import sys, os, json, math, time
import numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ffdriver import run, stats
from ff_lib import random_irreducibles, irreducibles

HERE = "/home/user/qrh/11-energy-wall/function-field"


def marks_for(q, lam, nm, seed=1):
    if lam == 0:
        return [(1,)]
    if lam <= 2 or q ** lam <= 3000:
        irr = irreducibles(lam, q)
        rng = np.random.default_rng(seed * 101 + q + 7 * lam)
        idx = rng.choice(len(irr), size=min(nm, len(irr)), replace=False)
        return [irr[i] for i in sorted(idx)]
    return random_irreducibles(lam, q, nm, seed=seed)


def measure(q, kappa, lam, nu, nm=4, budget=3e9, minrows=300, seed=1, threads=3, control=True, allmarks=False):
    marks = marks_for(q, lam, nm if not allmarks else 10 ** 9, seed)
    nm = len(marks)
    rows = None
    full_cost = q ** (kappa + nu) * nm
    if full_cost > budget:
        nsub = max(minrows, int(budget / (q ** nu * nm)))
        if nsub < q ** kappa:
            rng = np.random.default_rng(seed + 1000 * kappa + nu)
            rows = sorted(rng.choice(q ** kappa, size=nsub, replace=False).tolist())
    t0 = time.time()
    R = run(q, kappa, nu, lam, marks, threads=threads, rows=rows, tag=f"c{q}_{kappa}_{lam}_{nu}")
    st = stats(R)
    out = dict(q=q, kappa=kappa, lam=lam, nu=nu, nmarks=nm, subset=rows is not None,
               D=2 * kappa + lam - 2, window=(kappa < nu < kappa + lam - 2), theta=st)
    if control:
        Rc = run(q, kappa, nu, lam, marks, threads=threads, rows=rows, control=12345 + seed, tag=f"cc{q}_{kappa}_{lam}_{nu}")
        out["control"] = stats(Rc)
    out["secs"] = time.time() - t0
    return out


SETS = {
    # large q, fixed degrees
    "largeq": [(q, k, l, n) for (k, l, n) in [(2, 1, 3), (2, 2, 3), (2, 2, 4), (3, 2, 4), (3, 3, 4), (3, 1, 4)]
               for q in (7, 13, 19, 31, 37, 43)],
    # fixed q = 7, growing degrees, lam = 2..6
    "q7": [(7, k, l, n) for (k, l, n) in [
        (3, 2, 4), (4, 2, 5), (5, 2, 6), (6, 2, 7),
        (3, 3, 4), (4, 3, 5), (5, 3, 6), (6, 3, 7),
        (4, 4, 5), (5, 4, 6), (6, 4, 7), (6, 4, 8), (7, 4, 8),
        (4, 5, 5), (4, 5, 6), (5, 5, 6), (5, 5, 7), (6, 5, 7), (6, 5, 8),
        (5, 6, 6), (5, 6, 7), (6, 6, 7), (6, 6, 8)]],
    # fixed q = 13
    "q13": [(13, k, l, n) for (k, l, n) in [
        (3, 2, 4), (3, 3, 4), (4, 3, 5), (4, 4, 5), (4, 5, 5), (4, 5, 6), (5, 4, 6), (5, 5, 6)]],
}

if __name__ == "__main__":
    name = sys.argv[1]
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 3e9
    outf = os.path.join(HERE, f"campaign_{name}.jsonl")
    for (q, k, l, n) in SETS[name]:
        try:
            res = measure(q, k, l, n, budget=budget)
        except Exception as e:  # keep going
            print("FAILED", q, k, l, n, e, flush=True)
            continue
        th, co = res["theta"], res["control"]
        print(f"q={q:2d} k={k} l={l} nu={n} win={int(res['window'])} rows={th['rows']:6d}{'s' if res['subset'] else ' '} "
              f"rho_th mean={th['rho_mean']:+.4f} rms={th['rho_rms']:.4f} | ctrl mean={np.mean(co['rho']):+.4f} rms={co['rho_rms']:.4f} | "
              f"m4 th={np.nanmean(th['m4']):.3f} ct={np.nanmean(co['m4']):.3f} | cross th={th['cross_rms'] or 0:.4f} ct={co['cross_rms'] or 0:.4f} | "
              f"Ncoh={th['N_coherent']:+.4f} op={th['opnorm2_over_diag']:.3f}/edge {th['rmt_edge']:.3f} | {res['secs']:.0f}s", flush=True)
        with open(outf, "a") as f:
            f.write(json.dumps(res) + "\n")
