"""
Driver for ffenergy.c: choose marks, run, parse exact output, compute statistics.
"""
import sys, os, subprocess, random, math, cmath, json, time
import numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ff_lib import irreducibles, monic_polys, resultant, deriv, Fq, trim

HERE = "/home/user/qrh/11-energy-wall/function-field"
SCR = os.environ.get("FFSCR", "/tmp/claude-0/-home-user-qrh/f453a333-8650-5a2b-a022-02397b1fb4ba/scratchpad")
OMEGA = complex(-0.5, math.sqrt(3) / 2)


def irreducibles_fast(d, q):
    """monic irreducibles of degree d via Rabin-type test using ff_lib arithmetic (d <= 4)."""
    if d <= 3 and q ** d <= 4000:
        return irreducibles(d, q)
    # generic: degree d irreducible iff no factor of degree <= d/2  (trial division by small irreducibles)
    from ff_lib import pdivmod
    small = []
    for e in range(1, d // 2 + 1):
        small += irreducibles(e, q)
    out = []
    for p in monic_polys(d, q):
        ok = True
        for P in small:
            if not pdivmod(p, P, q)[1]:
                ok = False
                break
        if ok:
            out.append(p)
    return out


def choose_marks(q, lam, nm, seed=1):
    if lam == 0:
        return [(1,)]
    irr = irreducibles_fast(lam, q)
    rnd = random.Random(seed * 1000 + q * 10 + lam)
    if nm >= len(irr):
        return irr
    return rnd.sample(irr, nm)


def run(q, kappa, nu, lam, marks, threads=3, tag=None, control=0, rows=None):
    tag = tag or f"q{q}_k{kappa}_n{nu}_l{lam}_m{len(marks)}"
    fin = os.path.join(SCR, f"in_{tag}.txt")
    fout = os.path.join(SCR, f"out_{tag}.bin")
    with open(fin, "w") as f:
        f.write(f"{q} {kappa} {nu} {len(marks)} {lam}\n")
        for P in marks:
            f.write(" ".join(str(c) for c in P) + "\n")
    t0 = time.time()
    env = dict(os.environ)
    env.pop("FFRANDOM", None)
    if control:
        env["FFRANDOM"] = str(control)
    env.pop("FFROWS", None)
    if rows is not None:
        env["FFROWS"] = ",".join(str(int(r)) for r in rows)
    subprocess.run([os.path.join(HERE, "ffenergy"), fin, fout, str(threads)], check=True, env=env)
    dt = time.time() - t0
    raw = np.fromfile(fout, dtype=np.int64)
    q_, k_, n_, nm_, l_, nk = raw[:6]
    rec = raw[6:].reshape(nk, 1 + 4 * nm_)
    kidx = rec[:, 0]
    A = rec[:, 1::4].astype(np.float64); B = rec[:, 2::4].astype(np.float64)
    S = A + B * OMEGA
    d = rec[:, 3::4].astype(np.float64)
    cop = rec[:, 4::4].astype(bool)
    absS2 = (rec[:, 1::4] ** 2 - rec[:, 1::4] * rec[:, 2::4] + rec[:, 2::4] ** 2)   # exact integers
    os.remove(fout)
    return dict(q=q, kappa=kappa, nu=nu, lam=lam, marks=marks, kidx=kidx, S=S, d=d, cop=cop,
                absS2=absS2, secs=dt)


def stats(R, coeff_sets=("coherent", "random")):
    S, d, cop, a2 = R["S"], R["d"], R["cop"], R["absS2"]
    nk, nm = S.shape
    out = {}
    E = np.array([a2[cop[:, j], j].sum() for j in range(nm)], dtype=np.float64)
    D = np.array([d[cop[:, j], j].sum() for j in range(nm)])
    rho = E / D - 1
    m4 = np.array([np.mean(a2[cop[:, j], j].astype(float) ** 2) / np.mean(a2[cop[:, j], j].astype(float)) ** 2 for j in range(nm)])
    out["rows"] = int(cop[:, 0].sum())
    out["cols_per_row"] = float(D[0] / cop[:, 0].sum())
    out["rho"] = rho.tolist()
    out["rho_mean"] = float(rho.mean()); out["rho_rms"] = float(np.sqrt(np.mean(rho ** 2)))
    out["m4"] = m4.tolist()
    # normalised cross-mark Gram
    Sc = np.where(cop, S, 0)
    G = Sc.conj().T @ Sc
    dg = np.sqrt(np.real(np.diag(G)))
    C = G / np.outer(dg, dg)
    off = C[~np.eye(nm, dtype=bool)]
    out["cross_rms"] = float(np.sqrt(np.mean(np.abs(off) ** 2))) if nm > 1 else None
    out["cross_max"] = float(np.max(np.abs(off))) if nm > 1 else None
    # marked energy N(c) = sum_k |sum_P c_P 1_{(P,k)=1} S_{k,P}|^2, c_P = a(P) q^{-lam/2}
    lam, q = R["lam"], R["q"]
    rng = np.random.default_rng(7)
    for name in coeff_sets:
        if name == "coherent":
            a = np.ones(nm)
        elif name == "random":
            a = np.exp(2j * np.pi * rng.random(nm))
        c = a * q ** (-lam / 2)
        Nval = np.sum(np.abs(Sc @ c) ** 2)
        Ndiag = np.sum(np.abs(c) ** 2 * D)
        out[f"N_{name}"] = float(Nval / Ndiag - 1)
    # operator norm of [S_{k,P}] relative to mean diagonal (W uniformly in the marks)
    sv = np.linalg.svd(Sc, compute_uv=False)
    out["opnorm2_over_diag"] = float(sv[0] ** 2 / np.mean(D))
    out["rmt_edge"] = float((1 + math.sqrt(nm / out["rows"])) ** 2)
    return out


if __name__ == "__main__":
    q, kappa, nu, lam, nm = map(int, sys.argv[1:6])
    marks = choose_marks(q, lam, nm)
    R = run(q, kappa, nu, lam, marks)
    st = stats(R)
    print(json.dumps({k: v for k, v in st.items() if k not in ("rho", "m4")}, indent=1))
    print("rho per mark:", np.round(st["rho"], 4))
    print("secs", R["secs"])
