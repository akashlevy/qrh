"""
numerics / step 4: the rank-one pairing T = sum_m phi(q_m/Q) xi(m) rho(m) B_m (wall4 s.7) -- tested with the
actual sextic Gauss-sum probe vector, since rho itself (the residue at w = 2/3) is not computable per row.

  A_{m,sigma}(Y) = sum_{s sf primary odd, s in sigma, (s,m)=1} gamma_1(s) conj chi_s(-m) W(N(s)/Y),
  W = sin^2 dyadic bump on (Y/2, Y],  sigma = ray class of s mod 4 (12 odd classes).
By Corollary D (wall4 s.7.1) A = A^sig + A^rest with A^sig = Y^{-1/3} W^(1/6+iv) rho(m): the theta vector rho is the
pole part (w = 2/3) of exactly this Dirichlet series.  gamma_1(s) = mu(s) G(s) conj beta(s), G = conj chi_s(4) Gamma(s)
(P78 eq. signal; checked in validate.py).  B_m = S_m are the reflected row sums of the coherent marked designs
(run_energy.py, beyond the wall), so T(u) = sum_k V_k u(k) S_k is the reflected form <Psi* u, Pi_V b> (wall4 s.7.2).

Statistic (row ray classes tau mod 4 absorb unknown ray characters xi, psi_1):
  R(u)   = sum_tau |sum_{k in tau} V u S|^2 / sum_k V^2 |u|^2 |S|^2     (random model: mean 1)
  eta(u) = sum_tau |...|^2 / ((sum V|u|^2)(sum V D))                     (W_theta needs eta <~ Z^eps; random ~ C/K)
"""
import sys, os, json, math, time, ctypes
import numpy as np
from nlib import Arith, W
from analyze import load

HERE = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(HERE, "libprobe.so"))
i64 = np.ctypeslib.ndpointer(dtype=np.int64, flags="C_CONTIGUOUS")
i8 = np.ctypeslib.ndpointer(dtype=np.int8, flags="C_CONTIGUOUS")
f64 = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
lib.probe.argtypes = [ctypes.c_int64, i8, ctypes.c_int64, i64, i64, i64, f64, f64, ctypes.c_int64, f64, f64, ctypes.c_int]


def Gamma(a, b):
    return (1 + 1j ** ((-b) % 4) + 1j ** (a % 4) + 1j ** ((b - a) % 4)) / 2


def cls4(A, B):
    """index of the odd class of A + B w mod 4 (12 classes)."""
    key = (np.asarray(A) % 4) * 4 + (np.asarray(B) % 4)
    odd = [a * 4 + b for a in range(4) for b in range(4) if not (a % 2 == 0 and b % 2 == 0)]
    lut = -np.ones(16, np.int64)
    for i, k in enumerate(odd): lut[k] = i
    return lut[key]


def svec(ar, Y):
    """s with Y/2 < N(s) <= Y: factor lists, gamma_1(s) W, class."""
    A, B, Nn, mu, Gt, fac = ar.elements(Y // 2, Y)
    lam = np.array([sum(int(ar.lam3[i]) for i in f) % 3 for f in fac])
    beta = np.exp(2j * np.pi * lam / 3) * Gt * np.conj((A + B * W) / np.sqrt(Nn))
    k4 = np.array([sum(int(ar.k4[i]) for i in f) % 6 for f in fac])
    G = np.exp(-2j * np.pi * k4 / 6) * np.array([Gamma(int(a), int(b)) for a, b in zip(A, B)])
    g1 = mu * G * np.conj(beta)
    w = np.sin(2 * np.pi * (Nn / Y - 0.5)) ** 2
    return A, B, Nn, fac, g1 * w, cls4(A, B)


def build_A(ar, rows, Ys, nth=2):
    t0 = time.time()
    ar.k4 = np.array([int(p.sextic_exp(4, 0)[0]) for p in ar.P])
    S = {Y: svec(ar, Y) for Y in Ys}
    need = sorted({i for Y in Ys for f in S[Y][3] for i in f})
    pidx = {p: j for j, p in enumerate(need)}
    nr = len(rows["A"])
    E = np.empty((len(need), nr), np.int8)
    for p, j in pidx.items():
        E[j] = ar.P[p].sextic_exp(-rows["A"], -rows["B"])
    E0 = np.where(E < 0, -1, 0).astype(np.int8)
    print(f"   E table {E.shape} ({E.nbytes/1e6:.0f} MB)  {time.time()-t0:.1f}s", flush=True)
    out = {}
    for Y in Ys:
        A_, B_, Nn, fac, coef, cl = S[Y]
        ptr = np.zeros(len(fac) + 1, np.int64); ptr[1:] = np.cumsum([len(f) for f in fac])
        fp = np.array([pidx[i] for f in fac for i in f], np.int64)
        ore = np.zeros(12 * nr); oim = np.zeros(12 * nr)
        lib.probe(nr, np.ascontiguousarray(E.ravel()), len(fac), ptr, fp, cl.astype(np.int64),
                  np.ascontiguousarray(coef.real), np.ascontiguousarray(coef.imag), 12, ore, oim, nth)
        out[Y] = (ore + 1j * oim).reshape(12, nr)
        out[(Y, "ns")] = len(fac); out[(Y, "w2")] = float(np.sum(np.abs(coef) ** 2))
        # random-model size of |A_{m,sigma}|^2: sum_{s in sigma, (s,m)=1} |coef|^2 (same kernel, trivial characters)
        nre = np.zeros(12 * nr); nim = np.zeros(12 * nr)
        lib.probe(nr, np.ascontiguousarray(E0.ravel()), len(fac), ptr, fp, cl.astype(np.int64),
                  np.ascontiguousarray(np.abs(coef) ** 2), np.zeros(len(fac)), 12, nre, nim, nth)
        out[(Y, "noise")] = nre.reshape(12, nr)
        print(f"   A(Y={Y}): {len(fac)} moduli s  {time.time()-t0:.1f}s", flush=True)
    return out, S, E, pidx


def check_direct(ar, rows, S, Aout, Y, nrow_check=5):
    """direct evaluation of A_{m,sigma}(Y) for a few rows with python sextic symbols."""
    A_, B_, Nn, fac, coef, cl = S[Y]
    rng = np.random.default_rng(0)
    worst = 0
    for m in rng.choice(len(rows["A"]), nrow_check, replace=False):
        acc = np.zeros(12, complex)
        for j, f in enumerate(fac):
            e = 0; ok = True
            for i in f:
                k = int(ar.P[i].sextic_exp(-int(rows["A"][m]), -int(rows["B"][m]))[0])
                if k < 0: ok = False; break
                e += k
            if ok: acc[cl[j]] += coef[j] * np.exp(-2j * np.pi * (e % 6) / 6)
        worst = max(worst, np.abs(acc - Aout[Y][:, m]).max() / max(1e-12, np.abs(acc).max()))
    return worst


def stats_T(u, S, D, V, rcls, ncls=12):
    """u, S, D, V on the window rows.  R = |T|^2 / random-model size;  eta = |T|^2 / (Cauchy-Schwarz-with-W size)."""
    num = np.zeros(ncls, complex); nu = np.zeros(ncls); nd = np.zeros(ncls)
    np.add.at(num, rcls, V * u * S)
    np.add.at(nu, rcls, V * np.abs(u) ** 2); np.add.at(nd, rcls, V * D)
    tot = float(np.sum(np.abs(num) ** 2))
    rnd = float(np.sum(V ** 2 * np.abs(u) ** 2 * np.abs(S) ** 2))
    need = float(np.sum(nu * nd))
    return tot / rnd, tot / need


def run(N0=50000, Fs=(1, 2, 4, 8, 16, 32), designs=("X:s30:p78", "X:s100:p78", "X:s300:p78", "X:s1000:p78", "X:t31x61:p78x2", "beta",
                                                     "XR:s100:p78", "XR:s1000:p78", "bP97", "bP409"),
        nrand=40, nth=2, Ymax=200000):
    t0 = time.time()
    d = load(os.path.join(HERE, "data", f"energy_N{N0}.npz"))
    keep = d["rowN"] <= N0 / min(Fs)
    d["S"] = d["S"][keep]; d["D"] = d["D"][keep]
    rows = dict(A=d["rowA"][keep].astype(np.int64), B=d["rowB"][keep].astype(np.int64), N=d["rowN"][keep], mu=d["rowmu"][keep],
                g2=d["rowg2"][keep], alpha=d["rowalpha"][keep])
    Ks = [N0 / F for F in Fs]
    Ys = sorted({int(2 ** round(math.log2(K * c))) for K in Ks for c in (1, 4)} | {Ymax})
    Ys = [Y for Y in Ys if 64 <= Y <= Ymax]
    ar = Arith(max(410000, 2 * N0))
    print(f"N0={N0}: rows {len(rows['A'])}, Y = {Ys}", flush=True)
    Aout, Sdat, E, pidx = build_A(ar, rows, Ys, nth)
    chk = check_direct(ar, rows, Sdat, Aout, Ys[0])
    print(f"   direct check of A at Y={Ys[0]}: max rel diff {chk:.2e}", flush=True)
    zeta = rows["mu"] * np.conj(rows["alpha"]) * np.conj(rows["g2"])       # Lemma Z scalar up to a ray character
    beta_r = rows["g2"] * np.conj(rows["alpha"])
    rng = np.random.default_rng(5)
    names = d["names"]
    rcls = cls4(rows["A"], rows["B"])
    out = dict(N0=N0, Ys=Ys, check=chk, windows=[])
    # signal diagnostics: mean |A|^2 per row relative to the random size sum_s |coef|^2
    diag = {}
    for Y in Ys:
        a2 = np.mean(np.abs(Aout[Y]) ** 2, axis=1).sum()
        diag[Y] = float(a2 / Aout[(Y, "w2")])
    out["A_power_over_random"] = diag
    print("   |A|^2 / (sum |coef|^2) by Y (all rows): " + ", ".join(f"{Y}:{v:.2f}" for Y, v in diag.items()), flush=True)
    for F in Fs:
        K = N0 / F
        x = rows["N"] / K
        win = (x > 0.5) & (x <= 1) & (rows["N"] > 1)
        V = np.where(win, np.sin(2 * np.pi * (x - 0.5)) ** 2, 0.0)[win]
        rc = rcls[win]
        rec = dict(F=F, K=K, rows=int(win.sum()), tests={})
        for dn in designs:
            j = names.index(dn)
            Sk = d["S"][win, j]; Dk = d["D"][win, j]
            tests = {}
            for Y in Ys:
                if not (K / 1.5 <= Y <= 6 * K or Y == Ymax): continue
                lab = "big" if Y == Ymax else f"{Y/K:.2g}K"
                if dn == designs[0]:
                    a2 = np.abs(Aout[Y][:, win]) ** 2; nz = Aout[(Y, "noise")][:, win]
                    rec.setdefault("signal_fraction", {})[lab] = dict(Y=Y, f=float(np.sum(a2) / np.sum(nz) - 1))
                for sig in range(12):
                    u = Aout[Y][sig][win]
                    if np.sum(np.abs(u)) == 0: continue
                    tests.setdefault(f"A(Y={lab})", []).append(stats_T(u, Sk, Dk, V, rc))
                    tests.setdefault(f"A*zeta(Y={lab})", []).append(stats_T(u * zeta[win], Sk, Dk, V, rc))
            # cross-scale estimate of the Y-stable (rho) component: lambda = |T(rho)|^2 / random-model size
            bigY = [Y for Y in Ys if Y >= Ymax // 4]
            for zlab, zf in (("rho", np.ones(win.sum())), ("pi=rho*zeta", zeta[win])):
                Tc, Er, fr = {}, {}, {}
                for Y in bigY:
                    a2 = np.abs(Aout[Y][:, win]) ** 2; nz = Aout[(Y, "noise")][:, win]
                    fr[Y] = max(1e-3, float(np.sum(a2) / np.sum(nz) - 1))
                    T = np.zeros((12, 12), complex); er = 0.0
                    for sig in range(12):
                        u = Aout[Y][sig][win] * zf
                        np.add.at(T[sig], rc, V * u * Sk)
                        er += float(np.sum(V ** 2 * np.abs(u) ** 2 * np.abs(Sk) ** 2))
                    Tc[Y] = T; Er[Y] = er
                lam = []
                for i1 in range(len(bigY)):
                    for i2 in range(i1 + 1, len(bigY)):
                        Y1, Y2 = bigY[i1], bigY[i2]
                        x = float(np.real(np.sum(Tc[Y1] * np.conj(Tc[Y2]))))
                        lam.append(x / math.sqrt(fr[Y1] * fr[Y2] * Er[Y1] * Er[Y2]))
                tests[f"cross-scale lambda({zlab})"] = [(float(np.mean(lam)), float(np.mean(lam)) * 0 + float("nan"))]
            for lab, u in [("zeta", zeta[win]), ("beta", beta_r[win]), ("conj beta", np.conj(beta_r[win])),
                           ("mu", rows["mu"][win].astype(complex)), ("one", np.ones(win.sum(), complex)),
                           ("alpha", rows["alpha"][win]), ("alpha^-2", np.conj(rows["alpha"][win]) ** 2),
                           ("alpha^3", rows["alpha"][win] ** 3)]:
                tests[lab] = [stats_T(u, Sk, Dk, V, rc)]
            tests["random"] = [stats_T(np.exp(2j * np.pi * rng.random(win.sum())), Sk, Dk, V, rc) for _ in range(nrand)]
            # adversarial reference: u = conj S (Cauchy-Schwarz extremal) -> eta = 1 + rho
            tests["adversary conj(S)"] = [stats_T(np.conj(Sk), Sk, Dk, V, rc)]
            rec["tests"][dn] = {k: dict(R_mean=float(np.mean([t[0] for t in v])), R_max=float(np.max([t[0] for t in v])),
                                        eta_mean=float(np.nanmean([t[1] for t in v])) if not k.startswith("cross") else None,
                                        eta_max=float(np.nanmax([t[1] for t in v])) if not k.startswith("cross") else None, n=len(v))
                                for k, v in tests.items()}
        out["windows"].append(rec)
        tp = rec["tests"][designs[1]]
        print(f"   signal fraction of A on window rows: " + ", ".join(f"{k}(Y={v['Y']}): {v['f']:+.3f}" for k, v in rec["signal_fraction"].items()))
        print(f"F={F:<4g} K={K:<8.0f} rows={rec['rows']:<6d} [{designs[1]}] " +
              "  ".join(f"{k}: R={v['R_mean']:.2f}/{v['R_max']:.2f}" + (f" eta={v['eta_mean']:.1e}" if v['eta_mean'] is not None else "") for k, v in tp.items()), flush=True)
    json.dump(out, open(os.path.join(HERE, f"probe_N{N0}.json"), "w"), indent=1)
    print(f"done {time.time()-t0:.1f}s")


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    run(**args)
