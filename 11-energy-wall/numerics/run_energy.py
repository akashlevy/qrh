"""
numerics / steps 2-3: marked reflected energies beyond the wall, exact symbols in Z[w].

One pass of row sums S_k = sum_n c(n) psi_k(n) for many coefficient vectors c, rows k sf primary odd,
N(k) <= KMAX, columns n sf primary odd, N0 < N(n) <= 2 N0, w = sin^2(pi(N/N0 - 1)); exact n-diagonal
D_k = sum_n |c|^2 1_{(n,k)=1}.  Saved to data/energy_N{N0}.npz for analyze.py.

Coefficient families (names):
  beta                      beta w                                    (unmarked theta)
  rand{s}                   e^{i th} w                                (random-phase control, |c| = |beta w|)
  bP{q} / rP{q}             beta chi_P^{-2} w  /  e^{i th} w 1_{(n,P)=1}     single-mark energies E_P (+ control)
  g{j}                      beta chi_P^{-2} w for the Gram slot primes (per-mark sums S_{k,P}; (H2) Gram)
  X:{slot}:{design}         beta X_d w, X_d = sum_P c_{d,P} chi_P^{-2}, with the exact mask 1_{(P,k)=1}
  XR:{slot}:{design}        e^{i th} |beta X_d w|                     (random-phase control, no mask)
Designs (c_P, all times q_P^{-1/2}):
  p78   conj(alpha(P))^2 gamma_2(P)   [P78 Cor. marked-reflection, one active mark, eta = 1, c_F = 1;
                                       chi_P(lambda^3)^2 = 1 identically]
  gauss gamma_2(P) (= tau^+_{P,2});  ones 1;  kub conj gamma_1(P);  rnd random unit phases
  2-slot "p78x2": conj(alpha(P1P2))^2 gamma_2(P1) gamma_2(P2) chi_P1(P2)^2 chi_P2(P1)^2 (theta cocycle),
         "rndx2": random unit phase per tuple.
"""
import sys, os, json, time, math
import numpy as np
from nlib import Arith, Cols, Rows, rowsums, mark_chi_m2, gamma_prime, W

HERE = os.path.dirname(os.path.abspath(__file__))


def pick_prime(ar, target, used):
    c = [i for i, p in enumerate(ar.P) if p.kind == 0 and p.N >= target and i not in used]
    return min(c, key=lambda i: (ar.P[i].N, i))


def slot_primes(ar, L0):
    return [i for i, p in enumerate(ar.P) if p.kind == 0 and L0 <= p.N < 2 * L0]


def psi_rows(rows, cols, idx, nth=2):
    """Psi for an arbitrary list of row indices (one C call)."""
    from nlib import _lib
    idx = np.asarray(idx, np.int64)
    lens = rows.ptr[idx + 1] - rows.ptr[idx]
    ptr = np.zeros(len(idx) + 1, np.int64); ptr[1:] = np.cumsum(lens)
    sel = np.concatenate([np.arange(rows.ptr[i], rows.ptr[i + 1]) for i in idx]).astype(np.int64) if len(idx) else np.zeros(0, np.int64)
    out = np.empty((len(idx), cols.n), np.float64)
    _lib.build_psi(cols.n, cols.A, cols.B, cols.N, len(idx), ptr, np.ascontiguousarray(rows.fq[sel]),
                   np.ascontiguousarray(rows.fr[sel]), np.ascontiguousarray(rows.fqi[sel]), rows.bits,
                   rows.off[:-1].copy(), out, nth)
    return out


def run(N0=20000, Kfac=1.0, singles=(7, 13, 31, 61, 97, 151, 211, 307), gram_L0=None,
        slots=(100, 400), slots2=(), designs=("p78", "gauss", "ones", "kub", "rnd"), seeds=4,
        nth=2, tag=None):
    t0 = time.time()
    KMAX = int(Kfac * N0)
    nmax = max(2 * N0, KMAX) + 10
    nmax = int(math.ceil(nmax / 10000) * 10000)
    ar = Arith(nmax)
    cols = Cols(ar, N0); rows = Rows(ar, KMAX)
    print(f"N0={N0}: cols {cols.n}, rows {rows.n} (N(k) <= {KMAX})  {time.time()-t0:.1f}s", flush=True)
    rng = np.random.default_rng(1000 + N0)
    names, C, dnames, Dc, dmap = [], [], [], [], {}
    def addD(name, vec):
        dnames.append(name); Dc.append(vec); return len(dnames) - 1
    def add(name, vec, dname):
        names.append(name); C.append(vec); dmap[name] = dname
    w = cols.w
    iw2 = addD("w2", w ** 2)
    add("beta", cols.beta * w, "w2")
    for s in range(seeds):
        add(f"rand{s}", np.exp(2j * np.pi * rng.random(cols.n)) * w, "w2")
    # single marks
    used = set(); single_info = []
    for q in singles:
        i = pick_prime(ar, q, used); used.add(i); pr = ar.P[i]
        chi = mark_chi_m2(pr, cols)
        nz = (chi != 0).astype(float)
        addD(f"w2P{pr.N}", w ** 2 * nz)
        add(f"bP{pr.N}", cols.beta * chi * w, f"w2P{pr.N}")
        add(f"rP{pr.N}", np.exp(2j * np.pi * rng.random(cols.n)) * w * nz, f"w2P{pr.N}")
        single_info.append(dict(idx=i, N=pr.N, a=pr.a, b=pr.b))
    # Gram slot (per-mark sums)
    gram_info = []
    if gram_L0:
        for j, i in enumerate(slot_primes(ar, gram_L0)):
            pr = ar.P[i]; chi = mark_chi_m2(pr, cols)
            dn = f"w2P{pr.N}_{i}"
            addD(dn, w ** 2 * (chi != 0))
            add(f"g{j}", cols.beta * chi * w, dn)
            gram_info.append(dict(idx=i, N=pr.N))
    # designs
    slot_specs = [("s%d" % L0, [L0]) for L0 in slots] + [("t%dx%d" % (a, b), [a, b]) for a, b in slots2]
    design_data = {}
    for sname, Ls in slot_specs:
        if len(Ls) == 1:
            ps = slot_primes(ar, Ls[0])
            M = np.stack([mark_chi_m2(ar.P[i], cols) for i in ps], axis=1)             # cols x L'
            al = np.array([complex(ar.P[i].a, 0) + ar.P[i].b * W for i in ps]); al /= np.abs(al)
            g2 = np.array([gamma_prime(ar.P[i], 2) for i in ps]); g1 = np.array([gamma_prime(ar.P[i], 1) for i in ps])
            qn = np.array([ar.P[i].N for i in ps], float)
            cdict = {"p78": np.conj(al) ** 2 * g2, "gauss": g2, "ones": np.ones(len(ps), complex),
                     "kub": np.conj(g1), "rnd": np.exp(2j * np.pi * rng.random(len(ps)))}
            cs = {d: cdict[d] / np.sqrt(qn) for d in designs}
            design_data[sname] = dict(kind=1, primes=ps, M=M, c=cs)
            for d, c in cs.items():
                X = M @ c
                bx = cols.beta * X * w
                dn = addD(f"D:{sname}:{d}", np.abs(bx) ** 2)
                add(f"X:{sname}:{d}", bx, f"D:{sname}:{d}")
                add(f"XR:{sname}:{d}", np.exp(2j * np.pi * rng.random(cols.n)) * np.abs(bx), f"D:{sname}:{d}")
        else:
            p1, p2 = slot_primes(ar, Ls[0]), slot_primes(ar, Ls[1])
            p2 = [i for i in p2 if i not in p1]
            M1 = np.stack([mark_chi_m2(ar.P[i], cols) for i in p1], axis=1)
            M2 = np.stack([mark_chi_m2(ar.P[i], cols) for i in p2], axis=1)
            def al(i): z = complex(ar.P[i].a, 0) + ar.P[i].b * W; return z / abs(z)
            g2a = np.array([gamma_prime(ar.P[i], 2) for i in p1]); g2b = np.array([gamma_prime(ar.P[i], 2) for i in p2])
            # cocycle chi_P1(P2)^2 chi_P2(P1)^2 and conj(alpha(P1 P2))^2
            coc = np.empty((len(p1), len(p2)), complex)
            for a_, i in enumerate(p1):
                pi_ = ar.P[i]
                kA = ar.P[i].sextic_exp(np.array([ar.P[j].a for j in p2]), np.array([ar.P[j].b for j in p2])).astype(np.int64)
                for b_, j in enumerate(p2):
                    kB = int(ar.P[j].sextic_exp(pi_.a, pi_.b)[0])
                    coc[a_, b_] = np.exp(2j * np.pi * 2 * (int(kA[b_]) + kB) / 6)
            qq = np.outer([ar.P[i].N for i in p1], [ar.P[j].N for j in p2]).astype(float)
            alp = np.outer([np.conj(al(i)) ** 2 for i in p1], [np.conj(al(j)) ** 2 for j in p2])
            cs = {"p78x2": alp * np.outer(g2a, g2b) * coc / np.sqrt(qq),
                  "rndx2": np.exp(2j * np.pi * rng.random(qq.shape)) / np.sqrt(qq)}
            design_data[sname] = dict(kind=2, p1=p1, p2=p2, M1=M1, M2=M2, c=cs)
            for d, c in cs.items():
                X = np.sum((M1 @ c) * M2, axis=1)
                bx = cols.beta * X * w
                addD(f"D:{sname}:{d}", np.abs(bx) ** 2)
                add(f"X:{sname}:{d}", bx, f"D:{sname}:{d}")
                add(f"XR:{sname}:{d}", np.exp(2j * np.pi * rng.random(cols.n)) * np.abs(bx), f"D:{sname}:{d}")
    Cm = np.stack(C, axis=1); Dm = np.stack(Dc, axis=1)
    print(f"coefficients: {len(names)} complex, {len(dnames)} diagonal patterns  {time.time()-t0:.1f}s", flush=True)
    S, D = rowsums(rows, cols, Cm, Dm, nth=nth)
    dcol = np.array([dnames.index(dmap[n]) for n in names])
    Dfull = D[:, dcol]
    print(f"row sums done {time.time()-t0:.1f}s", flush=True)
    # exact mask 1_{(P,k)=1} for the designs: recompute rows that have slot-prime factors
    rowset = [set(f) for f in rows.fac]
    for sname, dd in design_data.items():
        allp = set(dd["primes"]) if dd["kind"] == 1 else set(dd["p1"]) | set(dd["p2"])
        groups = {}
        for k in range(rows.n):
            hit = tuple(sorted(rowset[k] & allp))
            if hit: groups.setdefault(hit, []).append(k)
        nfix = 0
        for hit, idx in groups.items():
            Psi = psi_rows(rows, cols, idx, nth)
            hs = set(hit)
            Xs = []
            for d, c in dd["c"].items():
                if dd["kind"] == 1:
                    keep = np.array([0.0 if p in hs else 1.0 for p in dd["primes"]])
                    X = dd["M"] @ (c * keep)
                else:
                    k1 = np.array([0.0 if p in hs else 1.0 for p in dd["p1"]])
                    k2 = np.array([0.0 if p in hs else 1.0 for p in dd["p2"]])
                    X = np.sum((dd["M1"] @ (c * k1[:, None] * k2[None, :])) * dd["M2"], axis=1)
                Xs.append((d, cols.beta * X * w))
            B = np.stack([x for _, x in Xs], axis=1)
            Sg = Psi @ B; Dg = np.abs(Psi) @ (np.abs(B) ** 2)
            for j, (d, _) in enumerate(Xs):
                c_ = names.index(f"X:{sname}:{d}")
                S[idx, c_] = Sg[:, j]; Dfull[idx, c_] = Dg[:, j]
            nfix += len(idx)
        print(f"mask {sname}: {nfix} rows with slot-prime factors recomputed ({len(groups)} groups)  {time.time()-t0:.1f}s", flush=True)
    tag = tag or f"N{N0}"
    out = os.path.join(HERE, "data", f"energy_{tag}.npz")
    rowmarks = {}
    for s in single_info: rowmarks[f"bP{s['N']}"] = s["idx"]
    np.savez(out, S=S.astype(np.complex128), D=Dfull, names=np.array(names), rowN=rows.N, rowA=rows.A, rowB=rows.B,
             rowarg=rows.arg, rowmu=rows.mu, rowg2=rows.gamma2, rowalpha=rows.alpha,
             rowfac=np.array([";".join(map(str, f)) for f in rows.fac]),
             meta=json.dumps(dict(N0=N0, KMAX=KMAX, cols=cols.n, rows=rows.n, singles=single_info, gram=gram_info,
                                  slots={k: (dict(primes=[int(x) for x in v["primes"]]) if v["kind"] == 1 else
                                             dict(p1=[int(x) for x in v["p1"]], p2=[int(x) for x in v["p2"]]))
                                         for k, v in design_data.items()},
                                  Ln={k: [ar.P[i].N for i in (v["primes"] if v["kind"] == 1 else v["p1"] + v["p2"])]
                                      for k, v in design_data.items()},
                                  colsum_w2=float(np.sum(w ** 2)), time=time.time() - t0)))
    print(f"saved {out}  total {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    run(**args)
