"""Like onset.py, but rows are (h, nu) for a list of ray-class characters nu mod 9 (nu_{u,v}), applied to
every family (B, K, K0).  Optionally multiplies the cubic family K by its Kubota normaliser
L(3s - 1/2, alpha^3 nu^3) (family 'Kc+'/'Kc-'), i.e. coefficients sum_{m^3 | n} alpha(m)^3 nu(m)^3 N(m)^{1/2} k(n/m^3).

usage: python3 onset2.py prefix Xmax norms out.npz fams [nus]
   norms: comma list of allowed N(h) for prime h (e.g. 4,7,13,19) or 'max:<H>'; composite rows 'c:<N>' pairs
   nus  : e.g. 00,01,02,10,11,12,20,21,22   (default all 9)
"""
import sys, time, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__file__))
import cbn_lib as L
import onset as O


def main():
    prefix, Xmax, normspec, out, fams = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5].split(",")
    nus = [(int(s[0]), int(s[1])) for s in (sys.argv[6].split(",") if len(sys.argv) > 6 else
           ["00", "01", "02", "10", "11", "12", "20", "21", "22"])]
    t0 = time.time()
    P, S = L.load(prefix)
    S = S[S["N"] <= Xmax]
    A, B, Nn = S["A"], S["B"], S["N"]
    G = S["gre"] + 1j * S["gim"]
    al = L.cval(A, B) / np.sqrt(Nn)
    binS = O.Binner(Nn)
    # rows h
    H = []
    for tok in normspec.split(","):
        if tok.startswith("max:"):
            Hm = int(tok[4:]); H += [(i,) for i in range(len(P)) if 4 <= P["N"][i] <= Hm]
        elif tok.startswith("c"):
            n1, n2 = map(int, tok[1:].split("x"))
            I1 = [i for i in range(len(P)) if P["N"][i] == n1]; I2 = [i for i in range(len(P)) if P["N"][i] == n2]
            H += [(i, j) for i in I1 for j in I2 if i < j or n1 != n2]
        elif tok == "1":
            H.append(())
        else:
            H += [(i,) for i in range(len(P)) if P["N"][i] == int(tok)]
    # cube-completion support for Kc: pairs (n', m) with N(n') N(m)^3 <= Xmax, m primary coprime to 3
    needc = any(f.startswith("Kc") for f in fams)
    if needc:
        Al = np.fromfile(prefix + "_all.bin", dtype=np.dtype([("A", "<i4"), ("B", "<i4"), ("N", "<i4"), ("tau", "<i4")]))
        Ms = Al[(Al["N"] > 1) & (Al["N"] <= int(round(Xmax ** (1 / 3))) + 1)]; Ms = Ms[Ms["N"].astype(np.int64) ** 3 <= Xmax]
        ext = []   # for each m: indices of n' and the product coordinates
        for m in Ms:
            Nm3 = int(m["N"]) ** 3
            sel = np.nonzero(Nn <= Xmax // Nm3)[0]
            a3, b3 = int(m["A"]), int(m["B"])
            a3, b3 = L.mulw(*L.mulw(a3, b3, a3, b3), a3, b3)
            pa, pb = L.mulw(A[sel].astype(np.int64), B[sel].astype(np.int64), a3, b3)
            ext.append((sel, Nn[sel].astype(np.int64) * Nm3, pa, pb, int(m["A"]), int(m["B"]), int(m["N"])))
        e_sel = np.concatenate([e[0] for e in ext]); e_N = np.concatenate([e[1] for e in ext])
        e_mA = np.concatenate([np.full(len(e[0]), e[4]) for e in ext]); e_mB = np.concatenate([np.full(len(e[0]), e[5]) for e in ext])
        e_mN = np.concatenate([np.full(len(e[0]), e[6]) for e in ext])
        e_al3 = (L.cval(e_mA, e_mB) / np.sqrt(e_mN)) ** 3 * np.sqrt(e_mN)
        binE = O.Binner(e_N); nb = max(binS.nb, binE.nb); binE.nb = nb; binS.nb = nb
        print(f"# cube completion: {len(Ms)} m, {len(e_sel)} extra terms", flush=True)
    nb = binS.nb
    Wn = ["bump", "lg35", "lg60"]
    Xs = np.exp(np.arange(np.log(3.0), np.log(Xmax / 2.0), 1 / 12))
    Wm = O.wmats(Xs, nb, Wn)
    valid = {nm: Xs * np.exp(O.UMAX[nm]) <= Xmax for nm in Wn}
    rows = [(h, nu) for h in H for nu in nus]
    print(f"# {len(S)} sf n, {len(H)} h x {len(nus)} nu = {len(rows)} rows, {len(Xs)} X, load {time.time()-t0:.1f}s", flush=True)
    res = {f: {nm: np.zeros((len(rows), len(Xs)), complex) for nm in Wn} for f in fams}
    dia = {f: {nm: np.zeros((len(rows), len(Xs))) for nm in Wn} for f in fams}
    Nh = np.zeros(len(rows)); nuarr = np.zeros((len(rows), 2), int); hid = np.zeros(len(rows), int)
    mods = {}
    nuv = {nu: L.nu_char(A, B, *nu) for nu in nus}
    if needc:
        nue = {nu: L.nu_char(e_mA, e_mB, *nu) ** 3 for nu in nus}   # = 1 (cubic), kept for clarity
    k = 0
    for ih, h in enumerate(H):
        for i in h:
            if i not in mods: mods[i] = L.Mod(P[i])
        j = L.chi_exp([mods[i] for i in h], A, B); ok = j >= 0
        ch = np.where(ok, L.W3 ** np.where(ok, j, 0), 0)
        for nu in nus:
            Nh[k] = np.prod([P["N"][i] for i in h]); nuarr[k] = nu; hid[k] = ih
            for f in fams:
                orient = f[-1]; tw = np.conj(ch) if orient == "+" else ch
                if f[:-1] == "B": c = S["mu"] * G * nuv[nu] * tw
                elif f[:-1] in ("K", "Kc"): c = al * G * nuv[nu] * tw
                elif f[:-1] == "K0": c = G * nuv[nu] * tw
                s0, s1, d0, d1 = binS.sums(c)
                if f.startswith("Kc"):
                    ce = c[e_sel] * e_al3 * nue[nu]
                    t0_, t1_, u0_, u1_ = binE.sums(ce)
                    s0 = s0 + t0_; s1 = s1 + t1_; d0 = d0 + u0_; d1 = d1 + u1_
                for nm in Wn:
                    w, wp, w2, w2p = Wm[nm]
                    res[f][nm][k] = w @ s0 + wp @ s1
                    dia[f][nm][k] = w2 @ d0 + w2p @ d1
            k += 1
        print(f"  h {ih+1}/{len(H)} N={Nh[k-1]:.0f} {time.time()-t0:.1f}s", flush=True)
    sv = dict(Xs=Xs, Nh=Nh, nu=nuarr, hid=hid, H=np.array([str(h) for h in H]), fams=np.array(fams),
              Wn=np.array(Wn), Xmax=Xmax, kind=np.array([10 * max(1, len(H[hid[r]])) for r in range(len(rows))]))
    for f in fams:
        for nm in Wn:
            sv[f"T_{f}_{nm}"] = res[f][nm]; sv[f"D_{f}_{nm}"] = dia[f][nm]; sv[f"valid_{nm}"] = valid[nm]
    np.savez_compressed(out, **sv)
    print(f"# done {time.time()-t0:.1f}s -> {out}")


if __name__ == "__main__":
    main()
