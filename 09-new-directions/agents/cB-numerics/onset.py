"""Smoothed sums T_f(h, X) = sum_n c_f(n) chi_h(n)^{-+} W(N(n)/X) and diagonals D = sum |c|^2 W^2
for many rows h and a log grid of X, for several test functions W.

Families (n over primary n, N(n) <= Xmax; squarefree families restricted to sf n, automatically
coprime to h through chi_h(n) = 0):
  B+  : mu nu G~(n) conj(chi_h(n))   = mu nu g~(1,n) chi-bar_n(h)   <- the B_h of the question
  B-  : mu nu G~(n) chi_h(n)           opposite orientation (sextic index h^2)
  K+  : alpha G~(n) conj(chi_h(n))   = alpha g~(h,n): cubic Kubota, index h   (control, theory c = 1)
  K-  : alpha G~(n) chi_h(n)           cubic Kubota, index h^2               (control, theory c = 2)
  K0+ : G~(n) conj(chi_h(n))           untwisted Kubota (theta pole at 5/6)
  H1  : chi_h(n) over ALL primary n     Hecke L(s, chi_h)       (Euler product, c = 1)
  H2  : tau(n) chi_h(n) over ALL n      L(s, chi_h)^2           (Euler product, c = 2)
where chi_h(n) = (n/h)_3 = (h/n)_3, alpha(n) = n/|n|, nu = ray class character mod 9 (nu_u, nu_v).
Test functions (in u = log x): 'bump' = exp(-1/(1-t^2)), t = (x-1.5)/0.5 on (1,2);
'lg35' = exp(-u^2/(2*0.35^2)) cut at |u| <= 1.4;  'lg60' = exp(-u^2/(2*0.6^2)) cut at |u| <= 2.4.
Binned evaluation in log N with first-order correction (bin width 0.002).

usage: python3 onset.py prefix Xmax Hmax out.npz [families] [nu_u nu_v] [comp_pairs]
"""
import sys, time, itertools, numpy as np
sys.path.insert(0, __import__("os").path.dirname(__file__))
import cbn_lib as L

DELTA = 0.002


def wfun(name, u):
    w = np.zeros_like(u); wp = np.zeros_like(u)
    if name == "bump":
        x = np.exp(np.clip(u, -50, 50)); t = (x - 1.5) / 0.5
        m = np.abs(t) < 1
        tt = t[m]; ww = np.exp(-1 / (1 - tt * tt))
        w[m] = ww
        wp[m] = ww * (-2 * tt / (1 - tt * tt) ** 2) * (x[m] / 0.5)
        return w, wp
    sig, cut = {"lg35": (0.35, 1.4), "lg60": (0.6, 2.4)}[name]
    m = np.abs(u) <= cut
    w[m] = np.exp(-u[m] ** 2 / (2 * sig * sig)); wp[m] = -u[m] / sig ** 2 * w[m]
    return w, wp


UMAX = {"bump": np.log(2.0), "lg35": 1.4, "lg60": 2.4}


class Binner:
    def __init__(self, Nn):
        u = np.log(Nn.astype(float))
        self.b = np.floor(u / DELTA).astype(np.int64)
        self.nb = int(self.b.max()) + 1
        self.du = u - (self.b + 0.5) * DELTA

    def sums(self, c):
        nb = self.nb; b = self.b
        s0 = np.bincount(b, c.real, nb) + 1j * np.bincount(b, c.imag, nb)
        cd = c * self.du
        s1 = np.bincount(b, cd.real, nb) + 1j * np.bincount(b, cd.imag, nb)
        a = np.abs(c) ** 2
        d0 = np.bincount(b, a, nb); d1 = np.bincount(b, a * self.du, nb)
        return s0, s1, d0, d1


def wmats(Xs, nb, names):
    ub = (np.arange(nb) + 0.5) * DELTA
    out = {}
    for nm in names:
        U = ub[None, :] - np.log(Xs)[:, None]
        w, wp = wfun(nm, U.ravel()); w = w.reshape(U.shape); wp = wp.reshape(U.shape)
        out[nm] = (w, wp, w * w, 2 * w * wp)
    return out


def hlist_build(P, Hmax, comp_pairs=True):
    """rows h: lists of prime indices. split primes (both pi, conj pi) and inert primes with N <= Hmax,
    plus a set of composites."""
    H = []
    for i in range(len(P)):
        if P["N"][i] <= Hmax and P["N"][i] >= 4:
            H.append((i,))
    if comp_pairs:
        small = [i for i in range(len(P)) if 4 <= P["N"][i] <= 200]
        for i, j in itertools.combinations(small, 2):
            Nh = int(P["N"][i]) * int(P["N"][j])
            if Nh <= min(Hmax, 3000):
                H.append((i, j))
    return H


def main():
    prefix = sys.argv[1]; Xmax = int(sys.argv[2]); Hmax = int(sys.argv[3]); out = sys.argv[4]
    fams = sys.argv[5].split(",") if len(sys.argv) > 5 else ["B+", "B-", "K+", "K-", "K0+", "H1", "H2"]
    nu_u, nu_v = (int(sys.argv[6]), int(sys.argv[7])) if len(sys.argv) > 7 else (0, 0)
    comp = (sys.argv[8] == "1") if len(sys.argv) > 8 else True
    t0 = time.time()
    P, S = L.load(prefix)
    S = S[S["N"] <= Xmax]
    A, B, Nn = S["A"], S["B"], S["N"]
    G = S["gre"] + 1j * S["gim"]
    al = L.cval(A, B) / np.sqrt(Nn)
    nu = L.nu_char(A, B, nu_u, nu_v)
    base = {"B": S["mu"] * G * nu, "K": al * G, "K0": G.copy()}
    binS = Binner(Nn)
    need_all = any(f.startswith("H") for f in fams)
    if need_all:
        Al = np.fromfile(prefix + "_all.bin", dtype=np.dtype([("A", "<i4"), ("B", "<i4"), ("N", "<i4"), ("tau", "<i4")]))
        Al = Al[Al["N"] <= Xmax]
        binA = Binner(Al["N"])
        tau = Al["tau"].astype(float)
    nb = max(binS.nb, binA.nb if need_all else 0)
    binS.nb = nb
    if need_all: binA.nb = nb
    Wn = ["bump", "lg35", "lg60"]
    Xs = np.exp(np.arange(np.log(3.0), np.log(Xmax / 2.0), 1 / 12))
    Wm = wmats(Xs, nb, Wn)
    valid = {nm: Xs * np.exp(UMAX[nm]) <= Xmax for nm in Wn}
    H = hlist_build(P, Hmax, comp)
    print(f"# {len(S)} sf n, {len(H)} rows h, {len(Xs)} X, load {time.time()-t0:.1f}s", flush=True)
    mods = {}
    res = {f: {nm: np.zeros((len(H), len(Xs)), complex) for nm in Wn} for f in fams}
    dia = {f: {nm: np.zeros((len(H), len(Xs))) for nm in Wn} for f in fams}
    Nh = np.zeros(len(H)); kind = np.zeros(len(H), int)
    for k, h in enumerate(H):
        for i in h:
            if i not in mods: mods[i] = L.Mod(P[i])
        mh = [mods[i] for i in h]
        Nh[k] = np.prod([P["N"][i] for i in h]); kind[k] = len(h) * 10 + int(P["kind"][h[0]])
        j = L.chi_exp(mh, A, B)
        ok = j >= 0
        ch = np.where(ok, L.W3 ** np.where(ok, j, 0), 0)
        if need_all:
            ja = L.chi_exp(mh, Al["A"], Al["B"]); oka = ja >= 0
            cha = np.where(oka, L.W3 ** np.where(oka, ja, 0), 0)
        for f in fams:
            if f.startswith("H"):
                c = cha if f == "H1" else tau * cha
                bn = binA
            else:
                bs = base[f[:-1]]
                c = bs * (np.conj(ch) if f.endswith("+") else ch)
                bn = binS
            s0, s1, d0, d1 = bn.sums(c)
            for nm in Wn:
                w, wp, w2, w2p = Wm[nm]
                res[f][nm][k] = w @ s0 + wp @ s1
                dia[f][nm][k] = w2 @ d0 + w2p @ d1
        if k % 50 == 0:
            print(f"  h {k}/{len(H)} N={Nh[k]:.0f} {time.time()-t0:.1f}s", flush=True)
    sv = dict(Xs=Xs, Nh=Nh, kind=kind, H=np.array([str(h) for h in H]),
              hab=np.array([[(int(P['a'][i]), int(P['b'][i])) for i in h] + [(0, 0)] * (2 - len(h)) for h in H]),
              fams=np.array(fams), Wn=np.array(Wn), Xmax=Xmax, nu=np.array([nu_u, nu_v]))
    for f in fams:
        for nm in Wn:
            sv[f"T_{f}_{nm}"] = res[f][nm]; sv[f"D_{f}_{nm}"] = dia[f][nm]
            sv[f"valid_{nm}"] = valid[nm]
    np.savez_compressed(out, **sv)
    print(f"# done {time.time()-t0:.1f}s -> {out}")


if __name__ == "__main__":
    main()
