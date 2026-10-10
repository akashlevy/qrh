"""xmark: original-side numerics for the cross-mark terms of lem:marked (exact Z[w] arithmetic).

Rows u: all nonzero u in Z[w] with N(u) <= U, weight Phi = (1 - N(u)/U)^2.
Inverse witness M_u = R^{-1/2} sum_{n sf primary,(n,6)=1, R<N(n)<=2R} mu(n) W(N(n)/R) (u/n)_6.
Slot primes p: split primary, P < N(p) <= 2P.
Cross terms  T(p,p') = sum_u Phi |M_u|^2 (u/p)_6 conj (u/p')_6   (exact: sum_u|M Q|^2 = P^{-1} sum a a' T).
Reports:  m2 = sum Phi|M|^2 / sum Phi,  m4 = sum Phi|M|^4 / sum Phi   (pair-large-sieve input),
          rms of O = T(p,p')/sqrt(T(p,p)T(p',p')) over p != p' vs. random-phase control,
          coherent ratio  sum_{p!=p'} T / sum_p T(p,p)  (a = 1 slot),
          correlations of O with pair phases (p/p')_3, gamma_2(p) conj gamma_2(p').
Reuses the 09 Z[w] helpers (read only)."""
import sys, math, json
import numpy as np
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/energy-wall")
sys.path.insert(0, "/home/user/qrh/09-new-directions/agents/edge-theta")
from eisenstein import eisenstein_primes, ResidueField, mul, Z6   # noqa
import dual_energy as de                                            # noqa
from edge_lib import gauss_vec                                      # noqa

def rows(U):
    R = int(2 * math.sqrt(U)) + 3
    a = np.arange(-R, R + 1); A, B = np.meshgrid(a, a, indexing="ij"); A = A.ravel(); B = B.ravel()
    Nn = A * A - A * B + B * B
    s = (Nn > 0) & (Nn <= U)
    return A[s].astype(np.int64), B[s].astype(np.int64), Nn[s].astype(np.int64)

def main(U=40000, r=0.7, P=200, seed=1):
    rng = np.random.default_rng(seed)
    R = int(round(U ** r))
    A, B, Nu = rows(U); Phi = (1 - Nu / U) ** 2
    primes = eisenstein_primes(2 * max(R, P) + 10)
    fields = {p[1]: ResidueField(*p) for p in primes}
    pnorm = {p[1]: (p[2] if p[0] == "split" else p[2] ** 2) for p in primes}
    plist = sorted(fields, key=lambda p: pnorm[p])
    # columns
    cA, cB, cN, fac = de.primary_sf_elements(R, 2 * R, fields, plist, pnorm)
    keep = [i for i, f in enumerate(fac) if all(pnorm[p] != 4 for p in f)]
    fac = [fac[i] for i in keep]; cN = cN[keep]
    mu = np.array([(-1) ** len(f) for f in fac], float)
    Wc = np.sin(np.pi * (cN / R - 1)) ** 2
    need = sorted({p for f in fac for p in f}, key=lambda p: pnorm[p])
    expo = {}
    for p in need:
        expo[p] = de.sextic_exp_vec(fields[p], A, B).astype(np.int16)
    z6 = np.exp(2j * np.pi * np.arange(6) / 6)
    M = np.zeros(len(A), complex)
    for i, f in enumerate(fac):
        e = np.zeros(len(A), np.int16); zero = np.zeros(len(A), bool)
        for p in f:
            ep = expo[p]; zero |= ep < 0; e += np.where(ep < 0, 0, ep)
        v = z6[e % 6]; v[zero] = 0
        M += mu[i] * Wc[i] * v
    M /= math.sqrt(R)
    w2 = Phi * np.abs(M) ** 2
    m2 = w2.sum() / Phi.sum(); m4 = (Phi * np.abs(M) ** 4).sum() / Phi.sum()
    # slot primes
    sp = [p for p in plist if fields[p].kind == "split" and P < pnorm[p] <= 2 * P]
    X = np.zeros((len(sp), len(A)), complex)
    for j, p in enumerate(sp):
        ep = de.sextic_exp_vec(fields[p], A, B)
        v = z6[np.where(ep < 0, 0, ep) % 6]; v[ep < 0] = 0
        X[j] = v
    T = (X * w2) @ X.conj().T           # T[p,p'] = sum_u w2 chi_p conj chi_p'
    d = np.real(np.diag(T)); O = T / np.sqrt(np.outer(d, d))
    off = ~np.eye(len(sp), dtype=bool)
    rms = float(np.sqrt(np.mean(np.abs(O[off]) ** 2)))
    # random-phase control: same |M|^2 weights, random sextic values for the slot characters
    Xr = z6[rng.integers(0, 6, size=X.shape)] * (np.abs(X) > 0)
    Tr = (Xr * w2) @ Xr.conj().T; dr = np.real(np.diag(Tr)); Or = Tr / np.sqrt(np.outer(dr, dr))
    rms_r = float(np.sqrt(np.mean(np.abs(Or[off]) ** 2)))
    # predicted rms of a random sum: sqrt(sum w2^2)/sum w2
    pred = float(np.sqrt((w2 ** 2).sum()) / w2.sum())
    coh = float(np.real(T[off].sum()) / d.sum())
    coh_r = float(np.real(Tr[off].sum()) / dr.sum())
    # structured pair phases
    g2 = np.array([gauss_vec(fields[p], 2) for p in sp])
    cub = np.zeros((len(sp), len(sp)), complex)
    for i, p in enumerate(sp):
        for j, q in enumerate(sp):
            if i != j:
                k = fields[q].sextic_exp(p)      # (p/q)_6 ; cubic symbol = square
                cub[i, j] = Z6[(2 * k) % 6] if k is not None else 0
    def corr(Y):
        y = Y[off]; o = O[off]
        return float(np.abs(np.vdot(y, o)) / np.sqrt(np.vdot(y, y).real * np.vdot(o, o).real))
    gg = np.outer(g2, g2.conj())
    res = dict(U=U, r=r, R=R, P=P, rows=int(len(A)), cols=len(fac), nslot=len(sp),
               m2=round(float(m2), 4), m4=round(float(m4), 4), m4_over_m2sq=round(float(m4 / m2 ** 2), 3),
               rms_O=round(rms, 5), rms_O_random=round(rms_r, 5), rms_pred=round(pred, 5),
               coherent_ratio=round(coh, 4), coherent_ratio_random=round(coh_r, 4),
               corr_cubic=round(corr(cub), 4), corr_gauss=round(corr(gg), 4),
               null_corr=round(1 / math.sqrt(off.sum()), 4))
    return res

if __name__ == "__main__":
    for U, r, P in [(20000, 0.7, 150), (40000, 0.7, 150), (80000, 0.7, 150), (80000, 0.7, 400)]:
        print(json.dumps(main(U, r, P)), flush=True)
