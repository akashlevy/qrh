"""
R3 numerical illustration of the coefficient-blindness obstruction.

For columns n squarefree primary with D < N(n) <= 2D (exact sextic symbols in Z[w]) and rows
0 < N(u) <= H = D^kappa, compute the per-row energy

    E_c(H) = (1/#rows) * sum_u |D^{-1/2} sum_n c(n) chi_n(u) w(N(n)/D)|^2      (normalised by diagonal)

for three coefficient vectors of identical modulus |c(n)| = 1:
    mu   : c(n) = mu(n)                 (the twisted-Mobius family A_u(D) of M_kappa)
    one  : c(n) = 1                     (plain; nu trivial)
    rand : c(n) = random +-1            (generic)
and split each sum into the sixth-power rows u = a^6 (the target copies) and all other rows.

Prediction (coefficient-blind obstruction): for c = 1 the H^{1/6} copies each carry |sum w|^2/D ~ D,
so the copy part alone is ~ H^{1/6} D, exceeding the diagonal budget H once H^{5/6} < D (kappa < 6/5).
Any inequality valid for all |c(n)| <= 1 must therefore allow the copies to be large; only the Mobius
structure keeps them small.  (This is a sanity check of a trivial construction, not evidence for M_kappa.)
"""
import math, sys
import numpy as np
sys.path.insert(0, "/Users/akashlevy/Documents/riemann/strengthening")
import wall_experiment as we

def run(D=3000, kappas=(0.6, 6/7, 1.0, 1.2), seed=1):
    HMAX = int(D ** max(kappas)) + 1
    primes = [p for p in we.eisenstein_primes(2 * D + 10) if p[2] != 2]
    fields = {pr[1]: we.ResidueField(*pr) for pr in primes}
    pnorm = {pr[1]: (pr[2] if pr[0] == "split" else pr[2] ** 2) for pr in primes}
    plist = sorted(fields, key=lambda p: pnorm[p])
    cols = []
    for e in we.elements_upto(2 * D):
        if we.norm(e) <= D: continue
        pe = we.primary(e)
        if pe is None or pe != e: continue
        rem = e; fac = []; ok = True; Nn = we.norm(e)
        for p in plist:
            if pnorm[p] > Nn: break
            F = fields[p]; red = F.red(rem)
            if (red == 0 if F.kind == "split" else red == (0, 0)):
                pc = (p[0] - p[1], -p[1]); prod = we.mul(rem, pc)
                rem = (prod[0] // pnorm[p], prod[1] // pnorm[p]); fac.append(p); Nn //= pnorm[p]
                red = F.red(rem)
                if (red == 0 if F.kind == "split" else red == (0, 0)): ok = False; break
                if Nn == 1: break
        if ok and Nn == 1: cols.append((e, fac))
    rng = np.random.default_rng(seed)
    coefs = {"mu": np.array([(-1) ** len(f) for _, f in cols], float),
             "one": np.ones(len(cols)),
             "rand": rng.choice([-1.0, 1.0], size=len(cols))}
    Ncol = np.array([we.norm(e) for e, _ in cols], float)
    wts = np.sin(math.pi * (Ncol / D - 1)) ** 2
    diag = float(np.sum(wts ** 2)) / D
    rows = we.elements_upto(HMAX)
    ra = np.array([r[0] for r in rows], np.int64); rb = np.array([r[1] for r in rows], np.int64)
    rN = ra * ra - ra * rb + rb * rb; o = np.argsort(rN); ra, rb, rN = ra[o], rb[o], rN[o]
    # sixth-power rows u = a^6 (a primary or not: the ideal (a) and unit 1 * a^6; any unit times a^6 with
    # unit^6 = 1 gives the SAME element a^6 up to the choice of a, so just collect {a^6})
    six = set()
    R6 = int(HMAX ** (1 / 6)) + 2
    for a in we.elements_upto(R6 * R6):
        if we.norm(a) ** 6 > HMAX: continue
        x = (1, 0)
        for _ in range(6): x = we.mul(x, a)
        six.add(x)
    is6 = np.array([(int(a), int(b)) in six for a, b in zip(ra, rb)])
    need = sorted({p for _, f in cols for p in f}, key=lambda p: pnorm[p])
    tab = {p: we.sextic_table(fields[p], ra, rb) for p in need}
    S = {k: np.zeros(len(ra), complex) for k in coefs}
    for j, (e, fac) in enumerate(cols):
        ex = np.zeros(len(ra), np.int16); mask = np.ones(len(ra), bool)
        for p in fac:
            t = tab[p]; mask &= (t >= 0); ex += t
        ph = np.where(mask, np.exp(1j * math.pi / 3 * ex), 0) * wts[j]
        for k in coefs: S[k] += coefs[k][j] * ph
    for k in S: S[k] /= math.sqrt(D)
    print(f"D={D}: {len(cols)} squarefree primary columns; diagonal per-row energy {diag:.4f};"
          f" rows up to {HMAX}; sixth-power rows found: {int(is6.sum())}")
    print(f"  {'kappa':>6s} {'H':>8s} {'rows':>7s} {'#a^6':>5s} | " +
          " | ".join(f"{k:>4s}: total/diag  copies/diag  copies/(H^(1/6)*D)" for k in S))
    for kap in kappas:
        H = D ** kap; i = np.searchsorted(rN, H, side="right")
        n6 = int(is6[:i].sum())
        parts = []
        for k in S:
            E = np.abs(S[k][:i]) ** 2
            tot = E.sum() / (i * diag)
            cop = E[is6[:i]].sum() / (i * diag)
            pred = E[is6[:i]].sum() / (max(n6, 1) * D)      # per-copy energy relative to D
            parts.append(f"{tot:10.3f}  {cop:10.3f}  {pred:14.3f}")
        print(f"  {kap:6.4f} {H:8.0f} {i:7d} {n6:5d} | " + " | ".join(parts))

if __name__ == "__main__":
    run(D=int(sys.argv[1]) if len(sys.argv) > 1 else 3000)
