"""Exact rationality test for D_{k,P}(u) = sum_nu S_{k,P}(nu) u^nu (squarefree n), S in Z[omega].
For each row: find the shortest linear recurrence (Berlekamp-Massey over Q(omega), exact) satisfied by
S(0..numax); report recurrence (denominator Q(u)) and numerator N(u) = Q(u) D(u) mod u^{numax+1}."""
import sys
from fractions import Fraction as Fr
from collections import Counter
import numpy as np
sys.path.insert(0, "/home/user/qrh/11-energy-wall/function-field")
from ffdriver import run, choose_marks

class W:   # element a + b*omega of Q(omega), omega^2 = -1 - omega
    __slots__ = ("a", "b")
    def __init__(s, a, b=0): s.a = Fr(a); s.b = Fr(b)
    def __add__(s, o): return W(s.a + o.a, s.b + o.b)
    def __sub__(s, o): return W(s.a - o.a, s.b - o.b)
    def __mul__(s, o): return W(s.a * o.a - s.b * o.b, s.a * o.b + s.b * o.a - s.b * o.b)
    def conj(s): return W(s.a - s.b, -s.b)            # omega -> omega^2 = -1-omega
    def norm(s): return s.a * s.a - s.a * s.b + s.b * s.b
    def inv(s): n = s.norm(); c = s.conj(); return W(c.a / n, c.b / n)
    def __truediv__(s, o): return s * o.inv()
    def iszero(s): return s.a == 0 and s.b == 0
    def __repr__(s): return f"({s.a}{'+' if s.b >= 0 else ''}{s.b}w)"

def berlekamp_massey(seq):
    """shortest C with C[0]=1, sum_{i=0}^L C[i] s[n-i] = 0 for n >= L."""
    Cpoly = [W(1)]; B = [W(1)]; L = 0; m = 1; b = W(1)
    for n in range(len(seq)):
        d = seq[n]
        for i in range(1, L + 1):
            d = d + Cpoly[i] * seq[n - i]
        if d.iszero():
            m += 1; continue
        T = Cpoly[:]
        coef = d / b
        need = len(B) + m
        while len(Cpoly) < need: Cpoly.append(W(0))
        for i in range(len(B)):
            Cpoly[i + m] = Cpoly[i + m] - coef * B[i]
        if 2 * L <= n:
            L = n + 1 - L; B = T; b = d; m = 1
        else:
            m += 1
    return Cpoly[:L + 1], L

def exact_seq(R):
    S = R["S"][:, 0]
    Bv = np.rint(S.imag * 2 / np.sqrt(3)).astype(np.int64)
    Av = np.rint(S.real + Bv / 2).astype(np.int64)
    return Av, Bv


if __name__ == "__main__":
    q, kappa, lam, numax, nrows = map(int, sys.argv[1:6])
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    marks = choose_marks(q, lam, 1)
    R0 = run(q, kappa, 0, lam, marks, tag=f"rt0_{q}_{kappa}_{lam}")
    allrows = [int(r) for r, c in zip(R0["kidx"], R0["cop"][:, 0]) if c]
    rng = np.random.default_rng(seed)
    rows = sorted(rng.choice(allrows, size=min(nrows, len(allrows)), replace=False).tolist())
    seqs = {r: [] for r in rows}
    for nu in range(numax + 1):
        R = run(q, kappa, nu, lam, marks, tag=f"rt{q}_{kappa}_{lam}_{nu}", rows=rows)
        Av, Bv = exact_seq(R)
        for i, r in enumerate(R["kidx"]):
            seqs[int(r)].append(W(int(Av[i]), int(Bv[i])))
    stats = Counter(); ex = {}
    for r in rows:
        s_ = seqs[r]
        Cp, L = berlekamp_massey(s_)
        N = []
        for n in range(numax + 1):
            acc = W(0)
            for j in range(min(n, len(Cp) - 1) + 1):
                acc = acc + Cp[j] * s_[n - j]
            N.append(acc)
        degN = max([n for n in range(numax + 1) if not N[n].iszero()], default=-1)
        reliable = 2 * L + 2 <= numax + 1
        den = tuple((j, round(float(Cp[j].norm()) ** 0.5, 4)) for j in range(len(Cp)) if not Cp[j].iszero())
        key = (L, degN, den if len(den) <= 5 else "long", reliable)
        stats[key] += 1
        if key not in ex:
            ex[key] = (r, [round(float(x.norm()) / q ** n, 4) for n, x in enumerate(N[:degN + 1])])
    print(f"q={q} kappa={kappa} lam={lam} P={marks[0]} numax={numax} rows={len(rows)}: (BM length L, deg numerator, (j,|C_j|), reliable) -> count")
    for k_, c in sorted(stats.items(), key=lambda t: -t[1]):
        print("  ", k_, c, " ex row", ex[k_][0], "|N(nu)|^2/q^nu:", ex[k_][1])
