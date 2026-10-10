"""Move calculus of the lem:plain induction (exponent level), solved by DP (vectorized).

State at a transform level: represented plain lengths (c1,c2) (inherited verbatim by children),
width M, live slot length z, allowance V (target Z^{M+V}).
Invariants S = c1+c2-M, D = c1-c2.  Moves (P78 sec:plain):
  reflect plain i at the current level, c_i -> M - c_i:
        (S,D) -> (-D,-S) [plain 1], (D,S) [plain 2], (-S,-D) [both]          (free)
  remove slot length d pointwise: z -> z-d, V -> V - kappa d                  (old-eq:3.5, 3.16)
  transform (ball g = h, amplifier ell = h/a, a = 6 for sixth powers):
        diagonal  S + z + g - ell <= V                                        (old-eq:2.12)
        child     (S+g, D, z) at width M-g with allowance V+ell               (old-eq:2.13, 2.14)
  terminal: z = 0 (zero-slot lemma, unrestricted lengths); at M = 0 remove the rest.
need[S,D,z] = minimal allowance V at entry (<= 0 means: provable with target Z^M).
The diagonal does not consume V (diagonal and children are separate terms of one sum)."""
import numpy as np


def solve(m, z0, kappa=0.75, a=6.0, h=0.01, M0=1.0, D0=0.0):
    N = int(round(M0 / h)); Zn = int(round(z0 / h)); ell = h / a
    R = N + 2; idx = np.arange(-R, R + 1); n = len(idx)
    need = np.zeros((n, n, Zn + 1))
    for zi in range(Zn + 1):
        need[:, :, zi] = kappa * zi * h          # width 0: remove everything
    sgrid = idx[:, None] * h
    for Mi in range(1, N + 1):
        new = np.empty_like(need)
        new[:, :, 0] = 0.0
        for zi in range(1, Zn + 1):
            child = need[:, :, zi]
            shift = np.full((n, n), 9.0); shift[:-1, :] = child[1:, :]      # value at (s+1, d)
            T = np.maximum(sgrid + zi * h + h - ell, shift - ell)           # transform from rep (s,d)
            # T indexed [s,d]; reachable reps from state (s,d): (s,d),(-d,-s),(d,s),(-s,-d)
            Tf = np.minimum.reduce([T, T[::-1, ::-1].T, T.T, T[::-1, ::-1]])
            new[:, :, zi] = np.minimum(Tf, kappa * h + new[:, :, zi - 1])
        need = new
    s0 = int(round((2 * m - M0) / h)); d0 = int(round(D0 / h))
    return need[s0 + R, d0 + R, Zn]


if __name__ == "__main__":
    for m in [0.3, 0.4]:
        for kappa in [0.75, 1.0]:
            out = []
            for z0 in [0.03, 0.05, 0.08, 0.12]:
                V = solve(m, z0, kappa=kappa)
                pred = kappa * z0 - (1 - 2 * m) / 6
                out.append(f"z={z0:.2f}: V={V:+.4f} (paper {pred:+.4f})")
            print(f"m={m} k={kappa}: " + "  ".join(out), flush=True)
    print("amplifier ratio a (paper a=6), m=0.4, kappa=0.75: threshold z where V crosses 0")
    for a in [6.0, 4.0, 3.0, 2.0]:
        zs = [0.02 * i for i in range(1, 11)]
        thr = None
        for z0 in zs:
            if solve(0.4, z0, a=a) > 1e-9:
                thr = z0; break
        print(f"  a={a}: first infeasible z on grid = {thr}  (formula (1-2m)/(a kappa) = {0.2/(a*0.75):.4f})", flush=True)
    print("asymmetric start D0 (not available for S_m^2 witness), m=0.4:")
    for D0 in [0.0, 0.1, 0.2]:
        out = []
        for z0 in [0.05, 0.08, 0.12, 0.16]:
            out.append(f"z={z0:.2f}: V={solve(0.4, z0, D0=D0):+.4f}")
        print(f"  D0={D0}: " + "  ".join(out), flush=True)
