"""
single-mark: exponent bookkeeping for the single-mark wall SM(w).

All quantities are log_Z of the corresponding sizes.  Inputs: (ell, M) or (ell, w) with
w = M + ell - 1 (wall excess).  Prints every threshold quoted in NOTES.md.

  K = Z^M rows, L = Z^ell mark norm, N = K^2 L / Z reflected column length,
  H = Z^2 / K Moebius-dual family size (twisted Prop U), Z primal length.
"""
import math, sys

CATALAN = 0.915965594177219


def table(ell, M):
    w = M + ell - 1
    K, L, N = M, ell, 2 * M + ell - 1
    H = 2 - M
    out = {}
    out["w = M+ell-1 (= log N/K)"] = w
    out["K"], out["L"], out["N"] = K, L, N
    out["KL (AFE length, conductor (KL)^2)"] = K + L
    out["K sqrt(L) (AFE point, U_P picture, conductor K^2 L)"] = K + L / 2
    # current bound (K+N)N : excess over KN
    out["excess of (K+N)N over KN"] = max(0.0, N - K)
    # reflection + sextic sieve branch: lossless iff N >= K^{3/2} L
    out["lossless reflection regime starts at N = K^{3/2}L"] = 1.5 * K + L
    N2 = 2 * K + L - N            # main Ramanujan piece length K^2 L / N  (= Z here)
    refl = max(0.0, N2 - K, (K + N2) * 2 / 3 - K, ((K + L + N2) * 2 / 3 - L) - K)
    out["excess of reflection+sextic-sieve branch at this N"] = refl
    # maximal loss over the gap K < N < K^{3/2} L : at N = K^{6/5} L^{2/5}
    out["max excess over the gap (K^{1/5} L^{2/5})"] = K / 5 + 2 * L / 5
    # (H1) proved: sum_k |D(1/2)|^2 << K^{6/5} L^{2/5}  (AFE split X = K^{6/5} L^{2/5})
    out["proved mean square of D(1/2): K^{6/5} L^{2/5}"] = 6 * K / 5 + 2 * L / 5
    # Moebius dual (twisted Prop U)
    out["H = Z^2/K"] = H
    out["natural boundary of Prop U_P: H = ZL"] = 1 + L
    out["Prop U_P error exponent 2 + max(0, 1+L-H)"] = 2 + max(0.0, 1 + L - H)
    out["needed error exponent"] = 2.0
    out["Diag exponent H+1"] = H + 1
    out["sextic-sieve theta term (HZ)^{2/3} Z"] = (H + 1) * 2 / 3 + 1
    # pointwise: subconvexity threshold with effective conductor C = K^2 L
    C = 2 * K + L
    out["individual subconvexity needed: delta >"] = 0.25 - w / (2 * C)
    # Poisson in P (AVG-SM): cubic GL(1) family of size H' = N^2/L at length N
    Hp = 2 * N - L
    out["H' = N^2/L (cubic GL(1) family)"] = Hp
    out["needed power saving in H'-units, delta >="] = 0.5 - L / (2 * Hp)
    out["conductor of twist xi_k = K, as power of H'"] = K / Hp
    # Phragmen-Lindelof interpolation between the lossless regimes
    Delta = 0.5 * K + L
    out["PL strip width Delta = log(K^{1/2} L)"] = Delta
    c = 2 * (4 * CATALAN / math.pi ** 2)   # energy cost exp(c sigma Delta)
    out["PL energy cost per unit type sigma: c*Delta, c = 8G/pi^2"] = c * Delta
    # kernel decay: dyads that matter
    out["only dual dyads N' > N Z^{-2w/3} matter; lower end"] = N - 2 * w / 3
    return out


if __name__ == "__main__":
    pts = [(0.19, 0.85), (0.18, 0.86), (0.20, 0.84), (0.191, 0.8492)]
    if len(sys.argv) == 3:
        pts = [(float(sys.argv[1]), float(sys.argv[2]))]
    for ell, M in pts:
        print(f"=== ell = {ell}, M = {M} ===")
        for k, v in table(ell, M).items():
            print(f"  {k:62s} {v: .4f}")
