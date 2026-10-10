# 11 — Energy-wall program: coordinator notes

See `BRIEF.md` for the objects (core sum 𝒩, single-mark energy E_P, Proposition U) and the plan.

## Step 1 — numerics at the weighted geometry (finished; `numerics/NOTES.md`)
Exact Z[ω] arithmetic, N₀ = 1.25·10⁴ … 2·10⁵ (≤ 56k columns, ≤ 28k rows per window), F = N/K ∈ [1, 181],
12 angular row sectors, random-phase controls with the same moduli. Library re-validated (β = γ₂ᾱ, Möbius
partner, identity (∗): 540/540; independent brute-force check to ≤ 1.4e-13).
- **Single-mark energies E_P** (13 primes, N(P) = 7 … 1009): ρ_rms ∝ rows^{−0.497±0.009}·F^{+0.001±0.009}·N(P)^{−0.019±0.005};
  theta/control rms 0.990 ± 0.012; inside the unproved window F < N(P): b = −0.500, c = −0.006 ± 0.010.
- **Coherent core sum 𝒩** with P78-type marks c_P = ᾱ(P)²γ₂(P)/√N(P) (single slots L ≈ 37…1459, two-slot design
  L ≈ 3600 with the theta cocycle, exact mask): rows^{−0.483±0.020}·F^{+0.006±0.020}·L^{+0.019±0.011}; theta/control
  1.003 ± 0.010; deterministic off-diagonal ≤ 0.2% (F ≤ 4), ≤ 1% (F ≤ 64). (H2) cross-mark correlations at the
  Marchenko–Pastur edge; coherent Rayleigh quotient 0.81–1.12; principal row |S₁|²/d₁ ≤ 5.9 (no Patterson bias seen).
- **Rank-one pairing T** (via the actual probe vector A_{m,σ}(Y) and a cross-scale estimator of its ρ-part, with a
  planted-signal positive control): |T(ρ)|² = 0.99 ± 0.03 and |T(ρζ)|² = 0.93 ± 0.03 times random; |T|/Q^{11/12} ≈ K^{−0.49}.
- New exact symmetry: |S_k̄| = |S_k| for β and conjugation-symmetric designs (K/2 independent members).
- **Extrapolated margin** at (M, ℓ, w) = (0.851, 0.1915, 0.043): 𝒩 0.407 ± 0.018 (0.47 ± 0.04 for L/F ≥ 10, L² ≤ N),
  single marks 0.427 ± 0.008, controls 0.415–0.423; over M ∈ [0.82, 0.86], ℓ ∈ [0.18, 0.20]: 0.39–0.43 (Z-exponent).
- Excluded at these sizes: F^f growth with f ≳ 0.05 (unless constant ≲ 10⁻³), dependence on L or the level P²
  (|d| ≤ 0.04), excess over random phases, deterministic secondary terms above 0.2–1%, an (NL)^{2/3}-type term on
  non-principal rows, rank-one alignment above 1.05× random. Not excluded: log factors, late-onset secondary terms,
  failure of the random-matrix law before K = Z^{0.85}, unmodelled P78 details (η ≠ 1, punctures, powerful rows, …).
- **Conclusion:** W(w) — including the single-mark case and the rank-one wall — is true numerically with margin
  ≈ Z^{0.41} at the weighted geometry. This is evidence of truth only, not of provability.
