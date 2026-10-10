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

## Step 2 — single-mark wall SM(w) (finished; `single-mark/NOTES.md`) — no proof, failing step pinned
- Statement: E_P = Σ_{k~K}|Σ_{n~N} β ψ_k χ_P^{−2}(n) W|² ≪ K N Z^ε for N = K²L/Z ∈ (K, K Z^w], L ≥ Z^w (weighted:
  K = Z^{0.85}, N = Z^{0.89}, KL = Z^{1.04}). Kernel decay in lem:reflected-energy narrows the need [B] to
  E_P(K, N′) ≪ K N′ (N/N′)^{1/2} on N′ ∈ [K Z^{w/3}, K Z^w], top rows, O < 2w, ℓ_b < w/3.
- Current bound [P]: (K+N)N Z^ε; the N² is the Goldmakher–Louvel column term (P78 l.2644 / the NG/C² term of
  lem:hybrid-energy l.7598), attained by flat vectors in span{ψ_k}. Lemma 2 [P]: reflection + P78's sextic sieve
  is lossless only for N ≥ K^{3/2}L; maximal loss in between K^{1/5}L^{2/5} = Z^{0.246}.
- **Prop U_P [P, same status as Prop U]:** for Z ≤ H ≤ Z^{3/2}, a = μḠ,
  F_P(H) = Σ_{h~H}|Σ_{n~Z} a(n)χ_P(n)²(h/n)₆|² = Diag + O(Z^{2+ε}(1 + ZL/H)).
  SM(w) ⟺ removing the factor ZL/H at H = ZL·Z^{−w} [B]. Exact failing step: Poisson in h, μ-absorption and the
  completed reflection are exact; the only loss is the quadratic sieve at reflected length K²L/Z > K.
- **Convention trap [P]:** the other cubic twist (χ_P⁴ in edge-theta's normalisation) is dual to the j = 2 theta,
  with its wall at H = ZL²; a "χ_P⁴ family" run tests a different family. (Step 1 measured E_P directly in primal
  form with χ_P^{−2}, so it tested the right object.)
- Routes: (a) Perron/short (H1): the needed short-range statement is SM itself; Lindelöf for one twist gives SM(∞)
  for that P, subconvexity helps only for δ > 0.2394; the T₀^{1/2} loss of height averaging sits in the high side
  (signal lives in one O(1) Mellin window; the hybrid sieve returns exactly T_eff^{1/2}); Phragmén–Lindelöf
  interpolation fails (Paley–Wiener, cost Z^{0.457σ}). (b) Hecke/Kazhdan–Patterson: the mark is invisible (Ψ³ = ψ_k);
  Hecke amplifiers lengthen the sum. (c) **Correction: W needs only the P-average AVG-SM**, which follows from SM for
  all but L·Z^{−w} primes. **New duality (Poisson in P, β·conj γ₂ = ᾱ):** AVG-SM ⟺ a Prop-U-type asymptotic for a
  cubic family of GL(1) Hecke L-functions at H′ = N²/L — finite functional-equation group, but it needs a power
  saving ≥ 0.44 uniformly in a twist of conductor H′^{0.53}: out of reach. (d) cubic first moment: beats theta terms
  ((QN)^{2/3}), not the column term; flat configurations satisfy all its constraints. (e) spectral: no Kloosterman
  sums arise (every Gauss sum is a character); θ⊗χ_P^{−2} is residual spectrum.
- **Minimal missing input:** ‖Π_V b_P‖² ≪ K Z^ε for b_P = βχ_P^{−2}W and V = span{ψ_k} — Lindelöf on average for K
  quadratic twists of one fixed theta of large level (theta analogue of an open GL(2) problem). True numerically (step 1).

## Step 3 — function-field analogue (finished; `function-field/NOTES.md`) — W true with random-matrix margin, no transferable mechanism
- **Dictionary [P]** (A = F_q[T], q ≡ 1 mod 6; brute-force checked at q = 7, 13, 19, 0 failures): for squarefree n,
  γ_j(n) = μ(n)(−τ_j/√q)^{deg n}(n′/n)₆^j (Hasse–Davenport + residue formula), hence **β = γ₂ᾱ = c_{deg n}·χ̄₆(disc n)**;
  ψ_k(n) = χ₂(Res(k,n)), χ_P(n)^{−2} = χ̄₃(Res(P,n)); identity (∗) holds with G, 𝔯 → degree constants; ᾱ is a degree
  character. E_P and 𝒩 are exact integers (exact Z[ω] engine `ffenergy.c`, validated against brute force); κ = Mg, λ = ℓg, ν = κ + wg.
- **Faithfulness [N/P]:** the twisted series has effective degree D = 2κ + λ − 2 (κ ≥ 3), i.e. conductor K²L — exactly P78's
  reflection N ↔ K²L/N. **Correction to BRIEF/wall-spectral:** θ⊗ψ_k⊗χ_P^{−2} has conductor K²L (self-dual point KL^{1/2}),
  not (q_k q_P)²; the hard window is N ∈ (K, KL) (both N and its reflection beyond the wall). Over F_q[T]: κ < ν < κ + λ − 2,
  nonempty iff λ ≥ 4; empty for λ = 0 (the FF form of Proposition U).
- **Numerics [N]:**
  - large q = 7…43, fixed degrees: ρ = E_P/Diag − 1 → 0 by a deterministic O(1/q) term (ρq = 0.917 … 0.999 at (κ,λ,ν) = (2,2,3);
    ρq ≈ −1.2 at (3,2,4)); ρ → −1 for ν > D; fourth moment → 2.
  - fixed q = 7, 13, 19 configurations **in the hard window** (λ = 4–6, κ = 4–7, ν = κ+1, κ+2): |ρ| ∝ rows^{−0.503}, median
    |ρ|√rows = 0.86 (random-phase control 1.0), i.e. |ρ| ≈ 0.9·q^{−κ/2}, independent of ν − κ and λ; m4 ≈ 2.03; max E_P/Diag 1.13;
    4-mark coherent 𝒩 within ±0.08.
  - all 588 quartic marks at q = 7: every ρ_P ∈ [−0.13, 0.12]; random-weight 𝒩 within ±0.03; structured weights (a ≡ 1,
    a = conj γ₂(P)) deviate by O(1) in both signs (−0.58 … +0.22) along Aff(F_q)-invariant directions; one Aff-isotypic
    operator-norm outlier, 3.53 vs control 2.45 (κ = 4) and 3.94 vs 3.38 (κ = 5, row subset) — bounded at the three sizes testable.
    Coordinator check: T ↦ aT + b splits the 588 marks into exactly 16 orbits (12 × 42 + 4 × 21), matching the 16 distinct ρ_P.
  - extrapolated to the weighted geometry: ρ ≈ Z^{−M/2} ≈ Z^{−0.43}, the same margin as step 1 over Q(ω).
  - side result: the Patterson bias (size q^{5ν/6}) survives on the principal row k = 1 only (ᾱ is a degree character in FF).
- **Large q [C on (T1) no classes of weight > ν, (T2) geometric irreducibility]:** Proposition A: E_P = τ_F q^{κ+ν} + O(q^{κ+ν−1/2}),
  τ_F = trace of Frobenius on the commutant of the geometric monodromy of F = Gr^W_ν R^ν f_!L (numerically τ_F = 1 for ν ≤ D).
  Mechanism: Deligne for the family sheaf End F + Schur ("Lindelöf on average" in the q-aspect). Constants e^{O(degree)}: nothing
  at fixed q; the full coherent mark family needs generic vanishing for F_{P₁}⊗F_{P₂}^∨ (open).
- **Fixed q [P/L]:** Σ_k x^{deg k}D_{k,P}(y₁)D̄_{k,P}(y₂) has the affine C̃₂ FE group with invariant I = ã + b̃₁ + b̃₂; all linear parts fix
  δ = (1,1,1) = the wall direction. W sits at I = ½; every Weyl-group FE preserves I, so FEs + Bochner stop at I ≥ 3/2 — over
  F_q[T] exactly as over Q(ω) (first moment: finite C₂, rational, no wall). The FF theorems that enter the Tits cone
  (Bucur–Diaconu; Diaconu–Paşol–Popa; Whitehead; Bergström–Diaconu–Petersen–Westerland + Miller–Patzt–Petersen–Randal-Williams)
  all use an input along δ (extra FE, rationality, homological stability); none covers this cubic-metaplectic C̃₂ series.
  [C]: a DPP-type continuation uniform in P would give SM(w) with a power saving for every w in the window.
- **Mechanism [H]:** control along δ from the geometry of configuration spaces (Deligne + big monodromy at large q; homological
  stability at fixed q ≥ q₀). Both need Spec F_q[T] to be a curve over a finite field; no number-field substitute is known
  (reciprocity, Rankin–Selberg, Kuznetsov and the move calculus act as finite FEs and preserve I). **Exponent payoff: 0.**

## Program status (steps 1–3 complete) — negative for θ
- W(w) and SM(w) are true numerically over Q(ω) (margin ≈ Z^{0.41}) and over F_q[T] (random-matrix law, margin q^{κ/2} ↔ Z^{0.43}).
- No proof. Single-mark failing step: the quadratic-sieve column term at reflected length K²L/Z > K; minimal missing input
  ‖Π_V b_P‖² ≪ K Z^ε (Lindelöf on average for K quadratic twists of one theta of large level). In multiple-Dirichlet-series
  terms: W needs I = ½, every finite functional equation preserves I, and the missing input is control along the imaginary
  root δ, which function fields get from geometry with no known number-field analogue.
- Payoff for θ: 0. θ < 41/47 remains open; the proved-input frontier is still 0.8748866.
