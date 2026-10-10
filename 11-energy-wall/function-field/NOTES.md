# function-field — the energy wall over F_q[T] (step 3 of the energy-wall program)

Agent `function-field`, 10 Oct 2026. Labels: **[P]** proved (here, or standard with a cited reference),
**[C]** conditional (hypotheses stated), **[N]** numerical (exact integer data), **[H]** heuristic,
**[L]** literature (abstracts / search snippets only — arXiv and most publisher sites are blocked; nothing by lixiang90 was used).

## 0. Verdict

1. **Exact dictionary [P].** Over A = F_q[T], q ≡ 1 (6), Hasse–Davenport plus the residue formula make every Gauss sum of
   squarefree modulus explicit: γ_j(n) = μ(n)(−τ_j/√q)^{deg n}·(n′/n)₆^j. Consequently
   **β(n) = γ₂(n)ᾱ(n) = c_{deg n}·χ̄₆(disc n)** (|c_d| = 1): the cubic-theta coefficient is a sextic character of the
   discriminant. ψ_k(n) = χ₂(Res(k,n)), χ_P(n)^{−2} = χ̄₃(Res(P,n)), identity (∗) holds with G, 𝔯 → degree constants,
   and ᾱ becomes a degree character (a vertical shift). E_P, 𝒩 are exact integers; K = Z^M becomes deg k = κ = Mg exactly.
2. **The model is faithful [N + P].** The twisted series D_{k,P}(u) = Σ_n β(n)ψ_kχ_P^{−2}(n)u^{deg n} has effective degree
   **D = 2κ + λ − 2** (κ ≥ 3), i.e. conductor of degree 2κ + λ: the FF form of P78's reflection N = K²L/Z.
   The second-moment series has the **affine C̃₂** functional-equation group with invariant I = ã + b̃₁ + b̃₂
   (wall-spectral's group, re-derived over F_q[T]); its linear parts all fix δ = (1,1,1) in exponent space, which is
   **the wall direction ν = κ** [P]. The genuine window (both ν and its reflection beyond the wall) is κ < ν < κ + λ − 2.
3. **Numerics: W and SM are true in F_q[T] with random-matrix margin [N].** In 19 window configurations (q = 7, 13;
   λ = 4–6; κ = 4–7; ν = κ+1, κ+2): E_P/Diag − 1 has median |ρ|·√rows = 0.86 (random-phase control 1.0), pooled fit
   |ρ| ∝ rows^{−0.503}, fourth moment m4 = 2.03 (range 1.87–2.10; control 2.0), cross-mark correlations ≈ rows^{−1/2}, coherent 𝒩/Diag − 1 within
   ±0.08 (median ±0.02) for 4 marks. With **all 588 quartic marks** (q = 7): every ρ_P ∈ [−0.13, +0.12], cross-mark rms
   = rows^{−1/2}, random-weight 𝒩 within ±0.03, structured weights (a ≡ 1, a = conj γ₂(P)) off by O(1) in both signs
   (−0.58 … +0.22), and ‖S‖²_op has one Aff-isotypic outlier: 3.53 vs control 2.45 (κ = 4), 3.94 vs 3.38 (κ = 5) — bounded,
   ratio not growing (§3.4).
   Translated to P78's geometry this is ρ ≈ (q^κ)^{−1/2} = Z^{−M/2} = Z^{−0.43}: the same margin edge-theta measured over Q(ω).
4. **Large q, fixed degrees [C on two topological inputs; N].** E_P/Diag → τ_F ∈ Z, the trace of Frobenius on the
   commutant of the geometric monodromy of the weight-ν part F of R^ν f_!L, with error O_{κ,λ,ν}(q^{−1/2}). So
   lossless SM at large q follows from geometric irreducibility of F (and fails if F splits into Frobenius-stable pieces). Numerically τ_F = 1 whenever ν ≤ D (ρ ≈ c/q with |c| ≤ 1.6, q = 7…43;
   no visible q^{−1/2} term), τ_F = 0 for ν > D, and m4 → 2 (big monodromy). **Mechanism:** Deligne's RH for the
   *family* sheaf End F on the row space (square-root cancellation in the family variable) + Schur — not equidistribution of
   individual Gauss-sum phases. It gives nothing at fixed q (constants e^{O(degree)}) and does not cover the full coherent
   mark family (≈ q^λ/λ marks; pairwise q^{−1/2} savings lose q^{λ−1/2}).
5. **Fixed q, growing degree [P/L/H].** The first moment is finite type C₂ (rational in FF); the second moment is affine
   C̃₂. W ⟺ continuation of the C̃₂ series to I = ½ (inside the Tits cone I > 0). Weyl-group FEs preserve I, so with
   Bochner they reach only I ≥ 3/2 — **in F_q[T] exactly as over Q(ω)**. The FF theorems that reach the Tits cone
   (Bucur–Diaconu, Diaconu–Paşol–Popa for quadratic D₄^(1); Whitehead's affine axioms; BDPW + Miller–Patzt–Petersen–
   Randal-Williams for all quadratic moments, q large) all use an input **along δ** (rationality/extra FE in the
   imaginary-root direction, or homological stability). None covers this cubic-metaplectic C̃₂ series. Conditional on a
   DPP-type continuation with polynomial growth uniform in P, SM(w) would follow for every w in the window, with power saving.
6. **Mechanism and number-field substitute [H].** What beats the large sieve beyond the wall in F_q[T] is control of the
   δ-direction (fixed excess ν − κ, growing family) by the geometry of the configuration-space tower: Deligne + big
   monodromy (large q) or uniform homological stability (fixed q ≥ q₀). Both rest on Spec F_q[T] being a curve over a finite
   field (q → ∞ deformation; comparison with configuration spaces over C). **No number-field substitute is known**;
   reciprocity formulas, Rankin–Selberg and the move calculus are finite FEs and preserve I. The FF analogue supports the
   *truth* of W (random-matrix law, as over Q(ω)), not a proof route. **Exponent payoff for P78: 0.**

## 1. Setting and exact dictionary

A = F_q[T], q prime, q ≡ 1 (mod 6) (q = 7, 13, 19, 31, 37, 43 used). Fix a generator g of F_q^× and
χ₆(c) := c^{(q−1)/6} read in μ₆(C) (χ₆(g) = e^{2πi/6}); χ₃ = χ₆², χ₂ = χ₆³. For monic m, (a/m)₆ := χ₆(Res(m,a)),
Res(m,a) = Π_{m(θ)=0} a(θ). Additive character e(f) = exp(2πi·res_∞ f / q) (coefficient of T^{−1}).

| P78 (Z[ω]) | F_q[T] | status |
|---|---|---|
| primary n, norm q_n | monic n, \|n\| = q^{deg n} | — |
| log-norm, Z, K = Z^M, L = Z^ℓ, N = K²L/Z, w = M+ℓ−1 | degrees: g, κ = Mg, λ = ℓg, ν = 2κ+λ−g, ν−κ = wg | exact |
| χ_m(a) = (a/m)₆; cubic χ_m², quadratic ψ_m = χ_m³ | χ₆(Res(m,a)); χ₃(Res(m,·)); χ₂(Res(m,·)) | [P] |
| reciprocity bicharacter 𝔯 (mod 4) | (a/b)₆ = (−1)^{((q−1)/6)deg a·deg b}(b/a)₆: constant on degree shells | [P], G3 checked |
| γ_j(n) = q_n^{−1/2}Σ χ_n^j(x)e(x/n) | **γ_j(n) = μ(n)(−τ_j/√q)^{deg n}(n′/n)₆^j**, τ_j = Σ_{x∈F_q}χ₆^j(x)ψ₀(x) | [P], G1 brute force |
| Grössencharacter α, γ₂³ = μα | α(n) = ε₀^{deg n}, ε₀ = −γ₂(T)³: **unramified (degree) character** | [P], G4 |
| β = γ₂ᾱ on squarefree n | **β(n) = c_{deg n}·χ̄₆(disc n)** | [P], G5 |
| G(n) = χ̄_n(4)Γ(n) (ray class mod 12) | a degree constant (Γ(n) = γ₃(n) ∝ (−1)^{deg n} by Stickelberger) | [P] |
| (∗) βψ_kχ_P^{−2} = μG𝔯·conj γ₁(n;k³P⁴) | same identity, G𝔯 ↦ const^{deg n} | [P] |
| E_P, 𝒩, Diag | Σ_{k sf, deg κ, (k,P)=1}\|Σ_{n sf, deg ν} χ̄₆(disc n)χ₂(Res(k,n))χ̄₃(Res(P,n))\|²; Diag = Σ_k #{n sf : (n,kP)=1} | exact integers |

*Proof of the Gauss-sum formula [P].* For monic irreducible P of degree d with root θ, a/P = Σ_θ a(θ)/(P′(θ)(T−θ)) gives
e(a/P) = ψ₀(tr_{F_{q^d}/F_q}(a(θ)/P′(θ))); substituting x = P′(θ)y gives g_j(P) = χ_P^j(P′)·G(χ₆^j∘N, ψ₀∘tr), and
Hasse–Davenport gives G = (−1)^{d−1}τ_j^d. Twisted multiplicativity g_j(mn) = (χ_m(n)χ_n(m))^j g_j(m)g_j(n) and
(n′/P_i) = (P_i′/P_i)Π_{j≠i}(P_j/P_i) assemble this into (n′/n)₆^j for squarefree n. Then (n′/n)₆ = χ₆(±disc n) and
Stickelberger μ(n) = (−1)^{deg n}χ₂(disc n) give β = c_d·χ₆^{3+2}(disc n)·(const) = c_d χ̄₆(disc n). ∎
`verify_gauss.py`, `verify_beta.py`: all identities checked by brute-force Gauss sums (q = 7, 13, 19; deg ≤ 4): 0 failures;
β(n)/χ̄₆(disc n) takes exactly one value per degree.

So the FF energy is a character sum over configuration spaces:
S_{k,P}(ν) = Σ_{n ∈ Conf_ν(A¹ − Z(kP))(F_q)} χ̄₆(disc n)·χ₂(Res(k,n))·χ̄₃(Res(P,n)),
the trace function of a tame rank-one local system L on Conf_ν(A¹ − Z(kP)) with half-twist monodromy ζ₆^{−1}, monodromy −1
around the roots of k and ζ₃^{−1} around the roots of P. The cubic "interaction" χ₃(disc) between points is the FF form of
twisted multiplicativity (no Euler product); the Möbius dual side (∗) is multiplicative: M_h ∝ [u^ν] 1/L(u, χ_h^{(6)}).

**Exact geometries.** Any (M, ℓ, w) with rational entries is realised exactly, e.g. (0.85, 0.20, 0.05) ↔ (g,κ,λ,ν) = (20,17,4,18),
but exhaustive computation costs q^{κ+ν}. Computable proxies: δ = ν − κ ∈ {1, 2}, λ ∈ {4,5,6}, κ ∈ {4,…,7}, matching
ℓ/w = λ/δ ≈ 4.75 but not M/w ≈ 21. The laws below (in q^κ and q) are what transfers.

## 2. Structure of the twisted series (q = 7 exact profiles)

* **Symmetries [P].** Aff(F_q) (T ↦ aT + b) preserves monicity after rescaling, squarefreeness, disc and Res up to
  degree-constants, so |S_{k,P}| is invariant under simultaneous action on (k, P). Consequences seen in the data: for κ = 1,
  λ = 0, S_k(ν) = 0 whenever the scaling character χ₆^{4ν−ν²} is nontrivial (ν = 1, 2, 3; the stabiliser G_m of k = T acts
  by it); all λ = 1 marks and all λ = 2 marks are equivalent
  (one Aff-orbit each), so their ρ_P coincide and their cross-correlations are inflated (FF artefact; use λ ≥ 3);
  selection-rule zeros S_k(ν) = 0 for rows with nontrivial stabiliser (e.g. ν ≡ 2 mod 3).
* **Effective degree [N].** mean_k|S_k(ν)|²/q^ν ≈ (squarefree density ≈ 0.75) for ν ≤ D and collapses beyond
  (`dk_profile_q7.out`, `profile_rows_q7k5.out`):

  | κ \ λ | 0 | 1 | 2 | 3 |
  |---|---|---|---|---|
  | 2 | 2 | 2 | 3 | – |
  | 3 | 4 | 5 | 6 | 7 |
  | 5 | 8 | 9 | ≥ 9 | – |

  So **D = 2κ + λ − 2** for κ ≥ 3 (κ = 2: D = λ + 1 for λ ≥ 1; κ = 4, λ = 0 is borderline between 5 and 6). Conductor degree 2κ + λ ↔ K²L, consistent with P78's
  reflection N = K²L/Z (dual of ν is D − ν). Note: this says the self-dual point is KL^{1/2} and KL is the dual of the wall;
  wall-spectral's "conductor (q_kq_P)², AFE length KL" should read "conductor K²L, window (K, KL)".
  Beyond D the squarefree series has only a small tail (|S| ~ q^{ν/6}): the squarefree restriction makes the FE approximate
  (`fe_test_q7.out`: exact row-by-row FE only in the smallest cases, κ = 2 and 42/294 rows at κ = 3).
* **Window [P given D].** Large sieve: lossless for ν ≤ κ. FE + large sieve on the dual: lossless for ν ≥ D − κ = κ + λ − 2.
  Genuine window: **κ < ν < κ + λ − 2**, nonempty iff λ ≥ 4. For λ = 0 (unmarked) it is empty: the FF analogue of
  Proposition U (no wall for the unmarked family).
* **Patterson bias on the principal row [N]** (`bias_check.out`): Σ_{deg n = ν}β(n)χ_P^{−2}(n) has size (0.75–1.2)·q^{5ν/6}
  (zero in one class of ν mod 3). In F_q[T], ᾱ is a degree character and cannot remove it (over Q(ω), β = γ₂ᾱ is bias-free).
  It lives on k = 1 only (κ = 0) — the FF form of wall-spectral's (NL)^{2/3} term — and is absent for κ ≥ 1.

## 3. Numerics (exact; `ffenergy.c`, validated against brute force in `test_ffenergy.out`)

ρ_P := E_P/Diag_P − 1; control = i.i.d. uniform sextic phases on the same support; m4 = mean|S|⁴/mean(|S|²)²;
cross = normalised Gram entries G_{PP′}/√(G_{PP}G_{P′P′}); 𝒩 normalised by the tuple diagonal Σ|c_P|²Diag_P;
op = ‖[S_{k,P}]‖²_op / mean Diag_P against the RMT edge (1 + √(#marks/#rows))²; 's' = uniform random row subset.

### 3.1 Large q, fixed degrees (`campaign_largeq.jsonl`, `analysis.txt`)

| (κ,λ,ν) | D | q = 7 | 13 | 19 | 31 | 37 | 43 | law |
|---|---|---|---|---|---|---|---|---|
| (2,2,3) | 3 | +0.131 | +0.076 | +0.052 | +0.032 | +0.027 | +0.023 | ρ·q = 0.917, 0.983, 0.993, 0.998, 0.998, 0.999 |
| (3,2,4) | 6 | −0.180 | −0.093 | −0.066s | −0.035s | −0.010s | −0.010s | ρ·q ≈ −1.2 (full rows) |
| (3,1,4) | 5 | −0.232 | −0.119 | −0.079s | −0.016s | −0.046s | −0.080s | ρ·q ≈ −1.5 (q ≤ 19) |
| (3,3,4) | 7 | −0.001 | +0.003 | +0.013s | −0.011s | −0.009s | +0.013s | ρ = 0 ± RMT |
| (2,1,3) | 2 | −0.797 | −0.902 | −0.937 | −0.964 | −0.970 | −0.975 | (1+ρ)·q = 1.42 … 1.09 |
| (2,2,4) | 3 | −0.878 | −0.925 | −0.948 | −0.968s | −0.973s | −0.977s | (1+ρ)·q ≈ 1.0 |

Fits (full-row runs): (2,2,3): |ρ| ∝ q^{−0.96}; ν > D: 1 + ρ ∝ q^{−1.15}, q^{−0.85}. Theta m4 → 2 (1.93–2.10 at q ≥ 31)
for ν ≤ D; cross-mark correlations equal the control. **Reading:** ρ → τ_F − 1 with τ_F = 1 (ν ≤ D) or 0 (ν > D); the
approach is a deterministic O(1/q) term (diagonal-density effects), no O(q^{−1/2}) term is visible.

### 3.2 Fixed q, growing degree (`campaign_q7.jsonl`, `campaign_q13.jsonl`)

Genuine window, theta ρ_P per mark (4 random irreducible marks of degree λ) and |ρ|·√rows (theta / control):

| q | (κ,λ,ν) | rows | ρ_P | \|ρ\|√rows th/ct | m4 | 𝒩_coh | op/edge |
|---|---|---|---|---|---|---|---|
| 7 | (4,4,5) | 2057 | −0.028, −0.016, −0.028, −0.013 | 1.02 / 1.18 | 1.99 | −0.034 | 0.95 |
| 7 | (5,4,6) | 14399 | +0.014, +0.019, +0.008, +0.008, +0.024, +0.014 (6 marks) | 1.88 / – | – | +0.025 | 1.00 |
| 7 | (6,4,7) | 526s | +0.018, +0.071, −0.069, −0.010 | 1.16 / 1.19 | 2.04 | +0.048 | 0.92 |
| 7 | (7,4,8) | 256s | +0.010, +0.016, +0.042, −0.041 | 0.49 / 1.22 | 2.06 | +0.026 | 0.92 |
| 7 | (4,5,5) | 2058 | −0.030, −0.003, +0.017, −0.022 | 0.93 / 1.18 | 1.99 | −0.006 | 0.96 |
| 7 | (4,5,6) | 2058 | −0.021, −0.004, +0.010, −0.028 | 0.82 / 1.14 | 2.04 | −0.024 | 0.95 |
| 7 | (5,5,6) | 3652s | +0.008, −0.019, +0.019, −0.005 | 0.85 / 1.49 | 2.01 | −0.006 | 0.97 |
| 7 | (5,5,7) | 513s | +0.012, +0.017, −0.004, −0.071 | 0.84 / 1.29 | 2.10 | −0.023 | 0.92 |
| 7 | (6,5,7) | 526s | +0.016, +0.006, +0.105, −0.000 | 1.22 / 1.20 | 2.08 | +0.034 | 0.98 |
| 7 | (6,5,8) | 256s | +0.089, −0.031, −0.108, +0.035 | 1.18 / 0.60 | 2.05 | +0.023 | 0.89 |
| 7 | (5,6,6) | 3652s | −0.003, −0.001, −0.006, −0.020 | 0.63 / 1.48 | 2.02 | −0.022 | 0.97 |
| 7 | (5,6,7) | 513s | −0.054, −0.001, −0.018, −0.023 | 0.69 / 1.30 | 2.05 | −0.005 | 0.89 |
| 7 | (6,6,7) | 526s | −0.023, +0.081, +0.008, −0.023 | 1.00 / 1.19 | 2.05 | +0.020 | 0.94 |
| 7 | (6,6,8) | 256s | +0.127, −0.036, +0.038, −0.077 | 1.26 / 0.60 | 1.87 | +0.021 | 0.95 |
| 13 | (4,4,5) | 1253s | −0.013, +0.016, +0.001, −0.013 | 0.42 / 0.62 | 1.97 | +0.027 | 0.96 |
| 13 | (4,5,5) | 1253s | −0.009, −0.047, +0.019, +0.002 | 0.91 / 0.63 | 1.98 | −0.016 | 0.97 |
| 13 | (4,5,6) | 271s | +0.008, +0.090, +0.063, +0.096 | 1.20 / 0.97 | 2.03 | +0.075 | 0.95 |
| 13 | (5,4,6) | 283s | +0.058, −0.034, +0.084, +0.062 | 1.04 / 0.99 | 2.06 | +0.068 | 0.94 |
| 13 | (5,5,6) | 283s | +0.025, +0.044, −0.004, −0.076 | 0.77 / 0.96 | 2.04 | +0.015 | 0.88 |

Pooled over all q = 7 configurations with ν ≤ D (window and FE-covered): **|ρ_P| ∝ rows^{−0.503}**, median |ρ|√rows = 0.86
(control 1.01); q = 13: median 0.86 (control 0.76). Since rows ≈ q^κ, the law is **|ρ| ≈ 0.9·q^{−κ/2}**: the off-diagonal of
the theta energy is random-matrix size in the family, with no dependence on the excess ν − κ or on λ. No configuration
in or beyond the window has E_P > 1.13·Diag (the max, +0.127, is on a 256-row subset with RMT size 0.0625).
FE-covered cases (λ ≤ 3) behave the same; the only O(1) deviations are the λ = 2 single-orbit artefacts and tiny κ.

### 3.3 Fits in q^{deg} and in q

* fixed degrees, q → ∞: ρ = (τ_F − 1) + c(κ,λ,ν)/q + (RMT fluctuation ≈ q^{−κ/2}), |c| ≤ 1.6, τ_F ∈ {0, 1}.
* fixed q, degrees growing: |ρ| ≈ 0.9·(#rows)^{−1/2} = 0.9·q^{−κ/2}, independent of δ = ν − κ ∈ {1, 2} and λ ∈ {4,5,6}.
* at P78's weighted geometry this extrapolates to ρ ≈ Z^{−M/2} ≈ Z^{−0.43} against a need of Z^{ε}: the same
  Z^{0.43} margin edge-theta found over Q(ω) [N/H].

### 3.4 The whole mark family (`allmarks.jsonl`)

q = 7, all 588 monic irreducible quartics (λ = 4):

| (κ,λ,ν) | rows | ρ_P over all marks | 𝒩_coh / 𝒩_rnd / 𝒩_Gauss | cross rms (1/√rows) | op: theta / control / RMT edge |
|---|---|---|---|---|---|
| (3,4,4) | 294 | mean −0.027, [−0.130, +0.077] | −0.46 / +0.02 / −0.49 | 0.063 (0.058) | 5.88 / 6.30 / 5.83 |
| (4,4,5) | 2057 | mean −0.034, [−0.086, +0.001] | −0.58 / −0.014 / −0.23 | 0.024 (0.022) | **3.53** / 2.45 / 2.35 |
| (5,4,6) | 847s | mean +0.010, [−0.091, +0.119] | +0.10 / −0.022 / +0.22 | 0.035 (0.034) | **3.94** / 3.38 / 3.36 |

SM holds for every mark (E_P ≤ 1.12·Diag). Random weights give 𝒩 = Diag(1 ± 0.03); the structured weights a ≡ 1 and
a(P) = conj γ₂(P) (the FF analogue of edge-theta's Kubota-aligned slots) deviate by O(1) in both signs (−0.58 … +0.22):
they project onto Aff(F_q)-invariant directions of the mark set (an FF artefact of using the whole Aff-stable set of
irreducibles), and the deviation is bounded, not growing from κ = 3 to 5.
The operator norm at (4,4,5) has one non-degenerate outlier (next five singular values 2.58, degenerate: the 6-dim irrep of
Aff(F₇)); its singular vectors are Aff-isotypic and correlate (0.37/0.33) with χ₆⁵(P(0)), χ₆⁵(k(0)) (`outlier.py`).
It is a bounded constant: theta/control op ratio 0.93 (κ = 3), 1.44 (κ = 4, full rows), 1.16 (κ = 5, 847-row subset), so
W uniformly in the marks holds with constant ≲ 4 and no growth is visible (only three sizes; a full κ = 5 run costs 10^{12}).

## 4. Theory (a): the large-q limit

Fix (κ, λ, ν) and P. Let U = {k ∈ Conf_κ : (k,P) = 1} (smooth, geometrically irreducible, dim κ),
X = {(k,n) : k ∈ U, n ∈ Conf_ν, Res(k,n)Res(P,n) ≠ 0}, f : X → U, and L the rank-one tame Kummer sheaf with trace function
χ̄₆(disc n)χ₂(Res(k,n))χ̄₃(Res(P,n)). Then S_k = Σ_i (−1)^i tr(Frob_k | R^if_!L) (Grothendieck–Lefschetz).

Standard facts [P]: R^if_!L is lisse on a dense open U₀ ⊂ U (generic base change; over C, f is a Fadell–Neuwirth fibration);
R^if_!L = 0 for i < ν (Artin vanishing on the affine fibres + duality); R^if_!L is mixed of weights ≤ i (Weil II).

**(T1) [hypothesis: no big classes]** all R^if_!L have weights ≤ ν. Equivalently H^j(Conf_ν(C − Z(kP)), L^∨) carries no
classes producing weight > ν in H_c^{2ν−j}; it fails for k = 1 (Patterson bias, §2). For κ ≥ 1 it is supported by the
large-q data: E_P/Diag stays bounded (→ 1 or 0) as q goes from 7 to 43, whereas a weight > ν class on a positive
proportion of rows would force E_P/Diag ≳ q. Put F := Gr^W_ν R^ν f_!L, lisse and pure of weight ν on U₀, geometrically semisimple:
F^geom = ⊕V_i^{m_i}.

**Proposition A [C on (T1)].** As q → ∞ (q ≡ 1 mod 6), for fixed (κ, λ, ν),
  E_P = τ_F·q^{κ+ν} + O_{κ,λ,ν}(q^{κ+ν−1/2}),  Diag = q^{κ+ν}(1 + O((κ+λ+1)/q)),
where τ_F = tr(Frob | End_{π₁^geom}(F)) ∈ Z, |τ_F| ≤ Σm_i². Hence ρ → τ_F − 1. τ_F = 1 if F is geometrically irreducible;
τ_F = Σm_i² ≥ 2 if F is reducible with Frobenius-stable isotypic parts; τ_F = 0 if F = 0. For two marks, Σ_k S_{k,P₁}conj S_{k,P₂} = q^{κ+ν}·tr(Frob | Hom_geom(F_{P₂},F_{P₁})) + O(q^{κ+ν−1/2}).

*Proof.* By (T1), S_k = (−1)^ν tr(Frob_k|F) + O(Bq^{(ν−1)/2}), B = total Betti number of the fibre. Pure weight gives
conj tr(F_k) = q^ν tr(F_k^∨), so Σ_{k∈U₀}|tr F_k|² = q^ν Σ_k tr(Frob_k|F⊗F^∨) = q^{κ+ν} tr(Frob|(F⊗F^∨)_{π₁^geom})
+ O(B′q^{κ+ν−1/2}) by Grothendieck–Lefschetz on U₀ and Deligne (H^i_c, i < 2κ, of weight ≤ i). Coinvariants equal the
commutant because F is geometrically semisimple. Rows in U∖U₀ and cross terms contribute O(q^{κ+ν−1/2}). The squarefree
density of n coprime to kP is 1 − O((κ+λ+1)/q). ∎

**(T2) [hypothesis: irreducibility]** F is geometrically irreducible for ν ≤ D. Under (T1)+(T2), SM holds as q → ∞ with
E_P/Diag = 1 + O(q^{−1/2}). Numerically [N] (§3.1): τ_F = 1 for ν ≤ D in all four tested triples, τ_F = 0 for ν > D,
theta m4 → 2 (Larsen: big monodromy), and the observed approach is O(1/q), suggesting H^{2κ−1}_c(U₀, End⁰F) = 0.

**What mechanism.** Square-root cancellation in the *family* variable: Deligne's bound for H^{<2κ}_c(U, End F), plus Schur's
lemma for the main term. It is "Lindelöf on average" in the large-q aspect (mean |S_k|² = #terms·(1 + o(1))). It is not
"equidistribution of the Gauss-sum phases": each S_k is a sum of q^ν unimodular terms whatever their phases; what is used
is the joint distribution of the S_k (Katz's monodromy). The pointwise Deligne bound |S_k| ≤ Bq^{ν/2} is lossless in q but
B ≈ binom(D, ν) is exponential in the degree. The proof does not see the wall at all: at large q every fixed ν works.

**Limits [P/B].** (i) The O(q^{−1/2}) constant is e^{O(κ+λ+ν)}: nothing at fixed q. (ii) The full coherent 𝒩 has ≈ q^λ/λ
marks; pairwise Hom-vanishing gives cross terms O(q^{κ+ν−1/2}) each, so 𝒩_off ≪ q^{λ−1/2}·𝒩_diag — useless for λ ≥ 1.
W for all coherent weights needs ‖Δ‖_op small for the Gram matrix G = q^{κ+ν}(I + Δ), e.g. via middle-degree concentration
of H^*_c(U, F_{P₁}⊗F_{P₂}^∨) (cross terms O(q^{κ/2+ν}), so ‖Δ‖_HS ≲ q^{λ−κ/2}: the random-model condition L ≤ K^{1/2}).
That is a generic-vanishing statement for the family sheaf [C, open].

## 5. Theory (b): fixed q, multiple Dirichlet series

**Functional equations [P, given the FEs of the complete series].** With |x| = q^{−a}, |y_i| = q^{−b_i},
ã = a − ½, b̃_i = b_i − ½, the generating function Z(x,y₁,y₂) = Σ_k x^{deg k}D_{k,P}(y₁)D̄_{k,P}(y₂) has
* σ₁ (theta FE in y₁, degree D = 2κ+λ−2): (ã, b̃₁, b̃₂) ↦ (ã + 2b̃₁, −b̃₁, b̃₂); σ₂ likewise in y₂;
* σ₃ (k-sum: Σ_kχ_k(n₁n₂)x^{deg k} = L(x,χ_{n₁n₂})·(sf factors), x ↦ 1/(qx)): (ã, b̃₁, b̃₂) ↦ (−ã, b̃₁ + ã, b̃₂ + ã).
Roots α₁ = (−1,1,0), α₂ = (−1,0,1), α₃ = (2,−1,−1); Cartan matrix [[2,0,−2],[0,2,−2],[−1,−1,2]], det 0: **affine C̃₂**,
with invariant **I = ã + b̃₁ + b̃₂** (α₁+α₂+α₃ = 0) and an infinite-order translation σ₃σ₁σ₃σ₂ (checked: its powers grow
linearly). Absolute convergence a, b_i > 1 gives I > 3/2. The first moment Σ_k x^{deg k}D_{k,P}(y) has Cartan matrix
[[2,−2],[−1,2]]: **finite type C₂** (σ₁σ₃ of order 4, group of order 8), so it is a rational function in F_q[T] (finite-type
MDS over F_q(T) are rational [L/H]): the first moment has an exact finite formula, no wall.

**Where W sits [P].** E_P(κ,ν) = [x^κ y₁^ν y₂^ν]Z; the bound E ≪ q^{κ+ν+ε} on rays ν = κ(1 + w′) corresponds to
(a, b₁, b₂) = (1, ½, ½), i.e. **I = ½**, strictly inside the Tits cone I > 0 — wall-spectral's number-field statement,
unchanged. In exponent space the linear parts act by κ ↦ ν₁+ν₂−κ, ν_i ↦ 2κ − ν_i (+λ−2), with common fixed vector
**δ = (1,1,1) = the wall direction** (fixed excess, growing family). Every Weyl-group FE preserves I, so Bochner hulls of the
convergence region stay in I ≥ 3/2: the "move calculus stops at the wall" in MDS language, and it is the same in F_q[T].

**What FF theory adds [L].**
* Bucur–Diaconu (Mosc. Math. J. 2010): quadratic 4th-moment series (affine D₄^(1)) continued to its natural boundary over
  F_q(T) — first infinite FE group treated.
* Diaconu–Paşol–Popa (AJM): the axiomatic quadratic D₄^(1) series over F_q(T) continues meromorphically into the interior of
  the complexified Tits cone, polar divisor and residues determined, exact formula for a weighted 4th moment with infinitely
  many terms, via "a new type of functional equation" (not in the Weyl group).
* Whitehead (thesis; Compositio 2016, type Ã): axiomatic affine series, "meromorphic continuation to the largest possible domain".
* Sawin (JNT 2024): coefficients of very general MDS (arbitrary characters, twisted multiplicativity) as trace functions of
  perverse sheaves, uniquely characterised by local-to-global axioms (FEs left open); Hase-Liu (2024): convergence region
  and FEs for these series.
* Bergström–Diaconu–Petersen–Westerland + Miller–Patzt–Petersen–Randal-Williams: stable homology of braid groups with
  symplectic coefficients + uniform twisted stability ⇒ CFKRS asymptotics for **all** moments of quadratic L-functions over
  F_q(t), q > q₀(r) (Fresán's Bourbaki survey quotes q > 224(r+1) as the survey's form).
* Diaconu (F_q(t)) and Diaconu–Whitehead (Q), third moment of L(½,χ_d), finite D₄: secondary terms q^{3D/4} resp. x^{3/4}.

None of these covers our series (quadratic twists of the cubic metaplectic theta, mixed quadratic/cubic/sextic characters,
affine C̃₂). **Conditional statement [C/H]:** if Z continues meromorphically into the Tits-cone interior with polynomial growth
uniform in P and polar divisor = Weyl orbit of {b₁ + b₂ = 1, a = 1}, then E_P = Diag(1 + o(1)) + Σ(secondary terms from the
polar orbit) + O(q^{(κ+ν)(1−η)}) on every ray of the window, i.e. SM(w) with power saving for all w < ℓ (κ < ν < κ+λ−2).
The numerics show no secondary main term above the RMT level in the window (§3.2).

## 6. (c) Mechanism, and a number-field substitute?

**The FF mechanism is control along δ.** Every finite FE (Weyl group of C̃₂; P78's reflection R; Poisson + μ-absorption)
fixes the wall direction δ and the invariant I. Both FF routes beyond the large sieve act along δ:
* large q (Prop. A): Deligne on the family sheaf + big monodromy;
* fixed q (conjectural here): homological stability of the tower U_κ ← (fibres Conf_{κ+δ}(C − (κ+λ) pts)) with coefficients
  End F_κ, stabilising by adding one row point and one column point (i.e. moving along δ); if H^i(U_κ; End F_κ) is
  stable for i ≤ cκ with stable part = trivial summand, then E_P = Diag(1 + O(B^κ q^{−cκ/2})), a power saving for
  q > B^{2/c} (BDPW template). Obstruction [H]: F_κ is a Lawrence-type representation at a sixth root of unity of rank
  ≈ binom(D, ν) (exponential in κ); its coefficient degree ≈ 2ν ≈ 2κ leaves no Randal-Williams–Wahl range, and it is not
  obviously pulled back from algebraic representations of a stabilising family of arithmetic groups, which MPPRW need.
  On the Möbius-dual side the coefficients are Sym^ν ⊗ Sym^ν∨ of the superelliptic (y⁶ = h) unitary monodromy — closer to
  the MPPRW setting; this is the most promising FF route to a proof of SM [H].
* DPP's "new type of functional equation" also acts along the imaginary root [H: interpretation].

**Number-field substitute: none known [H, honest].** Both routes use that Spec F_q[T] is a curve over a finite field:
(i) the q → ∞ deformation at fixed degree has no analogue for Z[ω]; (ii) stability comes from the topology of
configuration spaces over C (comparison theorem). Number-field tools that look similar act as finite FEs (Motohashi/
spectral reciprocity swaps families; Rankin–Selberg sees only the diagonal; Kuznetsov moves along W-orbits), preserving I,
so they cannot reach I = ½ uniformly in the level. A substitute would have to be a statement that moments of this family
"stabilise" as all conductors grow along δ — which for P78's family is W itself. The FF analogue therefore confirms the
*truth* side (random-matrix law with margin q^{κ/2} ↔ Z^{0.43}) and locates the missing input (control along δ), but gives
**no new provable input for P78: exponent payoff 0.**

## 7. Files

| file | content |
|---|---|
| `ff_lib.py` | F_q[T] arithmetic, resultants, residue symbols, Laurent/residue additive character, Ben-Or irreducibility |
| `verify_gauss.py`, `verify_beta.py` → `verify_gauss.out` | brute-force checks G1–G5 (Gauss-sum closed form, reciprocity, Stickelberger, β = c_d χ̄₆(disc)) |
| `ffenergy.c` (`ffenergy`) | exact Z[ω] computation of S_{k,P}(ν) for all/selected rows and many marks; FFRANDOM control, FFROWS subsets |
| `ffdriver.py`, `test_ffenergy.py` → `test_ffenergy.out` | driver/statistics; validation against brute force (0 error) |
| `scan_q7.py`, `dk_profile.py`, `profile_rows.py`, `fe_test.py`, `rational_test.py` | first scans, degree profiles, FE tests |
| `campaign.py` → `campaign_{largeq,q7,q13}.jsonl` | the campaigns of §3.1–3.2 |
| `allmarks.py`, `outlier.py` → `allmarks.jsonl`, `allmarks_q7*.out` | full mark family, operator norm, outlier direction |
| `bias_check.py` → `bias_check.out` | Patterson bias on the principal row |
| `analyze.py` → `analysis.txt` | tables and fits |
