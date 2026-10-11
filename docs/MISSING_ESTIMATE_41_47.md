# Below 41/47 with the cubic-theta method: the missing estimate

*Research note, 10 October 2026. Akash Levy, working with AI assistants (Claude Opus 5.5 sessions). Not peer reviewed.*

Labels:
- **[P]** proved, relative to the cited lemmas of OpenAI's 7/8 paper (P78);
- **[B]** proved as exponent bookkeeping;
- **[M]** an exact statement about our exponent model;
- **[C]** conditional;
- **[N]** numerical;
- **[H]** heuristic.

## Summary

- **Best bound with proved inputs:** this method gives at best θ = 0.8748866 [M], Collet's value. Our 10499/12000 ≈ 0.8749167 is the Lean-verified bound.
- **41/47:** the next natural target, 41/47 ≈ 0.8723404, would follow [M] from any one of about ten estimates (§2). Every one of them lies beyond the edge of a large sieve.
- **The cleanest of them is a lossless energy wall** (§3). It asks for a second-moment bound for quadratic twists of twisted cubic theta series, only about Z^0.04 past the large-sieve edge. It reduces to a statement about one mark and then to a single projection bound.
- **Evidence:** the estimate holds numerically, with a large margin, over Q(ω) and in a function-field model (§4).
- **Why it is hard:** known methods stop exactly at the edge, for a structural reason we can name (§5).
- **What we're asking:** whether any technique gives an estimate of this kind (§6).

## 1. Setting

θ means a zero-free half-plane: no finite-order Hecke L-function over F = Q(√−3) has a zero with Re s > θ. Base change then gives the same for ζ and every Dirichlet L-function.

| Bound | Source |
|---|---|
| 7/8 | OpenAI's paper (P78), with a probe built from the cubic theta function |
| 0.874957 | P78's lemmas re-optimised: exact LP floor 33216105/37963129; κ-feedback limit ≈ 0.874957019 |
| 10499/12000 ≈ 0.8749167 | The high-side numerator bounded by a weighted fourth moment. Lean-verified with standard axioms only; [tracker PR #6](https://github.com/GTimG/QuasiRiemannTracker/pull/6) |
| 0.8748866 | The same inputs, optimised exactly. Hailey Collet, [tracker PR #9](https://github.com/GTimG/QuasiRiemannTracker/pull/9) |

**The exponent model.** It is `10-below-41-47/agents/combo/combo_model.py`.
- It contains every proved input, and its count LP is solved exactly.
- It reproduces the 10499/12000 certificate margin (−1.29·10⁻⁴) and Collet's endpoint (to 10⁻⁹).
- With proved inputs only, its optimum is 0.8748866 [M]. No lemma proved in this work moves it.
- The gap to 41/47 is **0.0025462**.

**Where the argument binds.**
- **Geometry:** mark (slot) length ℓ, and row range M = l_x + l_y with M + ℓ ≤ 1 (the "wall").
- **Low side:** on the wall it is B = 11/12 − ℓ/4, so 41/47 needs ℓ ≥ 25/141 ≈ 0.1773.
- **High side:** it currently closes only at ℓ ≈ 0.167.
- **The deficit:** +0.011 in exponent, at zero bins δ ≈ 0.47 (zeros near Re s = (1+δ)/2) with the largest prime amplitude.

## 2. Single inputs that give θ < 41/47 [M]

Each row below, added to the proved inputs, gives θ < 41/47 in the model.
- **Model settings:** κ = 35/47, thresholds found by bisection (`combo/thresholds41.py`, `combo/singles.py`).
- **Caveat:** these are exponent-level statements. Turning one into a theorem means re-running the paper's argument with the hypothesis at the new geometry. We have not written that out.

| # | Input | Proved value | Needed for 41/47 |
|---|---|---|---|
| 1 | **Lossless energy wall W(w)** (§3) | w = 0 | w ≥ 0.0407. With w ≥ 0.0359 it suffices if the same saving also enters `lem:marked`. |
| 2 | Rank-one wall W_θ, plus W(w) (`wall4` §7; W_θ is conditional on 04's Corollary D) | — | w ≥ 0.014 |
| 3 | Restricted large sieve on the rows carrying zeros (see below) | — | The full statement gives −0.029. A saving U^0.0137 in Σ_C\|S_u\| is enough. |
| 4 | Sixth moment of the numerator, Σ_u\|S_u\|⁶ ≪ U^{1+η} | η = 0.34 | η ≤ 0.097 |
| 5 | Power saving at the edge of `lem:plain` | 0 | 0.060 |
| 6 | Saving beyond the diagonal in the amplified inverse moment | 0 | 0.059 |
| 7 | Cross-mark cancellation in `lem:marked` | 0 | 0.155 |
| 8 | Long-witness amplification slope α | 5/6 | ≤ 0.583 |
| 9 | Zero density for the row family, N(σ) ≪ U^{A(1−σ)} on σ ∈ [0.59, 0.80] | — | A ≤ 2.125, against the best known for ζ itself, 30/13 ≈ 2.31 |
| 10 | Gram exceptional factor g | 1 | ≤ 0.739. This is equivalent to a pairing problem that is circular in the signal. |

**Row 3 in full:** Σ_{u∈C}|S_u|² ≪ U^ε(#C + U^{l_y/d}).
- S_u is the plain sum of the numerator, of length U^{n₀} with n₀ ≈ 0.46.
- C is the set of rows that carry zeros; #C ≈ U^0.59 at the binding bin.

**Partial savings do not combine** (`combo` §5). The five power-saving levers (rows 1, 3, 5, 6, 7) add almost exactly: each lever's share of its own standalone requirement sums to about 1 along the frontier, to within 3%. So 41/47 would need a common saving of about 0.0099 in all five lemmas at once. A saving of 10⁻³ in any single lemma gains at most 9·10⁻⁵.

## 3. The energy wall

### 3.1 Statement

**Notation** (as in P78):
- Z[ω], primary elements, norms q_n;
- the sextic residue symbol χ_m(n) = (n/m)₆, and the quadratic character ψ_k = χ_k³;
- β = γ₂ᾱ, the cubic-theta coefficient on squarefree n (|β| ≤ 1; γ_j are normalised Gauss sums, α a Grössencharacter);
- the probe scale Z; rows k with q_k ≍ K = Z^M;
- marks P prime with q_P ≍ L = Z^ℓ, and coefficients c_P = a(P)q_P^{−1/2} with |a(P)| ≤ 1.

After Cauchy–Schwarz and reflection, the low side of the probe is controlled by

    𝒩 = Σ_{k sf, q_k ≍ K} | Σ_P c_P 1_{(P,k)=1} Σ_{n sf} β(n) ψ_k(n) χ_P(n)^{−2} W(q_n/N) |²,   N = K²L/Z.

**The known bound.** P78's `lem:hybrid-energy` gives 𝒩 ≪ (K + N)·N·Z^ε.
- N > K exactly when M + ℓ > 1, i.e. past the wall.
- The term N² is the column term of the quadratic large sieve (P78 eq. (quadratic-large-sieve)). For general coefficients it is attained by vectors in span{ψ_k}.

**W(w).** 𝒩 ≪ K·N·Z^ε whenever N ≤ K·Z^w, i.e. M + ℓ ≤ 1 + w.

**What the probe needs.** It needs w ≈ 0.04 at M ≈ 0.85 and ℓ ≈ 0.19. There, K = Z^0.85, N = Z^0.89, and the marks are long: L = Z^0.19 ≫ N/K = Z^0.04.

**Only one dyad matters** [B] (`10-below-41-47/agents/wall4` §2). Past the wall, only one branch of P78's `lem:reflected-energy` loses anything. So W(w) is equivalent to a lossless bound for one core sum with:
- all slots active;
- no dual cube part;
- rows in the top Z^{w/2} of their range.

### 3.2 One mark, and the exact failing step

**Splitting 𝒩.** It splits into a mark-diagonal part Σ_P |c_P|² E_P and cross-mark terms. The cross-mark terms need the twists for P ≠ P′ to decorrelate; they are of random-matrix size numerically. The diagonal part involves

    E_P = Σ_{k sf, q_k ~ K} | Σ_{n sf} β(n) ψ_k(n) χ_P(n)^{−2} W(q_n/N) |².

**SM(w)** (single mark). E_P ≪ K·N·Z^ε for N = K²L/Z ∈ (K, K·Z^w], uniformly in P with q_P = L.

W needs only the average of E_P over P [B]. Each twist θ ⊗ ψ_k ⊗ χ_P^{−2} has conductor of size K²L. So the hard range is N ∈ (K, KL): there both the sum and its reflection lie past the wall.

**Proposition U_P** [P, given P78's reflected and hybrid energy lemmas] (`11-energy-wall/single-mark` §3.2). For Z ≤ H ≤ Z^{3/2} and a = μḠ,

    F_P(H) := Σ_{h~H} | Σ_{n~Z} a(n) χ_P(n)² (h/n)₆ W(q_n/Z) |² = Diag_P(H) + O(Z^{2+ε}(1 + ZL/H)).

**What SM(w) amounts to.** SM(w) is equivalent [B] to removing the factor ZL/H at H = ZL·Z^{−w}.
- Every step of the proof is an exact identity: Poisson in h, μ-absorption, and the completed reflection.
- The one exception is the quadratic large sieve (Goldmakher–Louvel), applied to the reflected sum of length K²L/Z > K. That single step loses exactly N/K.
- A related fact [P]: reflection plus P78's sextic sieve is lossless only when N ≥ K^{3/2}L.

**Minimal form of the missing input.** The single theta vector is b_P = (β(n)χ_P(n)^{−2}W(q_n/N))_n, with ‖b_P‖² ≍ N. We need

    ‖Π_V b_P‖² ≪ K·Z^ε,   V = span{ψ_k : q_k ~ K},

i.e. b_P has no more mass on the K-dimensional space of quadratic characters than a random vector would. This is the theta analogue of a GL(2) statement: Lindelöf on average for K quadratic twists of one form of large level, at lengths between K and the approximate-functional-equation length. We know no result of that kind.

## 4. Evidence that W(w) is true [N]

**Over Q(ω)** (`11-energy-wall/numerics`).
- **Setup:** exact Z[ω] arithmetic at N₀ = 1.25·10⁴ … 2·10⁵, with N/K up to 181; random-phase controls with the same moduli.
- **Single marks:** the relative off-diagonal ρ = E_P/Diag − 1 scales as rows^{−0.497±0.009}·(N/K)^{+0.001±0.009}. That is random-matrix size, with no growth past the wall.
- **Coherent 𝒩** with P78-type marks: theta/control ratio 1.003 ± 0.010.
- **Extrapolated margin** at the weighted geometry: about Z^0.41, against a need of Z^ε.

**Over F_q[T]**, q ≡ 1 mod 6 (`11-energy-wall/function-field`).
- **An exact dictionary** [P]: for squarefree monic n, β(n) = c_{deg n}·χ̄₆(disc n), ψ_k(n) = χ₂(Res(k,n)) and χ_P(n)^{−2} = χ̄₃(Res(P,n)). So E_P and 𝒩 are exact integers. The twisted series has effective degree 2κ + λ − 2, which is conductor K²L, as over Q(ω).
- **The hard window** (both the sum and its reflection past the wall): 19 configurations at q = 7, 13, with marks of degree 4–6. There |ρ| ≈ 0.9·q^{−κ/2}, the random-matrix size, independent of the excess and of the mark degree. The largest E_P/Diag is 1.13.
- **All 588 quartic marks at q = 7:** every E_P is within 13% of the diagonal.
- **Proposition A** [C, on two stated hypotheses about the cohomology of the family]. As q → ∞ with degrees fixed,

      E_P = τ_F·q^{κ+ν} + O(q^{κ+ν−1/2}),

  where τ_F is the trace of Frobenius on the commutant of the geometric monodromy (numerically τ_F = 1). It is proved by applying Deligne's theorem to the family sheaf End F. The constants grow exponentially in the degree, so nothing follows at fixed q.

## 5. Why known methods stop at the wall

1. **The loss is a sharp step.** The column term N² is sharp for general coefficients. So a proof must use the theta structure of b_P, not just |b_P| ≤ 1. Poisson, reflection and μ-absorption are exact, and E_P is self-dual under them (`wall4` §3). Every combination of these moves with large sieves preserves the excess. This is the "move calculus" of `03-six-sevenths` R1.
2. **A functional-equation invariant** [P, given the functional equations of the complete series; `function-field` §5].
   - **The functional equations:** the generating series of the second moment in (row degree, column degree, column degree) has an affine C̃₂ group of functional equations, with invariant I = ã + b̃₁ + b̃₂.
   - **Where W sits:** W corresponds to I = ½. Absolute convergence gives I > 3/2.
   - **Why that blocks the standard route:** every Weyl-group functional equation preserves I, so functional equations plus Bochner's tube theorem never get below I = 3/2.
   - **Over Q(ω):** the same group appears there [H] (`09-new-directions`, wall-spectral). It is heuristic there because μ₁₂ ⊄ Q(ω), so the standard multiple-Dirichlet-series theory does not apply directly. In multiple-Dirichlet-series terms, the wall is the imaginary-root direction δ = (1,1,1): fixed excess, growing family.
   - **By contrast:** the first moment has a finite C₂ group and no wall.
3. **How function fields get past it.** The function-field results that do enter this region all add an input along δ:
   - Bucur–Diaconu;
   - Diaconu–Paşol–Popa;
   - Whitehead;
   - Bergström–Diaconu–Petersen–Westerland with Miller–Patzt–Petersen–Randal-Williams.

   That input is an extra functional equation, rationality, or homological stability. All of them rest on F_q[T] being a curve over a finite field, and we know no number-field substitute.
4. **Other routes, each failing at a stated step** (`single-mark` §4–§8):
   - **Perron and short sums:** the hybrid (k, t) sieve gains exactly what the high-side signal loses, an uncertainty principle worth T₀^{1/2}.
   - **Individual subconvexity:** it helps only for δ > 0.2394.
   - **Hecke operators:** the mark is invisible to them.
   - **Spectral methods:** θ ⊗ χ_P^{−2} is residual spectrum, and no Kloosterman sums arise.
   - **The cubic first-moment technique:** it beats the theta term of the cubic sieve, not the column term of the quadratic sieve.
   - **Poisson in P:** it turns the P-average into a cubic family of GL(1) Hecke L-functions. That family has a finite group of functional equations, but it needs a power saving ≥ 0.44 uniformly in twists of conductor H′^0.53, where H′ = N²/L.

## 6. Questions

1. **Second moments just past the edge.** Is there a number-field mechanism that gives Lindelöf on average for the second moment of K quadratic twists of a fixed automorphic object (here a twisted cubic theta series of level about L²) at sum lengths slightly beyond K? An excess of w ≈ 0.04 with L ≈ K^0.22 would suffice for 41/47.
2. **The Poisson-in-P dual.** Can it be taken to a power saving ≥ 0.44, uniformly in twists of conductor H′^0.53? This is the cubic GL(1) family with a finite group of functional equations from §5.
3. **A restricted large sieve.** C is a set of about U^0.59 sextic characters that carry zeros near Re s = (1+δ)/2. Is any saving over the large sieve possible on C for Dirichlet polynomials of length U^0.46? A saving of U^0.0137 would be enough.
4. **Stabilisation along δ.** Is there a number-field statement that these moments stabilise as all conductors grow at fixed excess, i.e. an analogue of the geometric input in §5.3?

## 7. Where the details are

| Topic | Files |
|---|---|
| Exponent model, thresholds, additivity | `10-below-41-47/agents/combo/` (`combo_model.py`, `thresholds41.py`, `singles.py`, `NOTES.md`) |
| Reduction of W to one core sum; rank-one wall | `10-below-41-47/agents/wall4/NOTES.md` |
| Other session-10 agents (numerator barrier, Lemma D, slot powers, Theorem G, cross-marks) | `10-below-41-47/NOTES.md` §8–§17 |
| Numerics over Q(ω) | `11-energy-wall/numerics/NOTES.md` |
| SM(w), Proposition U_P, failing step, routes (a)–(e) | `11-energy-wall/single-mark/NOTES.md` (`exponents.py` reproduces every exponent) |
| Function-field model, Proposition A, functional-equation invariant | `11-energy-wall/function-field/NOTES.md` (exact engine `ffenergy.c`) |
| Earlier numerics, Proposition U, (H1)/(H2), the conditional 19/22 | `09-new-directions/` (agents `edge-theta`, `wall-spectral`, `h2`, `hls`) |

**Sources.**
- P78 is OpenAI, *The Quasi-Riemann Hypothesis: a zero-free half-plane Re s > 7/8* (30 Sep 2026), in [openai/math](https://github.com/openai/math); its lemma labels are quoted above.
- The 10499/12000 manuscript is in `08-weighted-numerator-10499-12000/paper/`.
- The function-field literature is listed with its status in `11-energy-wall/function-field/NOTES.md` §5.
