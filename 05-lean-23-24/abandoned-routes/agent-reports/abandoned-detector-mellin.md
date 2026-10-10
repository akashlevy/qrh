# Abandoned routes: 7/8 detector/energy chain and the mean-square → zero-free step

Scope: dead code (reach.tsv `reachable=0`) in `Detector/`, `Energy/`, `Mellin/`, `PrimeRows/`, `Arithmetic/`, `RowCompletion/`, `PrimeCounting/` and the top-level modules, at upstream commit `fd4aeeb2e`. All paths below are relative to `lean/OAI/NumberTheory/DirichletL/`, and all line numbers refer to `fd4aeeb2e`.

## Overview

The upstream git history has only 3 commits, so I used the import graph as a proxy for chronology. "Depth" below means the longest import chain from the root of the library (2926 modules, maximum depth 404). I also diffed each dead module against its live counterpart.

Five patterns account for most of the meaningful dead code:

1. **A ladder of six successively weaker "Mellin step" interfaces** for 23/24. Each one is a complete reduction of `ShortDraft.DirichletTarget` to an analytic hypothesis that was never proved.
2. **Paper2's original prime-row extraction.** It was built on prime rows `p⁶` and is conditional on `FixedTestPrimeMS`, which is never constructed. The integer rows `b⁶` (with a totient weight) replaced it. These give exactly the same exponent.
3. **Quantifier and order-of-choices refactors** in the 7/8 final assembly.
4. **A systematic `*Power` → `*Bounded` rewrite** of the energy-induction stage lemmas.
5. **Over-uniform contour estimates** in the detector.

Caveat: some code is dead only because the pruning roots are the new `QRH/*` spine. This includes upstream's own working, unconditional 23/24 assembly (cluster 3). It is not an abandoned route. The rest of the dead Detector/PrimeRows code is mostly scattered leaf lemmas of 2–30 lines with no route-level meaning.

---

## 1. The `DirichletTarget` ladder: six Mellin-step interfaces (about 140 declarations, about 2.5k lines)

**Main declarations**

Each step reduces `ShortDraft.DirichletTarget` (`Arithmetic/EisensteinCoordinates.lean:253`, itself dead) to a different hypothesis:

| # | Hypothesis | Declarations | Depth |
|---|---|---|---|
| 1 | Exponentially smoothed rational Möbius sum `∑ χ(n)μ(n)e^{-n/D} ≪ D^σ` | `ShortDraftMellin.dirichlet_exp_zero_free_of_power_bound` (`Mellin/ExponentialSmoothing.lean:614`); `ExponentialMobiusPowerBound` / `dirichletTarget_of_exponentialMobiusPowerBound` (`Mellin/LogProfiles.lean:145,151`); `TargetLocalReduction.dirichletTarget_of_localBound` (`Mellin/LocalNonvanishing.lean:38`) | 2–3 |
| 2 | Same sum over ideals of `O = ℤ[ω]` (base change, `L_K(s,χ∘N) = L(s,χ)L(s,χχ₋₃)`) | `pair_zero_free_of_normFiber_exp_power_bound{,_only,_principal_left,_principal_right}` (`Arithmetic/NormFiberMellin.lean:608–714`); `dirichletTarget_of_ideal_base_change` (`Hecke/NormFiberCharacters.lean:450`); `dirichletTarget_of_normFiber_exp_bound` (`Hecke/LocalEulerIdentities.lean:500`) | 22–24 |
| 3 | Sharp-cutoff prefix sums `≪ n^σ` | `dirichletTarget_of_normFiber_prefix_power` (`Hecke/PrimeExtraction.lean:38`) | 25 |
| 4 | Compactly supported smooth weight on `[1,2]` | `CompactMellinBridge.BaseChangeCompactTest`, `pair_zero_free_of_compact_weight_bound{,_principal_left,_principal_right}`, `dirichletTarget_of_compact_normFiber_tests` (`Mellin/FiniteEulerRestoration.lean:593, 389, 507, 551, 623`) | 26 |
| 5 | Step 4, plus the Euler factors at a finite set `S` deleted | `BaseChangeOutsideCompactTest`, `baseChange_outsideCompact_zero_free`, `dirichletTarget_of_compact_outside_tests` (`Mellin/CompactScale.lean:779, 716, 793`) | 28 |
| 6 | Prime-row mean square | See cluster 2 | 29 |

Generic helpers: `ShortDraft.zero_free_of_analytic_identity` / `..._regularized_...` (`EisensteinCoordinates.lean:354, 377`) and `ShortDraftMellin.dirichlet_zero_free_of_mellin_bound` (`ExponentialSmoothing.lean:407`).

**What it attempted.** paper2's Mellin step (Prop 3.3, `prop:mellin`), equivalently paper78's Prop `lem:continuation-criterion`: a power saving in a smoothed Möbius sum continues `1/L` and excludes zeros.

**How far it got.** Every rung is complete and sorry-free but conditional. Each new rung removes obligations from the previous interface. For example, `BaseChangeCompactTest` still asks the caller for `term_integrable`, `fubini` and `locally_integrable`; `BaseChangeOutsideCompactTest` keeps only `power_bound` and `mellin_nonzero`. The surviving descendants are:
- `CompactMellinBridge.pair_zero_free_of_identity_gt_two` (`CompactScale.lean:613`), which takes the single hypothesis `χ ≠ 1 ∨ ψ ≠ 1` and replaces three separate principal-left/principal-right/nonprincipal lemmas;
- `ShortDraftMellin.inverse_mellin_analytic` (the only live declaration of the 39 in `ExponentialSmoothing`);
- `corrected_weighted_mellin_identity` (`RowCompletion/ThetaNonvanishing.lean:30`).

**Obstruction.**
- *Evidence:* the only analytic input ever proved is the mean square `outside_mean_square_of_theta_models` (`MeanSquare/ThetaDensity.lean:269`, depth 152). It bounds finite sums over ideals prime to `S`, with weights `W(N𝔫/Z)` where `W` is smooth and `supp W ⊆ [a,b]`, `a > 0`. Rungs 1–4 ask for sums with weights that are not compactly supported (`e^{-n/D}`), over the wrong ring, or with all Euler factors present. None of these can be fed by that theorem.
- *Conjecture:* the agents wrote the Mellin bridge top-down, before the mean square existed, guessing the interface. They re-guessed each time the input's true shape became clearer.

**Assessment.** Nothing new mathematically, but tidy. The exponential-smoothing criterion, "Möbius sum with `e^{-n/D}` weight `≪ D^σ` ⇒ `L(χ, s) ≠ 0` for `Re s > σ`", is a clean, general lemma that could be upstreamed. It would not improve any result.

---

## 2. Prime-row extraction (paper2's `p⁶` rows) and the role of `PrimeCounting/` (about 75 + 20 declarations, about 1.8k lines)

**Main declarations**
- `CompactScaleBridge.FixedTestPrimeMS` (`Hecke/FiniteEulerFactors.lean:609`). It bundles three things:
  - a set of generators of prime ideals with norms in `[D^{11/60}/2, D^{11/60}]`;
  - a count `#selected ≥ D^{11/60-ε}`;
  - a mean square of the row sums over all rows of norm `≤ D^{11/10}`, bounded by `D^{21/10+ε}`.
- `FixedTestPrimeMS.bound` (`:635`) turns these into a bound `384·D^{23/24+ε}` on the plain Möbius sum.
- The arithmetic core is `baseChange_prime_extraction_from_MS` (`PrimeRows/NormFiberExtraction.lean:327`), which calls `finite_actual_prime_extraction_dyadic` (`GaussSum/QuadraticNormalization.lean:417`). Related declarations: `outsideCoeff_prime_extraction_{from_MS,with_constants}` (`NormFiberExtraction.lean:521, 578`).
- Mellin side:
  - `realCompactOutsideSum_bound_of_fixedTestPrimeMS` and `FixedTestPrimeMSConstants(.absorb)` (`Mellin/CompactScale.lean:56, 148, 174`);
  - `zero_free_of_fixedTestPrimeMSConstants(_concrete)` and `dirichletTarget_of_fixedTestPrimeMSConstants_concrete` (`Mellin/CompactWeights.lean:128, 361, 385`).
- An earlier variant used rational-prime rows: `RationalPrimeExtraction.actual_rational_prime_extraction` (`PrimeRows/RationalApproximation.lean:332`). It uses rows `r⁶` for rational primes `r ≤ D^{4/45}`, with the count from Mathlib's Chebyshev bounds (`rational_prime_three_interval_count`, `GaussSum/GaussianDualLimit.lean:428`, also dead). It yields `≪ D^{89/90+ε}`.

**What it attempted.** paper2's `eq:prime-extract`. Since `χ_𝔫(p⁶) = 1` unless `p | 𝔫`, averaging `|A_{p⁶}|²` over primes of norm about `Y = D^{(1+ϑ)/6}` isolates the Möbius sum.

**How far it got.** Complete conditional chain from `FixedTestPrimeMS` to `χ.LFunction ρ ≠ 0` for `Re ρ > 23/24`.
- *Evidence:* `FixedTestPrimeMS` and `FixedTestPrimeMSConstants` appear only as hypotheses. No term of either type is ever produced; I grepped the whole repository.
- The rational-prime 89/90 version is conditional on a mean square over rows up to `D^{16/15}`. That mean square was never stated anywhere.

**Obstruction.**
- *Evidence:*
  - The bridge sits at depth 26–29. The mean square it needs appears only at depth 152.
  - When the mean square landed, the very next module, `Mellin/SixthPowerAverage.lean` (depth 153), extracted with integer rows `b⁶`. It uses a smooth bump average, a totient density `φ(𝔫)/N𝔫`, and an Euler correction `E_φ` (`outside_totient_moebius_bound`, `:317`). This needs only lattice-point counting, so no prime ideal theorem (PIT) is required. MAP.md says the same: "integer rows b⁶ … no prime counting".
  - The rational-prime variant loses because `N(r) = r²`: only `Y^{1/2}` rows are available, which gives 89/90.
- *Conjecture:* the prime-ideal count was the blocker at the time, because no PIT for `ℚ(√-3)` existed when the extraction was written.

**Was `PrimeCounting/` an attempt to discharge it? No.**
- *Evidence:*
  - Its namespace is `SevenEighths.PNT`.
  - Its only importers are `Hecke/RayPrimes.lean`, `Hecke/PrimePowers.lean` and `RayPrimeNormalizer.lean`, all in the 7/8 chain.
  - Its depth is 226–227, far above the prime-row bridge at 26–29.
  - It is absent from the 23/24 closure (`audit/upstream/closure-report-after.md`).
- The dead part is a superseded warm-up:
  - `Progressions.residueClass_ratio_tendsto` (`:30`) proves `ψ(x;q,a) ~ x/φ(q)` via PNT+'s Wiener–Ikehara;
  - `PrimeMass.prime_residueClass_ratio_tendsto` (`:66`);
  - the full-ray-class `RayAsymptotic.classCoeff_ratio_tendsto` (`:184`), replaced by the live quotient version `quotientClassCoeff_ratio_tendsto` (`:286`).
- It nonetheless got far enough to discharge the count. The live `AnnularPrimeMass.primeAnnularSum_log_div_tendsto` (`:335`) and `weightedPrimeSumDeleted_eventually_lower_bound` (`:584`) give `≫ x/log x` prime ideals with norm in an annulus. Even the dead rational `Progressions.lean` suffices: primes `p ≡ 1 (mod 3)` split into ideals of norm `p`.

**External packages** (quick grep of `.lake/packages`):
- Mathlib has only Dirichlet's infinitude: `Nat.infinite_setOfPred_prime_and_eq_mod`, `NumberTheory/LSeries/PrimesInAP.lean:430`.
- PNT+ has the rational `WeakPNT_AP` (`PrimeNumberTheoremAnd/Erdos970/Wiener.lean:3543`).
- "Dedekind-PNT" (a PIT for number fields) exists only as blueprint prose (`Consequences.lean:2557`).
- So for `ℚ(√-3)` the in-repo `PrimeCounting` is the only PIT available.

**Assessment.** Easy to complete now, but not useful. Write `ϑ` for the mean-square exponent: rows reach norm `D^{1+ϑ}` and the square sum is `≪ D^{2+ϑ+ε}`. Then both prime and integer rows give `D^{1+ϑ/2}·Y^{-1/2}` with `Y = D^{(1+ϑ)/6}`, i.e. the exponent `11/12 + 5ϑ/12`. Both have `ϑ = 1/10` hard-coded, giving 23/24. Completing the prime route therefore improves nothing, and it would pull the Wiener–Ikehara/PIT machinery into the short proof. The only reason to complete it is fidelity to paper2's text: the plain `μ`-sum with no `E_φ`.

---

## 3. Upstream's own 23/24 assembly (dead only because of the new roots; about 17 declarations)

**Main declarations**
- `CanonicalRowCompletion.dirichletTarget_unconditional`, `ShortDraft.dirichlet_LFunction_ne_zero`, `ShortDraft.riemannZeta_ne_zero_of_re_gt` (`RowCompletion/ZeroFreeRegion.lean:197, 204, 210`)
- `dirichlet_open_strip_of_theta_models` (`RowCompletion/ThetaNonvanishing.lean:105`)
- `SixthPowerAverage.outside_totient_moebius_bound`, `weightedCompactMellin_analytic`, `ideal_totient_moebius_bound` (`Mellin/SixthPowerAverage.lean:317, 682, 701`)

**Status.** These are complete and unconditional, and they were upstream's working route. `QRH/Extraction.lean` and `QRH/ZeroFree.lean` re-prove them from the same lemmas (`average_tsum_cauchy`, `average_remainder_small_power`, `pair_zero_free_of_identity_gt_two`). Do not count these as abandoned.

Side effect of the pruning: the 23/24 statement for all Dirichlet L-functions is gone, which is harmless because it is subsumed by the retained 7/8 Dirichlet theorem.

---

## 4. 7/8 top-level assembly: wrong order of choices (about 25 declarations plus near-duplicate modules, about 700 lines)

### (a) Supremum assemblies

**Declarations.** `Supremum.seven_eighths_of_uniform_margin(_with_sentinel)` (`Supremum.lean:68, 79`) and `HeckeZeroSupremum.beta_le_seven_eighths_of_uniform_margin` (`Hecke/ZeroSupremum.lean:64`).

**What they need.** A margin bound for *every* zero of *every* character, which is a pure order-theoretic sup lemma.

**Replacement.**
- The live `HeckeCommonProbe.beta_le_seven_eighths` (`Hecke/CommonProbe.lean:27`) does three things:
  - it picks a single zero near `β` of a **primitive** character (`HeckePrimitiveSupremum.exists_primitive_zero_near_beta`);
  - it shows `β` equals the supremum over primitive characters (`beta_eq_primitive_supremum`, via `HeckeFiniteDeletion.exists_primitive_zero`);
  - it requires a `PrimitiveContract` for primitive characters only.
- This matches paper78 Prop `lem:continuation-criterion`, where `β_*` is defined over primitive characters (paper78 lines 375–390).
- Obstruction (*evidence*): the probe estimates need primitivity, since `PrimitiveContract` requires `IsPrimitiveOnIdeals`, so a bound over all characters was never available.
- The deleted `TraceConductor.lean` (24 declarations, `PrimitivePresentation`) looks like an earlier attempt to handle imprimitive characters directly (*conjecture*).

### (b) `RawMomentInput`

**Declarations.** `ProbeFinalAssembly.common_probe_of_raw_moments` (`Detector/FinalAssembly.lean:11`) and `RawMomentInput` (`Detector/FinalAssemblyMomentInput.lean:34`). `RawMomentInput` demands moment bounds for **all** `D : HighData`.

**Replacement.** The live chain uses `FineMomentInput` (∀ `D` with `D.ell j ≤ mesh(D.t)/200`) and `ChosenMomentInput` (∃ `D`) (`FinalAssemblyChosenData.lean:13, 21, 30`), plus `exists_high_data_fine`.

**Evidence.** This is exactly paper78 Prop `prop:parameter-order` (line 16196): the mesh from Lemma `lem:plain` is chosen *before* the slot count `K`, so the moments hold only for fine slot lengths. The "for all D" version is presumably unprovable.

### (c) Existential count constants

**Declarations.** `ProbeHighRowFamily.actual_high_probe_from_raw_moments` (`PrimeRows/HighFromMoments.lean:18`) and `actual_nonfloor_rows_saving` (`PrimeRows/NonfloorCollected.lean:16`), plus the `NonfloorIntegral` and `NonfloorNormalized` modules.

**Evidence.** These are line-for-line twins of the live `Detector/FinalAssemblyHigh.lean:19` and `FinalAssemblyCollected.lean:17`. A `diff` shows a single change: the row-count constants `∃ cB κB cH κH` were chosen inside the theorem, and the live versions take them as an input `(counts : CountParameters …)`. The moment hypothesis `SourceMomentBound` is indexed by `counts`, so constants chosen internally could never be matched against it.

### (d) Other

`Parameters.lean` (4 trivial numeric lemmas, deleted).

### Assessment

None of this is a mathematical dead end. These are complete theorems with the wrong quantifier shape, and there is nothing to salvage.

---

## 5. Energy induction: `*Power` → `*Bounded` (about 32 declarations, about 720 lines, 11 modules deleted)

**Main declarations** (all at depth 373–383, just below the final assembly at 403):
- `original_reference_power_robust` (`Energy/ReferenceRobustPower.lean:26`)
- `original_reference_power_live` (`ReferenceLivePower.lean:24`)
- `original_high_power` (`OriginalHighReflectionPower.lean:26`)
- `balanced_reference_from_zero_low` (`ZeroReferencePower.lean:20`)
- `balanced_reference_live` (`ZeroReferenceLive.lean:17`)
- `zero_from_same_width_growth` (`ZeroGrowthConclusion.lean:80`)
- `saturated_successor` (`CappedBandTransport.lean:21`)
- `actual_reference_reserve_previous` (`BackwardSourceRanges.lean:82`)

**What they attempted.** The stage bounds of the energy induction (paper78 §9348 "Finite propagation", §12477 "Fourth moments"). They have the form energy `≤ C·…·Z^{width+ε}` for the reference, zero and high branches.

**How far it got.** Complete. Each `*Power` module has a live `*Bounded` twin (`ReferenceRobustBounded`, `OriginalHighReflectionBounded`, `ZeroReferenceBounded`, `ReferenceLiveBounded`, `ZeroReferenceLiveBounded`). Even the theorem names coincide.

**Obstruction.** *Evidence* from a `diff` of the twins; there are exactly three changes:
1. A new upper bound `L ≤ Mcap+Bmask+rho/100` on the level parameter.
2. The positive/zero-low input is required at an arbitrary slack `∀ e ≥ 0`, not only at the internally chosen `d`. The conclusion becomes `Z^{width+e+ε}`.
3. The reserve lemma changes from `exists_reference_reserve` to `exists_bounded_separated_reference_reserve` (`StageReserveBounded`).

*Interpretation (conjecture, well supported):* `CappedWidthInduction` → `CertifiedExistence` iterates these bounds across stages. A fixed loss `d` and an unbounded `L` do not compose, so the induction hypothesis had to be made loss-parametric and capped.

**Assessment.** No value, but it shows a clear failure mode: proving a stage lemma before the induction that consumes it is fixed.

---

## 6. Detector: over-uniform contour estimates and orphaned leaves (about 10 large declarations, about 420 lines, plus about 300 lines of Gaussian seminorms)

**Main declarations**
- `ProbePrincipalContours.source_uniform_high_slices`, `residue_uniform_moments`, `residue_uniform_tails`, `source_uniform_joint_tails` (`Detector/PrincipalContours.lean:965, 1003, 1080, 844`)
- `ProbePrincipalTransport.uniform_high_{x,z}_joins` and `strip_integrable_of_uniform_slices` (`PrincipalTransport.lean:435, 487, 407`)
- The Schwartz-seminorm control of the Gaussian dyadic profile: `GaussianSeminorm.lean:96–150`, `GaussianDerivative.lean:59–117`, `GaussianSummation.lean:43`, `GaussianLow.lean:13`

**Evidence.** Compare `source_uniform_high_slices` (`:965`) with the live `continued_source_slices` (`:430`):
- the dead lemma gives one constant `K` uniform over all `η, S, J, T, b` and over the real part in a whole interval `σ ∈ [a,Bs]`;
- the live lemma fixes these and works on one line `Re = a`.

paper78 says explicitly that "only the positive power margins must be uniform in the target", so these uniform versions were unnecessary strengthenings.

**Deleted leaf modules:**
- `HighRowsFiniteBounds` (`rowClosedMarked_second_region`, `:79`)
- `InitialFubini` (`:12`)
- `LowGramHeight` (`:61`)
- `LowOriginalPools` (`:20`)
- `LowPrimeFamily` (`:70`)
- `GaussianLow`, `GaussianSummation`

Each is 16–77 lines: individual lemmas whose consumer was rewritten.

**Assessment.** Low value. These are complete but unneeded; at most they are a source of reusable integrability lemmas.

---

## 7. Minor: early arithmetic

- **`Arithmetic/EisensteinIdeals.lean`** (depth 2; 32 of 38 declarations dead, 573 lines): the products of local quadratic characters on `O` and the "does not factor through a smaller modulus" (exact-conductor) lemmas (`finiteQuadratic_detects_prime :240`, `correctedQuadraticProduct_not_factor_without_prime :834`). *Conjecture:* groundwork for the quadratic sieve, superseded when `QuadraticSieve/LogarithmicLoss.lean` (depth 104) proved the sieve by Heath-Brown iteration.
- **`ThetaRegularity.lean`** (24 of 26 dead): appears superseded by `PeriodicThetaRegularity`.
- **`RowCompletion/PolynomialHeight.lean`**: `polynomial_smooth_quadratic` and `polynomial_smooth_series_rowPhase` (`:504, :621`, about 210 lines), unused alternatives in the height bookkeeping.

---

## Worth exploring further

1. **A 15/16 short proof, by making the extraction parametric in `ϑ`.**
   - *Evidence:* the mean square is proved for every `σ > 1/20` (`MeanSquare/ThetaDensity.lean:269`, live variant `:138`). The extraction exponent is `11/12 + 5ϑ/12`, which tends to 15/16 as `ϑ → 1/20`.
   - *Next step:*
     - Replace the literals `11/60`, `11/10`, `21/10`, `137/120` and `23/24` in `QRH/Setup.lean`, `QRH/Extraction.lean` and `QRH/ZeroFree.lean` by `(1+ϑ)/6`, `1+ϑ`, `2+ϑ`, …, `11/12+5ϑ/12`.
     - Then take `ϑ ↓ 1/20` for each fixed `ρ` with `Re ρ > 15/16`.
   - MAP.md estimates 100–150 lines. The dead `CompactScaleBridge` code shows exactly where such constants get hard-coded.
   - This is the cheapest exponent improvement available for the short proof. It is still much weaker than the 7/8 already in the repository.
2. **Optionally, complete the prime-row route as a paper-faithful cross-check. Do not pursue it for exponent or length.**
   - Construct `FixedTestPrimeMSConstants` (`Mellin/CompactScale.lean:148`) from:
     - `QRH.meanSquare`, applied to `y ↦ W(2y)` (support `[1/2,1]`) at scale `2m`, so the rows reach `(2m)^{11/10}`;
     - a prime count, from `PrimeCounting/Progressions.residueClass_ratio_tendsto` with `q = 3, a = 1` (restore 2 small dead modules) or from the live `AnnularPrimeMass` (`:335`, `:584`).
   - Then apply `zero_free_of_fixedTestPrimeMSConstants_concrete` (`Mellin/CompactWeights.lean:361`).
   - Estimate: 150–300 lines. Result: 23/24 again, with paper2's plain `μ`-sum.
3. **Upstream the generic Mellin criteria** to Mathlib or PNT+:
   - `ShortDraftMellin.dirichlet_zero_free_of_mellin_bound` (`ExponentialSmoothing.lean:407`);
   - `inverse_mellin_analytic`;
   - `mellin_of_exponential_smoothing` (`:445`);
   - `CompactMellinBridge.pair_zero_free_of_identity_gt_two`.

   They are self-contained and depend only on Mathlib. This is cleanup, not research.

## Lessons / dead ends to avoid

- **Do not fix an interface before its input exists.** The prime-row/Mellin bridge (depth 26–29) was designed about 120 import levels before the mean square (depth 152). Its hypotheses never matched what was eventually proved. Six `DirichletTarget` rungs were written the same way.
- **Prime rows versus integer rows makes no difference to the exponent** (`11/12 + 5ϑ/12` either way). Integer rows avoid any prime ideal theorem. Rational-prime rows are strictly worse (89/90) because `N(r) = r²`.
- **Match the test-function class to the mean square:** smooth weights compactly supported in `(0,∞)`, over ideals prime to `S`. Exponential or sharp-cutoff weights cannot be fed by it.
- **Encode the paper's "order of choices"** (paper78 Prop `prop:parameter-order`) as Lean structures early:
  - `∀ D` became `∃ fine D`;
  - constants chosen inside a theorem became input structures (`CountParameters`);
  - a fixed stage loss `d` became `∀ e ≥ 0` with a capped `L`.

  Most dead 7/8 top-level and energy code is the "before" version of one of these changes.
- **Prove only the uniformity the paper needs.** Uniform-in-target contour bounds (cluster 6) and the "every zero" supremum assembly (cluster 4a) were unnecessary strengthenings. Restrict to primitive characters, as paper78 does.
- **Avoid case-split duplication.** The principal-left / principal-right / nonprincipal triples collapsed into one lemma with the hypothesis `χ ≠ 1 ∨ ψ ≠ 1`.
- **Reading `reachable=0`:** with the `QRH/*` spine as root, upstream's own unconditional 23/24 chain is dead but was never abandoned. Check for a live twin under a different namespace before calling something a failed route. In clusters 4 and 5, most dead modules have a near-identical live twin.
