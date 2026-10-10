# Abandoned branches: large sieve, descent and Gauss sums

**Scope.** `CubicSieve/`, `GaussSum/`, `QuadraticSieve/`, `MeanSquare/`, `Poisson/`, `Fourier/`, `ChineseRemainder/`, `Reciprocity/`. Paths are relative to `lean/OAI/NumberTheory/DirichletL/` at `fd4aeeb2e`, with line numbers from `audit/upstream/reach.tsv`.

## Overview

These directories hold about 1,200 dead declarations, roughly 20k lines.

- **Almost none of it is conditional.** Only one chain rests on a hypothesis that nothing ever constructs: `HasInitialCanonicalDensity`. Everything else is fully proved and simply unused.
- **The dominant pattern is not a different mathematical route.** Most dead code is an earlier, *less uniform* or *differently normalised* version of a theorem whose final form is live. The agents repeatedly had to move an existential constant (derivative order, polynomial degree, implied constant) in front of a parameter that the descent instantiates later. The parameters were the prime pool, the height twist θ, the window family, and the support radius.
- **The four questions have mostly negative answers.** There is no Goldmakher–Louvel attempt, no band-limited Φ, no direct analogue of Lemma remove-exclusions, and no descent with a general gap or cube cutoff.
- **Two dead pieces could matter for the next step.**
  - Two complete, abstract descent lemmas (`ShortDraftDescent`) and a σ-general initial invariant are natural scaffolding for a parameterised descent toward 11/12.
  - The paper-faithful prime-row extraction at ϑ = 1/10 is complete.

**How "dead" was checked.**
- I cross-checked against `audit/reach2324.tsv`, which records reachability from the *original* upstream 23/24 theorem. Only one of my dead declarations was on that path: `InitialMeanSquare.outside_mean_square_of_theta_models`, which the new QRH spine merely re-chains.
- The same holds outside my scope for `SixthPowerAverage.outside_totient_moebius_bound` and `dirichlet_open_strip_of_theta_models`.
- Everything else here was unused by both the 7/8 and the 23/24 proofs.
- "Top-of-chain" means unreferenced by any other declaration. I computed it with a token cross-reference over the original sources (`scratchpad/xref.py`).

---

## 1. Initial dual-mean-square density, version 1 (the only conditional chain)

**Declarations.**
- `InitialMeanSquare.HasInitialCanonicalDensity` (MeanSquare/CanonicalDensity.lean:619–633);
- `.retained_density`, `.retained_source` and `.nonzero_source_bound` (MeanSquare/RetainedDensity.lean:287, 467, 617);
- the target `outside_mean_square_of_canonical_density` (MeanSquare/CanonicalTailBounds.lean:81–148), about 500 lines in its dead closure.

**What it attempted.** This is paper2's Prop thm:ms obtained from the descent (Prop canonical): a bound for the initial dual energies feeds Poisson to give
∑_{z∈T} |rowSum(z)|² ≪ Z·Z^{1+σ}·Z^ε for σ > 1/20.

**How far it got.** The chain is complete *modulo* `hcan : HasInitialCanonicalDensity χ S hbad`, which no declaration ever produces. Its quantifier shape is

```
∀ δ>0, ∀ P, ∃ J, ∀ V₁ V₂ (smooth, compact support), ∃ C ... ≤ C*Z^δ*(X*lengthScale)^2
```

**Obstruction (evidence).**
- The live replacement `HasInitialCanonicalDensityAtRadius` (Eisenstein/SeedHeight.lean:588) differs in exactly two places:
  - it inserts `∀ A, 0≤A` *before* `∃ J`, together with `(∀s, V s≠0 → |s|≤A)`;
  - it adds `(∀z∈T, IsCoprime z.1 (∏P∈S,P))`.
- The live descent enlarges supports at every level:
  - `canonicalLog_actual_step` (RowCompletion/DensityTransfer.lean:149) produces children with `|s| ≤ 2*A+11`;
  - `outsideLogEnergy_finite_rank` (:256) recurses with `ih (2*A+11)`.
- So the derivative order J must be allowed to depend on the support radius. Version 1 demanded one J for all windows, which the recursion cannot supply.
- The live chain is `outside_mean_square_of_canonical_density_at_radius` (MeanSquare/ThetaDensity.lean:192), discharged by `hasInitialCanonicalDensityAtRadius_of_theta` (RowCompletion/InitialDensity.lean:91).

**Assessment.** This is a quantifier-order dead end, not a mathematical one. It is not worth reviving.

## 2. Non-uniform versions of Prop R and of the transfer

**Declarations: Prop R ladder.**
1. `CompletedGauss.completed_mean_square_of_exact_reflection` (GaussSum/CompletedDyadicRows.lean:183–246). Fixed W. MAP.md still cites it as Prop R.
2. `CompletedHeight.polynomial_mean_square_of_exact_reflection` (MeanSquare/PolynomialHeight.lean:293–360), with wrappers `completed_mean_square_norm_height` (:439), `completed_mean_square_twist` (:578) and `short_completed_average_twist` (:595).
3. `completed_mean_square_uniform_height` and `completed_mean_square_uniform_twist` (MeanSquare/UniformPolynomialHeight.lean:458, 497).
4. Live: `uniform_polynomial_mean_square_of_exact_reflection` (same file, :309).

**Declarations: transfer.**
- `FirstPassCubeLabels.first_passage_transfer` (CubicSieve/FirstPassage.lean:182–276). Also cited in MAP.md.
- The untwisted `full_uniform_{tupleSource, globalSecondSource, globalSecondRawSource, globalArithmeticBin}_transfer` (CubicSieve/UniformTransfer.lean:273–483).

**What they attempted.** paper2's Prop R (eq:R) and Lemma first-transfer / second-transfer, in progressively more uniform forms.

**How far they got.** All are complete and unconditional. Each is superseded by a strictly more uniform live sibling:
- *Prop R.* `HasPolynomialKernel` (RowCompletion/PolynomialHeight.lean:412) is `∃n ∃C ∀t`. `HasUniformPolynomialKernel` (CubicSieve/RankMargins.lean:294) is `∃n ∀w ∃C ∀t`: the twist degree must be uniform over the whole window family.
- *First passage.* `first_passage_transfer` chooses `∃ K` *after* the prime pool p, F, B, v, X₁, X₂. The live `first_passage_full_uniform_transfer` (CubicSieve/SquarefreeTransfer.lean:442) has an identical conclusion with `∃ K` *before* `∀ {ι} (p : ι → O)`.
- *Second passage.* The untwisted transfers give way to the `*_twisted_transfer` versions (UniformTransfer.lean:523–680). These add `frequencyTwist A θ` and pay `(1+‖θ₁‖)^(J+2)(1+‖θ₂‖)^(J+2)`. This is the "height twists instead of C^J norms" device in QRH/Milestones.lean.

**Obstruction (evidence).** The descent step fixes `d` before `J` (`∃d, ∀J, d≤2*J → ∃child …`), and expands test functions in height twists θ. Every lemma it calls therefore needs constants uniform in the pool, θ and the window. Each dead version fails exactly one of these uniformities.

**Assessment.** These are pure engineering dead ends. For a cleaner proof the lesson matters more than the code: state every lemma in its fully uniform form from the start.

## 3. Monolithic second-transfer variants (about 6.9k lines, the largest cluster)

**Declarations.** All live in `SecondPassArithmetic`, `SecondPassIntegration` and `FirstPassCubeLabels`. The largest top-of-chain theorems, with dead-closure sizes:

| Theorem | Location | Dead closure |
|---|---|---|
| `first_passage_squarefree_child_transfer` | CubicSieve/SquarefreeTransfer.lean:285–375 | about 2.2k lines |
| `firstCoreInputRow_binned_integrated_transfer` | QuadraticSieve/DualMiddleExpansion.lean:609–731 | about 1.7k lines |
| `firstCoreInputRow_complete_finite_transfer` | CubicSieve/SupportNorms.lean:453–555 | about 1.5k lines |
| `normalized_truncated_second_transfer` | QuadraticSieve/TruncatedSecondTransfer.lean:334–423 | — |
| `normalized_twisted_second_sector_transfer` | CubicSieve/FirstCoreCutoff.lean:148–277 | — |
| `normalized_second_sector_transfer` | CubicSieve/FrequencyCutoffs.lean:207–320 | — |
| `actual_second_sector_uniform_transfer` | CubicSieve/FrequencyBlocks.lean:614–702 | — |
| `firstCoreInputRow_squarefree_binned_integrated_transfer` | QuadraticSieve/LogarithmicSectors.lean:563–685 | — |
| `actual_raw_twisted_sector_transfer` | CubicSieve/TwistedProfiles.lean:139 | — |
| `globalPositiveBin_transfer` | QuadraticSieve/GlobalBinLoss.lean:376 | — |

**What they attempted.** paper2 §7 (Lemmas first-transfer, arithmetic-poisson, second-transfer) as *one* end-to-end inequality. Each goes from the first-passage kernel sum straight to a sum of child `densityChildEnergy` terms. The bespoke budget functionals include `firstCoreSquarefreeBudget` and the structure `SquarefreeBudgetData` (SquarefreeTransfer.lean:155, 260). The variants differ in how child targets are indexed (squarefree log targets, log bins, integrated vs finite, truncated vs complete).

**How far they got.** All are proved, and none is used.

**Obstruction (evidence plus conjecture).**
- The surviving architecture is *modular*:
  - the first passage outputs an abstract `firstOutputEnergy` with a uniform constant (`first_passage_full_uniform_transfer`);
  - the second passage acts on that energy (`full_uniform_*_twisted_transfer`);
  - `reopenedCanonicalRow_closed_step` (CubicSieve/ClosedRecursion.lean:414) closes the step.
- The monolithic versions fix the column to `multiplicativeCoreColumn Ψ m g` with Schwartz `g`. They emit children indexed by ad hoc targets that do not match the recursive `CanonicalStateCondition` (MeanSquare/UnitRows.lean:168).
- **Conjecture:** they were abandoned because their outputs could not be fed back into the recursion without re-binning.

**Assessment.** No mathematical content beyond the live route; this is bookkeeping churn. Avoid fused, bespoke-budget transfer statements.

## 4. Orientation and normalisation variants of the two-passage block

**Declarations.**
- `orientedCubeFirstBlock_input_transfer` (CubicSieve/RadialWindows.lean:261–388);
- `canonicalCubeFirstBlock_input_transfer` (CubicSieve/QuantitativePassage.lean:181–300);
- `globalCubeFirstBlock_input_transfer` (GaussSum/PowerfulParts.lean:446–563);
- `{global,canonical,oriented,rawCanonical}TwoPassage_quantitative` (QuantitativePassage.lean:81, :308; RadialWindows.lean:396; MeanSquare/HarmonicAggregation.lean:413);
- `globalTwoPassage_transfer` (PowerfulParts.lean:571);
- `globalTwoPassage_original_ranges` (CubicSieve/RowCaps.lean:145).

**What they attempted.** Cauchy–Schwarz over cube blocks into the product of two input-family energies.

**How far they got.** The three `*_input_transfer` proofs are near copies. They differ in two choices:
- the bit order (`leftBit/rightBit` vs `rightBit/leftBit`);
- where the 1/N(d) normalisation sits: outside the sum, or inside the kernel as `K/primeProductNorm`.

The live `padded_orientedCubeFirstBlock_input_transfer` (MeanSquare/HarmonicAggregation.lean:498) is the oriented version with a source scale `Ksrc` decoupled from the cap `K`. Its bound gains a factor `(Ksrc/K)`. The live `padded_*TwoPassage_quantitative` (CubicSieve/PaddedPassage.lean:84, 154) follow the same pattern.

**Obstruction (evidence).** Padding was needed so that row caps nest across levels. The rank-margin lemmas show the padding: `paddedCanonicalRatio`, and `canonical_retained_scale_closure` (live) vs `canonical_retained_unpadded_closure` (dead, CubicSieve/DensityBudgets.lean:666).

**Assessment.** No value; fix conventions early.

## 5. Quadratic large sieve: alternative recursions and unused corollaries (no Goldmakher–Louvel)

**Declarations.**
- *Earlier forms of the large-gcd recursion.*
  - `rayThirdNorm`/`thirdNorm` and `highKernelNorm_le_thirdNorm` (QuadraticSieve/GcdReduction.lean:352–428). The operator norm of each gcd-D, ray-c block is an `sSup` of Rayleigh quotients.
  - The `smooth*` analogue (QuadraticSieve/AnnularRecursion.lean:119–187; GcdReduction.lean:625–672).
  - The `annular*` analogue, ending in `sieveNorm_le_annular_recursion` (AnnularRecursion.lean:448–634).
- *Corollaries of the live `sieveNorm_sharp`.*
  - `finite_family_quadratic_large_sieve` (LogarithmicLoss.lean:479);
  - `unrestricted_squarefree_energy_bound` (UnrestrictedCharacters.lean:360), which splits rows by bad-prime sector at a cost of a factor 4;
  - `bounded_supported_quadratic_norm` (MeanSquare/TransformedColumns.lean:585–684), for non-squarefree columns;
  - `product_quadratic_norm` (QuadraticSieve/ProductColumns.lean:123), for product columns n·b;
  - `completed_quadratic_dyad` together with the `CompletedGauss.exceptionalBranch` local decomposition (CubicSieve/RowCaps.lean:219–436).
- *Sharp-cutoff main-term analysis.*
  - `QuadraticMainBoundary.main_boundary_localization` (CubicSieve/CubeLabels.lean:466);
  - `QuadraticMainOperatorBound.actual_quadratic_main_difference_bound` (CubicSieve/CrossFactorization.lean:710);
  - `QuadraticMainTermDifference.actual_main_terms_difference` (QuadraticSieve/MainTermDifference.lean:108).

**What they attempted.** paper2 Lemma quadratic (eq:Q), and the use made of it in Lemma squarefree-completed.

**How far they got.**
- The recursion variants stop at a *recursive inequality*: `sieveNorm M N ≤ sieveNorm K N + 16 ∑_j (∑_D annularThirdNorm … + highKernelNorm … )`. No lemma ever bounds `thirdNorm`, `smoothThirdNorm` or `annularThirdNorm`.
- The live chain bounds weighted gcd blocks directly (`highKernel_weighted_energy_ray_gcd_reduction`, GcdReduction.lean:589) and bootstraps through `HasSieveExponent.improve : α ↦ 2−1/α` (LogarithmicLoss.lean:461; PoissonComparison.lean:133).
- The corollaries are all complete.

**Obstruction.**
- *Evidence:* the sSup-defined block norms were never estimated.
- *Conjecture:* an exponent recursion on a single scalar predicate (`HasSieveExponent`) was easier to close than one carried on block norms.
- No Goldmakher, Louvel or GL name, citation or hypothesis appears anywhere in `DirichletL`. Every variant is Heath-Brown-style.

**Assessment.**
- `product_quadratic_norm` and `completed_quadratic_dyad` with `exceptionalBranch` are a direct, paper-shaped route to the last step of Lemma squarefree-completed. They handle n·b product columns and the non-quadratic local factors at p | (h, f) that MAP.md item 3 mentions.
- They are possibly shorter than the live exact-models route (see "Worth exploring" item 3).

## 6. Prime-row extraction at ϑ = 1/10 (paper2 eq:prime-extract)

**Declarations.**
- `PowerExtraction.from_mean_square_and_prime_rows` (GaussSum/SexticNormalization.lean:567–615), with the variants `_with_loss` (GaussSum/GaussianFiberMass.lean:89) and `_with_constant` (GaussSum/QuadraticNormalization.lean:51);
- `ActualEisensteinCubic.finite_actual_prime_extraction` (GaussSum/QuadraticTrace.lean:439–493), `…_dyadic_ball` (GaussSum/GaussianDualLimit.lean:279) and `…_dyadic_ball_with_loss` (GaussianFiberMass.lean:375–458);
- the rational-prime version `RationalPrimeExtraction.lower_count_and_extract` (Poisson/GaussianMoments.lean:278), built on `rational_prime_three_interval_count` (GaussianDualLimit.lean:428). That lemma is a full Chebyshev lower bound for primes in (Y/3, Y], proved from Mathlib.

Each of these closures is about 1k lines.

**What it attempted.** Extract the Möbius sum from rows u = p⁶. The numerics are visible in the statement:
- `hJ : D^(11/60−ε) ≤ J`, `hmean : J*A² ≤ 2D^(21/10+ε) + 2J·D^(49/30)`, and the conclusion `A ≤ 2D^(23/24+ε)`;
- that is, H = D^{11/10} and A² ≲ D^{1+ϑ}/D^{(1+ϑ)/6}, which gives 11/12 + 5ϑ/12.

**How far it got.** Complete, but stated for the older finite-prime-pool row model (`finiteSquarefreeRow P hgood (support n) u`). The Gaussian-prime count `hcount` stays a hypothesis; the rational-prime variant discharges it.

**Obstruction (evidence).**
- It was replaced by `SixthPowerAverage.outside_totient_moebius_bound`, which uses integer rows b⁶ with a bump average and the identity `pair_inverse_euler_identity`. That route needs no prime counting and no error term for n sharing a factor with r, but it changes the target to a totient-weighted sum.
- **Conjecture:** the switch came when the mean square was restated over `idealRowSum`/`outsideIdealsUpTo`; the prime version was never re-targeted.

**Assessment.** This is the only place ϑ enters numerically. It is a complete template for any ϑ: the exponents are linear in ϑ.

## 7. Descent scaffolding and parameterisation pieces (none were used)

**Declarations.**
- `ShortDraftDescent.binary_energy_descent` (ChineseRemainder/AdditiveCharacters.lean:344–395). If `rank` strictly drops, `energy² ≤ (cL·E(left))(cR·E(right))` and `cL·mass(left) ≤ mass`, `cR·mass(right) ≤ mass`, then `energy ≤ mass`.
- `ShortDraftDescent.integral_geometric_descent` (QuadraticSieve/MaskedPrincipalBounds.lean:297–336). This is the integral version: `energy x ≤ terminal x + ∫ w(t)·√E(left x t)·√E(right x t)` together with a budget inequality for `mass`. It matches exactly the shape of a transfer step with a continuous separation parameter t.
- `contraction_terminates_after_3001` (AdditiveCharacters.lean:397).
- `InitialMeanSquare.initial_invariant_power_general` (CubicSieve/PooledRows.lean:611). For row scale H = Z^{1+σ} the initial ratio is ≤ c·Z^{−σ}, for *every* σ: paper2's eq:initial-positive-gap. Its siblings hard-code σ = 1/10 and reserve Z^{1/20}: `initial_invariant_power` (:621) and `initial_invariant_reserved_margin` (:632).
- `canonicalRankMargin_padding` and `canonical_retained_unpadded_closure` (CubicSieve/DensityBudgets.lean:574, 666). These are σ-general, but only within the fixed depth 3000.

**How far they got.** Complete and unreferenced. No dead declaration generalises the cube cutoff 1/1000 (`fixedDepthRank`, CubicSieve/ArithmeticBins.lean:246), the depth 3000, or the reserve 1/40 (`canonicalPadding`, CubicSieve/RankMargins.lean:79). The live invariant is

```
K*‖m‖² ≤ (X*F)*Z^(-canonicalRankMargin (1/20) η r)
```

(`CanonicalStateCondition.invariant`, MeanSquare/UnitRows.lean:177).

**The four specific questions.**
- *Band-limited Φ:* none. No statement in these directories asserts compact support of a Fourier transform. The row weight `rowMajorant` (CubicSieve/RowMajorants.lean:182) has compact support in *physical* space.
- *Lemma remove-exclusions:* no analogue. The closest pieces are the inclusion–exclusion identity `CoprimeMobiusExtension.ideal_mobius_coprime_reindex` (QuadraticSieve/CrossSeparation.lean:613) and the bad-sector splitting in item 5. The live state instead carries a mask `m`, which grows down the tree (`fixedTripleMask m q = m·q₁q₂q₃`, CubicSieve/PrincipalNormalization.lean:273) and is charged against the gap in the invariant above.
- *Conjecture:* paper2's lemma keeps H and Σ = XF unchanged, so a Lean version would remove `‖m‖²` from the invariant and simplify `CanonicalStateCondition`. Whether this also frees exponent is unverified.
- *Parameterised descent* and *Goldmakher–Louvel* are answered negatively above (here and in item 5).

**Assessment.** This is the most promising dead material for 11/12: abstract and constant-free, not tied to 1/20 or 3000.

## 8. Infrastructure prototypes and duplicate analytic lemmas

**Declarations.**
- *Coupled log separation.* This is the analogue of Appendix B: F(R·e^{∑ a_j y_j}) is written as a t-integral of product phases, uniform in R. The dead variants are:
  - `LocalLogFourier.coupled_double_positive_log_separation` (Fourier/PositiveLogProfiles.lean:852–981), which separates two coupled factors at once;
  - `_uniform` (:745) and `_two_regimes` (:819);
  - `FourierBridge.coupled_double_log_separation(_schwartz)` (PositiveLogProfiles.lean:124; GaussSum/LocalProductAssembly.lean:54);
  - `AnalyticBridge.double_weighted_fourier_minkowski` (LocalProductAssembly.lean:130).
  - Live: the single-factor `coupled_positive_log_separation` (PositiveLogProfiles.lean:686) and its `_envelope` variant.
- *Rational-model local identities.* `PoissonFirstLocal.first_local_identity` and `first_prime_product_identity` over `ZMod n` with ℕ moduli (GaussSum/FiniteFourier.lean:227–366).
- *Rank-two Poisson.* `RankTwoPoisson.tensor_poisson` and `eisenstein_gaussian_poisson` (GaussSum/QuadraticNormalization.lean:529–861).
- *Gauss-sum CRT.* `GaussCRT.gauss_sum_crt_standard` (GaussSum/ChineseRemainder.lean:19–248).
- *Ramified phases.* `quadraticGamma(Neg)TraceLambda_eq_(neg_)I` (GaussSum/RamifiedTracePhase.lean:96, 253).
- *Diagonal term of Prop poisson-reduction.* `ConcretePrimeRowBridge.idealRowSum_smoothed_poisson_zero_split` (GaussSum/ActiveConductor.lean:483), which isolates the zero frequency Z.
- *Cusp coercivity.* Hardy-type cusp coercivity and Green identities: `cusp_chart_dirichlet_coercivity` (QuadraticSieve/CuspCoercivity.lean:420) and `actual_eisenstein_local_green` (QuadraticSieve/GreenIntegrals.lean:445). These belong to the theta side despite their location.

**Assessment.** All are complete. Most are superseded by O-native or more general live versions. The double separation is a reasonable tool if the second transfer is ever rewritten with two Fourier factors, but it is not needed now.

---

## Worth exploring further

1. **Rebuild the descent on `integral_geometric_descent` with κ, cutoff and depth as parameters.**
   - *Next step:* instantiate it with `rank = ⌈log K / (c·log Z)⌉`, `mass = C_r·Z^{ε_r}(XF)²` and `energy = outsideLogEnergy`.
   - Prove `hstep`/`hbudget` from `canonicalLog_actual_step` by threading κ through `canonicalRankMargin`, `fixedDepthRank_step` and `canonicalPadding_budget`. These are the only places 1/20, 1/40, 1/1000 and 3000 enter the step.
   - Seed it with `initial_invariant_power_general` and generalise the extraction (item 6) in ϑ.
   - This gives 15/16 cheaply if only the seed changes, and 11/12 once κ → 0 is allowed.
2. **Prototype Lemma remove-exclusions to drop the mask from `CanonicalStateCondition`.**
   - *Next step:* prove, for `outsideLogEnergy`, E_m ≤ τ(m)·∑_{d|m} E(K, X/N(d), F·N(d)), starting from `ideal_mobius_coprime_reindex`.
   - Then check whether `reopenedCanonicalRow_closed_step` still closes with `invariant : K ≤ XF·Z^{−margin}` (no `‖m‖²`).
   - If it does, the state and the 66 constant-bearing lemmas shrink. This is a conjecture with a cheap test.
3. **Check whether `product_quadratic_norm` with `completed_quadratic_dyad`/`exceptionalBranch` can replace the live exact-models machinery behind Prop R.**
   - *Next step:* compare the live closure of `uniform_polynomial_mean_square_of_exact_reflection` against these complete corollaries of `sieveNorm_sharp`.
   - Measure the line count before investing any further.

## Lessons and dead ends to avoid

- **State lemmas fully uniform from the start.** Constants should come before the prime pool, the twist θ, the window family and the support radius, and the polynomial degree should be uniform in the window. Most dead declarations here are non-uniform drafts of live theorems: items 1, 2 and 4.
- **Do not fuse the first and second transfers into one inequality with a bespoke budget structure.** About 7k lines were stranded this way (item 3). The modular "first passage → abstract output energy → twisted children" split survived.
- **Fix conventions once.** This covers the bit orientation, where the 1/N(d) normalisation sits, and source-scale padding (item 4).
- **Large-sieve recursions on sSup-defined block norms never closed.** The scalar exponent bootstrap α ↦ 2−1/α did (item 5).
- **Prototypes over `ZMod n` and ℤ² were rewritten over O and then abandoned.** Write O-native lemmas from the start (item 8).
- **Dead does not always mean abandoned.**
  - `outside_mean_square_of_theta_models` was live for the upstream 23/24 proof and is dead only because the QRH spine re-chains it.
  - MAP.md's counterparts for Lemma first-transfer (`first_passage_transfer`) and Prop R (`completed_mean_square_of_exact_reflection`) are dead. The live ones are `first_passage_full_uniform_transfer` and `uniform_polynomial_mean_square_of_exact_reflection`.
