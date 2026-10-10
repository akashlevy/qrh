# Abandoned proof attempts in `Descent/`, `Inversion/`, `Dictionary/`

All paths below are relative to `lean/OAI/NumberTheory/DirichletL/` at commit `fd4aeeb2e`. Line numbers come from `audit/upstream/reach.tsv`.

## Overview

These three directories formalize paper78 §9 (`sec:inverse`):

- `Descent/` proves the canonical marked estimate (`lem:canonical-moment`) by a finite descent that uses two masked Poisson transforms.
- `Inversion/` proves the one-Poisson initialization that deduces `lem:marked` (paper78 l.11601 ff.).
- `Dictionary/` connects the detector's Hecke fibres to that initialization.

The 23/24 descent (cutoff D^{1/1000}, at most 3000 levels) is in `RowCompletion/`, not in these directories.

There are 914 dead declarations, about 17.8k lines in total:

| Directory | Dead declarations | Dead lines |
|---|---|---|
| `Descent/` | 376 | 10.1k |
| `Inversion/` | 477 | 7.0k |
| `Dictionary/` | 61 | 0.8k |

None of them contains `sorry`. Every abandoned branch is a finished proof that was replaced by a more general restatement. In one case (cluster 2) the replaced theorem depends on an induction-hypothesis Prop that nothing ever proves.

The main mechanism is "generalize by copying". 39 of the 76 deleted modules (4.7k of 7.7k lines) have a surviving successor that shares at least 70% of its lines (difflib line ratio). The copies differ only in 2–6 binder lines, and the new copy often still imported the old one. No dead branch implements paper2's adaptive cube cutoff H_c, and none gives a better exponent.

---

## 1. Ladders of quantifier order and type uniformity (about 2.7k lines)

**Declarations (dead, then live):**

- `Descent/PrioritySourceEnergy.lean:16` `actual_priority_source_recursive`
  → `PriorityMovingRadius.lean:16`
  → `PriorityMovingRadiusWindowOrder.lean:17`
  → live `PriorityMovingRadiusWindowOrderUniform`.
- The same four-step ladder exists for:
  - `PriorityRetainedBound.lean:16` (125 lines)
  - `PriorityCellEnergy.lean:15`
  - `SecondBlockEnergy.lean:14`
  - `SecondDeletedMovingEnergy.lean:15` and `SecondDeletedMovingMass.lean:16`
- The GlobalPriority ladder runs from `GlobalPriorityRays.lean:16` ("retained") to `GlobalPriorityFilteredRays.lean:16` and then to `…Uniform`. Its other steps are:
  - `GlobalPriorityEnergy:18`
  - `GlobalPriorityAggregate.lean:32` (141 lines)
  - `GlobalPriorityFiltered{Energy:19, Aggregate:18, Step:18}`
  - `GlobalPriorityFixedChildStep.lean:18` (182 lines)
- The terminal (short-completion) energies form the same kind of ladder:
  - `Descent/ShortReflectedEnergy.lean:38` and `ShortReflectedPower.lean:38`
  - `Inversion/ShortCanonicalNormalization.lean:80,167`
  - `Inversion/CanonicalShortAttachment.lean:154` (95 lines)
  - all are replaced by the live `Reflected*UniformDegree` modules.

**What it attempted.** The second-Poisson energy for one dyadic cell, summed over the priority branches of `eq:inverse-slot-priority` and over rays. This is paper78 l.10654–11400, from "The second Poisson transformation" through "The positive canonical children".

**How far it got.** Every version is complete. Successive copies have similarity 0.86–0.97. The diffs show four generalizations:

- **(a) Fixed radius to per-cell radius.** `nonzeroChildFrequencyBall … R` becomes `… (R (index p x))` (the "moving radius").
- **(b) Retained to filtered.** The filtered version adds `let source:=InverseMomentGlobalRetainedGates.geometrySource p originalSource b X; hi≤b →`.
- **(c) Window order.** `(Aker J K:ℕ) : ∃ ω₁ ω₂ …, ∃ C` becomes `(Aker K:ℕ) : ∃ ω₁ ω₂ …, ∀ J:ℕ, ∃ C`.
- **(d) Uniformity.**
  - `variable {ι σ:Type}` at file level becomes `∀ J, ∃ C, ∀ {ι σ : Type}`.
  - For the terminal energies, `(q) … ∃ (degree) C Z₀` becomes `∃ degree, ∀ q≠0, ∃ C Z₀`.

**Obstruction (evidence).** The final induction predicate forces this shape.

- `RankEnergyExists n` (`CanonicalRankExistenceData.lean:35`) reads `∃degree, ∃B, ∀q≠0, ∃C Z₀, …, ∀{σ:Type} …`.
- The live successor step (`CanonicalRankExistenceSuccessor.lean`) first gets the four child windows `w11…w22` from `actual_complete_normalized_rank_step`.
- Only then does it call the induction hypothesis on each window to get degrees `J11…J22`, and it sets `J := J11+J12+J21+J22`.
- If the windows depended on `J`, as in version (c)'s predecessor, the argument would be circular. If `C` depended on the index types, `∀σ` could not come after `∃C`.

This is exactly paper78 `lem:finite-seminorm-propagation` (l.9367): the indices are "selected backward … fixed before any restriction". The paper states this abstractly. The agents rediscovered it only when the top-level induction failed.

**Assessment.** Not mathematically interesting. It is the costliest process failure in these directories.

## 2. Termination and parameter budget (about 0.6k lines)

**Declarations:**

- `Descent/Basic.lean` with the paper's constants:
  - `total_step_loss:157` (`D*eta ≤ eps/160 → D*(40η+τ+π) ≤ 3ε/4`)
  - `depth_margin:163`, `canonical_transition:52`, `childF_upper:93` (`childF ≤ F+11η`)
  - `initial_*:145–155`
- `SourceDepth.lean:32` `actual_source_path_lengths`. It follows an explicit sequence `x : (k:ℕ) → MarkedSecondSource …` and proves `F n ≤ F 0+11nη ∧ M n+3nd/2 ≤ M 0`.
- `ActualChildStateDepth.lean:12`: the same argument with 15nη.
- `ActualDepthBudget.lean:22` `actual_depth_parameters` (separate η, τ, π) and `:82` `actual_depth_row_negative`.
- `CanonicalLowerMoments.lean:19` (a Prop that contains `M<Mmax`) and `:32`.
- `FirstTwoPassLowerRank.lean:24` `original_two_pass_from_lower_moments` (86 lines).
- `Inversion/TerminalWidths.lean:534` `canonical_terminal_saving`. It reproduces paper78's terminal case exactly: `d ≤ cstar/200`, `η,τ,π ≤ cstar/1000` imply exponent `+π ≤ F₀ − cstar/4`.

**What it attempted.** The paper's backward induction on depth h ≤ D, with:

- c_h = c_* − 7hη
- D = ⌈(M_max+2)/d⌉+1
- F' ≤ F+11η
- a per-step loss of 40η+τ+π (`eq:canonical-reduction-contract`, l.10020; "Order of choices and termination", l.11474)

**How far it got.**

- The scalar lemmas and the path-length lemmas are complete.
- `CanonicalLowerMoments` is an induction-hypothesis Prop that nothing ever proves.
- A strong induction over a real `M<Mmax` is not well-founded as stated. So `FirstTwoPassLowerRank` is a conditional theorem whose hypothesis was never discharged.

**What replaced it (evidence).** A plain `induction n` on an ℕ rank (`CanonicalRankExistenceInduction.lean`), with:

- `rankRowCap n cutoff = 3(n−1)cutoff/4`
- `rankTotalCap = L−15nη`
- `rankMargin = c*/2+7nη`
- `rankLoss = 62nη`

`CanonicalRankMoments` is `CanonicalLowerMoments` with `M≤Mmax` (similarity 0.93).

`actual_recursion_parameters` fixes `cutoff := c/1000` and a single `η := min(…)`. The successor step passes η again in the `tau` and `eps` slots. As a result, the paper's separate tolerances τ, π, τ_ref, π_ref collapse into one η. The formal step loses more than the paper's step (11η becomes 15η, and 40η+τ+π becomes 62η), so the lemmas carrying the paper's constants became dead.

**Was there an H_c attempt? No.**

- Every Descent variant, dead or live, uses a fixed terminal test `cutoff ≤ ell+V` with `η ≤ cutoff/16` or `/32`.
- paper78's own scheme also uses a fixed d (`d ≤ c_*/200`).
- The dead "paper-faithful" code here follows paper78, not paper2.
- In this architecture the number of levels only enters the η budget. An adaptive cutoff would change neither the proof length nor the exponent.

**Assessment.** The live rank formulation is cleaner than the dead form that follows the paper.

## 3. Per-source data replaced by data that is constant on each cell (about 1.5k lines)

**(a) Exact cutoff replaced by a dyadic one.**

- `FirstFullRapidTail.lean:13` `firstExactPhysicalCutoff` uses the radius `Z^(firstPhysicalHeight … (Real.logb Z (primeProductNorm p D)) …)`.
- `:21` `full_first_whole_marked_rapid_tail` is 105 lines.
- Its live successor is `FirstDyadicTail.full_first_dyadic_marked_rapid_tail` (similarity 0.83). It uses `firstDyadicExponent` and adds `2≤Z^eta`.

**(b) Correlated radius replaced by a uniform one.**

- `SecondCutoffScale.lean:13,44`: `correlatedSecondRadius` uses ‖d‖² of the actual element.
- `Inversion/SecondUniformCutoff.uniform_radius_bound` replaces it with Z^(δ+η). The proof is otherwise identical.

**(c) Fixed slot lists replaced by deleted or varying lists.**

- Dead modules: `SecondAllModes.lean:58`, `SecondAllIntegral:17`, `SecondModeEnergy:17`, `SecondOriginalBound:17`, `SecondPhysicalEnergy:19`, `SecondPhysicalRecursive:19`.
- Their live successors are the `SecondDeleted*` modules (similarity 0.82–0.90). The live versions add
  `hd₁ : ∀ x∈source, ∀ i∈deleted₁ x, i∈x.cube.support∪x.firstCommon ∨ (Ideal.span {p i})∣x.quotient`
  and use `lists₁ i \ deleted₁ x`.
- Likewise `SecondBalancedPhysical.lean:19`, `SecondDensityAssembly`, and `SecondInheritedEnergy` became the `SecondVarying*` modules, with `lists : MarkedSecondSource → σ → Finset ι`.

**Paper references.**

- `eq:inverse-ratio-bounds` (l.10042): dyadic η-ratio localization.
- `eq:inverse-slot-priority` (l.9950), together with "restores its original slot lists separately by the child zeros".

**Obstruction (evidence).** A child is a canonical datum only if:

- its slot lists are fixed and disjoint, and
- one coefficient measure serves the whole cell.

Exact per-source radii or lists break both conditions. The dead versions are correct only in the special case "deleted = ∅, one radius".

**Assessment.** These were natural first drafts and hide no improvement.

## 4. Whole-cube priority vs. global priority (about 1.2k lines)

**Declarations.**

- `Descent/WholePriorityEnergy.lean:17` `whole_priority_retained_fixed_ray` (127 lines)
- `WholePriorityRays.lean:16` (111)
- `WholePriorityPrincipal.lean:42` (90)
- `WholePriorityMass.lean:13` (80)
- `WholePriorityCollapse`
- `WholePriorityBranches.lean:39,78`
- `Inversion/WholePriorityRapidTail.lean:92` `whole_priority_tail_weighted` (109) and `:249`
- `WholePriorityValidSource.lean:150–253`

**What it attempted.** These theorems bound the retained, principal and tail first-pass contributions for one fixed cube. Their signature takes a single `(cube : CubeCoordinates ι) (firstCommon firstDivisor : Finset ι)` and sums `∑ D∈pool.powerset, ‖selector D‖ …`.

**What replaced it.** The live `GlobalPriority*` modules take `(original : Finset (InverseFirstPriorityParents.Source ι Jo))`. There, `x.cube` varies with the source, and `∑ x∈original` sits inside the bound.

- `GlobalPriorityMass.global_priority_weighted_branch` replaced `whole_priority_weighted_branch`.
- `GlobalPriorityTail.global_parent_tail_rapid` replaced the whole-cube tails.
- Some whole-cube infrastructure survives as components of the global source (`wholeAssignedParents`, `wholeExtractedSupport`).

**Obstruction.** *Conjecture, supported by the paper.* Block (A) (`eq:canonical-block`, l.9984) keeps the sum over cube ideals b inside the square. A per-cube energy must instead be summed outside the Hilbert norm, which costs a factor ≍ Z^ℓ unless the cubes are reassembled. There is no Lean evidence of a specific inequality that failed.

**Assessment.** A dead end. The global form is the paper's form.

## 5. Initialization: choosing the column pool (Inversion, about 1.3k lines)

**Declarations.**

- `Inversion/InitialCommonPool.lean:38` `residual_eq_common_input`
- `InitialCommonForcing.lean:40`
- `InitialCommonTuples.lean:50`
- `InitialCommonTransform.lean:45` `original_tuples_fourier` (86 lines)
- `InitialLargePool.lean:98,139`
- `InitialLargeFourier.lean:20` `original_tuples_large_fourier` (86)
- `InitialCompletePool.lean:56`
- dead smoothed-Poisson encodings of the overlap, e.g. `InitialPoissonBridge.lean:299,363`

**What it attempted.** To write the residual polynomial R₀ (`eq:initial-residual-polynomial`, l.11717) as a row of the canonical family over a fixed pool of primes. The diffs show three successive choices of pool:

1. **Common pool.** `F := tupleColumns S (I.pi L) P j`, which depends on the assigned tuples.
2. **Large pool.** Any `F ⊇ tupleColumns`, stated as a hypothesis `hsub`.
3. **Excluded pool (live, `InitialExcludedFourier`, similarity 0.91 with version 2).** F lies outside a set E of excluded primes, with `hη : ∀c, ¬outside E c → η c = 0`.

**Why the third choice won (evidence).** `RankEnergyExists` fixes the pool to `primePool (InitialMeanSquare.outsideSquarefreeIdeals (reflectionExcludedPrimes q) D)`. Only the excluded version matches it.

**Assessment.** A design iteration that converged on the canonical choice. QRH's Milestones file likewise says that excluded primes are "carried as a mask m".

## 6. Explicit analytic towers and root-kernel seminorms (Inversion, about 1.4k lines)

**Declarations.**

- `Inversion/AmbientProfileTower.lean`:
  - `:145,196` `conditional_tower_integrable/bound`
  - `:448` `fullProfile_tower_uniform`
  - `:508` `first_second_profile_tower` (73 lines)
- `SecondProfileUniform.lean:393,448,568` (`first_second_tower_uniform`, 67 lines)
- `ClippingProfiles.lean`:
  - `:122` `rootSchwartz_seminorm` (85 lines)
  - `:266–322` `rooted_twist_*`
  - `:392` `first_root_budget` (`firstRootExponent y ≤ (9/2)*L`)
- `Descent/DescentProfiles.lean:75,132` and `DescentHeight.lean:37`

**What it attempted.**

- One joint Fourier density (a "tower") for the composite of the first Poisson profile (9 log coordinates) and the second (6 coordinates). The moment weight grows from `ambientWeight J` to `J+2K`. This is a concrete two-level case of `lem:finite-seminorm-propagation`.
- Explicit seminorm bounds for the exponentially rooted kernel (`eq:inverse-first-root-kernel`).

**What replaced it.** The live code uses common measures at each level (`second_profile_common_measure`) plus bounds polynomial in the height, carried inside the induction hypothesis: `CanonicalRankMoments` gives `≤ A*Z^(N+V+eps)*(1+‖s‖)^(2*degree)` for the test `childLogTest W s`.

**Obstruction (evidence).** The tower is complete but is fixed at two levels. The descent has n levels, and the polynomial form composes simply by adding degrees.

**Assessment.** Sound and fairly elegant, but less general than the induction needs. Not worth reviving.

## 7. Generic arithmetic tools that lost to specialized ones (about 1.5k lines)

**(a) Hybrid quadratic–cubic lemmas.**

- `Descent/Hybrid.lean`: `:237` `hybridRow_mobius`, `:256` `hybrid_first_mask_energy`, `:317` `hybrid_double_mask_energy`.
- `SieveComposition.lean:92` `quadratic_cubic_composition`. It is abstract: a quadratic sieve norm Q and a cubic sieve norm S give the bound Q·|C|·S.
- `Quadratic.lean:13,119`, `SlotProducts:87`, `CommonMeasure:39`.

These are the two mask resolutions of `lem:hybrid-energy` (l.7510), where the paper says it "resolves the two such masks". They were replaced by the live MaskedGrid, ExtractedGrid and HybridReduction modules, and by `quadratic_squarefree_product_column_energy`. *Likely reason:* the quadratic sieve needs squarefree products; `quadratic_product_column_energy` was replaced by its squarefree version.

**(b) Redundant punctures.** `CanonicalPunctureMask.lean:21–108` and `…Energy.lean:12–79`.

- Multiplying the puncture m by c leaves `normalizedColumnEnergy` unchanged when the primes of c lie outside the pool or already divide m.
- `reflection_pool_mask_energy` shows that `ActualFiber.maskElement q m` is redundant.
- It was unused because the mask costs only O(1) for fixed q.

**(c) Dictionary.**

- `InverseRawFiberSource.lean:18–109` requires `hlarge : modulus.absNorm < lower s*U^(widths s)`.
- `InverseMarkedReferenceCoefficients.lean` (normalized slot coefficients) and `InverseMarkedPadding.lean:51,123` (a padded, smoothed radial majorant) form a "padding" route.
- It was replaced by `InverseRawInitialEnergy.raw_rows_le_initial` and `deleted_rows_le_initial`.

**Assessment.** Item (a) is the most reusable dead code: short, general, and independent of the surrounding bookkeeping.

---

## Worth exploring further

1. **Bundle hypotheses in the surviving ladder (shrinks the 7/8 proof).**
   - The 16 live ladder modules total 2,191 lines: `GlobalPriority*Uniform`, `Priority*Uniform`, `SecondBlockEnergyWindowOrderUniform` and `SecondDeleted*`.
   - About 1,146 of those lines (52%) are theorem statements that restate the same ~40 hypotheses.
   - Next step: define a `structure` for the shared hypotheses of the second-pass cell. Thread it from `SecondDeletedMovingEnergyUniform` up to `GlobalPriorityFixedChildStepUniform`, then measure with `scripts/measure.py`.
2. **Try the generic hybrid lemmas in place of the grid machinery.**
   - The live hybrid modules total about 1,150 lines: `HybridReduction`, `MaskedGrid`, `ExtractedGrid`, `SourceHybrid`, `HybridBound`, `Sharp`, `Hybrid` and `Grid`.
   - Next step: restore `hybrid_double_mask_energy` and `quadratic_cubic_composition`, and try to derive the live `hybrid_fixed_block_sieve_norm` from them together with `sieveNorm_sharp`.
   - The expected saving is a few hundred lines, and only for the 7/8 proof.
3. **(Documentation only.)** The dead lemmas with the paper's constants (`Basic.total_step_loss`, `TerminalWidths.canonical_terminal_saving`) state paper78's numerics exactly. They could serve as a ready-made paper ↔ Lean crosswalk for a 7/8 "spine" file.
   - There is no exponent gain to be had. No dead branch is a cleaner descent or uses the adaptive H_c, so a paper2-faithful descent would have to be written from scratch. This matches the earlier conclusion that it would save at most ~15%.

## Lessons / dead ends to avoid

- **Fix the shape of the top-level induction predicate first.** The order should be windows, then degree, then level q, then constants, then index types. About 2.7k lines were proved twice or three times because this shape was settled late.
- **Do not use strong induction on a real parameter or explicit path sequences.** Discretize the row length into an ℕ rank.
- **Use data that is constant on each cell from the start.** Exact `logb` radii, correlated radii, fixed slot lists, per-cube sums and tuple-dependent pools all had to be rewritten.
- **Keep the cube sum inside the square.**
- **Carry bounds polynomial in the height through the induction hypothesis instead of building multi-level Fourier towers.**
- **Do not chase the paper's constants.** The formal proof needs 15η and 62η where the paper has 11η and 40η+τ+π. This does not affect the exponent.
- **Generalize by adding a parameter, not by copying.** About 60% of the deleted module lines are near-duplicates that the new copies still imported, so they were compiled with every build.
- **There is nothing hidden to mine here for 23/24.** The dead code has no `sorry`, no different exponent and no adaptive cutoff.
