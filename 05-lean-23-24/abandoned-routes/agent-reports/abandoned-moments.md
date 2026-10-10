# Abandoned proof attempts in `DirichletL/Moments/`

All paths are relative to `lean/OAI/NumberTheory/DirichletL/` at commit `fd4aeeb2e`. Line numbers come from `audit/upstream/reach.tsv`. "Live" means reachable from the final theorems.

## Overview

`Moments/` (namespace `SevenEighths.CenteredMoment*`) formalizes Lemma `lem:plain` of paper78 (§`sec:plain`, "Fourth moments with short prime factors", lines 12477–14985). That lemma is the centered fourth moment Σ_k |S(n₁)S(n₂)Q|² ≪ Z^{M+ε}. Its proof uses two Poisson transforms, the Θ-row rectangle cancellation, and a finite induction on the width M.

The pruning deleted 95 whole modules and 1,397 non-auxiliary declarations, about 23.4k lines in all: 6.4k in the deleted modules and 17.0k inside modules that survived. Of these, 597 are not referenced by name anywhere in the original tree, so they are the endpoints of dead chains.

Almost none of this is an unfinished attempt. There is no `sorry` or `axiom` in `Moments/`. The Prop-valued hypotheses found in dead files (`SmoothedBound`, `CanonicalBound`, `SharedBound`, …) are all proved in the same files. Most of the dead code consists of complete, superseded variants. A module was cloned, one hypothesis or quantifier was changed, and the whole chain was proved again. This leaves parallel families side by side: `Physical*` vs `FixedQ*`, `Tests` vs `Profiles`, unseeded vs seeded, unfrozen vs frozen.

The diffs between dead and live twins therefore show directly which hypothesis could not be discharged downstream. This machinery is used only by the 7/8 proof. The 23/24 proof (`audit/reach-qrh.tsv`) uses nothing from `Moments/`.

---

## 1. Θ-row (exceptional) block chain: three generations (≈2.4k dead lines, 39 modules)

**Goal.** This chain bounds the transformed rows whose inducing character lies in the fixed family Θ (paper §"Rows with inducing characters in Θ", l.14312). It counts such rows (`eq:exceptional-row-count`, l.14380) using the forcing ideal, then sums the blocks dyadically.

**Generations (evidence):**
- **Generation 1, fixed detector tests.** Example: `CommonExceptionalCanonical.lean:26 actual_common_canonical_bound`, together with `CommonExceptionalWindow`, `ExceptionalAsymmetricCanonical`, `ExceptionalCanonicalShell`, `ExceptionalHeight`, `SecondExceptionalOriginalBlock`/`CanonicalPair`/`Budget`. These are stated for `p q : Tests`. `Tests` is defined at `ExceptionalAmplitudePair.lean:16`: two `detectorSchwartz` functions with `index ≤ 2` and `sigma ∈ [0,1]`.
- **Generation 2, finite profile families.** `FiniteProfileExceptionalCommonCanonical.lean:26` and `FiniteProfileExceptionalPhysical{Block:53, Dyadic:38, Pair:35, Rows:36, Saved:42, Uniform:44}`. The same statements, now for `p q : Profiles wlo whi` with a seminorm set `Sprofile`.
- **Generation 3, live: `FiniteProfileExceptionalFixedQ{Block, Dyadic, Pair, Rows, Saved, Uniform, Window, Canonical}`.** A `diff` against generation 2 shows exactly one change in every theorem statement: the hypothesis `IsCoprime Q C →` is deleted. The `sourceBudget` definitions in `…PhysicalBlock`/`…PhysicalBudget` stay live and are reused. The canonical bound gains a constant factor `qcost := N(commonRadical Q Q)^(2/3)`.

The root cause is in the counting lemma. The dead `SecondExceptionalCount.lean:42 actual_canonical_exceptional_count` takes `(hQC : IsCoprime Q C)`. Its live replacement `FixedQForcingCount.lean:180 actual_canonical_exceptional_count_fixedQ` drops it. In its place it splits the forcing ideal into a part coprime to Q and a part dividing rad Q (`forcing_dvd_good_mul_radical`, `fixed_radical_power_loss`).

**Obstruction (evidence, with a conjectured mechanism).** `IsCoprime Q C` could not be discharged at the call site. In the live `Energy/` layer, Q is instantiated as `internalQ Q η₀ := Q ⊓ η₀.modulus` (`NaturalFixedRaySourceInternal.lean:15`, 267 references): the fixed modulus enlarged by the source character's modulus. C is the common divisor of the transformed columns, and nothing prevents it from sharing primes with η₀'s modulus.

**Assessment.** This is a hygiene fix, not a new idea, but the trick is clean: pay N(rad Q)^{2/3} = O(1) instead of assuming coprimality. About 2,000 lines were rebuilt to remove one hypothesis. The lesson: run the call-site discharge check on a hypothesis before cloning a whole chain that carries it.

## 2. Uniformity of the masked rectangle cancellation (≈1.5k lines; `NaturalFixedRaySource*`, `Cancellation`, `MaskedVolume`)

**Goal.** This is paper `lem:centered-lattice-cancellation` (l.14546). For Θ-induced characters, the main terms of S(X₁)S(X₂) − S(Y₁)S(Y₂) cancel when X₁X₂ = Y₁Y₂ = T, leaving a bound of size √T/L.

The cancellation itself was always proved. The four successive versions differ only in **how uniform the constant is**:
1. `MaskedVolume.lean:15` and `Cancellation.lean:52 masked_rectangle_saving`. The error is `#divisors(R) · Σ_{r mod c} |χ(r)| · pvControl W`, so it depends on the full modulus c of the character.
2. `HeckeCancellation.lean:29 actual_rectangle_saving` (dead). It uses `volumeControl Q χ W = #divisors(χ.modulus)·#(O/fixedPeriod Q)·pvControl W`, which depends only on the fixed Q. The definitions are still live; the theorem is not.
3. `NaturalFixedRaySourceDegree.lean:55 rectangle_uniform_degree` and `:158 slots_uniform_degree`. These give an explicit `Z^ε (1+‖t‖)^J`, but the quantifier order is `∃ J, ∀ Q W₁ W₂ (hs hW), ∃ C, …`: **C depends on the profiles W₁, W₂**. The workaround was to fix the profiles to the detector family: `NaturalFixedRaySourceDetector{Profile:66, Slots:37, Cancellation:36, Capped:14}`, `AllocatedDetectorCapped.lean:16`, and the `Tests` structure from cluster 1. Each is uniform in `(reverse, index ≤ 2, σ, t, h)`.
4. Live: `FiniteProfileExceptionalRectangle.lean:13 rectangle_source_control`, which is `∃ J S, ∀ Q, ∃ C, ∀ W₁ W₂, … ≤ C·(sourceControl S W₁·sourceControl S W₂)·Z^ε(1+‖t‖)^J·(√T/L)`. Here the profile dependence is an explicit finite seminorm product.

**Obstruction.** Evidence: generation 3 and the detector specialization were abandoned once the explicit seminorm form existed. The paper states that requirement: "the bound uses finitely many smooth seminorms and a fixed polynomial in the separated norm-twist heights" (lem:plain). Conjecture: the induction feeds the lemma profiles that are not detector functions, such as `logProfile`, `normPowerProfile` and `annularProfile ∘ paperRadialFourier`. So uniformity over a fixed parametric family cannot be closed recursively.

**Assessment.** This was the right refactor; the evolution matches the paper's uniformity clause step by step.

## 3. Decentered positive majorant for the source → children step (≈0.8k lines)

**Main declarations.** `ActiveCanonicalSource.lean:26 positiveChildren`, `:60 SmoothedBound`, `:81 uniform_smoothed_bound`, `:116 canonical_source_to_positive`. Its twin `CanonicalRetainedSource.lean` has the same names, using the "retained" allocation instead of the "active" one. Supporting modules: `ActiveSmoothedSource.lean:20`, `SmoothedRetainedSource.lean:20`, `ActiveSourceRetained.lean:18`, `SourceRetainedEnergy.lean:17`, `SecondSourceCost.lean:21`.

**What it did.** It bounds the restricted source energy by `D·N(L)^ε·Σ_a weight·Σ_z Φ(…)(‖allocatedPositiveRow … X₁ X₂‖² + ‖allocatedPositiveRow … Y₁ Y₂‖²)`. Each `SmoothedBound`/`CanonicalBound` Prop is proved, uniformly over subsets J. This is a finished chain.

**Obstruction (strong evidence).** The majorant is *decentered*: it replaces the rectangle difference by the sum of the X-square and the Y-square. Compare the live `RestrictedSourceNormalization.lean:30 normalized_restricted_child_energy`. It is an exact allocation decomposition whose children carry `liveProfile` (`SourceLiveColumn.lean:26`), and that is built from `profileCoefficient … X₁ X₂ Y₁ Y₂`, i.e. it keeps the `idealRectangle`. The children are then repackaged as bundled `Input`s via `AmplificationChildInput.lean:63 child`, which has 967 references. The paper is explicit that the comparison product's "common character, mask, and norm power" must be preserved through the transforms for the Θ main terms to cancel (sec:plain intro).

Conjecture: the decentered route can serve only the uncentered range A ≤ 5M/6.

**Assessment.** It is a valid shortcut for the uncentered range, and it is simpler than carrying the rectangle. But the final proof needed one route that works for both ranges.

## 4. Second-transform descent: hypothesis-shape iterations (≈1.0k lines)

**Main declarations (dead):** `SecondSourcePowerDescent.lean:67 original_subsets_power_descent`, `SecondFrozenPowerDescent.lean:47`, `FirstNestedSeededGaussianPower.lean:195 nested_power_descent`, `SecondUniformSubsetSource.lean:78`, `SecondLiveSourceDescent.lean:35`, `SecondNonexceptionalLiveSource.lean:32`, `SecondNonexceptionalCost.lean:51`, `SecondRetainedRatioCost.lean:52`, `SecondSourceSeededPowerDescent.lean:53`.

**What it did.** This is the "power descent" (paper `sec:plain-second-transform`, l.13596). The energy after the second transform is bounded by four terms, given an induction hypothesis on the children (`physicalFactors`, `SecondSourcePowerDescent.lean:43`):
- main: `Z^{2δ+ε}·√(E₁E₂)·…`;
- exceptional: `Z^{…−2r/3} K^{5/6}`;
- diagonal: `K·|Ŵ(0)|`;
- tail: `K·Z^{−saving}`.

**The variants form a grid, and only one corner is live.** All the changes below come from `diff`s of the theorem statements.
- **Unseeded → seeded.** The live version adds `N(seed) ≤ Z^{Bseed}` and `seededFactors := ![a0, a1, a2, a3/seedNorm]`, so the tail gains a factor 1/N(seed). This is paid for by asking for saving + Bseed (`SourceInputTailSeed.lean:18 paid_seed_power`). Without it, the tail does not decay in the Möbius seed. Conjecture: summing over seeds would then lose up to Z^{Bseed}.
- **Unfrozen → frozen.** In the dead versions, the child induction hypothesis is `childEnergy ≤ E·N(τχ.modulus)·dyadicScale(n 1)·(1+|v|)^{2J}`, required for *every* allocation b. In the live `SecondFrozen*` versions it becomes `E·childEnvelope η C D U n·max(1, 1/retainedRatio n)^{1/6}·…`, required only when `frozenCoefficient b … ≠ 0`. The dead hypothesis asked for a bound on children that never occur, measured in the child's own conductor, and it apparently could not be discharged by the induction. The live hypothesis uses an envelope that depends only on parent data, and accepts an explicit (n₂n₃/n₀n₁)^{1/6} penalty.
- A third live axis (`SecondFrozenReferenceSeededPowerDescent`) replaces η.modulus by a reference `qref`.

**Assessment.** This is a strong example of hypothesis tuning. The dead versions are correct theorems whose *premises* were unprovable at the recursive call. The `^{1/6}` ratio penalty is a real loss that the live proof accepts.

## 5. Zero frequency: the paper's "added row zero" estimate (≈1.5k lines)

**Main declarations.** `AddedZero.lean`: 33 dead declarations, for example `:409 idealAddedZero_source_energy`, `:400 added_zero_exponent : a/3 + 4θ/3 ≤ a − s + 2θ` (given `s ≤ (a+θ)/6`), and `:220 tuple_added_zero_energy` (`≤ C·B²·T⁻¹·Y^{4/3+ε}`). Also `AddedZeroUniform.lean:307`, `SourceZeroEnergy.lean:244`, `SourceSecondZeroEnergy.lean:300`, `ZeroMode.lean:345–456`, `ZeroMeanBridge.lean:353` and `NaturalCoreFloor.lean:15`.

**What it did.** `AddedZero` is a line-by-line formalization of paper78 l.13571–13590 (§`sec:plain-gauss-enlargement`). G(u,0) = 0 unless u is a sixth power. There are O(Y^{1/6}) sixth powers, and |G| ≤ Y^{1/2}. This gives Z^{a/3+4θ/3} ≤ Z^{a−s+2θ}. It is completely proved.

**Why dead (evidence plus conjecture).** In the paper, this estimate is needed only because "the row zero is added … at the final smooth ball". The formal proof never adds the row zero. The live radial weights carry a `keep : O → Prop` predicate (`RadialEligibleEnergy.lean:14`), and the floors assume `keep z → z ≠ 0`. The genuine first-transform zero frequency (paper l.13180–13190, the powerful-product count) is live, in `ZeroMode.lean:92 powerful_pair_count` and `:282 ideal_zero_energy_bound`. The other dead zero-energy theorems are older unbundled versions of `SourceInputZeroUniform.first_input_zero_bound`.

**Assessment.** This points to a possible simplification of the paper: drop the sixth-power added-zero paragraph by keeping z = 0 punctured, as the Lean proof does.

## 6. Exponent ledgers (≈0.85k lines)

**Main declarations.** `DescentLedger.lean`: `firstSaving` (live), with `:113 uncentered_exceptional_bound` and `exceptional_excess_identity` dead. `FinalStageLedger.lean:8 capped_loss : M/6 ≤ 2v/3 + max(M/4 − v) 0` and `:15 final_exponent`. `FirstSecondPaidLedger.lean:47 actual_stage_saving`. `FirstExceptionalPaidScales.lean:15 Branch`, a 31-field record of one branch's exponents, with `:87 actual_branch_stage`. Also `FirstExceptionalChildPaid` and `FirstSecondCapacityShift`.

**What it did.** It closes, in log_Z exponents, the per-stage inequality `A − 5M/6 − F₁ − F₂ − max(M/4−v,0) ≤ max(A−M,0) + 3σ + 5·loss/6`. This is exactly paper eq. (old-eq:2.19), l.14759–14767, with L = M/4 (old-eq:2.2). `uncentered_exceptional_bound` needs `hA : A ≤ 5M/6`, which is the paper's uncentered range (l.12594).

**Why dead.** It was superseded by the stage geometry in `Energy/*` (for example `Energy/ZeroReferenceGeometry`, `ReferenceLowBranchGeometry`, which carry `5*M/6` and `(6κ−1)ℓ ≤ M`).

**Numerical evidence.** `capped_loss` is *tight* at v = M/4: the worst case of the left side of (2.19) is exactly A − M. **There is zero slack** in the centered exceptional deficit.

**Assessment.** These are clean, linarith-checkable certificates of the paper's numerology, and the most readable record of where the 5/6 and M/4 come from.

## 7. Prime-slot bounds for general rows and live-slot removal (≈0.35k lines)

**Main declarations.** `PrimeSlotBound.lean:18`, `PrimeSlotUniform.lean:18`, `PrimeSlotExternal.lean:19`, `ActualLiveSlotRemoval.lean:19`, `PrimeSlot.lean:38–71`. Also `NaturalFixedRaySourcePrimeBound.lean:17` and `PrimeGate`.

**What it did.** It proves eq. (old-eq:3.5) (`sec:weighted`, l.12876): ‖normalized prime slot‖² ≤ C(1+|t|)^deg Z^loss P^κ under `hbeta : 51/100 ≤ beta` and `hκ : 2*beta-1 ≤ κ`. This is the paper's β_* ≤ (1+κ)/2. It is stated for arbitrary rows `(η m A z)` with `¬FixedInducingRow`.

`actual_live_slot_removal` is a neat trade-off lemma. It chooses R ⊆ J with `κ·Σ_{R} w ≤ F/6 + κ·mesh` and bounds the product over J by Z^{ε+F/6+κ·mesh} times the product over J∖R. That is, it removes slots by their pointwise prime bound, paid by capacity.

**Why dead (evidence).** The live slot bounds are `NaturalFixedRaySourceNaturalPrime.lean:35/73` and `NaturalFixedRaySourceRayPrime.lean:18`, stated for "natural" source characters with a relative character θ. The paper requires natural characters, with the redundant radical deleted at O(P^{-1/2} log Z) cost (l.12895). Conjecture: the general-row version did not match the shape of the actual source rows.

**Assessment.** The general-row removal lemma is reusable and simpler to state than the live capacity-removal modules.

## 8. Early element-level (CRT) Poisson identities and bottom stage (≈0.75k lines)

**Main declarations.** `Bridge.lean:80 triple_pair_poisson` and `:121 masked_triple_pair_poisson` (Poisson for ξ·χ_a·conj χ_b through an explicit `productCharacter` CRT). `Initial.lean:63 allocated_pair_poisson`. `FirstPoisson.lean:114 coprime_pair_poisson`. `FirstWholeAssembly.lean:16 finite_whole_children`. `BottomStage.lean` (the trivial band bound `≤ … Z^{M+σ}`).

**Status.** These are complete exact identities. The live route works with ideal-level correlations instead (`Support.lean:151 idealCorrelation`, `FullCorrelation`, `FirstWholeKernel.dilated_frequency_whole_kernel`), avoiding explicit generator normalizations such as λ² ∣ a−1.

**Assessment.** This is a representation change. It is of no research interest beyond showing that the ideal-level API was the better abstraction.

*Not clustered:* about 1.9k scattered dead lines of reflection and radial infrastructure (`ReflectedAnnuli`, `OriginalReflectionErrorMass`, `ReflectedUniformPair:61`, `UniformReflectionProfile`). These are leftover helper lemmas rather than a coherent alternative route.

---

## Worth exploring further

1. **Simplify the paper's zero-row step (cluster 5).** The Lean proof never uses paper78 l.13571–13590 (the sixth-power count for the added row zero), because rows stay punctured (`keep z → z ≠ 0`). *Next step:* confirm that no live declaration reaches `gaussRow_zero_*` or `sexticGauss u … 0`. Then rewrite §`sec:plain-gauss-enlargement` to keep the zero punctured, and delete the θ_N "added row" allowance from the C_*ξ budget.
2. **Slack analysis of the exponent ledger (cluster 6).** The ledgers are pure linear inequalities and already machine-checked. `capped_loss` shows the centered deficit is tight at v = L = M/4. *Next step:* extract the full set of stage inequalities. Use `DescentLedger`/`FinalStageLedger` plus the live `Energy/*Geometry` constraints: 5/6, M/4, 6κ−1 ≥ 7/2, the M/15 margins of eq:comparison-uncentered-margin, and the ^{1/6} ratio penalty from cluster 4. Run an LP or Fourier–Motzkin slack computation to find which constraints bind the 7/8 exponent. Any non-binding constraint is a candidate for simplification; a binding one is where a better exponent would have to come from. This is speculative as an exponent improvement, but cheap and informative.
3. **Reuse `ActualLiveSlotRemoval` and the decentered majorant for the uncentered range (clusters 3 and 7).** Both are finished, short and correct for A ≤ 5M/6. *Next step:* in a fork, check whether the live A ≤ 5M/6 band can be closed through `canonical_source_to_positive` plus `actual_live_slot_removal`, keeping the rectangle-preserving route only for A > 5M/6. If so, the uncentered half of the induction could become noticeably shorter.

## Lessons / dead ends to avoid

- **Check that a hypothesis can be discharged at its call site before cloning a chain that carries it.** `IsCoprime Q C` (cluster 1) and the unfrozen child hypothesis (cluster 4) each cost a full parallel family of 6–10 modules.
- **Quantifier order is the uniformity.** `∃ J, ∀ W, ∃ C` is useless for an induction over moving profiles. State constants as explicit seminorm products from the start (`sourceControl S W`); restricting to a parametric family such as `Tests` instead was a dead end.
- **Never decenter before the Θ cancellation.** Majorizing a rectangle difference by |S(X)|² + |S(Y)|² gives up the cancellation that the A > 5M/6 range needs (paper sec:plain intro).
- **Induction hypotheses should depend on parent data, not on the child's own conductor.** Restrict them to children that actually occur (`frozenCoefficient ≠ 0`). Buy decay in a Möbius seed explicitly, as `paid_seed_power` does with saving + Bseed.
- **Formalizing the paper verbatim can be wasted work.** `AddedZero` formalizes a paper paragraph that the chosen representation (punctured rows) makes unnecessary. Check the formal route's representation before formalizing a step.
- **No numerical slack was found in this area.** (2.19) is tight at M/4, and the ledgers close only with exact margins. Improvements will not come from re-tuning these ledgers.
