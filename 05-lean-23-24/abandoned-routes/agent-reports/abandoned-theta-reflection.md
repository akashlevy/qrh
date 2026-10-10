# Abandoned proof attempts: cubic theta, reflection, cusps, Hecke

All paths are relative to `lean/OAI/NumberTheory/DirichletL/` at commit `fd4aeeb2e`. Line numbers refer to the original files. "Dead" means `reachable=0` in `audit/upstream/reach.tsv`. The original tree has no `sorry`, `admit` or `axiom` under `DirichletL`, so every dead declaration is fully proved. A branch was abandoned because it was superseded or not needed, or because it rests on a hypothesis that nothing ever supplies.

## Overview

After notation and auto-generated declarations are removed, the four directories contain about 1,200 dead declarations totalling roughly 15k lines. Eisenstein accounts for 652 declarations and about 6.4k lines, Reflection for about 5.2k lines, Hecke for 2.3k and Cusp for 0.9k.

Most of this is not failed mathematics. It is proved work that a later, more general sibling superseded. In Reflection alone, about 2.8k dead lines have a live twin that has the same name plus a suffix (`_uniform_degree`, `_type_uniform`, `_dependent`, `_on`, `_of_controlled`).

Three clusters carry real mathematical content:

1. An early attempt to run cusp reflection with the Kubota group Γ(3) only, using the unaveraged theta. It handled only the principal cusp and was replaced by the SL2(Z)-averaged theta.
2. A classical spectral apparatus that the final residue construction never used: a self-adjoint Laplacian, meromorphic resolvents, and Green's identity on cusp boxes.
3. Alternative routes to the theta relation τ(hp³) = |p|τ(h).

No dead code tries to prove that Kubota's θ is SL2(Z)-invariant. Even so, the dead and live pieces together reduce that invariance to one explicit identity between Gauss-sum residues (Worth exploring, item 1).

---

## Cluster 1. Γ(3)-only theta and principal-cusp reflection

**Main declarations (all dead).**

The cusp-reflection chain for the Kubota group only:
- `Eisenstein/ConjugateCuspProfiles.lean:84-99` `cuspBarProfile_conjugate_reflection (G : CubicKubota.levelThree) …`
- `:101-122` `conjugateSource_principal_cusp_mellin_entire (G : levelThree) (hc : complexMatrix G 1 0≠0)`
- `:124-131` `conjugateSource_cusp_mellin_reflection`, the Mellin functional equation at a cusp
- `Eisenstein/ReflectedCuspMellin.lean:638-680` `cubicSourceResidualFunction_levelThree` and `…_cuspDerivative (G : CubicKubota.levelThree) …`

∂̄-derivatives and cusp decay for the *unaveraged* θ = `cubicResidualFunction`:
- `Eisenstein/TranslatedFourierResidues.lean:493-598`, ending in `cubicResidualFunction_wirtingerBar`
- `Eisenstein/ResidualCuspDecay.lean:92-213` `cubicResidualFunction_wirtingerBar_cusp_decay`

Principal-stratum arithmetic:
- `Eisenstein/FixedConductor.lean:144-162` `structure PrincipalStratumArithmetic`
- `:164-250` `exists_principalStratumArithmetic`, which requires `(hc3:(3:Eis)∣c0)`
- `:329-345` `exists_constructed_principal_reflection`, which proves `smoothedValue = reflectedValue`

Projection lemmas that were never applied to θ:
- `Eisenstein/SourceProjection.lean:612-640` `sourceProjection_fixed`. Its hypothesis is `hf : ∀g : levelTwo, ∀w, f (g•w) = levelTwoComplexCharacter g * f w`.
- `sourceProjection_idempotent` and `sourceEisenstein_automorphy`, in the same range
- `Eisenstein/SourceQuotient.lean:280-457` `kernelSourceProjection_fixed/idempotent/covariance`
- `Eisenstein/CuspAverage.lean:363` `kernelSourceProjection_cubicEisensteinResidue_ne_zero`

**What it attempted.** This is the analytic core of paper78 Prop. `prop:completed-reflection` (paper78.tex:1679ff) and of paper2's cusp reduction (paper2.tex:3076-3097). At a rational point a/c, the ∂̄-profile of θ̄ is rewritten as a Mellin-reflected ∂-profile at a cusp. The dead chain does this for g = G·H with G in the Kubota group Γ(3) (`levelThree`). Its arithmetic side builds the level-3 matrix only when 3 | c₀, which is the case H = I in paper2.

**How far it got.** It is complete for the principal cusp. The dead files contain no version for the cases `v_λ(c)=1` or `(c,λ)=1`.

**Why it was dropped (evidence).** Paper2 handles the case `(c,λ)=1` with the Weyl-type matrix `H = (u₀ −1; 1 0)`. It then invokes "θ(γw)=θ(w) (γ∈SL₂(ℤ))", asserted without proof (paper2.tex:3089-3093), to reduce to the three cusps γ₀, γ₊, γ₋.

With Γ(3)-automorphy only, those strata cannot be reached. The live code avoids the problem by averaging. It defines `cubicSourceResidualFunction := sourceProjection cubicResidualFunction` (`ResidualFourierExpansion.lean:499`). This is an average over the 24 cosets `levelTwo ⧸ levelThreeInTwo ≅ SL₂(F₃)`, and it equals ¼·(E at ∞ + three S·Tᵇ translates) (`RationalTranslations.lean:138-147` `sourceEisenstein_four_cusps`). The averaged function is automorphic for `levelTwo = Γ(3)·SL₂(ℤ)` by construction (`cubicSourceResidualFunction_automorphy`).

Every dead `levelThree` statement has a live `levelTwo` twin:
- `GaussSum/MellinSupport.lean:166-201` `cubicSourceConjugateFunction_cuspDerivative_levelTwo` and `cuspBarProfile_conjugate_reflection_levelTwo`
- `conjugateSource_levelTwo_three_cusp_mellin_entire`, over `cuspRepresentative j` (j : Fin 3). These are paper2's γ₀, γ₊, γ₋ (`Cusp/ResiduesModThree.lean:212`).
- `ControlledStratumArithmetic … (mode : Bool)` (`FixedFrequencyExpansion.lean:273-308`), which covers both "a₀ primary" and "c₀ primary" and replaces the principal-only structure

Two pieces of negative evidence:
- No declaration ever states θ's invariance under SL₂(ℤ), and no lemma of the form `kernelSourcePullback M cubicEisensteinResidue = cubicEisensteinResidue` exists.
- `sourceProjection_fixed` takes invariance as a hypothesis, and nothing applies it to θ.

The dead `kernelSourceProjection_residue_fourier_mask` (`TranslatedFourierResidues.lean:286-296`) records what the live proof actually uses: the h-th coefficient of the averaged function is `if e(h)=1 then ¼(θ_h + 3·θ^{opp}_h) else 0`. So the formal theta is a mixture of the ∞-expansion and the opposite-cusp expansion, not Kubota's θ itself.

**Assessment.** The averaging trick is clever and sound. It turns an unproved global symmetry into a 24-term average whose cusp data come from an explicit Bruhat decomposition. The dead Γ(3)-only chain is a genuine dead end for the full reflection, because the Weyl-case cusps need SL₂(ℤ).

---

## Cluster 2. Classical spectral and Green-identity machinery for the Eisenstein residue

**Main declarations (dead).**

The unbounded self-adjoint Laplacian:
- `Eisenstein/VariationalLaplacian.lean:442-474`: the weak gradient is closable and densely defined
- `:513-576`: the variational resolvent is self-adjoint, positive and dense
- `:610-691`, ending in `kernelEnergyLaplacian_selfAdjoint : IsSelfAdjoint kernelEnergyLaplacian`

Global meromorphy:
- `Eisenstein/MeromorphicResolvent.lean:250-262, 329-367, 412-441, 475-563`, including `kernelEisensteinResolvent_meromorphicAt` and `kernelEisensteinL2Correction_meromorphicAt/graph`
- `:119` `kernelThetaNull_finiteDimensional`
- `Eisenstein/ScatteringCoefficients.lean:94-167`: `cubicEisensteinResidue_ne_zero`, `…_limit`, and `cubicEisensteinCorrection_not_continuousAt`, which says the pole is genuine
- `Eisenstein/SeedDefects.lean:151-183`: seeds are entire in s
- `Eisenstein/WhittakerContinuation.lean:100-279, 506-631`: `continuedNonzeroScattering`, the meromorphic continuation of each nonzero Fourier coefficient

Green's identity:
- `Eisenstein/GreenIdentity.lean:77-179` and `Eisenstein/CuspIntegration.lean:77-455`, ending in `actual_eisenstein_cylinder_green`:
  `(s*(s-2)-star (t*(t-2)))*cuspBoxIntegral a b (eisensteinHermitianDensity s t) = ∫∫ greenBoundaryFlux … b - … a` for `2<s.re`, `2<t.re`

**What it attempted.** These are the two textbook routes to continuing the Eisenstein series to s = 4/3 and taking the residue θ:
- spectral theory of the self-adjoint Laplacian (meromorphic resolvent);
- the Green/Maass–Selberg route.

**How far it got.** The spectral pieces are complete but unused. The Green route stops at a single cusp box, inside the convergent range Re s, t > 2. It has no truncated fundamental domain and no Maass–Selberg relation.

**Why it was dropped (evidence).** The live construction is a Colin de Verdière-style seed (`heightPoincareField (cuspSeedProfile …)`, `SeedDefects.lean:87`) combined with a Fredholm pencil on the energy space:
- `FredholmAlternative.lean:139-262` `kernel_theta_spectral_fredholm`
- the explicit Laurent split `kernelSpectralResolvent z F = ((17/9)/(8/9-z))•kernelThetaMassProjection F + kernelRegularL2Resolvent z F` (`MeromorphicResolvent.lean:588-595`)
- an analytic numerator (`:597-663`)

The residue is then just a `Tendsto` limit (`kernelEisensteinResidueVector_limit`). Nothing downstream consumes `MeromorphicAt`, `IsSelfAdjoint`, finite-dimensionality of the theta space, or non-vanishing of θ; non-vanishing of the *averaged* theta comes from its cusp average instead (`CuspAverage.lean:378`). The Fredholm route also avoids flux bookkeeping on a three-dimensional fundamental domain, which is where the Green route stalled.

**Assessment.** The pruned code isolates a minimal analytic interface: the residue exists, lies in the 8/9-eigengraph, and has an explicit cusp average. That interface is reusable and a candidate for an upstream library lemma, an abstract "Fredholm pencil gives a simple-pole residue" result. The Green route is not worth reviving.

---

## Cluster 3. Routes to the cube relation τ(hp³) = |p|·τ(h)

**Main declarations (dead).**

The Euler-factor recurrence:
- `Eisenstein/BesselDerivatives.lean:398-411` `unramifiedCubicGaussSeries_cube_recurrence`
- `:534-551`: the same recurrence at the residue
- `:570-584` `unramifiedGaussResidue_cube_closed`:
  `(1-q)·r_k = (r₁ - q r₀) + (r₀ - r₁)·q^k`, where q = N(p)⁻¹ and r_k = τ(hp^{3k})

Hecke-correspondence infrastructure:
- `Eisenstein/InverseCorrespondence.lean:589-648` (`cubeCoverMatrix`, `cubeInverseL2Correspondence`)
- `Eisenstein/CorrespondenceAdjoint.lean:127-138` `kernelCorrespondenceL2_adjoint`

A constant-term variant:
- `Eisenstein/CubeEigenrelation.lean:152-221` (`cubeSourceDifference_cusp_period`, `…_cusp_average_zero`)

Other:
- `PrimeRecurrence.lean:335-387` (`exists_nonzero_sourceResidualFourierCoefficient`)

**What it attempted.** These prove the theta coefficient relation used in paper2's "transfer"/cube-factor descent and in paper78's coefficient support bound `|d(μ)| ≤ 27·3^{k/6}|b|`.

**How far it got.** All of it is proved.

**Why it was dropped.** The recurrence solves to A + B·q^k. That leaves one free constant, so it cannot by itself yield r₁ = r₀. *Conjecture:* this is why it was abandoned. *Evidence:* the dead closed form displays exactly this two-parameter family.

The live route is spectral:
1. D = cubeAverage θ − N(p)⁻¹θ is an 8/9-eigenfunction with cubic decay at every cusp (`ResidualOrthogonality.lean:293-427`).
2. Hence D is orthogonal to all levelTwo-translates of the residue (`:499` `cubicEisensteinResidue_source_translate_orthogonal`).
3. So D = 0 (`PrimeRecurrence.lean:291` `cubeSourceDifference_eq_zero`).
4. This feeds `CorrespondenceEnergy.lean:584` `sourceResidualFourierCoefficient_cube_of_eigenrelation`.

This needs neither adjoints nor multiplicity one, which is why the adjoint and finite-dimensionality lemmas died.

**Assessment.** The "decay + orthogonality ⇒ zero" lemma is the most reusable idea in the directory (see Worth exploring, item 1).

---

## Cluster 4. Completed-reflection identity in paper form, and cusp-coefficient variants

**Main declarations (dead).**

The paper form of the identity:
- `Eisenstein/SourceCuspKernels.lean:479-539` `completedT_smoothed_reflection`. It is conditional on a single-cusp hypothesis,
  `hreflect : ∀v>0, completedBesselProfile Ψ q v = A*((v:ℂ)^(-2:ℂ)*radialBesselProfile r a d ((Q*v)⁻¹))`,
  and it outputs the paper's `paperKernel (Vstar W) (27*q^2*d^2*X*r i/Q^2)` series, which matches paper78's K = (2π)⁴/27 and V♯.
- `Eisenstein/FixedCuspArrays.lean:162-185` `completedT_finite_smoothed_full_index`, the fully expanded variant
- `Eisenstein/RamifiedValuation.lean:542-577` `thetaTwistedCuspProfile_reflection`

Cusp-coefficient variants in `Cusp/`:
- `Cusp/ResiduesModThree.lean:331-555` `unitLowerEisenstein_fourier`, the Fourier expansion at the unit cusps (1 0; u 1)
- `Eisenstein/CuspDilation.lean:145-285`
- `Cusp/WeylRepresentatives.lean:140-202` `A3_principal/ramified/unramified_character`, which factor the Kubota character as (c₀/a)(a/r) for paper2's Appendix A.3 cases
- `Cusp/MatrixPhases.lean:193-281` `ShortDraftGaussJacobi.cubic_gauss_cube`, `gaussSum χ ψ ^ 3 = q * jacobiSum χ χ`
- E-primary square classes and Gaussian phase tables (`MatrixPhases.lean:154-434, 505-701`)

**What it attempted.** A literal formalization of paper78 eq. `(eq:reflection)` and the paper2 A.3 case analysis.

**How far it got.** It is complete, but the single-cusp hypothesis in `hreflect` is never met. The twisted profile is a *finite sum over h mod c* of cusp terms (`thetaTwistedCuspProfile_reflection`), not one radial profile.

**Why it was dropped (evidence).** The live replacements are `FullFrequencyExpansion.lean:454-543` `completedT_finite_smoothed_reflection` and `completedT_canonical_smoothed_reflection`. These sum over h mod c and keep the kernel abstract as `(d h).smoothedKernel W X`. The unit-cusp expansion was generalized by the live `InverseCorrespondence.lean:515` `unitCusp_eisenstein_fourier_base`, which covers diag(u)·lower(t). The factored character formulas were replaced by the live "literal" versions (`A3_ramified_literal_character`, `A3_unramified_literal_character`). The `ShortDraft*` namespaces in `Cusp/` look like remnants of an earlier draft.

**Assessment.** This is engineering supersession. It is not promising for further compute.

---

## Cluster 5. Reflection/: reflected-branch energy for 7/8 (42 modules deleted)

This is paper78 §`sec:marked-energy` (paper78.tex:7249), in particular `lem:hybrid-energy` (7510) and `lem:reflected-energy` (7747). It is not theta reflection.

**Main declarations (dead).**

(a) Quantifier-order clones. For example, `Reflection/CanonicalChoiceEnergy.lean:15-161` `canonical_original_choice_energy` and its live twin `CanonicalChoiceEnergyUniformDegree.lean:13-163` have identical bodies except that `∃ (degree : ℕ)` moves in front of `∀ {Nlevel} a c₀ mode … B`. Seventeen such pairs (1.8k lines) plus thirteen `_type_uniform` pairs (0.8k lines) follow the same pattern: CanonicalDyads, CanonicalEnergy, CanonicalFrozenEnergy, CompletedCanonicalEnergy, OriginalFull, GlobalBudget, MemberBudget, and others.

(b) A Cauchy–Schwarz loss of 3^{|φ|}. `Reflection/SumEnergy.lean:13-67` `summed_actual_weighted_reflected_energy` and `PhysicalEnergy.lean:14-132` `finite_physical_source_energy` sum over *all* branch assignments `e : φ→Fin 3` (one cusp per prime), and they pay `(Fintype.card (φ→Fin 3):ℝ) * ∑ e : φ→Fin 3, …`.

(c) Generalization clones:
- `Matrices.lean:12-73` versus the live `Family.lean:12-73` `…_dependent` (dependent index types)
- `Partition.lean:9-66` versus `Cross.lean` `ordered_cross_partition_on`
- `CanonicalPool`/`CanonicalWidth` with `choice : pool→Fin 6` versus the live `OptionalWidth.lean:12` with arbitrary `A : Finset pool`
- `MarkedSource`/`MarkedDualSource` versus `…_of_controlled`

**How far it got.** Everything is complete.

**Why it was dropped (evidence).** For (b), the live `FilteredEnergy`/`SurvivingEnergy` differ *only* in replacing `univ` with `E`. In the physical chain, `E = survivingFrozenBranches F jF (actualCuspColumn …)` (`SurvivingEnergy.lean:14-134`), so the loss becomes `branches.card`. Even the `_type_uniform` versions of the summed chain are dead (`UniformSieve.lean:143-286`, `UniformPhysical.lean:18-136`).

*Conjecture:* 3^{|φ|} is not Z^{o(1)} once the pool of marked primes can grow, so the margin `cstar/8` cannot absorb it. For (a) and (c), downstream callers needed constants uniform in the level and the cusp data, and the agents cloned theorems instead of refactoring them.

**Assessment.** The only mathematical lesson is to restrict to the branches that survive before applying Cauchy–Schwarz. The rest is churn.

---

## Cluster 6. Hecke/: detector and prime-extraction variants

These findings come from a delegated sub-audit; I spot-checked two diffs.

- **Witness-first detector** (`Hecke/DetectorSimultaneous.lean:14-130`, `DetectorAdaptiveCutoff`, `DetectorDynamicSelection`, `DetectorBatchCount`). Verified: `DetectorSimultaneous` is the live `DetectorSimultaneousSupport` minus the conjunct `tstar-r-ε≤m`. The dead proof even computes `hsupport.2` and discards it. That conjunct is the support lower bound that the detector note in paper78 §`sec:detector` needs. Per-row adaptive cutoffs became unnecessary once the witness theorem held uniformly in t; the live code uses amplitude classes with one cutoff per class instead.
- **Per-row ramification certificates** (`DetectorPrimeFamily.lean:92-118` `structure Ramification`; `DetectorPhysicalBins.lean:15-50` `physical_bins` takes a `cert` that nothing constructs). These were superseded by the live contrapositive `RowNonprincipal.principal_row_multiplicity` (`χ.residue = 1 → 6 ∣ count P`). The live conductor bound is also better.
- **23/24 for Dirichlet L via prime extraction** (`LocalEulerIdentities.lean:500-751`, `PrimeExtraction.lean:38-548`, `FiniteEulerFactors.lean:169-653`). The final results are conditional on two inputs that are never discharged:
  - the Mertens-type hypothesis `hprefix`, which is close to assuming the zero-free region;
  - `structure FixedTestPrimeMS` (`FiniteEulerFactors.lean:609-653`). It encodes paper2's prime extraction and needs `hcount : (2m)^{11/60-ε} ≤ selected.card`, a lower bound on the number of primes, and nothing builds it.

  The live `QRH.saving_of_meanSquare` averages over all ideals b with N(b) ≤ D^{11/60} using totient weights. This removes the need for a prime count.
- **Hand-built pole subtraction for −L′/L** (`LogDerivative.lean:19-81`, `RayQuotient`, `RayPrimes.lean:38-58`) duplicates the generic `PNT.RayAsymptotic` machinery.
- **`Hecke/Theta.lean`** is a weight-1 Gaussian lattice theta used for Hecke's continuation. It is *not* Kubota's cubic theta. The dead `completed_functional_equation` (221-226) is a three-line corollary of Mathlib.

**Assessment.** Swapping prime extraction for totient averaging is the one real simplification here, and it already sits in the live 23/24 proof. Nothing else in Hecke looks promising.

---

## Worth exploring further

1. **Prove that Kubota's θ equals its SL₂(ℤ) average, closing paper2's unproved assertion.** The pieces already exist:
   - The constant terms already agree (`kernelSourceProjection_residue_cusp_average`, live).
   - The coefficient at a non-horizontal h vanishes for θ (live `kernelCuspFourier_residue_eq_zero_of_not_horizontal`) and for the average.
   - The dead `kernelSourceProjection_residue_fourier_mask` gives the averaged coefficient as ¼(θ_h + 3θ^{opp}_h).

   So θ = Pθ reduces to θ_h = θ^{opp}_h for horizontal h. By the live `kernelCuspFourier_principal_residue` and `translatedCuspFourier_unramified_residue`, this means proving

   `((9*√3/2):ℂ) * principalArithmeticResidue h = unramifiedGaussResidue (9*h)` for horizontal h.

   That is a finite identity between two explicit cubic Gauss-sum residues. Equality of the functions then follows from the two Fourier–Bessel expansions (`cubicResidualFunction_fourier_full` and `cubicSourceResidualFunction_eq_bessel`). An alternative is to reuse the "decay at all cusps + orthogonality ⇒ 0" template from `cubeSourceDifference_eq_zero`.

   *Caveat (unverified):* invariance is asserted in paper2 and in Dunn–Radziwiłł, but nothing in the repo checks it. If the identity fails, that would itself be a finding about paper2.

2. **Extract the minimal residue interface.** The dead code shows that meromorphy, self-adjointness, non-vanishing of the pole, and dim(theta space) < ∞ are all unused. A next step is to package `FredholmAlternative` with the Laurent split in `MeromorphicResolvent.lean:588-663` as one generic lemma, then re-run the reachability audit to measure how much of the about 3.2k live Eisenstein declarations remains. This aims at a shorter proof, not a better exponent.

3. **Push the "surviving branches" idea further in the 7/8 energy bound.** The step from all 3^{|φ|} branches to `survivingFrozenBranches` was decisive enough that both earlier chains died. Next, check whether the live bound in `SurvivingEnergy`/`UniformPhysical*` still has slack in `branches.card`, for example whether the count can be bounded by a power of 2 rather than 3. If so, see whether that margin propagates to the `cstar/8` saving in `CanonicalChoiceEnergyUniformDegree`. *Speculative.*

## Lessons / dead ends to avoid

- Do not attempt reflection with Γ(3)-automorphy alone. The Weyl-case cusps of paper2 (`(c,λ)=1`) need SL₂(ℤ). Average first, or prove invariance first.
- Do not build global spectral theory. The proof needs only a local Laurent split at s = 4/3, not `MeromorphicAt`, an unbounded `IsSelfAdjoint`, or Green/Maass–Selberg on cusp boxes in H³.
- Euler-factor recurrences do not determine τ(hp³). They leave one free constant, so use the spectral argument "difference is a decaying eigenfunction ⇒ 0".
- Do not apply Cauchy–Schwarz over all 3^{|φ|} cusp branches. Restrict to the surviving branches first.
- Fix the needed uniformity (polynomial degree in θ before the level and cusp data, universe and type uniformity, dependent index types) before proving long energy chains. About half of Reflection's dead lines are clones that differ only in quantifier order.
- A hypothesis-carrying structure such as `FixedTestPrimeMS.hcount`, `Ramification`, or a single-cusp `hreflect` is a warning sign. In every case here, the live route won by eliminating the hypothesis, not by discharging it.
