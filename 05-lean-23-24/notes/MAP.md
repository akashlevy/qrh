# MAP: paper2 ("The Quasi-Riemann Hypothesis") against lean-base

*8 Oct 2026; pointer updated 9 Oct.* Lean paths below are relative to `lean-base/OAI/NumberTheory/DirichletL/` unless they start with `QRH`.

This file is the M1 dependency analysis of paper2 against lean-base. The current deliverables are elsewhere: the human proof is `paper/qrh.pdf`, its Lean counterpart is `QRH/Spine/*.lean`, and paper ↔ Lean names are in Appendix C of the paper. See `REPORT.md`.

## 1. Headline finding

lean-base already contains an unconditional, fully proved instance of paper2's argument at the fixed value ϑ = 1/10:

```lean
theorem OAI.ShortDraft.riemannZeta_ne_zero_of_re_gt (s : ℂ)
    (hs : (23/24 : ℝ) < s.re) (hs1 : s ≠ 1) : riemannZeta s ≠ 0   -- RowCompletion/ZeroFreeRegion.lean:210
```

- 23/24 = 11/12 + 5ϑ/12 at ϑ = 1/10.
- `#print axioms` gives only `propext`, `Classical.choice` and `Quot.sound`. The theorem was built locally with Lean 4.34.1.
- Its import closure has 286 modules (190k lines). It does not touch any module of the 7/8 probe, detector or moment chain: the 7/8 proof (2924 modules) is built *on top of* it.

After the user's decision (8 Oct: "shortest possible QRH proof, any constant is fine"), the project target changed:
- **new target:** a minimal standalone proof of ζ(s) ≠ 0 for Re s > 23/24, by tree-shaking and shrinking this closure;
- **not pursued:** 11/12.

See REPORT.md for that work. This file records the M1 dependency analysis.

## 2. Import boundary

- **Denied (`audit/policy.json`):**
  - `Nonvanishing` (the 7/8 end theorem) and anything that depends on it, or on `Detector.FinalAssemblyUnconditional`;
  - the 7/8 argument chain: `Detector/`, `Moments/`, `Descent/`, `Energy/`, `Inversion/`, `Dictionary/`, `Hecke.Detector*`, `Parameters*`, `GammaZeroDetector`, `GammaDetectorExpansion`, `Principal{SignalComparison,SlotEstimate,MellinGrowth,MellinResidues}`, `Supremum`, `Endpoint`.
- **Allowed:** the 23/24 closure.
  - It contains no denied module (`scripts/audit_imports.py` passes).
  - 23/24 > 11/12, so importing it would not trivialize an 11/12 proof.
  - Under the new target it *is* the proof being shrunk.
- `audit/modules.tsv` classifies all 2924 modules:
  - lines;
  - closure size;
  - in-23/24-closure, or depends-on-23/24;
  - which 7/8-chain directories appear in its closure.

## 3. paper2 step → Lean (all inside the 23/24 closure)

| paper2 | Lean (main declaration) | Fit |
|---|---|---|
| Setup: O, primary, χ_n = (·/n)₆ | `ActualEisensteinCubic.O = 𝓞 (CyclotomicField 3 ℚ)` (GaussSum/FiniteFourier.lean:487); `sexticChar` (:801); `primaryGenerator` (Arithmetic/IdealMobius.lean:232) | exact |
| Prop thm:ms (mean square) | `InitialMeanSquare.outside_mean_square_of_theta_models` (MeanSquare/ThetaDensity.lean:269) | **holds for every σ > 1/20**, not every ϑ > 0 |
| Extraction (eq:prime-extract) | `SixthPowerAverage.outside_totient_moebius_bound` (Mellin/SixthPowerAverage.lean:317) | **integer rows b⁶** with a bump average, no prime counting; σ = 1/10 is hard-coded here only |
| Mellin step, base change, ζ | `dirichlet_open_strip_of_theta_models` (RowCompletion/ThetaNonvanishing.lean:105); `pair_inverse_euler_identity` (Mellin/LocalNonvanishing.lean:172) | uses L(χ)·L(χχ₋₃) from Mathlib; no Hecke L-function is built |
| Lemma arithmetic (Gauss–Jacobi, convert, crt, recip) | `FiniteGaussPhase.*` (PrimeRows/ActiveGaussFactorization.lean:291–798); `SquarePhaseFactorization`; `ConjugateLevelAction` | exact, but over finite prime pools rather than squarefree n |
| Lemma poisson | `canonical_masked_radial_poisson_collected` (PrimeRows/IdealGenerators.lean:251) | adapter |
| eq:energy, Lemma remove-exclusions, Prop poisson-reduction | `outside_mean_square_source(_with_errors)` (MeanSquare/PoissonSource.lean:80, DiagonalDensity.lean:323) | different architecture: exclusion mask m carried in the energy; physical-space majorant plus tails instead of band-limited Φ |
| eq:T (completed sums) | `CompletedGauss.completedT` (Arithmetic/IdealMobius.lean:471) | exact |
| Prop canonical (finite descent) | `outsideLogEnergy_finite_rank` (RowCompletion/DensityTransfer.lean:255) | same skeleton; gap 1/20, cube cutoff Z^{1/1000}, depth ≤ 3000 all hard-coded |
| Lemma cube-reduction | `completed_cube_inverse` (GaussSum/CompletedExpansion.lean:218); `canonicalLog_actual_step` | exact (inversion), partial (split) |
| Prop transfer, Lemmas first/second-transfer, arithmetic-poisson | `reopenedCanonicalRow_closed_step` (CubicSieve/ClosedRecursion.lean:414); `first_passage_full_uniform_transfer` (CubicSieve/SquarefreeTransfer.lean:442; the earlier `first_passage_transfer` is dead); `masked_pair_second_poisson` (CubicSieve/SecondPoisson.lean:83) | partial: same two-Poisson idea, different bookkeeping |
| Prop R | `CompletedHeight.uniform_polynomial_mean_square_of_exact_reflection` (MeanSquare/UniformPolynomialHeight.lean:309; the earlier `completed_mean_square_of_exact_reflection` is dead) + `actual_fixedShellModels` (RowCompletion/ZeroFreeRegion.lean:174) | adapter; unconditional, every ε, not tied to ϑ |
| θ, Patterson coefficients, cusps, automorphy | `cubicSourceResidualFunction` (Eisenstein/ResidualFourierExpansion.lean:499), built as the residue of the metaplectic Eisenstein series with a Fredholm resolvent; `infinityCoefficient_squarefree_cube` (Eisenstein/SquarefreeCoefficients.lean:134); `three_cusp_decomposition` (Cusp/ResiduesModThree.lean:251) | adapter: the SL₂(ℤ)-averaged θ_src replaces θ |
| Prop lem:reflection, reflection-uniformity, theta-bounds | `completedT_eq_reflectedFiber_idealValue` (RowCompletion/ZeroFreeRegion.lean:139); `FixedReflectionRay` (Eisenstein/FixedFrequencyExpansion.lean:408); `three_cusp_coefficients_squarefree_cube_bound` (Eisenstein/RamifiedValuation.lean:329) | adapter |
| Lemma quadratic (Goldmakher–Louvel) | `CanonicalQuadraticSieve.sieveNorm_sharp` (QuadraticSieve/LogarithmicLoss.lean:474) | proved directly by Heath-Brown's α ↦ 2 − 1/α iteration; GL is never used |
| Appendix B (smooth separation) | log-profile separation (Mellin/LogProfiles.lean, Fourier/PositiveLogProfiles.lean) | different technique; Lemma smooth-mean-square has no general analogue |

### What 11/12 would have needed (for the record)
- **15/16 (inferred, not compiled):** generalize the extraction step only, using Y = Z^{(1+σ)/6} (about 100–150 lines). This is possible because thm:ms already holds for σ > 1/20.
- **11/12:**
  - make the descent constants into parameters: 66 live declarations in 20 files contain the literals 1/20, 1/40, 1/1000, 3000, e^{9000} and similar;
  - or write paper2's descent fresh, estimated at 4.5–7.5k new lines, reusing Prop R and the leaf identities.

## 4. Gaps and issues found in paper2

No mathematical errors were found. The arithmetic identities of Lemma arithmetic were checked numerically on 74 split primes and their products. Inert primes were not tested.

### Definite (minor)
1. **Line 747.** For principal ν, L_K^S has a pole at s = 1, so "the identity holds on Re s > 11/12" should exclude s = 1. This is harmless because ϱ ≠ 1.
2. **Lines 749–750.** "The omitted Euler factors are nonzero" is not needed for the contradiction.
3. **Outline, eq:intro-quadratic-twist (lines 408–411).** The twist is quadratic only when j_p = 1. At p | (h, f), or when v_p(h) ≥ 2, B_{p,j} is not quadratic. The body (lines 1975–2075, 2189–2197) handles this correctly.
4. **Line 2851.** R(c,c) = χ_c(−1) is used without proof. eq:quotient also needs G(1) = 1 and |G| = 1 on every ray class.
5. **Line 805.** γ_j for j < 0 is undefined at non-units. Read χ_n^{-1} as χ_n^5.
6. **eq:g, composite case.** The text does not mention the CRT step it needs.
7. **Lemma theta-bounds.** The support and the bound for d_± are claimed for all m ≥ −4 but proved only for m = −4. The claim is true vacuously, because d_± is supported at v_λ = −4 (Lean proves this, Eisenstein/RamifiedValuation.lean:181).
8. **Cube-reduction (lines 1349–1355).** When H_c < 1, b = 1 lies in the supremum, so the lemma is vacuous there. The argument still works through the transfer, but a formalization must not assume b ≠ 1.
9. **Notation clashes.** j is used for the induction index, the derivative order and a divisor. 𝓗 is used for both the Poisson scale and the dual scale.

### Worth checking
10. **Kubota character convention (line 3153).** It must be κ = (c/a)₃, as in DR (5.4). With (c/d)₃, B_{p,1} would be sextic rather than quadratic, and Lemma quadratic would not apply. Lean uses (c/a)₃, and the Lean proof machine-checks that this is consistent for its θ.
11. **Sketched or cited rather than proved:**
    - the multiplier formula eq:ray-multiplier (line 3175);
    - eq:quadratic-gaussian (r_η, line 2799);
    - SL₂(ℤ)-invariance of θ (line 3091; Lean uses an averaged θ_src instead);
    - the Phragmén–Lindelöf contour shift (Lean uses a Gaussian-regularized shift instead);
    - the hypotheses of Goldmakher–Louvel (Lean proves its own quadratic sieve).
12. **Ray class modulus.** It must contain a power of λ: G(ωn) ≠ G(n), because (ω/−2)₃ = ω. The phrase "sufficiently large" covers this only implicitly.
13. **Band-limited Φ (lines 1025–1028, 1470–1474).** It is asserted, not constructed: squaring doubles the support of Φ̂.
14. **‖W‖_{C^J} uniformity.** Prop R states an explicit ‖W‖_{C^J} dependence. In Lean the constants depend on W, with uniformity only along N^{it} twists.
15. **Appendix B.** ‖·‖_{C^q} is never defined. Lemma smooth-mean-square needs complex-valued test functions U.
