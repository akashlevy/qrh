# 10 — Toward θ < 41/47: working notes (10 Oct 2026, afternoon session)

Goal: θ < 41/47 = 0.8723404 for the Hecke family over Q(√−3), formalized in Lean.
**Status: no proof yet.** The goal needs a new analytic input; this file records what the
weighted-numerator model says that input must be. See `BRIEF.md` for the precise model and sources.

## 1. State of the art (checked 10 Oct ~14:00 UTC)
- Tracker main unchanged (best framework-verified ≈ 0.874957). Open PRs: #6 (ours, 10499/12000),
  #7 (Gao, retuning), #8 (Collet, 29/33 non-record), #9 (Collet, 0.87488664704 "quartic boundary"),
  #10 (website). Nothing public below 0.8748; no 41/47 claim anywhere (web search).
- PR #9 = our weighted numerator + lem:plain extended to κ ∈ [7/10, 3/4] ("variable κ") + the
  completed zero-slot moment. A re-tuning; the structure is unchanged.

## 2. The weighted model (new: the 09 session only modelled the pointwise-numerator high side)
Scripts in `model/`: `wmodel.py` (reproduces the 08 certificate: margin −1.3·10⁻⁴ at ℓ = 0.167),
`gmodel.py` (count as an LP over witness splits, with knobs), `deficit.py`, `scan_knobs.py`, `wall.py`,
`line.py`, `ideal.py`, `sigma.py`.

1. **Low side:** B = 5/6 + M/12 − ℓ/6 (general), = 11/12 − ℓ/4 on the wall M + ℓ = 1.
   41/47 needs ℓ ≥ 25/141 ≈ 0.1773 on the wall. Along the line B = 41/47 (M = 0.468 + 2ℓ) the
   wall point is optimal for the high side (`line.py`).
2. **Deficit at ℓ = 0.1773:** worst min(E_old, E_wt) ≈ +0.011 (best b ≈ 0.075). It lives entirely at
   large prime amplitude x ≥ 0.2, δ ∈ [0.2, 0.65] (`deficit.py`); every x = 0 bin closes.
3. **The numerator is the bottleneck, not the counts** (`ideal.py`): with the paper's counts but a
   numerator of size O(U^ε) on the bad rows, the high side would close up to ℓ ≈ 0.2 (θ → 13/15).
   At δ = 0.47, x = 1/2 the allowed average numerator exponent is 0.0857 vs 0.0995 (Hölder/4th moment)
   and 0.1037 (pointwise δn₀/2).
4. **Single inputs that would reach 41/47** (`gmodel.py`, `scan_knobs.py`, `wall.py`):

   | input | worst exponent at ℓ = 0.1773 |
   |---|---|
   | none (paper + weighted) | +0.0117 |
   | numerator sixth moment Σ|S_u(n₀)|⁶ ≪ U^{1+0.10} | −0.012 |
   | DIAG (marked r+z ≤ 1, plain n₁+n₂+z ≤ 1) | −0.008 |
   | marked capacity 0.8(1−r) + plain capacity 0.4(1−2m) | −0.0001 |
   | marked capacity 0.75(1−r) + long-witness slope α ≤ 0.75 | −0.0001 |
   | lossless wall W(0.045) | ≈ −0.001 (W(0.038): +0.0014) |
   | DH counts | −0.028 |
   | marked capacity (1−r) alone / plain capacity 1 alone | +0.003 / +0.0044 |
   | α = 0.7 alone | +0.0057 |

   Each is a "beyond the edge" or coefficient-sensitive statement listed as open in docs/OBSTRUCTIONS.md.

## 3. Quick checks that do NOT help (this session)
- **Numerator × plain witness (lem:plain on S_u·S_m):** exactly log-linear; when the count comes from
  the plain moment at the same m it reproduces the 4th-moment Hölder bound to the digit.
- **Copies for the numerator:** the plain mean square is already diagonal at row range U (reflection),
  so amplification over u·a⁶ gives nothing (copies are sparse: A copies in range UA⁶).
- **Sextic large sieve** is diagonal only up to length 1/2 ((KD)^{2/3} term), weaker than lem:marked.
- **β*-numerator:** zero-free half-plane Re > β* gives |L(w,χ_u)| ≪ U^ε only on Re w ≥ β* + ε; the
  resulting numerator exponent −(1−β*)l_y is far worse than Hölder.
- **Halász–Montgomery with the a priori zero-free region:** off-diagonal prime/Möbius sums are only
  ≪ P^{β*}, and β* ≈ 0.87 > δ, so no count gain (needs GRH-strength P^{1/2}).
- **Separate numerator/denominator bins:** the numerator character χ_u and the denominator ηχ̄_u lie in
  the same finite family 𝒳_u; the adversary puts both zeros at a, recovering the current bound. A
  joint count via products of witnesses for the two twists gives coefficients μ(k)·Σ_{d|k}η(d), whose
  divisor variable d enters the theta reflection like a mark (conductor d²), so it is not free.
- **Marks in long-witness amplification:** extend the effective length with slope 5/6 − δ > 0: worse.

## 4. Research agents
Five agents (numerator, marked-cap, amplitude, redesign, lit2) were launched at ~14:15 UTC and all
terminated immediately on an API spend limit; they produced no results.

## 5. Further checks (after the agents failed)
- **Numerator contour Re w = c with moment convexity.** Gabriel/Hölder interpolation between the
  4th moment on Re w = 1/2 and the bin's zero-free pointwise bound on Re w = a gives the exponent
  (c−1)l_y + d(R + θ_c(1−R)/4), θ_c = (a−c)/(a−1/2): affine in c, so an endpoint is optimal
  (c = 1/2, i.e. the weighted bound, at the binding bins: slope l_y − d(1−R)/(2δ) ≈ +0.10 > 0).
  Pointwise convexity alone gives the P78 choice c = 1−a (σ = δn₀/2). No interior gain.
- **Restricting the low-side Cauchy–Schwarz to the bad rows** (rows factor as S_u·M_u): bounding
  Σ_{u∈C}|S_u M_u| needs Σ_C|M_u|², i.e. a large-values bound for the length-Z Möbius sums with
  q_u < Z (below the diagonal, M_κ type). Same obstruction.
- **Several Y-scales / dyadic splits of the numerator:** the numerator is single-scale after the
  compact-support weight; multi-scale probes are sums of probes (additive bounds, no gain).
- **Heath-Brown splitting of the long inverse witness:** already covered by the 09 result that
  detector redesigns only give the adversary more freedom; products of several Möbius pieces have
  no moment lemma.

**Conclusion so far.** In the weighted model every route to 41/47 found here requires one of the
open inputs in §2.4. The cheapest (closest to existing technology) look like a lossless energy
wall W(0.045) or a marked capacity ≥ 0.8(1−r) with plain capacity ≥ 0.4(1−2m); both are versions of
the beyond-the-edge problems the 09 session documented. Nothing here is a proof, so there is nothing
yet to formalize; the Lean environment build (Mathlib + OAI + WeightedQRH) was started in
/home/user/leanbuild for when a proof exists.

## 6. Required counts and the Gram lever (`model/rneed.py`, `model/gram.py`)
- **Counts needed at x = 1/2** (ℓ = 0.1773, b = 0.075, weighted branch): R(0.3) 0.7468 → 0.7359,
  R(0.4) 0.6561 → 0.6391, R(0.47) 0.5897 → 0.5713, R(0.55) 0.5100 → 0.4938, R(0.6) 0.4575 → 0.4454.
  A drop of 0.011–0.018, still well above DH (1−δ). At x = 0.3 the gap is ≤ 0.006.
- **Gram exceptional factor g** (Gram term P_a^{g/6}, low side B = 5/6 + M/12 − ℓ/6 − b(1−g)/12):
  g = 0.9 → 0.8741, 0.85 → 0.8734, 0.8 → 0.8728, **0.75 → 0.87221 < 41/47** (coarse grid; b ≈ 0.19).
  The 04/07 analyses show the Gram term's Mellin coefficient Ŵ₁(1)Ŵ₀(5/6) equals the signal's, so
  reducing it is equivalent to the (V′) pairing problem, which is signal-circular.

## 7. Literature agent `lit2` (finished; `agents/lit2/LIT.md`)
- Tracker: **our PR #6 (10499/12000) was merged at 14:08 UTC; the catalogue on main lists
  `akashlevy-20261009-weighted-numerator` as framework-verified, the best verified record** (checked).
  No bound below 0.8748 claimed anywhere; no 41/47 claim.
- Calibrations at ℓ = 0.1773: a family zero-density N(σ) ≪ U^{A(1−σ)} closes 41/47 only with A ≲ 2.125 on
  σ ∈ [0.59, 0.80] — better than the best known exponent for ζ itself (2.308 at σ = 0.7); P78's counts
  already beat Ingham and Guth–Maynard at the binding bin. Numerator sixth moment needs excess λ ≲ 0.09;
  available routes give 0.34 (perfect sextic sieve), 0.45 (de Faveri), 0.37 (4th moment × Burgess).
- No published or preprint tool supplies any input of the BRIEF table (details and citations in LIT.md).

## 8. Agent `redesign` (finished; `agents/redesign/NOTES.md`) — negative
- PROVED (relative to P78's lemmas): high identity for ANY s-weight λ; the numerator N_u[λ] is linear
  in λ with a multiplicative, w-independent local factor π_u. Low side: any |λ| ≤ 1 costs (Q/Y)(1+P_a);
  the Gram bound (1 + P_a^{1/6}) holds only for fixed ray-character combinations (b > 0 otherwise costs
  up to Z^{5b/12}); at b = 0 the weight is free but the high side is worse (+0.0132).
- (a) mollified/Selberg/Λ/square/cube/AP/sector weights: all worse (+0.0145 … +0.0228) or break the
  principal-row structure. (b) Y = Y₁Y₂ product numerators: only a second moment for products (+0.094);
  even granting diagonal moments of every real order the envelope is θ ≥ 0.87273. (c) target-dependent
  annihilation: dimension count fails at the worst bin (#C ≈ Z^{0.486} vs Y ≈ Z^{0.441}); worst stays
  +0.0118. (d) coupled probes reduce to multi-scale weights; per-bin oracle geometry still +0.0113.
- New sufficient condition (open): a restricted large sieve Σ_{u∈C}|S̃_u|² ≪ U^ε(#C + U^{l_y/d}) on the
  bad set closes 41/47 (θ ≈ 0.8663); minimal need at the worst bin U^{0.770} vs proved U^{0.805}.

## 9. Agent `numerator` (finished; `agents/numerator/NOTES.md`) — negative, with a barrier
- PROVED (exponent level): for every witness class there is a flat configuration satisfying ALL
  available inputs at once (lem:plain at any lengths/capacities, the numerator × plain-witness moment,
  lem:marked, amplification, the bin's pointwise bounds, the sextic large sieve, and every Hölder or
  level-set combination of them) with the numerator exactly at σ_wt = min((1−R)/4 − q z_N/2, δn₀/2)
  (`barrier_check.py`). So the weighted numerator bound is optimal for this black-box lemma set.
- **Correction to BRIEF §2.4:** the N6 knob assumed excess η = 0. A sixth moment ≪ U^{1+η} closes 41/47
  only for η ≲ 0.088 (b = 0, n₀ = 0.463); η = 0.10 leaves +0.0012. Best provable excess ≈ 0.34; the
  general-coefficient truth is 3n₀ − 1 ≈ 0.39, so the d₃ structure (sixth-moment Lindelöf on average)
  would be required. Mixed moments: |S|⁴|S_m|² needs η′ ≤ 0.086; numerator prime capacity rate ≈ 1.4.
- Halász–Montgomery over C, even with Lindelöf for every pair character, only gives σ ≤ n₀/4 ≈ 0.11.

## 10. Agent `marked-cap` (finished; `agents/marked-cap/NOTES.md`) — negative, structural reduction
- Fine grid (δ step 0.01): saturating lem:marked (z = L for all r < 1) alone leaves +0.0019; the BRIEF pair
  "marked 0.8(1−r) + plain 0.4(1−2m)" is a coarse-grid artefact (+0.0003). Closing pairs: plain 0.45(1−2m)
  with saturated marks (−0.0002); 0.75(1−r) with α = 0.75 (−0.0009); saturated marks with α = 0.78.
- A shared lossless wall saving Z^w feeds lem:marked only as r + 2z ≤ 1 + w/2; 41/47 needs w ≈ 0.037.
- PROVED (Lemma D): for r + z ≤ m − c the tuple-diagonal part of lem:marked's initial canonical energy is
  ≪ Z^{F+ε} with no mark cost (factor n = τn′, drop |A_*(τ)χ_τ(k)| ≤ 1, absorb the tuple into a fourth-power
  label, apply lem:canonical-moment with z₀ = 0). Hence lem:marked at (r, z) reduces to bounding the
  cross-tuple terms τ ≠ τ′ for a slot part of length r + 2z − 1. P78's range r + 2z ≤ m is exactly
  "no cross-mark cancellation". Going further needs ≈ 60% of square-root cancellation in the marks' Gram
  matrix — an (H2)-type cross-mark decorrelation for twisted cubic theta, true numerically, no tool known.
- Every suggested route (μ-absorption of slot primes, Möbius/Mellin decompositions, Vaughan/Heath-Brown,
  primes kept on the original side, marks → labels, β*, slot-only witnesses, copies) lands on a precise
  step and gives nothing (details in the agent notes).
