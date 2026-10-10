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
