# Strengthening the OpenAI 7/8 argument: what actually moves θ

*Companion to `QRH-analysis.md`. Revised after an adversarial audit (multi-agent workflow, 8 Oct 2026) found modeling errors in the first draft. The corrections are in §1. All θ values come from the Part II bookkeeping in `strengthening/conditional_certificate.py`. With nothing relaxed, that bookkeeping reproduces the paper's 7/8 exactly.*

## 0. Short answers

**Does W(1/5) give about 6/7? No.** Correctly modeled, relaxing the marked-energy wall by w = 1/5 gives θ ≈ **0.864**. This is certified in exact arithmetic: W(1/5) ⇒ Re s > 108/125. The low side alone is bounded below by **13/15 − w/60** for every choice of geometry, so reaching 6/7 that way would need w ≥ 4/7.

**Where 6/7 does come from.** It is exactly what you get from one clean lemma replacement: the twisted-Möbius mean square below the large-sieve diagonal,

> Σ_{0 < N(u) ≤ D^κ} |A_u(D)|² ≪ D^{1+κ+ε},  with κ = 6/7,
>
> where A_u(D) = Σ_n μ(n) ν(n) χ_n(u) W(N n / D).

This is the companion paper's Prop. `thm:ms`, and the 7/8 paper's eq. `raw-moment`, pushed from rows H ≥ D^{1+ϑ} down to H = D^{6/7}. Two extraction routes then give exactly 6/7 (§2):

- **The companion's prime sixth-power extraction.** θ = max((6+5κ)/12, 1 − κ/6). The minimum over κ is attained at κ = 6/7, where θ = 6/7.
- **The 7/8 paper's sixth-power amplification with one fixed row (U ≍ 1).** This gives θ = (6+5κ)/12, again 6/7 at κ = 6/7.

If your friend's result is exactly 6/7, the most likely replaced lemma is this mean square, or a statement equivalent to it.

## 1. Corrections to the first draft (retracted claims)

An independent audit agent re-derived the bookkeeping from `paper.tex`. It found three errors in my earlier model of W(w); I have checked each one and confirmed it.

- **(a) The low-side energy loss was wrong.** I used (1 + 2ℓ′ − 2M′)/3 as the "exact maximum" of eq. 4.6. That silently assumes the rows have no powerful part (O = 0).
  - With O free, the maximum over admissible block data is the paper's own (1 + 3ℓ′ − 2M′)/4.
  - It is attained at O = M′ + ℓ′/2 − 1/2, H = (1 − ℓ′)/2, z_a = ℓ′, with T_d = 0. Then E − M′ = M′/2 + 3ℓ′/4 + 1/4 − M′ = (1 + 3ℓ′ − 2M′)/4.
  - This is **case A**, which W(w) does not touch. It forces ℓ ≤ (1 + 2w)/5 when M + ℓ = 1 + w, and hence low side ≥ 13/15 − w/60.
- **(b) W(w) does not relax Lemma `marked` to r + 2z < 1 + w.** Its puncture invariant F − M − z₀ > 0 contains the canonical diagonal "rows ≤ columns", which is r + z < 1. W(w) cannot relax that. Only the case-B part relaxes, and only by w/2. The correct capacity is z ≤ min((1 + w/2 − r)/2, 1 − r, (3 − 2r)/8), for r < 1.
- **(c) W(w) does not relax the amplification slope α = 5/6.** The condition H ≥ D^{1+c} is that same canonical diagonal with no slots, so α stays 5/6.

**Retracted:**
- the claims W(1/10) ⇒ 0.865 and W(1/5) ⇒ 0.858 / 6/7;
- the table "one wall controls three lemmas, θ(w) = 0.869, 0.864, 0.854, …";
- the payoff rows "energy wall alone +0.05/+0.10/+0.25 ⇒ 0.872/0.869/0.859".

The flawed script is kept as `strengthening/conditional_certificate_SUPERSEDED.py`, with a header explaining the errors.

**Still valid:**
- θ = 7/8 − 1/40000 (`certificate.py`, which uses the paper's formula);
- the lower bound ≈ 0.874957 for the paper's lemma set, and the GRH-for-rows barrier 13/15 (`barriers.py`). The corrected constraints are stronger, so those LP lower bounds still hold;
- the quartic absorption identity;
- the floor-lowering result, confirmed legitimate by its own audit;
- the numerical energy experiments, which are facts about case-B energies, now with a much smaller payoff.

## 2. The 6/7 lemma in detail

**Statement (M_κ).** For a fixed finite-order Hecke character ν over ℚ(√−3), every fixed W, every ε > 0, and H = D^κ:

> Σ_{0 < N(u) ≤ H} |A_u(D)|² ≪ D^{1+ε} H.

Paper2 proves this for κ = 1 + ϑ, ϑ > 0, which is the large-sieve diagonal regime.

**Why it gives 6/7.**
- *Prime extraction (paper2 eq. `prime-extract`).* The rows u = p⁶ with N(p) ≈ Y = H^{1/6} give |A₁|² ≪ D·H^{5/6} + D²·H^{−1/3}. So θ(κ) = max((6+5κ)/12, 1 − κ/6), which equals 6/7 at κ = 6/7. This is the *absolute floor* of that extraction: even κ → 0 cannot push it below 6/7.
- *Integer extraction (7/8 paper, Lemma `inverse-amplification` at U ≍ 1).* The exact identity M_u(D) = Σ_{d | rad a} μ(d) ψ_u(d) q_d^{−1/2} M_{ua⁶}(D/q_d) removes the D²/Y² term. Then θ(κ) = (6 + 5κ)/12, which is 6/7 at κ = 6/7, 3/4 at κ = 3/5, and tends to 1/2 as κ → 0.

**Relation to W(w).** M_κ is the *unmarked, rows-beyond-columns* case: raw moment with rows ≥ D^κ. It is not implied by the energy wall W(w); the audit's point (c) is exactly this distinction. In the probe architecture M_κ would also lift the amplification slope. The simple extraction above, however, needs *no* probe machinery at all.

**How plausible.** M_κ for κ < 1 beats the large sieve's column term for Möbius sums over sextic twists. Heuristically it is at least as strong as a density-hypothesis-type statement for the row family; as κ → 0 it becomes GRH. So it is a genuinely strong input, not a routine strengthening.

**Numerics** (`strengthening/mobius_meansquare.py`; exact Möbius and sextic symbols in ℤ[ω], all rows including sixth powers). The per-row energy divided by its diagonal value:

| κ | 0.5 | 0.6 | 0.7 | **6/7** | 0.9 | 1.0 | 1.1 |
|---|---|---|---|---|---|---|---|
| D = 3000 | 1.17 | 1.07 | 0.98 | **0.91** | 0.91 | 0.94 | 0.95 |
| D = 10⁴ | 0.98 | 0.96 | 0.97 | **0.94** | 0.94 | 0.94 | 0.94 |

So M_{6/7} holds comfortably at these sizes. **Caveat:** small-scale Möbius sums show square-root cancellation anyway, consistent with GRH. These numerics are therefore consistent with M_κ for any κ, but do not show that it is *provable*.

## 3. What the marked-energy wall W(w) is worth, corrected

The corrected exact certificate is `strengthening/conditional_certificate.py`, run via `certify.py`.

| Hypothesis | Proven zero-free half-plane | Geometry (l_x, l_y, ℓ) |
|---|---|---|
| W(0), the paper's lemmas | Re s > 7/8 (reproduced) | 17/48, 23/48, 1/6 |
| W(1/5) | Re s > 108/125 = **0.864** | 91/200, 93/200, 7/25 |
| W(w), any geometry | low side ≥ 13/15 − w/60 (`energy_lp.py`) | — |

The audit's float model gives W(1/10) ⇒ ≈ 0.868.

**Numerical evidence.** The case-B energies really do stay diagonal far past the wall (`wall_experiment*.py`):
- unmarked sums, to 8× the wall;
- single-mark sums, to 128×;
- full slot sums with cross terms, to 100×, at X = 3000 and X = 10⁴.

So W(w) is very plausibly *true*, but in this architecture it is worth little. The binding obstructions are:
- the case-A powerful-row term (1 + 3ℓ′ − 2M′)/4;
- the canonical diagonal "rows ≤ columns".

The canonical diagonal is the same diagonal that M_κ crosses.

## 4. Other levers (unchanged, small)

- **Row-count lemmas** (plain capacity c_P, inverse capacity c_M, α): each at most ≈ 0.002–0.008. Perfect counts cap at 13/15.
- **The Gram P_a^{1/6} term** is structural: it carries the signal's sixth-power frequencies.
- **Detector floor 51/100 → 1/2 + η** is legitimate (audited) but buys almost nothing alone.
- **Quartic setup over ℚ(i)** (`order_n_engine.py`, `quartic_routes.py`): ≈ 0.81 via the theta route, which needs Patterson's quartic conjecture. The μ-absorption identity was verified on 538 primes. The Eisenstein route needs a better quartic large sieve.

## 5. Recommendation

To reach 6/7 or below, the lemma to replace is the **twisted-Möbius mean square below the diagonal (M_κ)**, not the marked-energy wall:
- κ = 6/7 ⇒ θ = 6/7;
- κ = 3/5 ⇒ θ = 3/4;
- κ → 0 ⇒ GRH.

Through the simple extraction it needs nothing else from either paper. The open problem is proving M_κ for some κ < 1: square-root cancellation on average for Möbius sums over a family of sextic twists smaller than the sum length. Theta reflection doesn't help here: after Poisson in u, the dual cubic-theta mean square has rows D^{2−κ} > D, which is the regime where reflection loses.
