# Pushing below Liu's 0.874957: results on the four ideas

> **CORRECTION (9 Oct 2026, after the audit workflow `wf/`): the conditional improvement below does NOT hold.**
> Step (ii), bounding the theta pairing via (V) and Heath-Brown, applies to the wrong sum.
> - The reflected row B_m carries the paper's reflection scalar ζ_m (eq. reflection-row-phase, l.1832; l.3017 "ζ_R depends on R").
> - For squarefree m, ζ_m = κ₀ψ₁(m)·conj(γ₂(m))⁴, a cubic Gauss sum times an angle character (Lemma Z, `wf/V_gaps.md`). Three agents confirmed this independently, numerically to 10⁻¹⁵.
> - So the pairing is Σ_m ρ(m)ζ_m(μ/m)₂, not Σ_m ρ(m)(μ/m)₂. Lemma (V) covers only the latter.
> - Without a bound on the true pairing the LP gives exactly 0.874957 (Liu), even with the P²/Y term removed. The slot-free variant gives 0.8947.
> - **Hence 0.8665, 0.871 and 0.864 are all unsupported.** The certificates verified only the bookkeeping.
>
> What survives (verified by independent referees):
> - the theta-pole decomposition of the Gram excess;
> - the signal-cancellation obstruction (the excess cannot be cancelled, only bounded);
> - (FE)/(LoA6) and Corollary D (`wf/FE_LoA6_proof.md`), via Kubota 1969, not Brubaker–Bump, whose hypothesis μ₁₂ ⊂ F fails here. This gives ‖A − A^sig‖² ≪ (Q/Y)Z^ε, with no P_a^{1/6} and no P²/Y term;
> - the slot-free Part II machinery and zero-slot counts (`wf/PARTII_slotfree.md`).
>
> **The single missing input is now (V′):** a non-Cauchy–Schwarz bound for Σ_m ρ(m)ζ_m(λ⁴μ/m)₂ W(N(m)/Q), uniform in μ ≲ Q²/Z. Any uniform power saving s₀ beats Liu:
>
> | s₀ | θ |
> |---|---|
> | 0.002 | 0.8733 |
> | 0.01 | 0.8668 |
> | the large-sieve shape | 0.8662 |
>
> Note: ρ·ζ = ρ·(the row's own sextic Gauss sum)·(Hecke)·(ray). On squarefree parts this is a residue of the sextic A₂ Weyl-group multiple Dirichlet series (BBFH), a purely sextic object. Its conductor under the quadratic twist decides (V′):
> - dual length L²/Q gives a saving;
> - L³/Q or L⁴/Q gives none.

> **Follow-up on (V′) (workflow `wf2/`, both reports refereed): no saving with current tools.**
> - On squarefree rows, π(m) = ρ(m)ζ_m is the theta residue of the sextic twisted A₂ multiple Dirichlet series. The quadratic twist enters as the twisting parameter μ³, not as a ramified character.
> - The surviving functional equation is z ↔ 5/6 − z with factor N(μ)^{3(1/2−z)}. That gives a heuristic conductor N(μ)³ and dual length L³/Q = Z^{1.13}, which exceeds Q.
> - The referee found a stronger obstruction: the dual coefficients ρ(sμ^k) couple s and μ, so no large sieve applies. Even k = 2 gives no saving without a μ-independent character kernel; a saving would need k < 1.53 at M = 0.826. With the dual route the LP gives min(0.874957, 0.934) = Liu's value.
> - Numerics (Q up to 6000, exact ζ_m on 8280 rows):
>   - π is not simpler than ρ;
>   - the cubic-type degeneracy is excluded;
>   - its quadratic twists are at the random level, so (V′) is plausibly *true*;
>   - the test cannot resolve the conductor beyond degree 1.
> - The literature WMDS theory assumes μ₁₂ ⊂ F, which fails here. Those identifications are [H] or [R-].
>
> **Net status:** proving (V′) needs a "beyond the edge" large sieve with few moduli: family size L, conductor ~L³, against a sequence with no known automorphic kernel. That is the same type of input that blocks Ideas 2 and 3. **No unconditional improvement on 0.874957 is available along Idea 1 with current technology.**

> **Direct attack on (V′) (workflow `vls/`: 5 attacks, each refereed; target in `TARGET_VLS.md`): not proved.**
>
> | Attack | Result | Key finding |
> |---|---|---|
> | A, functional-equation geometry | no saving | Complete search of the exact 3-node model U = (s,m,A): 24 chambers, 4 residue-adapted types, none beats Cauchy–Schwarz at M ∈ [0.80, 0.86]. Also corrects an orientation error in `wf2/A2_theory.md` §4.2: the faithful 4-node object is of infinite type. |
> | B, Gauss-sum structure | no saving | β's structure is already used up by the paper's reflection, so any argument using only β is capped near 0.89·Q^{11/12}. A proof must use the automorphy of π. |
> | C, moments and spectral | no saving | The second-moment Coxeter group is hyperbolic for k = 3, and its convexity hull gives only the trivial bound. No Kuznetsov formula or spectral large sieve exists on the needed cover. T contains the full zero-detecting signal, so signal-blind bounds are circular. |
> | D, probe redesign | impossible | **Rigidity:** Möbius absorption needs the row twist j = 1, but ζ-freeness needs j ∈ {2,5}. So F1 is intrinsic to every admissible probe. |
> | E, numerics | not a proof | Q up to 9.6·10⁴: (V′_LS) holds at the random rate, (L/Q)^{1.008±0.007}, with saturation excluded by ×30–50. (V′_β) and the Type I sums are at the random level. The conductor exponent is still undeterminable. |
>
> **Conclusion.** (V′) is very likely *true*, but proving it needs a beyond-the-edge mean value. Equivalently, it needs Lindelöf-on-average for the twisted family Π(z; μ) with exponent θ_π < about 1.31 (the large sieve gives 3/2, full LoA gives 1). No current technique supplies this. If it were proved at the random rate, θ would be about 0.852 (LP, conditional [H]).


*9 Oct 2026. Working files are in `gram_split/`. The detailed notes are `IDEA1_NOTES.md`, `V_derivation.md`, `D_derivation.md` and `LIT_ideas234.md`.*

## 0. Original bottom line (SUPERSEDED by the banner above; kept for the record)

Idea 1 works, conditionally. The two missing inputs are (FE)/(LoA6) and (V), both about the sextic theta function and Kubota series. Each has a derivation, cross-checked numerically, but neither is a finished proof. Write

> A = A^sig + A^rest,

where A^sig is the theta-pole (sextic theta) component of the additive polynomial. Then:

- **(a) No Cauchy–Schwarz on the theta part.** The P_a^{1/6} Gram term is exactly ‖A^sig‖². Its pairing with the reflected row B is bounded directly. A Voronoi step turns the quadratic twist of θ_6 into a cubic-character dual sum, and Heath-Brown's cubic large sieve then applies.
- **(b) The rest has a clean Gram.** ‖A^rest‖² ≪ (Q/Y) Z^ε. This follows from Lindelöf-on-average for the sextic Kubota series, and it also removes the paper's P²/Y Gram term.

Together these allow a slot-free geometry. Under hypotheses (V) and (LoA6) below, the Part II bookkeeping closes at **θ ≈ 0.8662**, against 0.874957 for Liu.

| assumptions | θ (LP optimum) | exact certificate (stdlib Fractions, r4 checker) |
|---|---|---|
| (V) + (LoA6), base probe only (no slots, ℓ = 0) | 0.86624 | **ALL PASS at σ0 = 1733/2000 = 0.8665**; geometry (989/5000, 3141/5000, 0) |
| (V) + (D) + Part II with slots and marks (P²/Y kept) | 0.87072 | **ALL PASS at σ0 = 871/1000**; geometry (2671/10000, 2671/5000, 107/2500) |
| (V) + (LoA6), slot-free, plus the provable detector floor 1/2+η | 0.86361 (δ0 = 1/200); 0.86330 (δ0 = 1/1000) | **ALL PASS at σ0 = 108/125 = 0.864**, δ0 = 1/200; geometry (1817/10000, 3213/5000, 0); 79,327 boxes, log `cert_floor_8640.log` |
| reference: paper's lemmas (Liu) | 0.874957 | — |

**Row counts at small or zero slot supply.** The slot-free variant is consistent with the paper's own lemmas; details in `IDEA1_NOTES.md` §6.
- The witness proposition holds for every t ∈ [1, 3/2] with no slots.
- Every witness case then has a slot-free count bound: no-slot inverse amplification, or zero-slot plain.
- This gives R_0(t) = max(L(t), 1 − 2δt/3), which is exactly the R4 model at S = 0.
- The paper's remark "t = 1 gives 1 − 2δ/3" is only one choice of t. Using t = 1 alone would give no improvement (0.8751).

The model also gives:
- 0.8631 with an "ideal" large sieve, which de Faveri–Dunn–Hoffstein show is not available in general;
- about 0.837 if the theta pairing were negligible. The (V) derivation rules that case out: the large-sieve coefficients are genuinely oscillating Gauss-sum products.

**Status.** This is a conditional improvement. Both hypotheses are statements about the sextic theta function and Kubota series that standard technology should prove. Neither is a below-diagonal mean value of the M_κ or W(w) type. Neither has been proved here, and the (V) derivation still has the gaps G3–G4 listed in §3.

**Two negative results** (rigorous, from the Mellin structure of the paper's high side):
- The P_a^{1/6} excess and the principal signal have the same scaling, Y^{-1/3} at fixed Q.
- Precisely, the signal coefficient is K̂(5/6), the Mellin transform at 5/6 of the theta term's weight K(μ) = ∫ W_1(t) t^{1/6} W_0(μ/t) dt/t.
- Hence no combination of probes, weights or scales can cancel the excess without also cancelling the signal. The theta pairing has to be *bounded*, not removed.

## 1. Idea 1 in detail

The full working is in `IDEA1_NOTES.md` §1–7, `V_derivation.md` and `LoA6.md`.

**Conditional theorem (slot-free form).** Assume:
- **(LoA6)** Σ_{N(m)≤Q} |D_c(1/2+ε+it; m)|² ≪ Q^{1+ε}(1+|t|)^R. This is to hold for the sextic Kubota series D_c(w;m) = Σ_s c(s) g̃(s,−m) N(s)^{-w}, for each of the finitely many ray-class twists c used by the probe, and also at the finitely many other cusps.
- **(V)** The Voronoi formula of `V_derivation.md` §5.2 for Σ_m ρ_c(m)(μ/m)_2 W(N(m)/Q), with the dual-coefficient bounds G3–G4. Here ρ_c is the residue of D_c at w = 2/3.

Then the base probe of the 7/8 paper, with no slots, satisfies the low bound |I_η(X,Y,Z)| ≪ Z^{l_x/2+ε} when b/6 ≤ 2s(M). Here s(M) = −½·max(M−1, 7M/3−2, −M/6), M = l_x+l_y, b = l_y−l_x. Combined with the paper's own high-side lemmas, this gives ζ(s) ≠ 0 for Re s > 0.8665 (exactly certified), and 0.8633 with the provable detector floor.

**Proof skeleton.**
1. **Split.** A = A^sig + A^rest, where A^sig = Y^{-1/3} Ŵ_1(1/6+iv) ρ_c (the residue at the theta pole), and ‖A^rest‖² ≪ (Q/Y)Z^ε by (LoA6) (contour shift; see `LoA6.md`).
2. **Rest.** Apply Cauchy–Schwarz to ⟨A^rest, B⟩, with the paper's row norm ‖B‖² ≪ Z^{M+ε}. This gives (X)^{1/2}Z^ε.
3. **Theta part.** Use the reflected form of B (paper eq. reflection): B_m = Σ_{μ≲L} β_μ (μ/m)_2 + …, with L = Q²/Z. By (V), ⟨ρ_c, B⟩ is a bilinear form in conjugate cubic symbols, plus main terms at square μ totalling O(Q^{5/6+ε}). Heath-Brown's cubic large sieve then gives |⟨ρ_c,B⟩| ≪ Q^{5/12}(L + L²Q^{-2/3})^{1/2}Q^ε + Q^{5/6+ε}. This is below ‖ρ‖‖B‖·Z^{-b/12} exactly when b/6 ≤ 2s(M).
4. **Counts.** The paper's high side is used unchanged. The zero-slot row counts are R_0(t) = max(L(t), 1−2δt/3), from the witness proposition and the slot-free count lemmas (§6 of the notes).
5. **Close.** The exact Fraction certificate `exact_cert_theta.py` checks all low, floor, small-row and high-side conditions at the stated geometry.

## 2. Ideas 2–4

Full citations are in `LIT_ideas234.md`, with each fact labelled [V] (verified in the source), [S] (from a secondary summary) or [I] (inference).

**Idea 2: spectral (Kuznetsov) treatment of the dual form.** No route found.
- A Bruggeman–Kuznetsov formula on the cubic cover over Q(ω) exists: Louvel, Forum Math. 2014, Thm 4.2. He places the exceptional spectrum in (0, ≤1/6] ∪ {1/3}.
- There is no spectral large sieve for metaplectic Maass forms. Watt's inequality is non-metaplectic.
- There is no Kuznetsov-based second moment for theta coefficients in the regime with fewer rows than columns. The affine "beyond the edge" second moment that the dual form needs is unavailable for any family.
- Verdict: a long shot, unchanged.

**Idea 3: proving the energy wall W(w) directly.** No route found.
- The only Lindelöf-on-average results for Kubota or theta series are at or inside the edge:
  - DdFDS Lemma 6.6, cubic;
  - DFDH, quartic;
  - DDHL (1.11), which gets A^{3/2}, not A^{1+ε}, for a family of A members with conductor A².
- There is no mean-value theorem for L(½, θ_3 ⊗ χ_d).
- I also checked whether the reflected and unreflected forms of the marked row help. They do not. With marks, the reflected length Q²Z^{ℓ'−1} and the unreflected column count Z^{1+ℓ'} both exceed Q in the wall region. Factoring out χ_D(m) fails because the marked primes are summed inside the row. So the wall really is "beyond the edge" on both sides.
- Note: the Voronoi trick of Idea 1 does not transfer. For B, Voronoi only toggles between the two forms that the paper already uses.
- Verdict: a long shot, unchanged. W(w) is not needed for the Idea-1 result (the energy wall is inactive at the new optimum).

**Idea 4: other families.** Sextic is confirmed unique.
- Davenport–Hasse with p = 2 gives g(χ)g(χ^{m+1}) = χ̄(4) g(χ²) g(χ^m) for χ of order 2m.
- It linearises Möbius × Gauss sum into a single Gauss sum of order m only when m is odd:
  - n = 6 → cubic, theta known (Patterson; Kazhdan–Patterson);
  - n = 10 → quintic, theta unknown;
  - n = 4, 8, 12 give only "square-absorption".
- Larger fields containing μ_6 have the same tools with the same exponents in norms. I see no gain.
- Quartic over Q(i) remains dead: the conductor is N(h)³.

## 3. What remains to make this unconditional

0. **What the (D) agent established** (`D_derivation.md`).
   - **Lemma D1.** The cross-Gram is G(Y₁,Y₂) = diagonal + Q^{5/6}(Y₁Y₂)^{-1/3} Ŵ Ŵ̄ K₀ + O((Q/Y)(1+P²/Y)Z^ε). It is rank one in the scale variable, with an explicit convergent constant K₀ (6.594·10⁻³ in the numerical model).
     - Status: [R]/[R-] from the paper's Poisson identity.
     - It refines the paper's frequency classification: only k = ε₀C⁵a⁶ carries a main term.
   - **The strong (D) follows from the Kubota functional equation (FE) alone**, i.e. ‖A − A^sig‖² ≪ (Q/Y)Q^ε with no P²/Y term:
     - the large-sieve step can use the paper's own planar additive sieve;
     - the residue term in the approximate functional equation is handled by bootstrap.
   - **Numerics** (Q = 1500, 3000, 6000; Y up to 5·10⁵):
     - the residual ‖A − A^sig‖² equals 0.95–0.97 × the exact diagonal, flat out to P²/Y ~ 10⁸;
     - the fitted surviving excess is −0.2 ± 8%, −3 ± 6% and +4 ± 5% of the theta excess.
   - **Caveat.** D1 alone does not identify the rank-one vector with ρ; that needs (FE). A finite-scale alternative (§6.5(c)) avoids (FE) for the Gram part, but moves the problem to a pairing bound for A(X).
1. **(LoA6): Lindelöf on average for sextic Kubota series in the m-aspect.**
   - Statement: Σ_{N(m)≤Q} |D(1/2+ε+it; m)|² ≪ Q^{1+ε}(1+|t|)^R.
   - Expected route: the n = 6 analogue of David–de Faveri–Dunn–Stucky Lemma 6.6 (cubic). It uses the approximate functional equation of Kubota series (Brubaker–Bump 2006) and the sextic large sieve (BGL; de Faveri Thm 1.1 gives Ξ_6(Q, Q^{1/2}) ≪ Q^{1+ε}).
   - No GRH is needed. DFDH say the n ≥ 5 analogue is "expected".
   - Consequence: ‖A − A^sig‖² ≪ (Q/Y) Z^ε.
2. **(V): Voronoi for quadratic twists of the sextic theta, uniform in the modulus.**
   - Statement: Σ_m ρ(m) χ_μ(m) W(N m/Q) has dual length ≍ N(μ)²/Q. The dual is twisted by the conjugate cubic character (n/μ)_3‾, with root number ε(μ) = γ₂(μ)/γ₃(μ). There are no main terms for squarefree μ ≠ 1; main terms occur only for μ ∈ uλ^k·(squares), and total O(Q^{5/6+ε}).
   - Derivation: `V_derivation.md`. It uses Lemma M: the multiplier of θ_ρ at the cusp a/μ is χ_μ(a), proved from Kubota's character on Γ(72).
     - The orientation is fixed three independent ways.
     - The sextic numerics agree: only j = 1 is coherent, and F₁(c)·γ₅(c) shows a common phase across primes, ratio 0.686 ± 0.004.
     - For the exact cubic analogue, the same computation reproduces the paper's own proven local factor χ_p^{−j−2}.
   - Remaining gaps:
     - G3: mean squares of θ₆'s coefficients at the cusps 1/μ₀. This needs a Zagier-regularised Rankin–Selberg, which de Faveri–Dunn–Hoffstein do for n = 3, 4.
     - G4: non-squarefree dual indices, which need the Kazhdan–Patterson local theory of θ₆.
     - G1, G2: bookkeeping.
   - The cubic-cover analogue, Dunn–Radziwiłł Lemma 5.2 / Prop 5.1, is already uniform in the modulus.
   - Afterwards one applies Heath-Brown's cubic large sieve, which is sharp in general by DFDH. The coefficients x_μ = β_μ ε(μ) are oscillating Gauss-sum products, so there is no extra cancellation for free.
3. **(P): bookkeeping inside Part II.** The slot-free variant needs only the base probe, with no marks and no compensation. What remains:
   - non-squarefree rows, where the reflected twist is conj χ_{m·rad(m)²}(μ);
   - the v-integral (harmless, since Ŵ_1(1/6+iv) decays);
   - the ray-class restriction s ∈ σ (θ with a character);
   - a check that the paper's high-side lemmas hold at ℓ = 0 and l_x ≈ 0.2. The R4 model and the exact checker treat S = 0 as admissible, but the paper itself only works at ℓ > 0.
