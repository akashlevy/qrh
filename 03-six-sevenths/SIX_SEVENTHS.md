# Toward θ = 6/7: what is proved, what is conditional, and what blocks it

*8 Oct 2026. Companion to `QRH-analysis.md` and `STRENGTHENING.md`.*

**Method.** Two multi-agent workflows ran:
- five independent research routes;
- a focused check of the one surviving lead (the quartic transplant).

Every positive claim was then adversarially refereed and the exact certificates were re-run. Route notes and scripts are in `six_sevenths/routes/`; literature text is in `six_sevenths/lit/`.

## 0. Bottom line

**I have not proved that Hecke L-functions over ℚ(√−3) (or ζ, or Dirichlet L-functions) are zero-free in Re s > 6/7.** No route reached a complete, unconditional proof.
- Every argument that reaches 6/7 needs a mean-value estimate *below the large-sieve diagonal*. Nobody can currently prove one, and §3 shows why the existing machinery cannot.
- So there is no 6/7 theorem to formalize in Lean. I have not written a Lean proof, and I will not present a conditional or heuristic argument as one.

What **is** established:

| # | Statement | Status |
|---|---|---|
| 1 | **Theorem C1.** The D-averaged mean square M_{6/7} implies Re s > 6/7. The hypothesis is Σ_{N(u)≤D^{6/7}}\|A_u(D)\|² ≪ D^{13/7+ε} (or its dyadic D-average). | **Proved** (§1). The hypothesis is open. |
| 2 | Inside the 7/8 paper's Part II, the energy wall relaxed by w = 4/7 (W(4/7)) gives θ = 6/7 exactly. | **Exact certificate.** The hypothesis is open. |
| 3 | **Refined powerful-row count.** In Lemma reflected-energy, the case-A loss (1+3ℓ′−2M′)/4 becomes (1+2ℓ′−2M′)/3. | **Proved, new, elementary.** It does not move θ on its own (still 0.874957). |
| 4 | Item 3 plus W(2/7) gives θ ≤ 0.857094 < 6/7. | **Exact certificate** (80 235 boxes, worst margin −4.7·10⁻⁹). W(2/7) is open. |
| 5 | The detector floor 51/100 can be lowered to 1/2 + η. | **Provable.** It does not move θ on its own. |
| 6 | **Obstructions** (§3): the dual target is false below the diagonal; a move-calculus invariant; coefficient-blind impossibility; a 13/15 floor for any row-count input; a consistency barrier for bootstrapping. | Proved; the dual-target falsity is GRH-conditional in general and unconditional for κ < 1/4. |

The best *unconditional* exponent for this lemma set is **0.874957069799 = (1507 − 2√921)/1653**. Baiying Liu posted it on 8 Oct 2026, and it matches our exact LP floor 0.8749570. This is the "already known" re-parametrization. Liu also gets the 34999/40000 = 7/8 − 1/40000 from `certificate.py`.

## 1. Theorem C1 (conditional, complete proof)

**Hypothesis M_{6/7}.** For every finite-order Hecke character ν of K = ℚ(√−3), every W ∈ C_c^∞((0,∞)), and every ε > 0,

  ∫_D^{2D} Σ_{0<N(u)≤D^{6/7}} |A_u(y)|² dy/y ≪_{ν,W,ε} D^{13/7+ε},  where A_u(y) = Σ_{(n,S)=1} μ(n)ν(n)χ_n(u)W(Nn/y).

**Conclusion.** Every finite-order Hecke L-function over K is zero-free in Re s > 6/7, and so is every Dirichlet L-function, including ζ (with the principal-character pole at s = 1 allowed).

**Proof.**
1. *Setup.* Take H = D^{6/7} and Y = H^{1/6} = D^{1/7}. There are J ≍ Y/log Y primes p with Y/2 < N(p) ≤ Y and (p,S) = 1.
2. *Close to the target row.* Since χ_n(p⁶) = 1_{p∤n}, we have |A_1(y) − A_{p⁶}(y)| ≪_W y/Y for every such p. The rows p⁶ are distinct and have norm ≤ H. (This is paper2 §3.)
3. *Average.* So |A_1(y)|² ≤ (2/J)Σ_u|A_u(y)|² + O(D²/Y²). Integrating dy/y over [D,2D] gives
   ∫|A_1|² ≪ D^{13/7+ε}·D^{−1/7} + D^{12/7} = D^{12/7+ε}.
4. *Mellin transform.* Cauchy–Schwarz on dyadic blocks makes 𝓜_W(s) = ∫A_1(D)D^{−s}dD/D holomorphic on Re s > 6/7. There it equals Ŵ(s)/L_K^S(s,ν).
5. *Contradiction.* At a putative zero ϱ, the choice W = y^{−ϱ}φ(y) from paper2 §3 gives 0 = Ŵ(ϱ) > 0.
6. *Dirichlet case.* Base change L_K(s,χ∘N) = L(s,χ)L(s,χχ_{−3}) handles Dirichlet L-functions, exactly as in paper2. ∎

Both extraction terms equal D^{12/7} exactly at κ = 6/7. That makes 6/7 the floor of paper2's prime extraction: θ(κ) = max((6+5κ)/12, 1 − κ/6) is minimized at κ = 6/7. The 7/8 paper's integer extraction (Lemma inverse-amplification at U ≍ 1) gives (6+5κ)/12, again 6/7 at κ = 6/7.

**How strong the hypothesis is.** M_{6/7} is *strictly stronger than the theorem it proves.*
- It contains the row u = 1, so by itself it gives zero-freeness in Re s > 13/14.
- It is equivalent to a zero-density estimate for the sextic twist family, N(σ;Q) ≪ Q^{(13−14σ)/6}. That is stronger than the density hypothesis Q^{2(1−σ)} for every σ > 1/2 (`routes/R1_direct_meansquare.md` §7, C2).
- The 6/7 theorem for all Hecke characters gives only Σ_u|A_u|² ≪ H·D^{12/7}, so it does not imply M_{6/7}.

## 2. Inside the 7/8 paper's Part II (route R4, refereed)

The exact Part II model is `routes/r4/model.py`. It reproduces the paper's lemma-set optimum 0.8749568, gives 13/15 under GRH for rows, and matches R_* to 10⁻⁹.

| Replacement | θ | Provable? |
|---|---|---|
| none | 0.874957 | — |
| refined powerful-row count (item 3) | 0.874957 | **yes** |
| detector floor 1/2+η | 0.874957 | **yes** |
| Gram bound without the P_a^{1/6} term | 0.859796 | no — that term is a genuine main term at the sixth-power frequencies |
| amplification with M_{6/7} | 0.869792 = 167/192 | no |
| GRH for every row L-function | 13/15 | no |
| W(4/7) | **6/7 exactly** (certified at (4/7, 4/7, 3/7) and (22/49, 34/49, 3/7)) | no |
| refined count + W(2/7) | **0.857094 < 6/7** (certified) | W(2/7) no |
| refined count + W(w), threshold | w* ≈ 0.28493 | — |

**What W(w) asks for.** W(w) says the reflected cubic-theta energy in Lemma probe-row-norm loses nothing while T_d ≤ H + w. The excess over the wall is the column term NB of the quadratic large sieve (Goldmakher–Louvel / Heath-Brown), which is sharp for general coefficients. So W(w) is again a "Lindelöf on average with fewer rows than columns" statement for the specific theta coefficients.
- Numerically it looks true: per-row energies stay flat out to 8–128× past the wall (`strengthening/wall_experiment*.py`).
- No known technique proves it.

**Retracted R4 side claim.** R4 claimed that θ = 6/7 ⟺ the argument closes with row range h = 6/7. The referee refuted this. Only θ ≥ 1 − h/6 holds, and the 6/7 optima lie on faces of the LP, not only at h = 6/7. So the "fingerprint" reading of exact-6/7 vertices is not evidence about anyone's proof.

## 3. Why the existing machinery cannot reach 6/7 (proved obstructions)

1. **Paper2's dual target is false below the diagonal** (R1, Prop. O1).
   - Setup: paper2 proves the mean square by Poisson in u, μ-absorption, a theta reflection and the quadratic large sieve. For H = D^κ with κ < 1, its intermediate hypothesis E(𝓗,X) ≪ D^ε X (eq. `auxiliary-target`, with 𝓗 = X^{2−κ}) is false.
   - Under GRH, E ≫ 𝓗. Unconditionally the same holds for κ < 1/4 (via the 7/8 theorem).
   - The loss is the absolute value on μ(f) in eq. `weighted`. In the exact identity, the "fake diagonals" cancel down to a Möbius sum over the gcd variable.
   - Numerics at X = 6000 confirm this. The positive form equals its diagonal (1.000 ± 0.007), while the true coprime-pair form is 10⁻⁴–10⁻² of it. *The cancellation M_κ needs exists, but it lives only in the signs μ(v)*, which every positivity/large-sieve step discards.
2. **Move calculus** (R1). Track (δ, ε) = (target − diagonal, target − column term).
   - Poisson swaps δ and ε; theta reflection maps (δ, ε) to (δ, 2δ − ε).
   - Both preserve |ε − δ| = 1 − κ, and from κ < 1 every reachable node has min(δ, ε) ≤ −(1 − κ).
   - So no arrangement of Poisson/reflection/large-sieve steps proves M_κ, even granting optimal large sieves at every level.
3. **Coefficient-blind impossibility** (R3, Theorem B). Any large-values bound valid for all coefficients |c_n| ≤ 1 must allow about H^{1/6} rows (the copies u = a⁶) to be large. This rules out the large sieve, Halász–Montgomery, Huxley subdivision and Guth–Maynard as sources of the needed input.
   - Theorem A: the multiplicity-counted large-values hypothesis is *equivalent* to the conclusion, so there is no slack.
4. **Part II floor** (R3, Theorem D; exact LP). With the paper's low-side lemmas, *any* row-count input, even GRH for every row, gives θ ≥ 13/15 > 6/7. The density hypothesis for rows already attains 13/15.
5. **Bootstrap barrier** (R2). Consider a counterfactual where one target has a zero at β ≤ 13/15 and every other primitive Hecke L-function satisfies GRH. At the exponent level it is consistent with the conclusions of every lemma in both papers. So no continuation/bootstrap argument built from those lemmas can reach 6/7.

**Literature** (R5).
- Nothing published supplies a usable input.
- de Faveri's optimal n-th order large sieve (arXiv 2610.04045, 2 Oct 2026) is genuinely new. It improves Part I's sextic row count but changes θ in neither paper.
- de Faveri–Dunn–Hoffstein (arXiv 2607.07911) show Heath-Brown's cubic large sieve is sharp.
- Gao–Zhao's quadratic analogue of M_κ needs GRH and works only above the diagonal.

## 4. What a proof of 6/7 would have to contain

At least one of the following, none of which is available:
1. a "column functional equation" for μ against sextic characters (GRH-type information);
2. a power-saving asymptotic for the second moment of quadratic twists of the cubic theta function at unbalanced length X = 𝓗^{7/8}, saving 𝓗^{−1/8} (an affine-type, "beyond the edge" moment);
3. a non-positive (sieve-theoretic or spectral) treatment of the coprime-pair dual form that keeps the signs μ(v);
4. a lossless marked moment / energy beyond the wall (W(w) with w ≳ 0.285, together with the new refined count);
5. a different engine altogether: see §5 (the quartic lead).

## 5. The quartic lead: refuted

The idea was to transplant the whole argument to ℚ(i) with quartic symbols, Kubota's quartic Eisenstein series, and de Faveri's new optimal quartic large sieve (arXiv 2610.04045). R5's heuristic model gave θ ≈ 0.854. A focused workflow derived the true geometry, and the result is negative. Its notes and scripts are in the scratchpad `wf67q/`; I stopped the referee stage, because only a reversal of the conductor claim could revive the route.

- **The row twist enters as a cube.**
  - After Poisson in u, absorbing μ switches the Gauss sum to the conjugate character. So the dual sums are partial sums of g̃₄(λ⁻²h³, c), not g̃₄(h, c).
  - Brute force for every primary z up to norm 65, including composite z, confirms this: with shift h³ the ratio is independent of h to about 10⁻¹³; with shift h it spreads by 2.0.
- **So the reflection conductor is N(h)³.**
  - The reflection is the Kubota–Diaconu functional equation, which carries the factor N(ν)^{1−s} with ν = λ⁻²h³.
  - The conductor is therefore neither N(h)² nor N(h), and the dual length is N(h)³/X.
  - This exponent rests on the functional equation together with a growth argument. It was not visible numerically, because a large fixed 2-adic factor hides it below X ≈ 10⁵.
- **The dual twist stays quartic, and there is no polar term.**
  - The sextic case drops to a quadratic twist through the automorphy of the cubic theta function. Quartic Gauss sums are Eisenstein coefficients, not theta coefficients, so that mechanism is unavailable.
  - The series is entire for angular frequency ℓ = −1, and Suzuki's theorem makes the residue vanish at cube shifts anyway. The polar-term absence was checked numerically up to norm 2·10⁶.
- **Outcome.**
  - With unconditional inputs (de Faveri n = 4 or BGL), the transplanted Part II gives **θ = 1, i.e. no zero-free half-plane**: the optimum sits on a plateau where the low bound equals the trivial bound.
  - Patterson's quartic conjecture changes nothing, since it enters none of the derived exponents.
  - In paper2's simpler architecture, the transfer recursion closes only when the row range satisfies H > D^{4/3}, while the fourth-power extraction needs H < D^{4/3}. These are exactly complementary.
  - The earlier 0.854 came entirely from the wrong N(h)² conductor.

## 6. Lean

The OpenAI formalization of the 7/8 theorem is 2 924 modules and about 486 000 lines on Lean 4.34.1. No Lean toolchain is installed on this machine. Since there is no unconditional 6/7 proof, there is nothing to formalize as a 6/7 theorem. Two things *can* be formalized honestly if wanted:
- Theorem C1 as a conditional statement (hypothesis M_{6/7}), on top of the existing Hecke/Mellin infrastructure.
- The exact rational certificates: the W(4/7) and refined count + W(2/7) certificates, the 13/15 floor, and the energy dual certificates. These are pure linear arithmetic over ℚ.
