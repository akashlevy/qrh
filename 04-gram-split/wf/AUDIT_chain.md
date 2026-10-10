# Adversarial audit of the Idea-1 chain (i)–(iii)

*9 Oct 2026. Scope: the slot-free claim of `REPORT.md` §1. The claim is that (LoA6) + (V) give the low bound for the base probe at σ0 = 1733/2000 = 0.8665, with geometry (l_x, l_y, ℓ) = (989/5000, 3141/5000, 0).*

*Labels:*
- *[R] rigorous;*
- *[R-] routine but not written;*
- *[H] heuristic;*
- *[G] gap.*

*Paper line numbers refer to `formalize_11_12/papers/paper78.tex`.*

*Supporting scripts, in this directory:*
- *`check_row_phase.py`: the numerical check in §1.3;*
- *`audit_lp_no_theta.py`: the LP with and without step (ii), §1.5.*

## Verdict: FLAWED

The chain has a fatal gap in step (ii).

The reflected completed row B_m contains a row scalar ζ_m. This scalar depends on m through Gauss sums: up to ray-class phases, it is a normalised cubic Gauss sum of m times an angle character of infinite order. Two notes treat it as "a fixed ray phase":
- `IDEA1_NOTES.md` line 49;
- `REPORT.md` §1 step 3.

That is false (§1 below). As a result, the theta pairing ⟨ρ, B⟩ is a bilinear form in the sequence m ↦ ρ(m)ζ_m, not in m ↦ ρ(m).

The Voronoi lemma (V) of `V_derivation.md` is a statement about the automorphic sequence ρ(m)(μ/m)₂. Even if it were fully proved, it would not bound the pairing that actually occurs. There is no known automorphic structure for ρ(m)ζ_m, so the saving s(M) has no support.

Without step (ii), the LP gives:
- θ = 0.894737 at zero slots;
- θ = 0.875092 at the best slot supply, i.e. the paper's value (Liu's 0.874957 up to grid resolution).

Step (i) alone gives nothing. So all three headline numbers fail with (ii): 0.8665, 0.871 and 0.8633.

The rest of the bookkeeping is essentially sound:
- the exponent translation;
- the absence of circularity;
- the zero-slot row counts;
- the exact certificate, which I re-ran and which passes.

There are also some secondary gaps, listed in §6.

---

## 1. Fatal issue F1: the reflected row carries an m-dependent root number ζ_m

### 1.1 What the paper says [R]

**The base-probe row.** B_{m,σ}(Z) = Σ_θ a_θ 𝒞_{V_G}(Z; m, ν_σθ) (around line 3440). Here the character is Ψ = ν_σθ·χ_•(m), so the moving local primes of the reflection are the primes of m, with exponents j_p = v_p(m) mod 6.

**The reflection (line 1763).** It has the form ζ·Σ_μ d(μ)α(μ)ϑ(λ⁴μ) q_μ^{-1/2} ∏_p B_p(λ⁴μ) V♯(q_μ Z/q_c²), with c = c_F·r, where r is the product of the active primes. The scalar ζ is given by eq. `reflection-row-phase` (line 1832):

  ζ = −(i/81) ᾱ(c)² φ̂(h₀) κ̄_F ∏_{inactive}(1−q_p^{-1}) ∏_{active} χ_p(σ_p)^{-2} ω_{p,j_p},

with ω_{p,1} = τ⁻_{p,1} τ⁺_{p,3} χ_p(ε_p)^{-3}. Here:
- α(c) = c/|c| is the angle character (line 575), a Hecke character of infinite order;
- τ⁻_{p,1} = q_p^{-1/2} Σ_y χ_p(y) e(−y/p) is the normalised **sextic** Gauss sum at p.

**The paper flags the R-dependence itself.** Lemma `unmarked-reflected-block` (line 2957ff.) says:
- line 3017: "The scalar ζ_R of the reflection depends on R but not on μ";
- line 3044: it absorbs this scalar into the bounded row coefficient ρ(R) = 81ζ_R χ_R(x_F)³.

For mean squares this is harmless, which is why the paper never needs the phase. For a pairing with a fixed row vector ρ_6(m), it is essential.

### 1.2 Identification of ζ_m [R]

Take a squarefree row m. All its primes are active with j_p = 1. By eq. `reflection-cross-prime-phase` (line 2271), χ_p(σ_p)^{-2}ω_{p,1} = ξ_{p,1} χ_p(c/p)⁴, where ξ_{p,1} ∝ τ⁻_{p,1}τ⁺_{p,3} times symbols of fixed elements.

- **Prime factor.** At a prime, the Davenport–Hasse identity of eq. `signal` (line 965), γ₁γ₂ = μαG, together with γ₂³ = μα, gives

  ᾱ(p)² τ⁻_{p,1} τ⁺_{p,3} = (ray-class phase of p) · \overline{γ₂(p)}⁴.

- **Composite m.** The cross-prime phases χ_p(c/p)⁴ combine with the twisted multiplicativity of γ₂, where (χ_pχ_{p'})^{±2} appears, to leave only ±1 quadratic-reciprocity signs.

Hence, for squarefree m,

  **ζ_m = Φ(m) · \overline{γ₂(m)}⁴ = Φ(m) μ(m) \overline{α(m)} \overline{γ₂(m)},**

where Φ is a finite-order ray-class function, twisted-multiplicative through a fixed ±1 bicharacter. So ζ_m is, up to ray phases, the conjugate of the phase of the *cubic theta coefficient* at m. It is not a ray-class phase.

### 1.3 Numerical confirmation (`check_row_phase.py`) [R as computation]

The check uses all 6028 split primary primes with N(p) ≤ 60000, and computes w(p) := ᾱ(p)²τ⁻_{p,1}τ⁺_{p,3}.
- **w·γ₂⁴ is a class function.** In all 432 classes mod 72, w(p)·γ₂(p)⁴ is constant within the class, with zero exceptions.
- **w itself is not.** Within a class, w(p) takes as many distinct values as there are primes, spread over the circle:
  - the median within-class max |w − w₀| is 1.99;
  - |mean w| ≈ 0.17 over 20 primes.

So the row phase is exactly (ray constant) × \overline{γ₂(p)}⁴, as derived. By Heath-Brown–Patterson equidistribution of cubic Gauss sums, including twists by Grössencharacters, it is not a ray-class function for any modulus.

### 1.4 Consequence for (ii) [R that (V) does not apply; G for any repair]

For squarefree m the exact low-side theta pairing is

  ⟨ρ, B⟩ = Σ_terms Σ_μ β_μ Σ_m Ω(q_m/Q) ξ(m) (q_m/Q)^{-iv} · **ρ(m) ζ_m** · (λ⁴μ/m)₂ · V♯(q_μ Z / q_{c_F m}²).

The m-coefficient is ρ(m)ζ_m = ρ(m)·Φ(m)μ(m)ᾱ(m)γ̄₂(m). The analysis in `V_derivation.md` covers only ρ(m)·(μ/m)₂:
- the cusp multiplier of Lemma M;
- the d-sum of §5.2 that produces the conjugate cubic dual character;
- the "no main term for squarefree μ ≠ 1" statement;
- the main terms at square μ.

Multiplying by ᾱ(m)γ̄₂(m) destroys the automorphy in m:
- it is not a congruence twist;
- it is not a finite-order character;
- it cannot be absorbed into translates of θ_ρ.

So none of these steps applies to the actual sum. `V_derivation.md` §5.2 shows that ε(μ) is genuinely oscillating in μ, where it is harmless. The analogous m-side phase ζ_m was never computed.

**Heuristic [H] on whether a repair is plausible.** ρ(m)·\overline{(\text{cubic theta coefficient})(m)} is a Rankin–Selberg-type product θ₆ × θ̄₃. Even in the best case, where it behaved like a GL₂×GL₂ object (it does not, since the cover multipliers κ₆ and κ̄₃ do not cancel), a quadratic twist of conductor μ would have conductor ≍ N(μ)⁴. Then:
- the dual length would be L⁴/Q;
- at the certified geometry, L = Q²/Z = Z^{0.652} and Q = Z^{0.826}, so L⁴/Q = Z^{1.78} ≫ Q;
- so even an ideal functional equation gives no saving.

I see no route to a nontrivial bound for ⟨ρ, B⟩ beyond Cauchy–Schwarz.

**The two sides cannot be swapped.** The other order, Voronoi in m on the unreflected row Σ_x a(x)χ_x(m), is the coherent j = 1 twist with conductor q_x ≍ Z. Its dual length is Z²/Q > Q, which is worse. So both orders fail.

### 1.5 Impact on the numbers (`audit_lp_no_theta.py`, model `theta_split_model.py`)

| low-side Gram model | S ≈ 0 (slot-free) | best over slot supply S |
|---|---|---|
| paper (g = max(0, b/6, 2b − l_y + d)) | 0.894737 | 0.875092 (S = 0.2) |
| (i) only: P²/Y dropped, theta part by CS (g = max(0, b/6)) | 0.894737 | 0.875092 |
| (i)+(ii) as claimed (LSmarked, P²/Y dropped) | 0.866197 | 0.866197 |

Every improvement comes from step (ii). Step (i) does nothing, because the P²/Y term is not binding at the paper's optimum.

The conditional theorem of `REPORT.md` §1, "(LoA6)+(V) ⇒ |I| ≪ Z^{l_x/2+ε} when b/6 ≤ 2s(M)", is **false as an implication**: (V) as stated concerns Σ_m ρ_c(m)(μ/m)₂W, which is not the pairing that occurs.

A correct conditional statement would need a hypothesis (V′): the claimed saving for Σ_m ρ_c(m)ζ_m(λ⁴μ/m)₂W(q_m/Q) V♯(…), uniformly in μ. That is not a statement about θ₆ alone, and it has no known route [G].

The same objection applies to the slotted 0.871 certificate, whose marked reflected rows carry the same ζ, and to the 0.8633 variant.

---

## 2. Check 1: the low-side inequality and its LP translation

**Setup [R-].** The exact separation is eq. `low-separated` (line 3470). Put Q = q_{b*}XY ≍ Z^M. A^sig is the residue term of the contour shift, and A^rest is the rest. Then

  |I| ≤ (2π)^{-1}Q^{-1/2} Σ_σ ∫|Ŵ₀(iv)| (‖A^rest_{σ,v}‖ ‖B_σ‖ + |⟨A^sig_{σ,v}, B_σ⟩|) dv.

- The v-dependence enters only through (q_m/Q)^{-iv} and Ŵ₁(1/6+iv)(1+|v|)^{R'}. Both are integrable against Ŵ₀(iv).
- The σ-sum is finite.

**Rest term [R].** ‖A^rest‖² ≪ (Q/Y)Z^ε and ‖B‖² ≪ Z^{M+ε}, from the unmarked completed moment at ℓ′ = 0 with energy term 0 at M = 0.826. So the rest term is Q^{-1/2}(Q/Y)^{1/2}Q^{1/2} = X^{1/2}, i.e. g = 0.

**Theta term [R, arithmetic given the claimed bound].**
- The Cauchy–Schwarz scale is Q^{-1/2}Y^{-1/3}‖ρ‖‖B‖ ≍ Q^{5/12}Y^{-1/3} = X^{1/2}P^{1/12}, i.e. g = b/6. This uses ‖ρ‖² ≪ Q^{5/6+ε}, which follows from the bootstrap in `D_derivation.md` §6.2 [R-].
- The claimed bound Q^{5/12}(L + L²Q^{-2/3})^{1/2} + Q^{5/6}, with L = Z^{2M−1}, gives the ratio Z^{-s} with −2s = max(M−1, 7M/3−2, −M/6).
- So g = max(0, b/6 − 2s(M)), and the [2b − l_y + d] bracket is dropped under LoA6.

**The cap s ≤ M/12** comes from the main terms O(Q^{5/6+ε}), i.e. Q^{-1/12} relative to Q^{11/12}.
- In probe units these main terms are X^{1/3}Z^ε = Z^{C(5/6)}.
- Since b/6 − M/6 = −l_x/3 < 0, the cap never binds once g is floored at 0. It is harmless.
- Its derivation (`V_derivation.md` §6, constants c_𝔞 at the cusps, G5) is moot after F1.

**At the certified geometry:**
- b/6 − 2s(M) = 1076/15000 − 1090/15000 = −14/15000 < 0, so g = 0;
- the binding piece is Heath-Brown's (UV)^{2/3} term, 7M/3 − 2;
- the low side is 5/6 + l_x/6 = 8663/10000, a margin of 2·10⁻⁴ below σ0.

**Verdict on Check 1.** The translation is correct [R] given the claimed bound. The claimed bound itself fails (F1).

## 3. Check 2: circularity

**No circularity in the logic [R].**
- The continuation criterion (line 402ff.) needs an unconditional upper bound |J| ≪ Z^{C(σ0)+ω}.
- The ingredients claimed for (ii) use no information about zeros of L(s,η): Voronoi, Heath-Brown's large sieve, LoA6 and the paper's Gram bound.
- The high side alone uses β* (Lemma plain with κ = 2β*−1, which is irrelevant at z = 0, and the bin ceiling).
- So if a valid bound for ⟨ρ, B⟩ existed, the argument would not be circular.

**One misstatement [R that it is wrong; cosmetic for the logic].** `V_derivation.md` line 37 and `REPORT.md` say the square-μ main terms are "exactly the signal". This cannot be right:
- those terms are unconditionally O(Q^{5/6+ε}), i.e. Z^{C(5/6)} in probe units;
- under the contradiction hypothesis the signal f_η is ≍ Z^{C(β*)} with β* > σ0 > 5/6.

The signal actually arises from the coherent j = 1 main terms of the *unreflected* pairing. Concretely, Σ_x a(x)V(x/Z)·[φ(x)/τ(χ̄_x)]·cQ^{5/6}N(x)^{-5/6}, where a(x)γ₁(x) ∝ μ(x)G(x), gives the 1/L sum.

In the (V) picture it would therefore have to sit inside the large-sieve part. Logically that is acceptable, because the bound is an upper bound. But it shows that the (V) bookkeeping of "where the signal lives" was not checked against the unreflected form. With F1 corrected, the reflected-side main-term analysis does not apply at all.

## 4. Check 3: the reflected form of B for the slot-free base probe

**Reflected length [R].** L = q_c²/Z with c = c_F·F_act·R. For squarefree m, q_c ≍ q_m ≍ Q, so L ≍ Q²/Z. For general m, L ≤ C·q_{rad m}²/Z. The marked-length formula is irrelevant at ℓ = 0.

**Twist for squarefree m [R].** The column factor is ∏_{p|m} χ_p(λ⁴μ)³ = (λ⁴μ/m)₂. **But it comes multiplied by ζ_m. This is F1.**

**Non-squarefree rows [G].** For p² | m, the column factors are χ_p^{-j-2} (j = 2, 3, 5), a Ramanujan factor (j = 4), or active/inactive j = 0 factors.
- The paper's own treatment (Definition `unmarked-reflected-block`) freezes F and sums over squarefree R coprime to F.
- Voronoi in R would need m ↦ ρ(FR)·1_{R sqfree, (R,F)=1}·ζ_{FR}. That brings in:
  - the sextic coefficients at non-squarefree arguments, which are not determined (Kazhdan–Patterson);
  - squarefree sieving with unbounded moduli;
  - ζ again.
- None of this is addressed. `V_derivation.md` lists it under "(P) items not addressed".

**Base probe versus marked row [R].** In the slot-free case the row is the unmarked B_{m,σ}(Z), so the object itself is correct. The error is only in its reflected description.

## 5. Check 4: the exact certificate `exact_cert_theta.py`

**Re-run [R as computation].** Command:

`python exact_cert_theta.py theta=1 dropq=1 lx=989/5000 ly=3141/5000 ell=0 sigma0=1733/2000`

Output (33 s):

```
h = 4011/5000, supply S = 0, low = 8663/10000 ~ 0.866300
PASS low side <= sigma0   [0.8663]
PASS geometry: M-2ell>0, ly-ell>0, lx-ell>0
PASS small rows exponent < 0   [-0.162452]
PASS floor bin E(h+zeta) < 0   [-0.000275667]
PASS high side E<0 on [D0,DMAX]x[0,1/2], all d in [d_min,h+zeta] (3055 boxes)   [-0.000503433]
ALL CHECKS PASS
```

**What it encodes.**

*`gram()`.* With theta = 1 and dropq = 1, it returns max(0, b/6 + max(M′+ℓ′−1, 7M′/3+2ℓ′−2, −M′/6)). At ℓ = 0 this is exactly g = max(0, b/6 − 2s(M)) [R].

*`energy()`.* Its value is 0 at M = 0.826, ℓ′ = 0. The default wE = 4/7 is a hypothetical relaxation, but it is inert here.

*High side.*
- It is the paper's eq. `stage-high-exponent` (line 6092) and its specialisation eq. `common-high-exponent` (line 15726), with general (h, l_y, ℓ, σ0), z₀ = 17/50, d_min = 1/100 and floor δ₀ = 1/50.
- It requires E < 0, which is conservative, since E < Δ = β* − σ0 suffices.
- R = R_t is evaluated exactly at the box corner of smallest δ and q, which is valid because R is nonincreasing in both. It is minimised over t ∈ {t_float, 1, 3/2}, which is valid because the count holds for each t.
- At S = 0 all capacities vanish. R_t reduces to max_r min(max(1, (1+5r)/6) − δr, 1 − δr [r < 1], 1 − 2δ(t−r)), which is R₀(t) = max(L(t), 1 − 2δt/3).
- Margins are tiny but strict:
  - floor bin: 2.8·10⁻⁴;
  - high bins: 5.0·10⁻⁴;
  - low side: 2·10⁻⁴.

**Zero-slot row counts [R-].** These close the concern in `IDEA1_NOTES.md` §6.
- Prop. `detector-witness` (line 4510) holds for every t ∈ [1, 3/2] and does not involve slots.
- The no-slot inverse count, eq. `no-slot-inverse-count` (line 15282) and the Remark after the proposition, holds for all r.
- Lemma `plain` case z = 0 (line 12531) has "no restriction on n₁, n₂", so 1 − 2δm holds for every m ≤ 1/2 + O(ε).
- The paper's own t = 1 computation with x = 0 (line ~15400), repeated at general t, gives R₀(t) = max(L(t), 1 − 2δt/3).
- No Δ/4 capacity-comparison loss arises at zero slots.

**What the certificate cannot see.**
1. F1 (garbage in, garbage out).
2. **Scope of the paper's high-side lemmas [G/R-].** Several lemmas are stated only for σ0 ∈ [7/8, 1):
   - `fixed-bin-contour` (line 5811);
   - `bin-contour-accounting` (line 6024);
   - `principal-residue-interface` (line 6178);
   - `outer-row-tails` (line 6345).

   The certificate uses σ0 = 0.8665 < 7/8. The exponent algebra in their proofs is general in (l_x, l_y, ℓ, σ0), via Definition `stage-high-data` (line 5583) with C(s) = s + l_x/2 − 1 + h/6. So this is probably routine, but nobody has checked it. The same applies to the ℓ = 0 principal-signal extraction for the uncompensated probe at l_x ≈ 0.2 (`REPORT.md` §3 item (P) already lists this).
3. The detector-floor variant δ₀ = 1/1000 (0.8633) changes the detector proof, which uses δ ≥ 1/50 and σ ≥ 51/100 (lines 4504, ~4580). It is unverified, and outside this audit's main scope.

## 6. Step (i): LoA6 ⇒ ‖A^rest‖² ≪ (Q/Y)Z^ε

**Contour shift [R].** Ŵ₁ is entire with rapid vertical decay. With polynomial vertical growth of D [R-], Cauchy–Schwarz in t gives Y^{-1+2ε}·Σ_m|D(1/2+ε+it;m)|². The residue at w = 2/3 is exactly Y^{-1/3}Ŵ₁(1/6+iv)ρ_c [R].

**LoA6 [R-, conditional on the Kubota functional equation with explicit cusp data, G-a].**
- **AFE.** The approximate functional equation has length (Q(1+|t|)^c)^{1/2}.
- **Large sieve.** `D_derivation.md` §6.2 uses the planar additive large sieve with points x/s. The fractions are distinct reduced fractions with separation ≥ X^{-1}, which gives (Q + X²)Σ|b_s|²φ(s)/q_s ≪ Q^{1+ε}. This is correct [R]; it needs no multiplicative sieve and removes G-b.
- **Residue term.** The residue term X^{1/6}ρ is handled by the bootstrap ‖ρ‖² ≪ Q^{5/6+ε}, which uses the paper's own Gram bound at Y = Q^{1/2+δ} [R-].

**Minor error in `LoA6.md` line 38 (harmless).** The normaliser written there is wrong.
- Written: ζ_K(6w−2)/ζ_K(6w−1).
- Correct: ζ_K(6w−3)/ζ_K(6w−2). Its pole is at 6w − 3 = 1, i.e. w = 2/3, and the zeros of the denominator lie in Re w < 1/2. This is consistent with "ζ(6s−5), s = w+1/2" in `IDEA1_NOTES.md` §5.
- The conclusion (only the w = 2/3 pole in Re w > 1/2) is unaffected.

**Ray twists.** The probe's c(s) is a finite combination of ray characters. These depend on the target via b_*, which is allowed: only exponents need to be uniform.

**Status.** (i) is plausible and close to standard [R- / G-a], but it is worthless without (ii): see the table in §1.5.

## 7. Summary of statuses

| step | claim | status |
|---|---|---|
| (i) | ‖A^rest‖² ≪ (Q/Y)Z^ε from the Kubota FE via LoA6 | [R-] modulo the FE with cusp data [G-a] |
| (i) | P²/Y and P^{1/6} both removed from A^rest | [R-] |
| (ii) | the reflected row is ζ_m·Σβ_μ(λ⁴μ/m)₂V♯, with ζ_m = Φ(m)\overline{γ₂(m)}⁴ | [R] (paper eqs. 1832, 2271; line 3017; numerics) |
| (ii) | (V) applies to ⟨ρ, B⟩ | **false**: (V) concerns ρ(m)(μ/m)₂, but the pairing involves ρ(m)ζ_m(μ/m)₂ |
| (ii) | the saving s(M) for the actual pairing | **[G], no known route**; heuristically a conductor ≍ N(μ)⁴ gives none |
| (ii) | main terms at square μ "are the signal" | wrong: they are O(Z^{C(5/6)}) unconditionally; harmless |
| (ii) | (V) internal gaps G3, G4, Lemma M for ray twists, non-squarefree rows | [G] / [R-] (moot after F1) |
| (iii) | LP translation g = max(0, b/6 − 2s(M)), with the cap and the v/σ bookkeeping | [R] arithmetic; [R-] analysis |
| (iii) | zero-slot counts R₀(t) = max(L(t), 1 − 2δt/3) | [R-] (paper's own argument at general t) |
| (iii) | high-side lemmas at σ0 < 7/8, ℓ = 0, l_x ≈ 0.198 | [G]/[R-] (paper states σ0 ∈ [7/8, 1)) |
| (iii) | exact certificate at σ0 = 1733/2000 | re-run: ALL CHECKS PASS; it encodes exactly the claimed exponents |
| — | without (ii) | θ = 0.8947 (S = 0), 0.8751 (best S): no improvement on Liu |

---

## Referee report

*Independent skeptical referee, 9 Oct 2026. Goal: refute the audit's main conclusions if possible. Labels as above. New supporting scripts in this directory: `referee_composite_phase.py`, `referee_fine_S.py`.*

### Verdict on the audit: NOT REFUTED (agree: the chain is flawed at (ii))

I could not find an error that changes any of the audit's conclusions. F1 is correct, and it is fatal for (ii) as stated. Every headline number (0.8665, 0.871, 0.8633) depends on (ii). Without (ii) the LP returns exactly the paper/Liu value. I found four minor inaccuracies in the audit (R4 below); none of them affects the verdict.

### R1. F1 re-derived independently from the paper [R]

**Paper text checked.**
- eq. `reflection` (line 1763) and eq. `reflection-row-phase` (line 1832): ζ contains ᾱ(c)² ∏_active χ_p(σ_p)^{-2} ω_{p,j_p}, with ω_{p,1} = τ⁻_{p,1}τ⁺_{p,3}χ_p(ε_p)^{-3}.
- eq. `reflection-cross-prime-phase` (line 2271): χ_p(σ_p)^{-2}ω_{p,j} = ξ_{p,j}χ_p(c/p)^{2j+2}.
- Line 3017 ("ζ_R … depends on R but not on μ") and line 3044 (ρ(R) = 81ζ_Rχ_R(x_F)³). Both quotes are accurate.

**The row really is the reflection input with moving primes = primes of m.** In the probe (lines 3398–3420), B_{m,σ} uses T(t+1/2, ν_σθ·χ_•(m)) with Ψ(a) = ν_σθ(a)χ_a(m). By sextic reciprocity χ_a(m) = R(a,m)χ_m(a), where R is a finite bicharacter, so 𝒫 = {p | m} and j_p = v_p(m) mod 6.

**Algebra for squarefree m, c = c_F·m.**
- From τ⁻_{p,1} = χ_p(−1)γ₁(p), τ⁺_{p,3} = γ₃(p) and eq. `signal` (γ₁γ₂ = μαG, G = χ̄_p(4)γ₃, γ₃² = χ_p(−1)), I get τ⁻_{p,1}τ⁺_{p,3} = −α(p)χ̄_p(4)γ̄₂(p).
- The cross-prime factors are ∏_{p|m} χ_p(m/p)⁴. The twisted multiplicativity γ₂(ab) = χ_a(b)²χ_b(a)²γ₂(a)γ₂(b) (line 1043) gives ∏_{p|m}γ̄₂(p) = γ̄₂(m)·∏_{pairs}χ_p(p′)²χ_{p′}(p)². The two pair products cancel **exactly**.
- Hence ζ_m = Φ(m)·μ(m)ᾱ(m)γ̄₂(m) = Φ(m)·γ̄₂(m)⁴.
- Here Φ(m) = (const)·κ̄_F(r mod M²)·χ_m(a₀) for a fixed element a₀ built from λ, 2, −1 and c_F. This is a ray-class function, and in fact multiplicative.

**Numerics.**
- I re-ran `check_row_phase.py` and reproduced its output exactly: 6028 split primes, 432 classes mod 72, 0 exceptions, median within-class spread of w equal to 1.993.
- New composite check, `referee_composite_phase.py`: squarefree m = p₁p₂ (6612 rows), with γ₂(m) computed directly as a Gauss sum mod m and the full prime factor including χ_p(σ_p)^{-2}χ_p(ε_p)^{-3}.
  - W(m)·γ₂(m)⁴ is constant on all 432 classes of m mod 72, with 0 exceptions.
  - W(m) itself has median within-class spread 1.995.
- So the cross-prime cancellation is confirmed.

**Possible rescue checked and closed [R].** Could another m-dependent unimodular factor in eq. `low-separated` cancel ζ_m? No.
- ξ(m) is a fixed residue character mod b_* (line 3350ff.).
- The A-side carries only g_{χ_s}(s,−m), whose modulus is s, not m.
- (q_m/Q)^{-iv} and Ω are norm functions.
- Nothing in the separation absorbs γ̄₂(m)ᾱ(m).

### R2. F1 is fatal for (ii), not a bookkeeping slip [R for "(V) does not apply"; H for "no repair"]

Without structure in m ↦ ρ(m)ζ_m, the best available bound is Cauchy–Schwarz.
- The quadratic large sieve in the column direction gives Σ_{μ≲L}|Σ_m c_mχ_μ(m)|² ≪ (L+Q)‖c‖².
- Since L = Q²/Z < Q, this is Q‖c‖², which reproduces ‖ρ‖‖B‖ exactly.
- Any saving would need a new hypothesis (V′) on correlations of θ₆ with θ̄₃ (a Grössen-twisted cubic theta) against χ_μ. I agree that (V′) is a genuine [G] with no known route.

**Structural reason [H].** The probe's signal comes from the Davenport–Hasse collapse γ₁γ₂ = μαG between the cubic-theta row and the sextic data. ζ_m is the reflected-side incarnation of the same cubic Gauss sum. A method that sees only ρ paired with quadratic characters therefore throws away exactly the θ₆ × θ̄₃ correlation in which the signal lives.

### R3. The audit's other checks, re-verified

- **LP without (ii).** I re-ran `audit_lp_no_theta.py` and reproduced its output. `referee_fine_S.py` refines the slot-supply optimisation for "(i) only" (P²/Y dropped, theta part by Cauchy–Schwarz):
  - θ = 0.874957 at S = 0.2055, (l_x, l_y, ℓ) = (0.3549, 0.4783, 0.1668);
  - this is identical to the paper/Liu value, with or without the P²/Y piece;
  - the P²/Y piece is slack by at least 0.065 at the optimum (2b − l_y + d = −0.231 at d = 0 and −0.065 at d = ℓ, against b/6 = 0.021).
  
  So "(i) alone gains nothing" holds exactly. The audit's 0.875092 is grid resolution, as it says.
- **Certificate.** I re-ran `exact_cert_theta.py theta=1 dropq=1 lx=989/5000 ly=3141/5000 ell=0 sigma0=1733/2000`: ALL CHECKS PASS, with the same margins (low 0.8663, floor −2.757e-4, high −5.034e-4, 3055 boxes).
  - gram() and energy() encode exactly what the audit says.
  - energy() = 0 at ℓ = 0 for both wE = 0 and wE = 4/7.
  - The audit's arithmetic b/6 − 2s(M) = −14/15000 is correct.
- **Zero-slot counts.** I checked the proof of Prop. `detector-counts` (lines 15226–15320). The zero-slot branches are eq. `no-slot-inverse-count` (for r ≥ 1 it gives L(t)), zero-slot plain for m ≥ 1/2, and the unweighted 1 − δr and 1 − 2δm otherwise. Combined with Prop. `detector-witness` (line 4510, valid for every t ∈ [1, 3/2] and slot-free), they give R₀(t) = max(L(t), 1 − 2δt/3). Agree [R-].
- **Scope of the high-side lemmas.** Lines 5811, 6024, 6178 and 6345 do state σ0 ∈ [7/8, 1). Agree [G/R-].

### R4. Minor corrections to the audit (none changes a conclusion)

1. **"False as an implication" is overstated.** What is shown is that the implication is *unsupported*: (V) concerns Σ_mρ(m)(μ/m)₂W, while the pairing needs ρ(m)ζ_m. Nothing shows that (LoA6)+(V) hold while the conclusion fails. The correct wording is "the stated conditional theorem does not follow; its proof uses the wrong sequence".
2. **Φ is exactly multiplicative.** The cross-prime symbols cancel identically (R1), so no ±1 bicharacter survives. This is harmless.
3. **The swapped-order remark (§1.4) is imprecise [H].** For the coherent j = 1 twist, the dual character is principal: the a-sum is a Ramanujan sum c_x(n), and by eq. `signal` a(x)/τ(χ̄_x) ∝ μ(x)·(ray). So Voronoi in m on the unreflected row does not produce an oscillatory dual of length Z²/Q. Instead it produces Möbius sums over x ≍ Z, which is the 1/L(s,η) structure of the paper's Poisson/high-side representation. It cannot give an unconditional low bound, so the conclusion "neither order rescues (ii)" stands; only the stated reason ("dual length Z²/Q > Q") is too naive.
4. **`check_row_phase.py` uses a simplified prime factor.** Its w(p) omits χ_p(σ_p)^{-2}χ_p(ε_p)^{-3}. For a prime row these are ray-class, so the omission is harmless. The composite check above includes them.

### Bottom line

The audit stands. The reflected base-probe row carries the m-dependent root number ζ_m = Φ(m)·γ̄₂(m)⁴ (paper, exact; numerically confirmed for prime and composite m). The (V)-based saving in (ii) therefore does not apply to the pairing that actually occurs.

Without (ii) the LP gives exactly 0.874957, so no improvement on Liu is established by this chain. The status of 0.8665, 0.871 and 0.8633 is: **not established**. They would require the new, unsupported hypothesis (V′).
