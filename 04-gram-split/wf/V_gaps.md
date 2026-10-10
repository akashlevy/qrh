# Closing the gaps in (V) and (P) for the theta pairing

*Working notes, 9 Oct 2026. They build on `REPORT.md`, `IDEA1_NOTES.md`, `V_derivation.md`, `D_derivation.md` and `LoA6.md`. Scripts are in this directory:*
- *`check_rowscalar.py`, `check_rowscalar2.py`: numerical checks for §1;*
- *`lp_after_F1.py`, `lp_const_saving.py`: LP consequences.*

Labels:
- **[R]** rigorous, given the cited results;
- **[R-]** routine but not written out;
- **[S]** taken from a secondary statement and not re-derived;
- **[H]** heuristic;
- **[G]** gap;
- **[N]** numerical evidence only.

---

## 0. Verdict

**F1: a new, fatal gap in step (ii).** Lemma (V) is applied to the wrong sum.
- In the paper's completed reflection, the reflected row carries a row scalar ζ_m (`paper78.tex` eq. `reflection-row-phase`).
- Earlier notes treated ζ_m as "a fixed ray phase". It is not.
- For a squarefree good row R, Lemma Z (§1.2) gives
  ζ_R = κ₀ · ψ₁(R) · \overline{γ₂(R)}⁴,
  where κ₀ depends only on finitely many classes and ψ₁ is a ray character. Here γ₂(R) is the normalised cubic Gauss sum.
- Because γ₂³ = μα, the factor \overline{γ₂}⁴ equals μᾱ·\overline{γ₂}, which is a Hecke character times a conjugate cubic Gauss sum. It is not a Hecke character.
- Status [R], and checked numerically to 10⁻¹⁵ on 2248 primes and on products of two primes.
- Consequence. The inner row sum of the theta pairing is
  F♮(μ) = Σ_R ρ(R) · \overline{γ₂(R)}⁴ · ψ(R) · (μ/R)₂ · w(R).
  It is not Σ_R ρ(R)(μ/R)₂ w(R), which is the sum that (V) treats.
- The sequence ρ·\overline{γ₂}⁴ is a product of a θ₆ coefficient and a θ₃ coefficient. No Voronoi formula is known for it, and a degree-4 functional equation would not help (§1.5).
- So claim (ii), and with it the certified value 0.8665, is **not supported** by (V) + Heath-Brown.

**What (i) alone gives.** If the theta pairing gets no power saving, the LP returns exactly **0.874957**, Liu's value. This holds with the P²/Y Gram term removed by (LoA6), because that term is inactive at Liu's optimum. The slot-free variant gives 0.8947. So (i) alone gives nothing.

**Any saving would help.** A uniform power saving Z^{−s₀} in ⟨ρ,B⟩ over Cauchy–Schwarz would give the following (`lp_const_saving.py`):

| s₀ | 0.002 | 0.005 | 0.01 | 0.02 | 0.04 |
|---|---|---|---|---|---|
| θ | 0.87333 | 0.87090 | 0.86683 | 0.86151 | 0.85323 |

The idea is therefore still alive, but it needs a new input: a non-trivial bound for the pairing of ρ·\overline{γ₂}⁴ with quadratic characters. Section 1.6 says what that input is and why the unreflected route cannot supply it.

**The requested items.** Each is done relative to (V) as stated, i.e. in the "ζ-free model". They would carry over to any repair that goes through θ₆-Voronoi.

| item | status |
|---|---|
| **G3** (mean square at the cusps 1/μ₀) | **Closed** relative to (LoA6 at all cusps) + (G1 at all cusps) + the paper's Gram bound for ray-twisted coefficients. Status [R-]. No Rankin–Selberg is needed. Section 2. |
| **G4** (non-squarefree dual indices) | **Closed** relative to G3, de Faveri Thm 1.1 [V, statement checked in arXiv:2610.04045] and the Kazhdan–Patterson local relations at good primes [S]. The needed bound is (G4*): Σ_{N(n)≤N, c³∣n} \|τ_𝔞(n)\|² ≪ N^{1+ε} N(c)^{−2}. The KP relations give exactly this, with equality only in the case v_p = 3. The normalisation of the KP relations is [S] and must be checked: a loss N(p)^δ in \|τ(mp³)\|/\|τ(mp)\| is tolerable only for δ ≤ 0.318 at the certified geometry. Section 3. |
| **Non-squarefree rows** | [R]: no column produces a coherent twist in the free row variable R. The R-twist is quadratic at column primes and has sextic exponent 2j_p ∈ {0,2,4} at the frozen active primes (a new cross term found here). The Ramanujan case j_p = 4 is purely column-side. Block counting: dyads with powerful norm Z^O ≥ Z^{4s(M)} are handled by Cauchy–Schwarz alone ([R] given G3 and the paper's moment lemma). The remaining dyads cost no exponent in the ζ-free model [R-]. Section 4. |
| **Non-squarefree columns** | [R-], with log losses only. Section 5. |
| **v-integral** | [R]. It disappears exactly: the theta part is Σ_m K(q_m/Q) ξ(m) ρ_σ(m) B_{m,σ}(Z), with K(x) = ∫ W₁(t) t^{1/6} W₀(x/t) dt/t. Section 6. |
| **Ray class s ∈ σ** | [R-]. It is a finite combination of Kubota series at the cusps of Γ(72f). All lemmas apply cusp by cusp. Section 7. |

---

## 1. F1: the row scalar of the completed reflection

### 1.1 The paper's formula

The low side pairs A_{m,σ,v}(Y) with B_{m,σ}(Z) = Σ_θ a_θ 𝒞_{V_G}(Z; m, ν_σθ) (`paper78.tex` ≈3440–3470).
- 𝒞 is the completed sum of Prop. `completed-reflection`, with Ψ(A) = ν(A) χ_A(m).
- Lemma `reflection-row-sectors` converts χ_A(m) into ∏_{p∣m} χ_p(A)^{j_p}, with j_p = v_p(m) mod 6.

The reflection eq. `reflection` then gives

  𝒞(m) = Σ_terms ζ · Σ_μ (d(μ) α(μ) ϑ / √q_μ) · ∏_{p active} B_p(λ⁴μ) · V♯(q_μ X / q_c²),   c = c_F r,

where r is the product of the active primes, and

  ζ = −(i/81) · ᾱ(c)² · φ̂(h₀) · κ̄_F · ∏_{p inactive} (1 − q_p^{−1}) · ∏_{p active} χ_p(σ_p)^{−2} ω_{p,j_p}

(eq. `reflection-row-phase`), with:
- σ_p = λ²c/p and ε_p = −((λ³c/p)σ_p)^{−1};
- ω_{p,1} = τ⁻_{p,1} τ⁺_{p,3} χ_p(ε_p)^{−3}.

The paper says explicitly that all moving-prime dependence outside B_p sits in ζ. It needs only |ζ| ≤ 1/81, because it uses only mean squares (Lemma `unmarked-reflected-block` absorbs ζ into a bounded row coefficient ρ(R)).

**In a pairing with a fixed row vector, ζ_m is part of the row sequence.**

### 1.2 Lemma Z: evaluation for squarefree rows [R]

**Statement.** Let R be a squarefree product of good primes (all active, j_p = 1), with no frozen part. Then

  ζ_R = κ₀ · ψ₁(R) · \overline{γ₂(R)}⁴ = κ₀ · ψ₁(R) · μ(R) · ᾱ(R) · \overline{γ₂(R)},

where:
- κ₀ = −(i/81) φ̂(h₀) κ̄_F ᾱ(c_F)². It depends only on h₀ and on the class of R mod M².
- ψ₁(R) = χ_R(λ)^{−1} χ_R(c_F)^{−2} \overline{χ_R(4)} χ_R(−1) is a finite-ray character of R with conductor supported on S.

**Proof.**
1. **Local factor at p ∣ R.**
   - Since c/p is a unit mod p and χ_p⁶ = 1 on units, χ_p(ε_p)^{−3} = χ_p(−λ⁵(c/p)²)³ = χ_p(−1) χ_p(λ)³.
   - χ_p(σ_p)^{−2} = χ_p(λ)^{−4} χ_p(c_F)^{−2} χ_p(R/p)^{−2}.
   - τ⁻_{p,1} = χ_p(−1) γ₁(p) and τ⁺_{p,3} = γ₃(p).
   - So χ_p(σ_p)^{−2} ω_{p,1} = χ_p(λ)^{−1} χ_p(c_F)^{−2} χ_p(R/p)^{−2} γ₁(p) γ₃(p).
2. **Twisted multiplicativity.** The CRT computation gives γ_j(ab) = (χ_a(b) χ_b(a))^j γ_j(a) γ_j(b), hence
   ∏_{p∣R} γ₁(p) γ₃(p) = γ₁(R) γ₃(R) ∏_{p≠p′} χ_p(p′)^{−4}.
   Together with ∏_p χ_p(R/p)^{−2} = ∏_{p≠p′} χ_p(p′)^{−2}, the cross terms combine to ∏ χ_p(p′)^{−6} = 1.
3. **Signal identity.** Eq. `signal` gives γ₁γ₂ = μαG with G = \overline{χ(4)} γ₃. Eq. `fixed-gauss-phases` gives γ₃(p)² = χ_p(−1), so γ₃(R)² = χ_R(−1). Hence
   γ₁(R) γ₃(R) = \overline{χ_R(4)} χ_R(−1) γ₂(R)².
4. **The α factor.** α is multiplicative on squarefree arguments, since γ₂(ab)³ = γ₂(a)³ γ₂(b)³. So ᾱ(c)² = ᾱ(c_F)² ᾱ(R)².
5. **Conclusion.** Finally ᾱ(R)² γ₂(R)² = \overline{γ₂(R)}⁶ γ₂(R)² = \overline{γ₂(R)}⁴. ∎

**Numerical checks [N].**
- `check_rowscalar.py`: γ₁(p)γ₃(p) = \overline{χ_p(4)} χ_p(−1) γ₂(p)² holds to 2.6·10⁻¹⁵ on all 2248 split primes of norm ≤ 2·10⁴.
- Same script: γ₁γ₃·α^k is *not* a function of p mod 72 for any k = 0,…,5. The median phase spread within a class is about 4 rad, against 0 for the ray control.
- `check_rowscalar2.py`: step 2 and step 3 hold to 3.5·10⁻¹⁵ for 60 products R = π₁π₂.

### 1.3 General rows

Rows have the form m = (frozen k) · R, with frozen active set F_act and exponents j_p. Repeat the computation, including the R-dependence of σ_p and ε_p at the frozen primes:
- p ∈ F_act with j_p ∉ {0, 4} contributes χ_p(R)^{−2+2(j_p+2)};
- for j_p = 4 the contribution is χ_p(R)^{−2};
- R's own primes contribute χ_R(F_act)^{−2} = ∏_{p∈F_act} χ_p(R)^{−2}, by cubic reciprocity for primary elements.

Hence

  ζ_{kR} = ζ^{fr}(k, class data) · ψ₁(R) · \overline{γ₂(R)}⁴ · ∏_{p∈F_act} χ_p(R)^{2j_p}.   [R]

The new cross factor ∏_p (R/p)₃^{j_p} is a cubic character of R modulo rad(F_act). It matters in §4.

### 1.4 Consequence for (V) and for claim (ii)

Use the notation of §6 for the v-integral. The theta pairing for the squarefree-row block is

  T = κ₀ Σ_μ β_μ Σ_R a♮_R (μ/R)₂ V♯(q_μ X / (q_{c_F}² q_R²)),
  a♮_R = K(q_R/Q) ξ(R) ρ_σ(−R) ψ₁(R) \overline{γ₂(R)}⁴.

Lemma (V) is about Σ_R ρ(R)(μ/R)₂ W(q_R/Q). It therefore bounds the different pairing ⟨ρ, ζ̄B⟩, not ⟨ρ, B⟩.

The sequence R ↦ \overline{γ₂(R)}⁴ = μᾱ·\overline{γ₂} has these properties:
- At primes it is equidistributed, by Heath-Brown–Patterson [S] (numerically, the within-class spread above).
- It is twisted-multiplicative only up to Gauss-sum cross terms.
- Multiplying θ₆ coefficients by it is not a Hecke twist. Lemma M, the Kubota multiplier computation and the Voronoi formula of `V_derivation.md` §5 therefore do not apply.

**F1 is fatal for the derivation of (ii) as written [R].** The value 0.8665 remains conditional on an input that has not been formulated, let alone derived.

### 1.5 Would a functional equation for a♮ help? [H]

Suppose R ↦ ρ(R)\overline{γ₂(R)} were the coefficient sequence of a degree-4 object, a θ₆ ⊗ \bar{θ}₃ Rankin–Selberg over K.
- Its quadratic twists would have conductor ≍ N(μ)⁴.
- The dual length would be L⁴/Q = Z^{7M−4}, with L = Q²/Z and Q = Z^M.
- At M = 0.826 this is Z^{1.78}, against Q = Z^{0.83}.

So even the most optimistic Voronoi picture gives no saving. Only a "random-sequence" model would give one, and that is unprovable here.

### 1.6 The unreflected viewpoint

Why ζ_m cannot be wished away [R-, Voronoi as in `V_derivation.md` §5].

Without reflecting B:

  T = Σ_A c_A q_A^{−1/2} V_G(q_A/Z) Σ_m ρ(m) χ_A(m) K(q_m/Q) ξ(m),   c_A = γ₂(A) ᾱ(A) Ḡ(A) ν(A).

**The m-sum is coherent.** The twist χ_A(m)¹ is the coherent case j = 1 (`V_derivation.md` §3). Voronoi gives a main term (a Ramanujan-type d-sum φ(A)) plus a dual sum with Ramanujan sums c_A(n), of length q_A²/Q ≍ Z^{2−M}.

**The A-side Gauss sum cancels exactly.** Since 1/γ₅(A) = χ_A(−1)γ₁(A),

  c_A / τ(χ̄_A) = q_A^{−1/2} χ_A(−1) γ₁γ₂ᾱḠν(A) = q_A^{−1/2} μ(A) χ_A(−1) ν(A).   [R]

**The two pieces of T.** So T = 𝒮 + T′, where:
- 𝒮 = Q^{5/6} Σ_A μ(A) ν′(A) (φ(A)/q_A) q_A^{−5/6} c_{A₀} V_G(q_A/Z) is exactly the zero-detecting signal (cf. `IDEA1_NOTES` §2);
- T′ = Q^{11/12} Σ_A ν′(A) q_A^{−2} V_G Σ_n τ_{A₀}(n) μ(A) c_A(n) Ω♯(·), with μ(A)c_A(n) = μ(A)² when (A,n) = 1.

**Why this route is useless.**
- The best unconditional bound I can get for T′ is worse than the Cauchy–Schwarz scale when M < 1. The dual is longer than Q, the A-sum does not oscillate, and the squarefree-counting error is Z^{1/2} q_{rad n}^{1/2}. This gives T′ ≪ Q^{−7/12} Z^{3/2+ε}, which exceeds Q^{11/12} when M < 1.
- More importantly, writing T = 𝒮 + T′ only re-expresses the signal. It gives no unconditional upper bound for I.

**The structural point [H, but suggestive].** In the unreflected form, the Gauss sum γ₂(A) is converted into μ(A) by the coherent Voronoi step. In the reflected form, the Möbius information can only be carried by the row sequence, and it reappears as \overline{γ₂(R)}⁴. Any unconditional bound for T with a power saving is therefore a genuinely arithmetic statement about ρ·\overline{γ₂}⁴ against quadratic characters.

### 1.7 Possible repairs, none established [H]

**(a) A twisted double Dirichlet series.** Consider

  Z(w, z; μ) = Σ_{s,m} g̃₆(s,−m) \overline{γ₂(m)}⁴ ψ(m) (μ/m)₂ q_s^{−w} q_m^{−z},

an A₂-type series mixing sextic and cubic Gauss sums.
- The closest precedent is Brubaker–Friedberg–Hoffstein (Invent. 2005), whose Z₆ and Z₃ link cubic twists, 6-fold Kubota series and π × cubic theta.
- Its residue at w = 2/3 is the L-series of a♮ ⊗ (μ/·)₂.
- Even if functional equations in z exist, §1.5 suggests the dual is long.

**(b) Change the probe so that A's theta component is ρ·ζ̄.** This changes the high side. It is not explored here.

**(c) What is actually needed.** The precise hypothesis is

> (V♮): |Σ_μ β_μ Σ_R a♮_R (μ/R)₂ V♯| ≪ ‖ρ‖ ‖B‖ Z^{−s₀}, uniformly over the blocks of §4.

The LP table in §0 shows that any s₀ > 0 helps.

---

## 2. G3: mean squares of θ₆ coefficients at the cusps (closed relative to LoA6) [R-]

**Lemma G3.** Let f be the conductor of the finitely many ray characters used (§7). For every cusp 𝔞 of Γ(72f) (in particular 𝔞 = 1/μ₀) and every N ≥ 1,

  Σ_{n∈Λ_𝔞, q_n≤N} |τ_𝔞(n)|² ≪_ε N^{1+ε},   τ_𝔞(n) = q_n^{1/12} ρ_𝔞(n).

This covers all n, squarefree or not, including S-parts.

**Inputs.**
1. **(LoA6)_𝔞 for every cusp, with the AFE residue term.** This is `LoA6.md` with gaps G-a and G-c. Since the Kubota functional equation mixes cusps, LoA6 has to be proved for all cusps at once anyway; it is not an extra assumption.
2. **(G1)_𝔞.** ρ_𝔞(m) is (a fixed constant times) Res_{w=2/3} D^𝔞(w;m).
3. **The paper's Prop. `probe-gram`** for coefficient functions that are finite combinations of ray characters, and for moduli that are fixed units times primary elements.

**2.1 The coefficients of E at 𝔞 = 1/μ₀** [R-].
- Put σ_𝔞 = L(μ₀). Since Γ(72) is normal in SL₂(O), σ_𝔞^{−1}Γ(72)σ_𝔞 = Γ(72), and the cusp has width 72.
- Bottom rows of δ ∈ Γ(72)σ_𝔞 satisfy (c′, d′) ≡ (μ₀, 1) mod 72.
- The Kubota symbol of δσ_𝔞^{−1} ∈ Γ(72) is ((c′ − d′μ₀)/d′)₆ = (c′/d′)₆.
- Since d′ ≡ 1 mod 72 is a local sixth power at 2 and λ (T1), the Hilbert symbols at v ∣ 6 are trivial. Reciprocity then gives (c′/d′)₆ = χ_{c′}(d′).
- So D^𝔞(w;m) is a Kubota series Σ_{c′≡μ₀ (72)} \overline{g(c′, ·m)} N(c′)^{−w}, with fixed phases.
- Write c′ = u s with u a unit and s primary. Then χ_{us} = χ_s and e(hx/(us)) = e(u^{−1}hx/s), so g(us, h) = g(s, u^{−1}h).
- Hence D^𝔞 is a finite combination of the ray-twisted series D_c(w; u′m) of `LoA6.md`, with rows rescaled by fixed units and λ-powers.

**2.2 Gram at 𝔞** [R-]. By 2.1, A^𝔞_m(Y) is a finite combination of the paper's A_{u′m}(Y) with ray-twisted coefficients c. Rescaling by u′ permutes the rows. So

  ‖A^𝔞(Y)‖²_{q_m≍N} ≪ (N/Y)(1 + P^{1/6} + P²/Y) N^ε,   P = Y²/N ≥ 1.

The proof of `probe-gram` uses only |c| ≤ 1 and the fixed-modulus structure of c in the frequency classification.

**2.3 Bootstrap and exponents** [R given the inputs].
- Choose W₁ ≥ 0, not identically zero, and v = 0, so that Ŵ₁(1/6) > 0.
- By (LoA6)_𝔞 with its AFE residue term (`D_derivation.md` §6.2, Addition 2),
  ‖A^{rest,𝔞}(Y)‖² ≪ Y^{−1}(N + N^{1/6} ‖ρ_𝔞‖²) N^ε.
- Hence
  ‖ρ_𝔞‖² ≤ 2 Ŵ₁(1/6)^{−2} Y^{2/3} (‖A^𝔞‖² + ‖A^{rest,𝔞}‖²)
  ≪ [ N^{5/6} + N Y^{−1/3} + Y^{8/3}/N + Y^{−1/3}(N + N^{1/6}‖ρ_𝔞‖²) ] N^ε.
- At Y = N^{3/5} the terms are N^{5/6}, N^{4/5}, N^{3/5}, N^{4/5} and N^{−1/30}‖ρ_𝔞‖². The last is absorbed, since ‖ρ_𝔞‖² is finite for each N.
- So Σ_{q_n≍N} |ρ_𝔞(n)|² ≪ N^{5/6+ε}, i.e. Σ_{q_n≍N} |τ_𝔞(n)|² ≪ N^{1+ε}.
- The constraints used are P ≥ 1 (Y ≥ N^{1/2}) and Y ≤ N^{11/16}.

**2.4 Small N.** For N ≤ Z^{ε}, use the pointwise bound |τ_𝔞(n)| ≪ q_n^{O(1)}, which follows from moderate growth of θ₆ [R].

**Remarks.**
- This confirms `IDEA1_NOTES` §7, and the arithmetic there is correct.
- Two corrections:
  1. The residue term of the AFE must be bootstrapped simultaneously. This is harmless.
  2. "LoA6 at 𝔞" is LoA6 for unit-rescaled rows and ray twists, by 2.1. It is not a new kind of statement.
- No Zagier regularisation is needed.
- The constant terms c_𝔞 (G5) play no role.

---

## 3. G4: non-squarefree dual indices

### 3.1 The input from de Faveri [V]

I read the statement in arXiv:2610.04045v1. Theorem 1.1 bounds Θ_n(A,B), the operator norm with **both** moduli a and arguments b running over n-th-power-free ideals prime to S:

  Θ_n(A,B) ≪_ε (AB)^ε (A + B + A^{1−1/n} B^{2/n} + A^{2/n} B^{1−1/n})   (n ≥ 3).

- For n = 3 this is Heath-Brown's A + B + (AB)^{2/3}, now for **cube-free** moduli and arguments. The cube-free modulus a = a₁a₂² gives χ_a = χ_{a₁} \overline{χ_{a₂}}.
- The family χ_a is the Friedberg–Hoffstein–Lieman extension of (·/a)₃.
- On primary elements it agrees with our cubic symbols up to finitely many ray characters, which are absorbed into the coefficients [R-].

### 3.2 Reduction [R]

After Voronoi (`V_derivation.md` §7) the form is Σ_μ Σ_n x_μ y_n \overline{(n/μ)₃}, with y_n = τ_𝔞(n) q_n^{−it} for q_n ≤ N ≍ L²/Q.

Write n = u n_S n₀ c³, with:
- u a unit;
- n_S the S-part;
- n₀ good and cube-free;
- c good.

Then \overline{(n/μ)₃} = \overline{(u n_S/μ)₃} \overline{(n₀/μ)₃} 1_{(c,μ)=1}. The first factor is a ray character of μ (absorbed into x_μ). The indicator restricts x_μ.

Group by (u, n_S, c) and apply Θ₃(L, N/q_{n_S c³}) to each group:

  |Σ| ≤ ‖x‖ Σ_{u,n_S,c} Θ₃(L, N/q_{n_S}q_c³)^{1/2} ‖y^{(u,n_S,c)}‖.

- **S-parts.** There are O(log² N) values of n_S with norm ≤ N. Using only ‖y^{(n_S)}‖ ≤ ‖y‖ ≪ N^{1/2+ε} (G3) costs O(log² N) in the L-term. In the other terms the sum over n_S converges.
- **Cube parts.** Put K(c) := ‖y^{(c)}‖² / N^{1+ε} (for fixed u, n_S). The three terms of Θ₃ give
  L^{1/2} N^{1/2} Σ_c K(c)^{1/2} + N Σ_c q_c^{−3/2} K(c)^{1/2} + L^{1/3} N^{5/6} Σ_c q_c^{−1} K(c)^{1/2}.

**The required bound.** Let K(c) ≤ q_c^{−κ}. Then:
- **Sufficient everywhere:** κ ≥ 2, which costs log N.
- **Needed by this argument to avoid a power loss, regime L ≥ N²:** κ ≥ 2 (the L-term dominates).
- **Needed by this argument, regime L < N²** (the (LN)^{2/3} term dominates): κ ≥ log L / log N.

At the certified geometry, L = Z^{0.652} and N = Z^{0.478}, so the requirement is κ ≥ 1.36. Hence:

> **(G4*)** Σ_{q_{n₀}≤N/q_c³, n₀ cube-free} |τ_𝔞(u n_S n₀ c³)|² ≪ N^{1+ε} q_c^{−2}, uniformly in good c and in the S-part.

### 3.3 Proof of (G4*) from the local relations at good primes [S] + [R]

**The relations.** In the unitary normalisation (Σ|τ|² ≍ X), Kazhdan–Patterson periodicity and the Hecke relations give, for a good prime p with p ∤ m:
- |τ(mp^{6b+r})| = N(p)^{b/2} |τ(mp^r)|;
- τ(mp⁵) = 0;
- |τ(mp⁴)| = |τ(m)|;
- |τ(mp³)| = |τ(mp)|.

These come from τ(mp^i) = G_{i+1}(m,p) τ(mp^{4−i}) with |G| = 1 (Bröker–Hoffstein §1, CFH (2.5), per `LIT_ideas234.md` [S]).

**Why they hold at every cusp** [R-]. These are local statements about the spherical vector at p. θ_𝔞 = θ∘σ_𝔞 has the same local component at every good p, because σ_𝔞 ∈ SL₂(O_p). So the relations hold at every cusp 𝔞.

**Consistency checks.** Periodicity N(p)^{1/2} matches:
- the classical theta (n = 2): τ(Np²) = p^{1/2}τ(N);
- Patterson's d(nb³) = |b| d(n) for n = 3;
- the triangular residue formula ρ(mp⁶) ≈ ρ(m), which I derived from the twisted multiplicativity of g(p^j s′, ·).

**Reduction map.** Fix c, and at p ∣ c put a = v_p(c) and e = 3a + r₀ with r₀ ∈ {0,1,2}. Write e = 6b + r, and let r′ = r for r ≤ 2, r′ = 1 for r = 3, r′ = 0 for r = 4; for r = 5 the coefficient vanishes. The map n₀ ↦ n♭ (p-exponent r′ at p ∣ c, unchanged elsewhere) is injective for fixed c. Its effect is

  |τ(n)|² = ∏_p N(p)^{b_p} |τ(n♭)|²,   q_{n♭} = q_n / ∏ N(p)^{e_p − r′_p}.

So G3 gives ‖y^{(c)}‖² ≪ (N/D)^{1+ε}, with D = ∏ N(p)^{e_p − r′_p − b_p}. The local exponent e − r′ − b satisfies:

| a | e | b | r → r′ | e − r′ − b | ≥ 2a? |
|---|---|---|---|---|---|
| 1 | 3, 4, 5 | 0 | 3→1, 4→0, 5: τ = 0 | 2, 4 | yes (equality at e = 3) |
| 2 | 6, 7, 8 | 1 | r′ = r | 5 | yes |
| 3 | 9, 10, 11 | 1 | 3→1, 4→0 | 7, 9 | yes |
| a ≥ 5 | — | ⌊e/6⌋ | — | ≥ 2.5a − 7/3 | yes |

(The case a = 4 also holds by direct check.)

Hence D ≥ q_c², and **(G4*) holds with κ = 2 exactly**. Status [R] given G3 and the [S] relations.

**Caveat** [S → must check]. The case e = 3 is tight. If the true relation were |τ(mp³)| = N(p)^δ |τ(mp)|, then κ = 2 − 2δ. That is fine at the certified geometry only for δ ≤ 0.318. Someone should confirm the normalisation of G_{i+1} in KP / Bröker–Hoffstein.

**A route avoiding KP [G].** One could bound the mean square on the row sublattice c³O directly. The paper's Gram bound restricted to m ∈ kO is the Gram for the twisted coefficient c(s)χ̄_s(k), with rows Q/q_k. Its theta part needs a restricted (LoA6). Using de Faveri's Θ₆, that restricted LoA6 has size Q^{1+ε} q_k^{−5/6} for the (s, km) = 1 part. With Y = Q^{3/5} the bootstrap closes and gives Σ_{kO} |ρ|² ≪ (Q/q_k)^{5/6} for q_k ≤ Q^{0.12}. But the terms with (s,k) > 1, which carry the periodicity, and [G-b] of `LoA6.md` reappear. Not done.

---

## 4. Non-squarefree rows in the theta pairing

### 4.1 Structure [R, from the paper]

Lemma `unmarked-completed-moment` organises the rows as follows:
- m = u · m_S-part · m_pow · m_sup · R, where m_pow (the maximal powerful factor) and m_sup (the valuation-one S-part) are frozen;
- R is the free squarefree good part;
- the frozen good primes of m_pow carry j_p = v_p(m) mod 6, and are active or inactive (j_p = 0 may be inactive);
- q_{F_act}² ≤ q_{m_pow} (eq. `unmarked-powerful-radical`).

The column factor is χ_R(λ⁴μ)³ · ∏_{p∈F_act} B_p(λ⁴μ), i.e. \overline{χ_{m rad(m)²}(μ)} away from j_p ∈ {0, 4}. The cusp denominator is c = c_F F_act R.

### 4.2 Coherence: no column gives an effective j = 1 twist in the free variable [R]

As a function of the summed variable R, the twist is:
- **Column side:** χ_R(λ⁴μ)³ = (fixed) · (μ/R)₂. This is sextic exponent 3 at the column primes, for every column and every frozen choice. The B_p factors at frozen primes depend on μ only.
- **Row scalar (§1.3):** ∏_{p∈F_act} χ_p(R)^{2j_p}, i.e. exponent 2j_p ∈ {0, 2, 4} mod 6. It is never ≡ 1. It is ≡ 0 (just coprimality) exactly when j_p ∈ {0, 3}.

So in the ζ-free model, Σ_R ρ(kR)(·) has a main term only from columns whose squarefree kernel μ′ is trivial (μ ∈ uλ^k □, as in `V_derivation.md` §6). The row scalar's cross factor has no such effect.

**A caution about phrasing.** Viewed as a twist of the *full* index m, a frozen prime with v_p(m) ≡ 3 mod 6 contributes χ_μ(p)^{−(j_p+2)} = χ_μ(p)¹ (a "j = 1" exponent). But that prime is not summed, so this is a constant, not a coherent twist.

### 4.3 The Ramanujan case j_p = 4 [R]

B_p(x) = q_p^{−1/2}(−1 + q_p 1_{p∣x}) acts on columns only. The paper splits it into two summands:
- an unrestricted column with weight q_p^{−1/2};
- columns in pO with weight q_p^{1/2}.

The column ℓ²-mass is Σ_μ |β_μ|² |B_p|² ≈ q_p^{−1} + O(1) = O(1), because Σ_{p∣μ} |β_μ|² ≍ q_p^{−1}. With ω(F_act) = O(log Z / log log Z), all such factors together cost Z^ε.

The case j_p = 0 active (B_p = q_p^{−1/2} χ_p^{−2}) is likewise column-only.

### 4.4 Block counting

Let s = s(M) be the target saving. At the certified point s = −½·max(M−1, 7M/3−2, −M/6) = 0.0363.

**(a) Large powerful part, O ≥ 4s: Cauchy–Schwarz suffices** [R given G3 and the paper's eq. `unmarked-completed-moment` at N = 1, M < 1].
- The paper's moment bound gives ‖B‖²_{O-dyad} ≪ Z^{O/2 + max(M−O, 2M−O−1) + ε} = Z^{M − O/2 + ε}.
- G3 gives the trivial bound ‖ρ‖²_{O-dyad} ≤ ‖ρ‖² ≪ Q^{5/6+ε}.
- So the O-dyad contributes ≤ Q^{5/12} Z^{M/2} Z^{−O/4}, i.e. a saving Z^{−O/4} over the global Cauchy–Schwarz scale. That is ≥ Z^{−s} once O ≥ 4s ≈ 0.145.

**(b) Small powerful part, O < 4s, in the ζ-free model** [R-].
- Expand 1_{k∣m} additively and absorb the frozen cubic character ∏(R/p)₃^{j_p}. The cusp denominators then divide k F′ μ′, with F′ ∣ rad(k).
- The frozen sextic Gauss / Kloosterman sums mod k depend on the dual index n only, so they are absorbed into y_n.
- Columns satisfy L_k ≤ L/q_k, since q_{F_act}² ≤ q_k.
- The dual length satisfies N_k ≤ q_k³ L_k² / Q ≤ q_k N.
- Per frozen k: |T_k| ≲ Q^{5/12} q_k^{−1/2} Δ(L_k, N_k)^{1/2}, where Δ(L,N) = L + N + (LN)^{2/3}.
- Summing over the ≍ Z^{O/2} powerful k in the dyad:
  Σ_k |T_k| ≲ Q^{5/12} (L Z^{−O} + Z^O N + (LN)^{2/3})^{1/2}.
- Here Z^O N ≤ L because Z^O ≤ Z^{1−M}, which holds for O < 4s ≈ 0.145 < 1 − M = 0.174.
- So every small-O dyad is bounded by the squarefree-block bound Q^{5/12} Δ(L,N)^{1/2}.
- There are O(log Z) dyads.

**Conclusion.** In the ζ-free model, powerful rows do not destroy s(M). With F1, the same block structure holds, but each block's free sequence carries \overline{γ₂(R)}⁴ and (V) does not apply to any block.

---

## 5. Non-squarefree columns [R-]

The true columns are μ = u λ^k n b³ with n squarefree and b arbitrary, not the squarefree μ of `V_derivation.md` §7.
- Write nb = μ′t² with μ′ squarefree. The row twist is then (μ′/R)₂ 1_{(R,t)=1}.
- Expanding 1_{δ∣R} (δ ∣ t) additively adds cusps x/δ. The multiplier there is χ_δ(x)^{+1}, and the dual picks up the sextic character χ_δ(n) (twist exponent 0 at δ), times τ(χ̄_δ).
- For fixed δ, χ_δ(n) is a function of n alone and is absorbed into y_n. Heath-Brown or de Faveri is then applied over μ′.
- The column mass with square part t is ≪ q_t^{−2} log L, from |β_μ|² ≤ 729/(q_n q_b²) with b ⊇ t. So Σ_{δ∣t} costs Σ_t τ(t) q_t^{−1} ≪ log² L.
- The dual length for the δ-cusps is q_{μ′}² q_δ² / Q ≤ L²/Q.

This is only bookkeeping. It is moot until F1 is repaired.

---

## 6. The v-integral [R]

In eq. `low-separated`, insert A_{m,σ,v} = Y^{−1/3} Ŵ₁(1/6+iv) ρ_σ(−m) + A^rest.
- ρ_σ does not depend on v: the residue at w = 2/3 of Y^{−1}(2πi)^{−1} ∫ Ŵ₁(w − 1/2 + iv) Y^w D_σ(w; −m) dw is Y^{−1/3} Ŵ₁(1/6+iv) ρ_σ(−m).
- The m-sum is finite and both Mellin transforms decay rapidly, so the v-integral can be done first:

  (1/2π) ∫ Ŵ₀(iv) Ŵ₁(1/6+iv) x^{−iv} dv = ∫₀^∞ W₁(t) t^{1/6} W₀(x/t) dt/t =: K(x).

- K is smooth and compactly supported in (0, ∞). Since Ω = 1 on the support,

  I^sig = Q^{−1/2} Y^{−1/3} Σ_σ Σ_m K(q_m/Q) ξ(m) ρ_σ(−m) B_{m,σ}(Z).

So there is no v-integral at all in the theta part. This K is the K of `IDEA1_NOTES` §2. The rest part keeps its v-integral, with the polynomial (1+|v|)^{R′} of LoA6 against the rapidly decaying Ŵ₀Ŵ₁.

---

## 7. The ray-class restriction s ∈ σ [R-]

c(s) = 1_{s∈σ} χ_s(b_*) / (τ ξ(s)) is a finite combination of finite-ray characters ψ of s: the indicator by orthogonality on T, χ_s(b_*) by Lemma `fixed-numerator-ray`, and 1/ξ(s).
- **Automorphy.** For each ψ, D_ψ is the Fourier coefficient of E_ψ = Σ_{Γ_∞\Γ(72)} ψ(class of the prime-to-6 part of c mod f) κ̄(γ) y(γv)^s. This weight is invariant under right multiplication by Γ(72f). So E_ψ, and its residue θ_ψ, are automorphic on Γ(72f) with multiplier κ [R].
- **Lemmas.** Lemma M holds verbatim at level 72f, with μ₀ classes taken mod 72f. G3 and G4 hold at each of the finitely many cusps of Γ(72f).
- **Remaining bookkeeping.** The varying S-part of c is the bookkeeping of G1.

---

## 8. Updated gap list for the theta pairing

| gap | status after this note |
|---|---|
| **F1: row scalar ζ_m ∝ \overline{γ₂(m)}⁴** | **[G], fatal for (ii) as derived.** Needs (V♮) of §1.7(c), or a different probe. |
| G1: identification of ρ with θ-coefficients (all cusps, S-parts) | [R-], unchanged |
| G2: Hankel kernel | [R-], unchanged; only decay is used |
| G3: cusp mean squares | **closed** relative to (LoA6 at all cusps) + G1 + the paper's Gram for ray twists (§2) |
| G4: non-squarefree dual indices | **closed** relative to G3 + de Faveri Thm 1.1 [V] + the KP local relations [S, normalisation to confirm] (§3) |
| G5: constant terms | irrelevant for upper bounds |
| (P) rows | coherence [R]; Ramanujan case [R]; blocks [R] for O ≥ 4s, [R-] otherwise (§4) |
| (P) columns | [R-] (§5) |
| (P) v-integral | [R], eliminated (§6) |
| (P) ray class | [R-] (§7) |
| (P) high side at ℓ = 0 | not addressed here (see `REPORT.md` §3.3) |

**Bottom line.** All the bookkeeping around (V) can be closed, relative to LoA6 and standard KP local theory. But (V) bounds ⟨ρ, ζ̄B⟩, not ⟨ρ, B⟩. The real pairing involves ρ(m)·\overline{γ₂(m)}⁴, and nothing currently gives it a power saving. Without such a saving the method returns Liu's 0.874957.

---

## Referee report

*Independent skeptical review, 9 Oct 2026. Sources checked: `paper78.tex` (Prop. `completed-reflection` statement and proof, eqs. `reflection-local-fourier`, `reflection-row-phase`, `reflection-cross-prime-phase`, `signal`, `fixed-gauss-phases`, `low-separated`, Lemma `unmarked-completed-moment`, Prop. `probe-gram`); de Faveri arXiv:2610.04045 (Thm 1.1, read from the arXiv HTML); Bröker–Hoffstein arXiv:1312.0568 §1 (read from the PDF). I re-ran `check_rowscalar.py`, `lp_after_F1.py` and `lp_const_saving.py`.*

**Overall.** I could not refute any main conclusion. F1 is correct, and its consequences are wider than the verdict says (item 2). The LP numbers reproduce exactly. The G4 normalisation caveat can now be removed in the authors' favour (item 5). The heuristic in §1.5 and the description in §1.7(a) overstate how hopeless a repair is (item 4), but that does not change the verdict.

### 1. F1 is correct [R, checked against the paper's proof]

**Where the Gauss sum comes from.** In the proof of Prop. `completed-reflection`, eq. `reflection-local-fourier` gives the finite Fourier coefficient of the row character χ_p^j:

  C_{p,j}(h) = q_p^{−1/2} τ⁻_{p,j} χ_p(h)^{−j}.

- For an active row prime with j_p = 1, this carries the normalised sextic Gauss sum τ⁻_{p,1} = χ_p(−1)γ₁(p) of the row prime itself.
- Summing h against the cusp multiplier then produces τ⁺_{p,3}.
- This is exactly the factor ω_{p,1} in eq. `reflection-row-phase`.

**Nothing else in ζ can cancel it.**
- ᾱ(c)² is the Hecke angle c̄²/|c|², coming from c^{−2} = ᾱ(c)²/q_c in eq. `reflection-reflected-mellin`.
- Every other factor is a residue symbol of fixed elements or of row elements.
- So ζ_R ∝ (ray)·ᾱ(R)²·γ₁(R)γ₃(R) is unavoidable. Since γ₃ = χ(4)G is a ray function and conj(γ₂) = μᾱḠγ₁ (eq. `signal`), this is the same as the authors' conj(γ₂(R))⁴ = ᾱ(R)²Ḡ(R)γ₁(R).

**Lemma Z, re-derived by hand.**
- Steps 1–5 are correct.
- The general-row exponent in §1.3 agrees with the paper's own eq. `reflection-cross-prime-phase`, χ_p(σ_p)^{−2}ω_{p,j} = ξ_{p,j}χ_p(c/p)^{2j+2}. Adding −2 from R's primes, via cubic reciprocity, gives 2j_p.
- This also holds for j_p = 4 (−2−2 ≡ 2 ≡ 8) and for an active j_p = 0 (4−2−2 = 0).

**Numerics.** I re-ran `check_rowscalar.py`:
- the identity holds to 2.6·10⁻¹⁵ on 2248 split primes;
- for every k = 0,…,5, the within-class phase spread of γ₁γ₃α^k is 3.9–4.3 rad;
- the ray control has spread 0.

`check_row_phase.py` (primes up to 6·10⁴) did not finish in 10 minutes. It is not needed.

**An independent structural confirmation [H, but strong].**
- In the ζ-free model, the s = 5/6 main terms of the reflected pairing come from columns μ ∈ uλ^k·□ (`V_derivation.md` §0(c)). Their total is ≍ Q^{5/6} and is dominated by μ = O(1). It carries no arithmetic weight over A ≍ Z.
- The exact unreflected computation (§1.6, c_A/τ(χ̄_A) = q_A^{−1/2}μ(A)χ_A(−1)ν(A), which I re-checked) shows that the s = 5/6 residue of Σ_m ρ(m)B_m q_m^{−s} is the Möbius sum Σ_{A≍Z} μ(A)ν′(A)(φ(A)/q_A)q_A^{−5/6}V_G(q_A/Z).
- These two cannot agree identically in Z. So the ζ-free reflected form never equalled ⟨ρ,B⟩: the Möbius information must be carried by the row sequence, which is what F1 says.
- It also follows that the statement in `V_derivation.md` §0(c), "this total is exactly the s = 5/6 pole, i.e. the signal", is wrong as stated.

### 2. F1 reaches further than the verdict says

All three conditional values in `REPORT.md` §0 rest on the (V) + Heath-Brown saving (`LS`/`LSmarked` in `theta_split_model.py`):
- 0.8665, slot-free;
- 0.871, slotted, with marks;
- 0.864, with the detector floor.

All three are therefore unsupported. The slotted version is no better off: its marked slot primes are active in the reflection, so they also bring their own τ⁻ Gauss sums into ζ.

The "conditional theorem" in `REPORT.md` §1, and `IDEA1_NOTES.md` §3–5, should be restated with (V♮) of §1.7(c) in place of (V). What survives:
- (LoA6) ⇒ (D), which alone gives no gain;
- Lemma D1;
- the negative results.

### 3. LP consequences reproduced exactly

**`lp_after_F1.py`.**
- 0.874957 for the paper's model.
- 0.874957 with gram_q = 0 (P²/Y dropped). This confirms that the P²/Y term is inactive at Liu's optimum, (lx, ly, ℓ) = (0.3549, 0.4783, 0.1668).
- 0.874957 with caseA refined.
- Slot-free (P²/Y dropped, theta part by CS): 0.894728 at S = 10⁻⁴, 0.893892 at S = 0.01, 0.890364 at S = 0.05.

**`lp_const_saving.py`.** The §0 table reproduces exactly:

| s₀ | 0 | 0.002 | 0.005 | 0.01 | 0.02 | 0.04 |
|---|---|---|---|---|---|---|
| θ | 0.874957 | 0.873332 | 0.870895 | 0.866832 | 0.861510 | 0.853229 |

The LS-shaped saving with P²/Y dropped gives 0.86624 at (0.1978, 0.6282, ℓ ≈ 0).

### 4. The heuristics in §1.5 and §1.7(a) are overstated [H]

These points do not affect the verdict, but the "fatal" language should be read as "an open gap with no known tool", not as a proven obstruction.

1. **The object is purely sextic.** Exactly, a♮(R) = ρ(R)·γ₁(R)·(Hecke ᾱ²·ray), where γ₁(R) = g̃₆(R,1). For coprime squarefree s and R,
   g̃(s,−R)g̃(R) = χ̄_s(−R)g̃(s)g̃(R),
   which up to a ray/sign factor is the coefficient H(s,R) of the sextic A₂ Weyl-group multiple Dirichlet series (Brubaker–Bump–Friedberg–Hoffstein, Annals 2007, for all n). So Σ_R ρ(R)γ₁(R)N(R)^{−z} is, on its squarefree coprime part, the w-residue of a pure 6-fold A₂ object. It is not a series "mixing 6-fold and 3-fold Gauss sums".
2. **Degree counting is unreliable for theta coefficients.**
   - In the cubic analogue, Patterson/Kazhdan–Patterson give τ₃(d) = G₁(1,d) at squarefree d (Bröker–Hoffstein §1). So θ₃(d)·conj(G₁(1,d)) ≡ 1, a degree-1 object, while naive "θ₃ ⊗ θ̄₃ Rankin–Selberg" counting would predict degree 4.
   - For n = 6 the plausible heuristics disagree:
     - Rankin–Selberg gives dual length L⁴/Q = Z^{1.78};
     - a GL(3)^{(6)} Eisenstein-series-induced-from-θ₆ picture gives L³/Q = Z^{5M−3} = Z^{1.13} at M = 0.826. This is still above Q, but it would fall below Q for M < 3/4;
     - the one-variable functional equation of the residue uses two simple reflections, which might suggest conductor N(μ)².
3. **The phase of a♮ at primes is not known.** Bröker–Hoffstein's data, Conj. 5.7, give τ₆(π)² as a cubic Gauss sum times a "patternless" algebraic factor (powers of 3 and of elements of norm 7). Depending on orientation, a♮(π)² is then either conj(γ₂)·Hecke·F or Hecke·F. So even whether a♮ oscillates like a cubic Gauss sum is open.

The sentence "even the most optimistic Voronoi picture gives no saving" should therefore read: "no tool is known, and the heuristics are inconclusive."

### 5. The G4 normalisation caveat is resolved in the authors' favour [V, referee-checked]

**Bröker–Hoffstein's normalisation matches.** They normalise
  θ^{(n)}(z) = τ_n(0)y^{1−1/n} + y Σ τ_n(m)K_{1/n}(4π|m|y)e(mx),
with
  τ_n(m) = N(m)^{1/(2n)} Res_{2s=1+1/n} A_m(s),   A_m(s) = Σ_d g_n(m,d)N(d)^{−2s}.
- With w = 2s − 1/2, A_m(s) is exactly D(w;m), with the pole at w = 2/3 for n = 6.
- So their τ₆ is exactly the authors' τ = q^{1/12}ρ.

**Their relations (from Kazhdan–Patterson).** For (m,p) = 1 and 0 ≤ j ≤ n−2:
- τ_n(mp^j) = G_{j+1}(m,p)τ_n(mp^{n−2−j}), with G_j normalised to absolute value 1;
- τ_n(mp^{n−1}) = 0;
- τ_n(mp^n) = N(p)^{1/2}τ_n(m).

**Consequence.**
- |τ(mp³)| = |τ(mp)| exactly, so δ = 0.
- |τ(mp⁴)| = |τ(m)|.
- (G4*) holds at ∞ with κ = 2, as the authors derived.
- Transfer to the other cusps is still [R-] (local sphericity at good p), as stated.

**de Faveri.** Thm 1.1 is as quoted: Θ_n is the supremum over ‖λ‖₂ = 1 of Σ^♭_{N𝔞≤A} |Σ^♭_{N𝔟≤B} λ_𝔟χ_𝔞(𝔟)|², where ♭ means n-th-power-free ideals in I(S), and the bound is (AB)^ε(A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}).

**Arithmetic.** I checked:
- the local-exponent table;
- the injectivity of n₀ ↦ n♭ for fixed c;
- the threshold κ ≥ log L/log N = 0.652/0.478 = 1.364.

### 6. Minor points (no effect on the conclusions)

- **§2.3.** "Y ≤ N^{11/16}" is not a hypothesis of Prop. `probe-gram`. It is only the condition Y^{8/3}/N ≤ N^{5/6}, which is irrelevant at Y = N^{3/5}. The bootstrap arithmetic (N^{5/6}, N^{4/5}, N^{3/5}, N^{4/5}, N^{−1/30}‖ρ‖²) is correct. The simultaneous all-cusp bootstrap is needed because the AFE's dual pieces carry residues at other cusps; the authors say this.
- **§2.1.** Cusps 1/μ₀ with λ | μ₀ or 2 | μ₀ give moduli c′ with an S-part, so "c′ = us, s primary" needs G1-type S-bookkeeping. The [R-] label is appropriate.
- **§1.2, step 4.** α(a) = a/|a| is completely multiplicative by definition (`paper78.tex` l. 575). The argument via γ₂³ is unnecessary but harmless.
- **§4.4(a), checked.**
  - The paper's moment bound at N = 1 with M < 1 gives Z^{M−O/2}.
  - Cauchy–Schwarz with ‖ρ‖²_O ≤ Q^{5/6} saves Z^{−O/4}.
  - 4s(M) = 0.145 at M = 0.826.
  - Z^O N ≤ L if and only if O ≤ 1 − M = 0.174.
- **§6, checked.** With the paper's convention Ŵ₀(iv) = ∫W₀(r)r^{iv}dr/r, one has K̂(s) = Ŵ₁(s+1/6)Ŵ₀(s), and the v-integral identity holds.
- **Labels.** "Lemma Z [R]" and "§4.2 coherence [R]" are justified relative to the paper's displayed formulas. The G4 label can be upgraded from [S] to [V] for the normalisation at ∞ (item 5).

### 7. Recommendation

Accept the verdict "flawed". F1 is a genuine, previously missed gap. It is not a proof that the theta-split idea cannot work.

Required edits elsewhere:
- mark all three certificates in `REPORT.md` §0 as resting on the unproved (V♮);
- correct `V_derivation.md` §0(c), which identifies the square-column main terms with the signal.

Most promising direction: formulate (V♮) as a statement about quadratic twists, in the row variable, of the residue of the sextic A₂ Weyl-group multiple Dirichlet series.
