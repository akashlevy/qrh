# Lemma V: Voronoi summation for quadratic twists of the sextic theta coefficients

Working notes, 9 Oct 2026. Companion scripts in this directory:
`cubic_theta.py`, `exp_cubic_voronoi.py`, `exp_cubic_scaling.py`, `exp_sextic_jtable.py`. The older scripts are `exp_twist.py` and `exp_voronoi.py`.

Status labels used below:
- **[R]** rigorous, given the cited theorems.
- **[R-]** routine but not written out in full here.
- **[H]** heuristic.
- **[G]** genuine gap.
- **[C]** needs the theta function at a cusp other than ∞, or bounds for its coefficients.

---

## 0. Conclusions

**(a) Dual length.** The dual length is C·N(μ)²/Q. C is a fixed constant, depending only on the level 72, the cusp width and the mod-4 class splitting. It is irrelevant for exponents but numerically large (the cubic numerics indicate C ≳ 10⁴). The dual weight decays like (1+x)^{-A} in x = C N(n)Q/N(μ)², and is O(x^δ) as x → 0. [R-; standard; **not** numerically verifiable at the accessible Q, see §8.2]

**(b) The dual character is cubic, never principal.** Write χ_μ(x) = (x/μ)₆ and (x/μ)₃ = χ_μ(x)². The quadratic twist (m/μ)₂ = χ_μ(m)³ produces the dual character

  χ_μ(n)^{-2} = \overline{(n/μ)_3},

with the unimodular factor

  ε(μ) = γ₂(μ)/γ₃(μ) = μ(μ) α(μ) \overline{χ_μ(4)} \overline{γ₁(μ)}.

This is up to fixed ray-class phases, which depend only on μ mod 72 and on the class of m mod 4 used to convert (μ/m)₂ into (m/μ)₂. The "Ramanujan / principal" alternative is impossible for the quadratic twist.
- The multiplier of θ₆ at the cusp a/μ is χ_μ(a)^{±1}, a primitive sextic character (§4).
- Hence a twist χ_μ^j combines to χ_μ^{±1−j}. This is principal only for j ≡ ±1, never for j = 3.
- The exponent is pinned down three independent ways (§3), for ρ defined by D(w;m) = Σ_s g̃(s,−m)N(s)^{-w}. The multiplier is χ_μ(a)^{+1}, and the coherent twist is j = 1.
- In general, the twist χ_μ^j (j ≢ 0, 1) has dual character χ_μ(−n)^{1−j} and unimodular factor τ(χ_μ^{j−1})/τ(χ_μ^{−j}). [R, given Kubota and the identification of ρ in §3.]

**(c) Main terms.** For squarefree primary μ ≠ 1 the d-sum at n = 0 vanishes identically. It is a complete sum of a primitive cubic character. So there is **no main term at all**. [R]
- In general, main terms occur exactly for columns at which the quadratic twist is principal: μ ∈ u λ^k · (squares). The notes say "cubes", which is not right.
- In the support of the reflected columns, μ = u λ^k n b³, this means nb = □, that is μ = u λ^k n⁴ c⁶.
- For such μ, F(μ) is a product of local factors times the untwisted main term, which is ≍ Q^{5/6}.
- The total is O(Q^{5/6+ε}), dominated by small μ. This total is exactly the s = 5/6 pole, i.e. the *signal* of §2 of IDEA1_NOTES. It is not an error term.

**(d) Bound for T.** Heath-Brown's cubic large sieve gives

  |T| ≪ Q^{5/12+ε} (L + L²/Q + L²/Q^{2/3})^{1/2} (μ squarefree, μ ≍ L, |β_μ| ≍ N(μ)^{-1/2}).

Compared with the Cauchy–Schwarz / Goldmakher–Louvel bound Q^{11/12+ε}, the ratio is exactly the claimed saving

  [L/Q + L²/Q² + L²/Q^{5/3}]^{1/2}.

Moreover, T = O(Q^{-A}) when L ≤ (Q/C)^{1/2}Q^{−ε}, because the dual sum is empty. [R modulo gaps G3, G4]

---

## 1. Set-up and conventions (paper78 §"Arithmetic and analytic preliminaries")

**Basic objects.**
- K = Q(ω), O = Z[ω], λ = 1+2ω.
- e(z) = exp(4πi Im z/√3). Primary means ≡ 1 mod 3.
- χ_c(x) = (x/c)₆, with the *denominator* c and zero on nonunits. The cubic symbol is (x/c)₃ = χ_c(x)², and the quadratic symbol is (x/c)₂ = χ_c(x)³.
- γ_j(c) = q_c^{-1/2} Σ_{v mod c} χ_c(v)^j e(v/c), and τ(χ_c^j) := √q_c · γ_j(c).
- Sextic reciprocity: χ_b(a) = R(a,b)χ_a(b).

**Gauss sums and the coefficients ρ.**
- g(s,h) = Σ_{d mod s} χ_s(d) e(hd/s) and g̃ = g/√N(s).
- D(w;m) = Σ_{s primary, (s,6)=1} g̃(s,−m) N(s)^{-w}, and ρ(m) := Res_{w=2/3} D(w;m). The numerics restrict D to squarefree s; see G1.
- The "unitary" coefficients are τ(m) := N(m)^{1/12} ρ(m).
- Numerically, Σ_{N m≍Q}|ρ|² ≍ Q^{5/6}, i.e. Σ_{N m≍Q}|τ|² ≍ Q. This is also what Rankin–Selberg predicts (§6).

**The twisted sum.**

  F(μ) = Σ_m ρ(m) (μ/m)₂ W(N(m)/Q),

with μ primary squarefree, (μ,6) = 1 and N(μ) = L < Q.

**Reduction to ψ_μ(m) = (m/μ)₂.** By quadratic reciprocity in K, (μ/m)₂ = R₂(μ,m)·(m/μ)₂. Here m is written as a unit times a primary element, and R₂ depends only on the classes of μ and m modulo a fixed modulus N₁ | 12^∞.

Split m into these classes. Each restriction m ≡ η (mod N₁) is a finite combination of additive twists e(hm/N₁), so θ is replaced by finitely many translates θ(z + h/N₁). Every argument below applies verbatim to each translate. The cusp denominators become c_F·μ with c_F | N₁ fixed, and the multiplier computation of §4 is unchanged (see the remark there).

From now on, write ψ_μ(m) = (m/μ)₂ = χ_μ(m)³. It is **primitive** mod μ for squarefree μ, since each local component is the nontrivial quadratic character.

---

## 2. Known theorems used

**T1. Kubota (J. reine angew. Math. 222 (1966); Kyoto lectures 1969).** Assume K ⊃ μ_n and n² | N₀. Then κ(γ) = (c/d)_n for c ≠ 0, and κ(γ) = 1 for c = 0, is a character of Γ(N₀) ⊂ SL₂(O_K).
- Here (c/d)_n is the n-th power residue symbol of the lower-left entry c modulo the lower-right entry d. This is "Kubota's symbol (c/d)".
- We use n = 6 and N₀ = 72 = 8λ⁴.
- At N₀ = 72 one also has 1+72O ⊂ (K_v^×)⁶ for v | 6:
  - at 2, 1+2O₂ consists of cubes, and 1+8O₂ consists of squares;
  - at λ, 1+λO_λ consists of squares, and 1+λ⁴O_λ consists of cubes (Hensel applied to (1+x)³ = 1+3x+…).
- Hence for γ ∈ Γ(72) the Hilbert symbols at 2 and λ are trivial. The product formula then gives

  (c/d)₆ = (d/c′)₆ := ∏_{p|c, p∤6} (d/p)₆^{v_p(c)} = χ_{c′}(d),

  where c′ is the prime-to-6 part of c. So on Γ(72) the two forms "(c/d)₆" and "(d/c)₆" coincide.

**T2. Kubota (1969) and Kazhdan–Patterson ("Metaplectic forms", Publ. IHES 59, 1984), for GL₂ and n = 6.**
- The Eisenstein series E(v,s) = Σ_{Γ_∞\Γ(N₀)} κ̄(γ) y(γv)^s on H³ continues meromorphically.
- It has a simple pole at s = 1 + 1/6. The residue θ₆ is a nonzero automorphic form with multiplier κ. It has moderate growth and lies in L²(Γ(N₀)\H³).
- At every cusp 𝔞, θ₆ has the expansion

  θ₆(σ_𝔞 v) = c_𝔞 y^{5/6} + Σ_n τ_𝔞(n) y K_{1/6}(4π|n|y) e(nz).

- The m-th coefficient of E is ∝ |m|^{s−1} Σ_c (Kubota Gauss sum) N(c)^{-s}. This is why the pole of D(w;m) sits at w = s − 1/2 = 2/3.

**T3. Kazhdan–Patterson (structure).**
- For the n-fold cover of GL₂, the theta representation has a *unique* Whittaker model only for n = 2, 3.
- For n = 6 the coefficients ρ(m) at squarefree m are not determined by Hecke theory, so there is no closed formula. This is why everything below uses only automorphy.
- The KP periodicity theorem and the local exceptional Whittaker functions control τ at prime powers. This is needed for G4 and is not carried out here.

**T4. Patterson ("A cubic analogue of the theta series", Crelle 296, 1977), in the form of Dunn–Radziwiłł (Annals 200, 2024, §5 and App. A) as recorded in paper78 eq. (reflection-infinity-coefficient).**
- The coefficient of θ̄₃ at λ^{-3}nb³ is 3^{5/2}|b| g̃₃(n), for n squarefree and b primary.
- This is used only for the cubic cross-check (§8).
- The cubic analogue of the present lemma is paper78's Prop. "completed reflection": twist exponent j_p ↦ local dual factor χ_p(x)^{−j_p−2}.

**T5. Heath-Brown (Israel J. Math. 120, 2000, Thm 2), cubic large sieve over Q(ω).** For squarefree primary a and b,

  Σ*_{q_a≤U} |Σ*_{q_b≤V} c_b (b/a)₃|² ≪ (UV)^ε (U + V + (UV)^{2/3}) Σ|c_b|².

This is paper78 eq. (cubic-large-sieve).

**T6. Goldmakher–Louvel (Math. Proc. Camb. Phil. Soc. 154, 2013), quadratic large sieve over number fields.** The bound is (M+N)(MN)^ε. It is used only for the comparison bound in (d).

**T7. Dunn–Radziwiłł (2024).**
- The cubic theta function and its coefficients at the three cusps. These are the source of the column weights β_μ through paper78's reflection.
- The "theta-Voronoi + cubic large sieve" template that is copied here.

**T8 (not available off the shelf, see G3).** Zagier's renormalised Rankin–Selberg method (1981). It would give Σ_{N(n)≤X}|τ_𝔞(n)|² ≍ X for θ₆ at each cusp 𝔞; θ₆ is L² but not cuspidal, since |θ₆|² ~ y^{5/3}.

---

## 3. Which theta function has coefficients ρ, and its multiplier

### 3.1 Eisenstein side [R-, bookkeeping G1]

Write γ = (a b; c d). Then γ(z) = a/c − 1/(c²(z + d/c)), so the m-th Fourier coefficient of E_κ := Σ κ̄(γ) y(γv)^s is

  ∝ Σ_c N(c)^{-s} Σ_{d mod c} κ̄(c,d) e(md/c).

On Γ(72), κ̄(c,d) = \overline{χ_{c′}(d)} (T1), and

  Σ_d \overline{χ_{c′}(d)} e(md/c′) = \overline{g(c′,−m)}.

The 2- and λ-parts of c give finitely many fixed phases (depending on m mod 72 and c′ mod 72²) and ray-class twists.

So E_κ carries \overline{g(c,−m)}. Conjugating gives Ē_κ, which has multiplier κ̄, and whose m-th coefficient carries g(c,−m) up to m ↦ −m. The map z ↦ −z, conjugation by diag(i,−i), changes the multiplier by (−1/d)₆ = 1 on Γ(72).

Hence the form θ_ρ whose coefficients are ρ (in the unitary normalisation τ) has

> **multiplier κ̄(γ) = (c/d)₆^{-1} = \overline{χ_{c′}(d)} on Γ(72).**

### 3.2 Independent check from the s-side (Poisson in m) [R as algebra; the polar-term reading is H]

Fix s primary with (s,μ) = 1, and the twist ψ = χ_μ^j. Poisson summation in m modulo sμ gives

  Σ_m g(s,−m) ψ(m) W(N m/Q) = (Q/N(sμ)) Σ_{k∈O} S(s,μ;k) Ŵ(√Q k/(sμ)).

Using the CRT splitting 1/(sμ) ≡ μ̄_s/s + s̄_μ/μ, one finds exactly

  S = N(s) τ(ψ) · χ_s(k) \overline{χ_s(μ)} · ψ(s) \overline{ψ(k)} = N(s) τ(ψ) R(μ,s) · χ_s(k) χ_μ(s)^{j−1} · χ_μ(k)^{−j}.

The polar (θ-) part comes from those k for which s ↦ χ_s(k)χ_μ(s)^{j−1} is principal, i.e. kμ^{j−1} ∈ K^{×6}·(fixed). For μ = 1 this reproduces exactly the pole at w = 2/3 with residue ≍ Q^{5/6}.
- For j = 1 the polar terms are k = c⁶. These survive the factor χ_μ(k)^{−1}, so there is a main term.
- For j ∈ {2,3,4,5} the polar terms are k ∈ μ^{1−j}K^{×6} ⊂ μO. These are killed by χ_μ(k)^{−j} = 0, so there is no main term.

So **the twist χ_μ^j of ρ is coherent iff j ∈ {0,1}.** This orientation is intrinsic to the definition of D and independent of Kubota-symbol conventions.

### 3.3 Numerical check [numerics, see §8]

`exp_twist.py` and `exp_sextic_jtable.py` (Q = 3000) report the cross-scale estimate of ⟨|F_j|²⟩ / Σ|What ρ|². For N(μ)²/Q ∈ [4,64) (56 primes):

| j | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| ratio (s.e.) | 6.0 (0.75) | 0.65 (0.31) | 0.65 (0.27) | 1.19 (0.39) | 1.71 (0.44) |

At N(μ)²/Q ≥ 64 (468 primes) the j = 1 ratio is still 2.15 (0.17); the others are 0.6–1.6.
- Only j = 1 is coherent, consistent with 3.1 and 3.2.
- j = 0 is unusable: the estimator is biased by a deterministic, scale-dependent component (ratio −44, see §8.3).

---

## 4. Lemma M: the multiplier at the cusp a/μ [R, given T1]

**Lemma M.** Let μ be primary with (μ,6) = 1, and let μ₀ ≡ μ (mod 72) be a fixed representative of its class. Put σ = L(μ₀) := (1 0; μ₀ 1). For a coprime to μ, choose the representative a ≡ 1 (mod 72), d ≡ a^{-1} (mod μ), d ≡ 1 (mod 72), and b = (ad−1)/μ. Then g_a = (a b; μ d) ∈ SL₂(O), g_a ≡ σ (mod 72), and

  κ(g_a σ^{-1}) = (d/μ)₆ = χ_μ(d) = \overline{χ_μ(a)}.

*Proof.*
1. Since ad ≡ 1 mod μ and mod 72, we get b ∈ 72O, so g_a ≡ σ mod 72.
2. For a = 1, g₁σ^{-1} = (1 0; μ−μ₀ 1) has κ = 1.
3. Put G = g_a g₁^{-1} = (a−bμ, b; μ(1−d), d) ∈ Γ(72). By T1,

  κ(G) = (μ(1−d)/d)₆ = (d/μ)₆ · (d/(1−d)′)₆ = χ_μ(d),

  because d ≡ 1 mod (1−d).
4. κ is a homomorphism on Γ(72), so κ(g_aσ^{-1}) = κ(G)κ(g₁σ^{-1}). ∎

*General version.* For any a, a′ with g_a ≡ g_{a′} mod 72, put G = g_a g_{a′}^{-1}. It has lower row (μ(d′−d), da′−μb′) with D := da′−μb′ ≡ 1 mod (d′−d), so κ(G) = χ_μ(da′). The same proof covers cusp denominators c_F·μ with c_F | 72^∞ fixed (needed for the translates of §1): the extra factor (D/c′_F)₆ is constant on the class.

**Consequence.** θ_ρ has multiplier κ̄, so

  θ_ρ(g_a v) = χ_μ(a) · θ_ρ(L(μ₀) v).

This exponent +1 matches the coherent twist j = 1 found in 3.2/3.3. With θ in place of θ_ρ (multiplier κ) one gets χ̄_μ(a): the conjugate convention, coherent twist j = 5.

**The dual function** is θ_ρ at the cusp 1/μ₀: θ_ρ^{(μ₀)}(v) := θ_ρ(L(μ₀)v), with coefficients τ_{μ₀}(n).
- Since (μ₀, 3) = 1, the cusp 1/μ₀ is *not* Γ(72)-equivalent to ∞. **[C]**: the dual coefficients are theta coefficients at other cusps. There are finitely many: one for each class μ₀ of primary units mod 72 (and of the h/N₁ translates).
- This is exactly the sextic counterpart of the cubic Θ_t = θ(L(t_b)w) in paper78's proof.

---

## 5. The Voronoi formula and the d-sum

### 5.1 Additive form [R-]

Write v = a/μ + w with w → 0 in H³. Then g_a^{-1}(a/μ + w) = −d/μ − 1/(μ²w), and Lemma M gives θ_ρ(a/μ + w) = χ_μ(a) θ_ρ^{(μ₀)}(−d/μ − 1/(μ²w)).

The usual Mellin/Hankel argument then gives, for W smooth with compact support in (0,∞):

  Σ_m τ(m) e(am/μ) W₁(N m/Q) = χ_μ(a) · (Q/N(μ)) · [ Σ_{0≠n∈Λ} τ_{μ₀}(n) e(−nd/μ) W₁^♯(C N(n) Q/N(μ)²) + M_{μ₀}(Q,μ) ].

Here:
- W₁(t) = t^{-1/12}W(t), so that the left side is Q^{1/12}·Σ ρ e(am/μ) W.
- Λ is the frequency lattice at the cusp (cusp width 72), and C is the corresponding constant.
- W^♯ is the Hankel transform of the GL₂(C) representation with Whittaker function yK_{1/6}. By analogy with paper78 eq. (reflection-kernel) it has kernel R(t) = ∏_± Γ(1+t±1/12)/Γ(1−t±1/12). The exact Γ-factors are **not verified [G2]**; they do not affect any exponent.
- W^♯(x) ≪_A (1+x)^{-A} and W^♯(x) ≪ x^{δ} (x → 0), for some δ > 0.
- M_{μ₀} is the constant-term contribution c_{μ₀}y^{5/6}. It is independent of d, and of size ≍ Q^{-1/12}·N(μ)^{1/6}.

That last size gives Σ_m ρ e(am/μ)W ⊃ χ_μ(a) c_{μ₀} Q^{5/6} N(μ)^{-5/6}. For μ = 1 this is the Q^{5/6} main term.

**(a) Dual length:** N(n) ≲ C N(μ)²/Q · Q^ε.

### 5.2 The quadratic twist [R]

ψ = χ_μ³ is primitive, so ψ(m) = τ(ψ̄)^{-1} Σ*_{a mod μ} ψ̄(a) e(am/μ), and τ(ψ̄) = τ(χ_μ³) = √N(μ)·γ₃(μ). The d-sum is

  Σ*_a ψ̄(a) χ_μ(a) e(−nā/μ) = Σ_d χ_μ(d)^{3−1} e(−nd/μ) = χ_μ(−n)^{-2} τ(χ_μ²) = \overline{χ_μ(n)²} · √N(μ) γ₂(μ),

using χ_μ(−1)² = 1. So

> **F(μ) = ε(μ) · Q^{11/12} N(μ)^{-1} · Σ_{0≠n∈Λ} τ_{μ₀}(n) \overline{(n/μ)₃} W₁^♯(C N(n)Q/N(μ)²),  ε(μ) = γ₂(μ)/γ₃(μ),**

up to fixed phases depending on μ₀ and on the N₁-classes of §1. There is **no** main term: the n = 0 sum Σ_d χ_μ(d)² vanishes. The terms with (n,μ) ≠ 1 vanish as well.

Using paper78 eq. (signal), γ₁γ₂ = μ(μ)α(μ)G(μ) with G = \overline{χ(4)}γ₃, so

  ε(μ) = μ(μ) α(μ) \overline{χ_μ(4)} \overline{γ₁(μ)}.

This is a genuinely oscillating unimodular function (a sextic Gauss sum), not a fixed ray phase. For the large sieve it is harmless, since it is absorbed into the coefficients.

This confirms the notes' heuristic in the form "(Q/q_μ) Q^{-1/12} Σ a′(n)(n/μ)₃^{±1}", with the sign −1 (conjugate cubic symbol) in the present orientation and a′ = τ_{μ₀}.

**General twist table.** For ψ = χ_μ^j with j ≢ 0, the dual character is χ_μ(−n)^{1−j}, times τ(χ_μ^{j−1})/τ(χ_μ^{−j}). The cases are:
- j = 1: the Ramanujan sum c_μ(n), with n = 0 term φ(μ). This is coherent, and the "alternative" scenario belongs here.
- j = 2: χ_μ^{-1}, sextic.
- j = 3: χ_μ^{-2}, conjugate cubic.
- j = 4: χ_μ^{-3}, quadratic.
- j = 5: χ_μ^{-4} = χ_μ², cubic.

Cubic analogue: Patterson's θ̄₃ has multiplier χ_μ(a)^{4}, measured in §8.1. The same computation gives dual χ_μ(−n)^{4−j} = χ^{−2−j}. This is *exactly* paper78's proven local factor B_p(x) = χ_p(x)^{−j_p−2}, including the Ramanujan case j_p = 4. That is a non-trivial consistency check of the whole mechanism and of the conventions.

---

## 6. Main terms (c)

- **Squarefree μ ≠ 1.** No main term, exactly. [R]
- **General columns.** Write the odd part of λ⁴μ as μ′t² with μ′ squarefree. Then (λ⁴μ/m)₂ = (fixed ray character)·(m/μ′)₂·1_{(m,t)=1}.
  - The coprimality factor expands into cusps b/δ′ with δ′ | t. There Lemma M gives the combined character χ_{μ′}^{−2}·χ_{δ′}^{+1}.
  - The n = 0 term therefore survives only for μ′ = 1 and δ′ = 1, i.e. **μ ∈ uλ^k·□**.
  - In the reflected support μ = uλ^k n b³ (n squarefree), this means nb = □, i.e. **μ = uλ^k n⁴ c⁶**.
  - "μ a cube" (the notes' wording) is *not* sufficient: (b³/m)₂ = (b/m)₂ is still a nonprincipal twist.
- **Size.** For such μ, F(μ) = ∏_{p|t}(1 − N(p)^{-1}) · c₁ Q^{5/6} + (dual sum of length N(t)²/Q-type). Here c₁ is the constant term of θ_ρ at the cusp 0/1 for the relevant classes. Its value, possibly zero, is not computed **[C]**.
- **Totals.**
  - Model weights |β_μ| ≍ N(μ)^{-1/2} on squares: O(Q^{5/6}) per dyadic L-block.
  - True reflected weights |β_μ| ≲ |b|/√(q_n q_b³) with b = nc²: O(Q^{5/6} L^{-1/8}) per block.
  - Either way the total is O(Q^{5/6+ε}), dominated by μ = 1.
- **Interpretation.** These main terms are the Q^{5/6} pole of Σ_m ρ(m)B_m, i.e. the signal of IDEA1_NOTES §2. They must be kept as signal, not bounded as error.
- **Caveat for the LP.** Relative to the Cauchy–Schwarz scale Q^{11/12} they are Q^{-1/12}. This exceeds the large-sieve term (L/Q)^{1/2} when L < Q^{5/6} (model weights) or L < Q^{2/3} (true weights). So any LP that subtracts the signal must treat this term separately; it is not covered by [L/Q + …]^{1/2}.

---

## 7. The bilinear form and the cubic large sieve (d)

Take μ squarefree primary, N(μ) ≍ L, |β_μ| ≍ L^{-1/2}, so that Σ|β_μ|² ≍ 1.

**Step 1: separate variables.** Split μ by its class μ₀ mod 72 (finitely many) and by the N₁-classes. Then separate W₁^♯(C N(n)Q/N(μ)²) with a Mellin integral on Re t = 0, after a smooth truncation at x ≤ Q^ε. This gives

  |T| ≤ Q^{11/12+ε} L^{-1} · max_t |Σ_μ Σ_n x_μ y_n \overline{(n/μ)₃}|,

where
- x_μ = β_μ ε(μ) N(μ)^{2it}·1_{μ≡μ₀} is arbitrary, with |x_μ| = |β_μ|;
- y_n = τ_{μ₀}(n) N(n)^{-it}, for N(n) ≤ N := C L²/Q · Q^ε;
- the fixed units, λ- and 2-parts of n are moved into ray characters of μ and absorbed into x_μ.

**Step 2: Cauchy–Schwarz and T5.**

  |Σ_μ Σ_n| ≤ ‖x‖ ‖y‖ Δ^{1/2},  Δ = (L + N + (LN)^{2/3})(LN)^ε.

**Step 3: the inputs.**
- ‖x‖² ≍ 1.
- ‖y‖² ≪ N^{1+ε} (G3, **[C]**).
- Non-squarefree n (G4).

**Result.** With N^{1/2} = L/Q^{1/2},

  |T| ≪ Q^{11/12} L^{-1} · (L/Q^{1/2}) · (L + L²/Q + L²/Q^{2/3})^{1/2} = **Q^{5/12}(L + L²/Q + L²/Q^{2/3})^{1/2} · Q^ε.**

**Comparison bound.** Apply Cauchy–Schwarz in m and T6:

  |T| ≤ ‖ρW‖₂ · ((Q+L)Σ|β|²)^{1/2} Q^ε ≍ Q^{5/12}·Q^{1/2} = Q^{11/12+ε}.

**Ratio:** [L/Q + L²/Q² + L²/Q^{5/3}]^{1/2}, as claimed in IDEA1_NOTES.
- Since L < Q, the middle term is redundant.
- The (UV)^{2/3} term of Heath-Brown is the bottleneck exactly when L > Q^{2/3}.
- An "ideal" large sieve (U+V) would give (L/Q)^{1/2}.
- If L ≤ (Q/C)^{1/2}Q^{−ε}, then N < 1 and T = O(Q^{-A}). This is strictly better than the formula.

**Remark.** x_μ = β_μ ε(μ) is itself a product of cubic and sextic Gauss sums: β_μ ∝ g̃₃ and ε ∝ α·\overline{γ₁}. Heath-Brown's sieve does not use this structure. Any further saving would have to come from it.

---

## 8. Numerical checks

### 8.1 Exact cubic analogue (Patterson's θ̄₃) — `cubic_theta.py`, `exp_cubic_voronoi.py`

**The sequence.**
- a(x) = |b|·g̃₃(n) for x = nb³ ≡ 1 mod 3 and x odd. The Gauss sum uses ě; it is computed exactly from prime Gauss sums by twisted multiplicativity.
- Twisted multiplicativity was verified directly on products of two primes.
- The masks x ≡ 1 (3) and x odd are congruence conditions, so a(x) is the coefficient sequence of a finite combination of translates of θ̄₃.
- Q = 2·10⁵ gives 28 411 values of x.

**Additive twists measure the multiplier directly.** For small μ, the additive twist A(a/μ) = Σ a(x)W e(xa/μ) is dominated by the constant-term main term. The fit of A(a/μ)/A(1/μ) to χ_μ(a)^k gives **k = 4** for every μ tested (N(μ) = 7, 13, 19, 31, 37). The mean deviation is 0.27 at N(μ) = 7, against ≥ 1.25 for every other k.
- The deviations are the size expected from the dual sum: |A(a)| ≈ 360 against a random scale of 93 at N(μ) = 7.
- |A(a)| decays like N(μ)^{-2/3} = N^{-(1−1/n)}, as predicted by the constant term y^{2/3}: 360, 262, 208, 140, 155 observed, against 360, 238, 185, 133, 119 predicted.

**Multiplicative twists.** R_j = |F_j(μ)|²/Σ|aW|², by bins of N(μ)²/Q from 2^{-12} to 2^{5}:
- j = 0 (coprimality): R ≈ 160–214.
- j = 4 (= k): R ≈ 76 → 16, decaying in N(μ) as predicted by φ(μ)τ^{-1}N^{-2/3}. At N = 7, R₄ = 76 matches |Σ_a χ̄⁴A(a)|²/|τ|² ≈ (6/√7·360)²/8660 = 77.
- j = 1, 2, 3, 5: R = O(1) at every μ. There is no main term. The quadratic twist j = 3 behaves like the other non-principal j.

### 8.2 Dual length — `exp_cubic_scaling.py`, `exp_cubic_dualspec.py` (inconclusive by design)

**Scaling test.** Gauss sums of all primes up to 4·10⁵ are cached in `cubic_gauss_cache_400000.npz`. We tracked R_j(μ) for Q = 2.5·10⁴, …, 4·10⁵, with fixed μ of norm 7 to 400.
- The main-term classes scale as predicted:
  - j = 0: R grows 55 → 320, against the prediction Q^{2/3}: ×6.3.
  - j = 4: R grows 16 → 142.
- For j = 1, 2, 3, 5 there is **no** systematic decrease. R_j stays O(1), fluctuating between 0.01 and 3, even at N(μ)²/Q = 10⁻⁴.

**Dual-spectrum test.** For split μ with N(μ) = 61 to 523, take f(b) = A(1/b)χ_μ(1/b)^{-4} and Fourier-transform it over b ∈ F_p. Voronoi says this transform is supported on the residues of dual vectors n with N(n) ≲ C N(μ)²/Q.
- To capture 90% of the total energy (the r = 0 term included), one needs residues of all lattice vectors up to norm X₉₀ ≈ 0.24 N(μ), for every lattice scale c ∈ {1,2,λ,3,6,3λ,9,18}.
- That is exactly the value for an unlocalised spectrum (3.63·X₉₀ ≈ 0.9p).

**Interpretation.** These results do not contradict Voronoi, but they show that the fixed constant C is large. Both tests are consistent with C ≳ 10⁴:
- the masks x ≡ 1 (3) and x odd, the λ^{-3} frequency normalisation and the level all enter C;
- in paper78's cubic reflection, C = 81·N(c_F)² with c_F | L and (18) | L.

The regime C N(μ)²/Q ≪ 1 would need Q ≳ 10⁴·N(μ)², i.e. Q ≳ 5·10⁶ even for N(μ) = 7. Computing all Gauss sums up to that range is O(Q²) work, so it is out of reach here.

**So (a) rests on the standard theory, not on numerics.** The sextic case (§8.3) has the same problem; its level is 72, so C is presumably larger still.

### 8.3 Sextic data (Q = 3000, cross-scale estimator)

The j-table of §3.3 identifies the coherent exponent j = 1, which is the orientation used here. The data cannot test (a) or (c) for the quadratic twist, for two reasons:
1. The noise per prime is about 1–3 in ratio units, and the cross-scale estimator is biased by deterministic, scale-dependent components. The j = 0 ratio is −44, coming from the k = 1 Poisson term Q Y^{-1/2}Ŵ(√(Q/Y)), which oscillates in Y.
2. C N(μ)²/Q ≥ 1 for every prime μ when Q = 3000.

The j = 3 values (0.6–4) are consistent with "no main term, random-size dual sum".

---

## 9. Gaps and heuristic steps

**G1 [R-]. Identification of ρ with theta coefficients.** ρ is defined by D over primary s with (s,6) = 1 (squarefree only, in the numerics). It must be identified with a finite linear combination of coefficients of finitely many translates and ray-twists of θ_ρ on Γ(72).
- The 2/λ-parts of the Kubota moduli give fixed phases and ray characters.
- Restricting to squarefree s is related to the full series by the triangular relation

  D_full(w;m) = Σ_{s₂|m^∞} g̃(s₂,−m) N(s₂)^{-w} · D_sf^{R(·,s₂)}(w; m s₂⁴),

  which comes from twisted multiplicativity and χ_{s₁}(s₂)² g(s₁,−m) = g(s₁,−m s₂⁴).
- Standard, but not written out. It does not affect exponents.

**G2 [R-]. Exact Hankel kernel.** The Hankel kernel for yK_{1/6} on H³ (Γ-factors, cusp-width constant C) is not verified. Only the decay of W^♯ is used.

**G3 [G], [C]. Mean square of the dual coefficients.** The bound Σ_{N(n)≤N}|τ_{μ₀}(n)|² ≪ N^{1+ε} is needed at the cusps 1/μ₀, which are not equivalent to ∞.
- Route: Zagier-renormalised Rankin–Selberg for θ₆ at each cusp. The extra poles from |c_𝔞|² y^{5/3} lie at s < 2, so the main pole gives ≍ X.
- Not found in the literature for the 6-fold cover over Q(ω).
- Numerically supported only at ∞ (Σ|ρ|² ≍ Q^{5/6}).

**G4 [G], [C]. Non-squarefree dual indices.** Heath-Brown requires squarefree b. Write n = n₁n₂²n₃³, so that (n/μ)₃ = (n₁/μ)₃·\overline{(n₂/μ)₃}·1_{(n₃,μ)=1}. Then group by (n₂,n₃) and apply T5 to n₁.
- This costs only (log N)^{O(1)}, *provided* Σ_{N(n₁)≤X}|τ_{μ₀}(n₁n₂²n₃³)|² ≪ X^{1+ε}(N(n₂n₃))^ε uniformly.
- That needs the local (KP) theory of θ₆ at prime powers, at every cusp. Not done.
- Alternatively, extend T5 to cube-free arguments with mixed exponents.

**G5 [C]. Constant terms.** The constants c_𝔞 of θ₆ at the cusps are only needed for the *size* of the main terms in (c), where an upper bound suffices. They are not computed and may vanish for some classes.

**G6 [H]. The rank-one claim.** The identification "Y^{1/3}A_m(Y) = What·ρ(m) + noise", i.e. that the measured component is ρ, belongs to (D), not V. The present lemma is about the true ρ.

**G7. What is NOT a gap.**
- The cubic, non-principal character in (b): the combined character χ_μ^{±1−3} cannot be principal for either sign convention.
- The vanishing of the n = 0 term for squarefree μ ≠ 1.
- The exponent bookkeeping in (d), given G3 and G4.

**(P) items not addressed.** Non-squarefree rows m, where the row twist is not quadratic: for those rows, check that no column produces an effective twist j ≡ 1 in m, which would be coherent. Also the marked rows B^J, the v-integral (harmless: W(t)t^{-iv} with polynomial seminorms), and the ray class s ∈ σ.
