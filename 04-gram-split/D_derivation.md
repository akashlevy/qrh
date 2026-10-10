# Lemma D: the exceptional-frequency asymptotic of the cross-Gram, and the clean Gram bound

*Working notes, 9 Oct 2026. These notes treat Prop. `probe-gram` of `formalize_11_12/papers/paper78.tex` (lines ~8340–8560) as the source of the P_a^{1/6} term. Companion code and data are in `gram_split/D_numerics/` (§7).*

Status labels follow `V_derivation.md`:
- **[R]** rigorous, given the paper's lemmas or the cited theorems;
- **[R-]** routine but not written out in full;
- **[H]** heuristic;
- **[G]** genuine gap;
- **[N]** numerical evidence only.

---

## 0. Summary

1. **Lemma D1 (exceptional asymptotic) [R / R-].** Put Y = √(Y₁Y₂) and P = Y₁Y₂/Q, with Y₁ ≍ Y₂ and P ≥ 1. The cross-Gram satisfies

   G(Y₁,Y₂) := Σ_m Ω(q_m/Q) A_m(Y₁;v₁) \overline{A_m(Y₂;v₂)}
            = 𝒟 + Q^{5/6}(Y₁Y₂)^{-1/3} Ŵ(1/6+iv₁) \overline{Ŵ(1/6+iv₂)} K₀ + O((1+|v|)^B (Q/Y)(1 + P²/Y) Z^ε).

   - 𝒟 ≪ Q/Y is the k = 0 diagonal. It vanishes for disjoint supports.
   - K₀ is an explicit, absolutely convergent constant: a sum over the common factor C of local densities (§3).
   - The exceptional frequencies contribute the main term plus O((Q/Y)Z^ε).
   - The P²/Y part of the error is the paper's non-exceptional bound, unchanged.

   A refinement of the frequency classification is needed. Only the frequencies k = ε₀ C⁵ a⁶ carry a main term, where C = gcd(s₁,s₂) and ε₀ is a unit times an S-monomial. The other "exceptional" frequencies of the paper, k = a⁶e with e | C^∞ and e ≢ C⁵, have mean zero.

   For the numerical model (c(s) = μ²(s), rows coprime to 6, Ω = bump on (½,1)), K₀ = 6.594·10⁻³.

2. **D1 does not imply (D) by itself.**
   - D1 fixes the inner products among the vectors A(Y) for √Q ≤ Y ≲ Q^{11/16}, and shows that their excess is a single vector (rank one in the scale variable).
   - It does not identify that vector with ρ, the Kubota residue. ρ = lim_{Y→∞} Y^{1/3}A(Y)/Ŵ is determined by scales where the paper's non-exceptional bound is trivial.
   - The missing input is analytic information on the sextic Kubota series. With the functional equation (FE), (D) follows in the *strong* form ‖A − A^sig‖² ≪ (Q/Y)Q^ε, with no P²/Y term. This is the route of `LoA6.md`, and I agree with it.
   - Two additions here: the paper's own planar *additive* large sieve suffices, which removes LoA6's gap [G-b]; and D1 together with (FE) gives ‖ρ‖²_Ω = K₀Q^{5/6}(1 + O(Q^{-1/36+ε})).
   - A Rankin–Selberg unfolding is not needed.
   - Scale combinations alone cannot replace ρ. One finite-scale variant, the signal S(Y) = λA(X) with X ≈ Q^{2/3}, does get the clean Gram from D1 alone. But it moves the whole difficulty into a Voronoi estimate for A(X) (§6.5).

3. **Numerics [N].** Q = 1500, 3000, 6000, twelve scales Y = 256…524288; the first 9 scales reproduce `A_Q3000.npy` bit-for-bit.
   - **Residual.** The ρ-free cross-scale estimator of ‖A(Y) − Y^{-1/3}Ŵρ‖² is flat in Y: it equals c·Q/Y with c = 0.0169–0.0174, independent of Q, and it is 95–97% of the exact k = 0 diagonal.
   - **No excess.** The fitted Y^{1/3} excess is (−0.2 ± 8)%, (−3 ± 6)% and (+4 ± 5)% of the full theta excess. The raw Gram diagonal has the full excess.
   - **Off-diagonals.** These are ≈ 0.8·K₀Q^{5/6} at large scales, with the Q^{5/6} scaling confirmed. A brute-force evaluation of the exceptional sum on the actual s-lists confirms the density formula for K₀.

---

## 1. Setting and the exact Poisson identity

**Notation.**
- K = Q(ω), O = Z[ω], q_x = N(x).
- e(z) = exp(4πi Im z/√3). O is self-dual for this pairing, and the self-dual measure is dμ = (2/√3) dA.
- s is primary (≡ 1 mod 3), χ_s(x) = (x/s)₆, g(s,h) = Σ_{x mod s} χ_s(x) e(hx/s), and g̃ = g/√q_s.
- Ŵ(z) = ∫ W(r) r^z dr/r.

**The additive polynomial.** This is the paper's A_{m,σ,v}, up to the coefficient c:

  A_m(Y;v) = Y^{-1} Σ_s c(s) W(r_s) r_s^{-1/2+iv} g̃(s,−m),   r_s = q_s/Y.

The two coefficient settings are:
- **(M) Numerical model.** c(s) = 1 on squarefree primary s with (s,6) = 1. Rows m range over (m,6) = 1 with weight Ω(q_m/Q). This is `exp_rank1.py`.
- **(P) Paper.** c(s) = 1_{s∈σ} χ_s(b_*)/(τξ(s)) on all primary s prime to S. The row weight Ω is a smooth radial majorant of the row ball.

Write H(t) := Ω̂(η) at |η|² = t, where Ω̂(η) = ∫ Ω(|x|²) e(−xη) dμ(x). H is Schwartz in η.

**Exact Poisson identity [R].** Poisson on O gives

  Σ_m Ω(q_m/Q) g(s₁,−m) \overline{g(s₂,−m)} = Q Σ_{j∈O} F(s₁,s₂;j) H(Q q_j/(q_{s₁}q_{s₂})),   (1.1)

with F the full correlation of the paper, eq. `old-eq:2.11`.

The mask 1_{(m,6)=1} in model (M) is handled by writing it as Σ_{d|6} μ(d) 1_{d|m}. Since g(s,−dm') = \overline{χ_s(d)} g(s,−m'), the d-term is (1.1) with Q replaced by Q/q_d and an extra factor \overline{χ_{s₁}(d)} χ_{s₂}(d).

Hence

  G(Y₁,Y₂) = (Y₁Y₂)^{-1} Σ_{s₁,s₂} c c̄ W W r^{…} (q_{s₁}q_{s₂})^{-1/2} Σ_{d|6} μ(d) (Q/q_d) \overline{χ_{s₁}(d)} χ_{s₂}(d) Σ_j F(s₁,s₂;j) H((Q/q_d) q_j/(q_{s₁}q_{s₂})).

The normalisation of (1.1) was checked numerically with a Gaussian test function (`checks1.py`).

---

## 2. The correlation and a refined classification of frequencies

Write s_i = C n_i with C = (s₁,s₂). Then j = Ck, and by Lemma `full-correlation`

  F(Cn₁,Cn₂;Ck) = χ_{n₁}(k) \overline{χ_{n₂}(−k)} R(n₁,n₂) L_C(n₁,n₂;k).

In model (M), C, n₁, n₂ are pairwise coprime and squarefree, so the local formula gives

  L_C = χ_C(n₁) \overline{χ_C(n₂)} c_C(k),   c_C(k) = Σ_{δ|(C,k)} μ(C/δ) q_δ   (the Ramanujan sum).

**Fixed-modulus facts [R, checked numerically].**
- R(x,y) = χ_y(x)/χ_x(y) = ±1 is a bicharacter of the classes of x, y mod 4: the 2-adic quadratic Hilbert symbol. The primary normalisation makes the cubic part of the reciprocity exact.
- The symbols n ↦ χ_n(ζ₆), χ_n(λ), χ_n(2), with λ = 1−ω, are characters of n mod 72 on primary n.
- Both facts were checked on all primes up to norm 3·10⁴: 0 inconsistencies at modulus 4 for R, and at modulus 72 for the symbols. Modulus 2 and modulus 36 fail.

**Classification of k ≠ 0 relative to C.** The n₁-dependence of F is χ_{n₁}(k)χ_C(n₁) times fixed-modulus factors. By reciprocity, χ_{n₁}(p)^{e}χ_p(n₁) = R(p,n₁)^{e} χ_p(n₁)^{e+1}. This character is principal in n₁ exactly when

- 6 | v_p(k) at every good prime p ∤ C, and
- v_p(k) ≡ 5 (mod 6) at every p | C.

That is, kC is a sixth power up to units and S-primes. The three classes are:

- **(E) exact-exceptional:** k = ε₀ C⁵ a⁶, with
  - a ∈ O∖{0} taken mod units;
  - ε₀ ∈ E₀ = {ζ₆^u λ^i 2^j : 0 ≤ u,i,j ≤ 5} in model (M), or units times S-monomials mod sixth powers in setting (P).
  
  This decomposition of k is unique.
- **(P) partial:** k = ε₀ e a⁶, where e | C⁵ is sixth-power-free and e ≠ C⁵. The n₁-character contains Π_{p|C} χ_p^{e_p+1}, which is non-principal at some p | C. So its complete mean is zero.
- **(N) non-exceptional:** some good p ∤ C has 6 ∤ v_p(k). These are exactly the paper's non-exceptional frequencies.

The paper's "exceptional" class is (E) ∪ (P), and its count (RP/C₀)^{1/6}·q_C^ε is the count for (P) with e = 1. Only (E) can carry a main term, and the number of (E)-frequencies with q_k ≤ P/q_C is ≍ P^{1/6}/q_C, not (P/q_C)^{1/6}.

**Lemma 2.1 [R].** On (E), in model (M),

  F(Cn₁,Cn₂;Ck) = φ(C) χ_{n₁}(ε₀) \overline{χ_{n₂}(−ε₀)} R(n₁,n₂) R(n₁,C) R(n₂,C) 1_{(a,n₁n₂)=1}.

*Proof.*
- χ_{n₁}(k) = χ_{n₁}(ε₀) \overline{χ_{n₁}(C)} 1_{(a,n₁)=1}, and χ_C(n₁)/χ_{n₁}(C) = R(n₁,C).
- The n₂ factor is the same computation, conjugated.
- C⁵ | k gives c_C(k) = φ(C). ∎

---

## 3. Statement of Lemma D1 and the constant K₀

**Lemma D1.** Fix W, Ω and the arithmetic data. Let Y₁ ≍ Y₂, put Y = √(Y₁Y₂) and P = Y₁Y₂/Q, and assume √Q ≲ Y ≤ Q, all in fixed polynomial ranges of Z. Then

  G(Y₁,Y₂) = 𝒟 + 𝓔 + 𝓟 + 𝓝,

with:
- 𝒟 = (k = 0 term) = (Q H(0)/2) (Y₁Y₂)^{-1} Σ_s |c(s)|² W(r₁)W(r₂) r₁^{-1/2+iv₁} r₂^{-1/2−iv₂} φ(s)/q_s ≪ Q/Y. It vanishes if the two annuli are disjoint.
- 𝓔 = Q^{5/6}(Y₁Y₂)^{-1/3} Ŵ(1/6+iv₁) \overline{Ŵ(1/6+iv₂)} K₀ + O((1+|v|)^B (Q/Y) Z^ε).
- 𝓟 ≪ (1+|v|)^B (Q/Y) Z^ε.
- 𝓝 ≪ (1+|v|)^B (Q/Y)(P²/Y) Z^ε. This is the paper's bound, verbatim.

**The constant, model (M).**

  K₀ = ĝ(0) Σ_{d∈{1,λ,2,2λ}} μ(d) q_d^{-5/6} Σ_{C sqfree primary, (C,6)=1} φ(C) q_C^{-3} Dens(C) Ẽ(C,d).

The ingredients are:
- **Archimedean factor.** ĝ(0) = κ_K · (2/√3)^{2/3} π^{2/3} (Γ(1/6)/Γ(5/6)) · ∫₀^∞ Ω(u) u^{-1/6} du, with κ_K = π/(3√3) = Res ζ_K.
- **Density of n-pairs.** Dens(C) = Dens(1) Π_{p|C} (1 + 2(q_p−1)/q_p²)^{-1}, where
  Dens(1) = κ_K² (2/3)² (3/4)² Π_{p∤6} (1 − 5q_p^{-2} + 6q_p^{-3} − 2q_p^{-4}).
- **Class average.** Ẽ(C,d) is an average over (x₁,x₂) ∈ G₇₂², where G₇₂ is the group of 432 primary classes mod 72 coprime to 2:

  Ẽ(C,d) = E_{x₁,x₂}[ 1{χ_{x₁}(ζ₆)=χ_{x₂}(ζ₆)} · \overline{χ_{x₂}(−1)} · \overline{χ_{x₁}(d)} χ_{x₂}(d) · R(x₁,x₂)R(x₁,C)R(x₂,C)
        · (2/3)/(1 − 3^{-1/6} χ_{x₁}(λ)\overline{χ_{x₂}(λ)}) · (3/4)/(1 − 4^{-1/6} χ_{x₁}(2)\overline{χ_{x₂}(2)}) ].

  The geometric factors are the ε₀-sums over λ^i 2^j, 0 ≤ i,j ≤ 5, weighted by q_{ε₀}^{-1/6}. The indicator is the sum over the six units, which cancels the 1/6 from counting a mod units.

**Values for Ω = bump(·;½,1)** (`D_numerics/k0.py`, `k0parts.py`):

| quantity | value |
|---|---|
| ∫Ω u^{-1/6} | 0.201503 |
| ĝ(0) | 1.418352 |
| Dens(1) | 0.067190 |
| Ẽ(1,d), for q_d = 1, 3, 4, 12 | 0.131446, 0.109453, 0.104329, 0.086873 (all real) |
| C = 1 inner sum | 0.065723 |
| share of K₀ from C = 1 | 95.0% |
| share from each of the two primes of norm 7 | 1.3% |
| **K₀** | **6.594·10⁻³** |

**Setting (P).** The structure is identical, but the local data change:
- the class averages run over the ray class σ with the characters χ_·(b_*)/ξ;
- ε₀ runs over the larger set of S-monomials;
- non-squarefree s add, at p | C with multiplicity c, the condition v_p(k) ≡ −c (mod 6) and the local factors of eq. `correlation-local`.

So K₀ changes by an absolutely convergent Euler-type correction. [R-, not computed.]

**Consistency.** K₀ is the theta mean-square constant: under (FE), Σ_m Ω|ρ(m)|² ~ K₀Q^{5/6} (§6.3).

---

## 4. Proof of D1

### 4.1 Frequencies (E): the sum over a first [R]

Fix C, ε₀, d and (n₁,n₂). The k-sum over (E) is

  (1/6) Σ_{a≠0, (a,n₁n₂)=1} H(|a|^{12}/T),   T = q_{s₁}q_{s₂} q_d/(Q q_C⁶ q_{ε₀}).

**Lattice lemma.** Put g(y) = H(|y|^{12}), a Schwartz function on C, and L = T^{1/12}. Then

  Σ_{a≠0, (a,n)=1} g(a/L) = ĝ(0) T^{1/6} φ(n)/q_n + O(Σ_{b|n} min(1, T^{1/6}/q_b)),

where ĝ(0) = ∫ g dμ.

*Proof.* Write the coprimality condition with Möbius: the sum is Σ_{b|n} μ(b) Σ_{a'≠0} g(ba'/L). Poisson on the lattice bO, whose μ-covolume is q_b, gives (L²/q_b) ĝ(0) − g(0) + O_A((L/|b|)^{-A}) when |b| ≤ L. When |b| > L, both the sum and L²/q_b are small. ∎

**The (q_{s₁}q_{s₂}/Q)^{1/6} main term.** The archimedean constant comes from

  ĝ(0) = (2/√3) π ∫₀^∞ H(u⁶) du = (π/(3√3)) ∫₀^∞ H(t) t^{-5/6} dt.

The 2-D Fourier transform of |ξ|^{-5/3} is π^{2/3}(Γ(1/6)/Γ(5/6))|x|^{-1/3}. Accounting for the self-dual pairing gives

  ∫₀^∞ H(t) t^{-5/6} dt = (2/√3)^{2/3} π^{2/3} (Γ(1/6)/Γ(5/6)) ∫₀^∞ Ω(u) u^{-1/6} du,

checked numerically to 6 digits (`checks1.py`). So the a-sum is

  (1/6) ĝ(0) (q_{s₁}q_{s₂}/Q)^{1/6} q_C^{-1} q_{ε₀}^{-1/6} q_d^{1/6} · φ(n₁n₂)/q_{n₁n₂} + O(τ(n₁n₂) min(1,T^{1/6})).

The m-row weight enters only through ∫Ω u^{-1/6}. This is the density N(m)^{-1/6} of |ρ(m)|².

*Numerical caveat.* For the Gevrey bump Ω, H(t) decays slowly: H(50) ≈ 4·10⁻³ and H(200) ≈ 3·10⁻⁴. So g is a chirp, and the lattice asymptotic is reached only for T^{1/12} far above the chirp frequency. At T ≤ 10⁹ the exact lattice sum differs from its main term by up to ±40%. At T = 10¹² the difference is 2%. This is irrelevant for the theorem, because it is an O(1) error per (C,n) with an Ω-dependent constant. But it shapes the finite-Q numerics (§7.3).

### 4.2 The n-sums [R-]

Insert Lemma 2.1 and the a-sum main term. The remaining (n₁,n₂)-sum has:
- smooth weights W(r_i) r_i^{-1/2±iv_i} (q_C q_{n_i})^{1/6−1/2};
- the arithmetic weight μ²μ² 1_{(n₁n₂,6C)=1} 1_{(n₁,n₂)=1} (φ(n₁)φ(n₂)/(q_{n₁}q_{n₂})) · ψ(x₁,x₂), where ψ is a function of the classes x_i of n_i mod 72.

Expand ψ in characters θ₁⊗θ₂ of G₇₂². The two-variable Dirichlet series is

  L(z₁,θ₁) L(z₂,θ₂) H_θ(z₁,z₂),

where H_θ is an absolutely convergent Euler product in Re z_i > 1/2. The local factor is 1 + (1−1/q)(θ₁(p)q^{-z₁} + θ₂(p)q^{-z₂}) for p ∤ 6C.

Only θ₁ = θ₂ = 1 has a pole at (1,1). Its residue is Dens(C)·E_{x₁,x₂}[ψ]. Shifting to Re z_i = 1/2+ε gives relative error O((q_C/Y)^{1/2−ε}).

The integrals of the smooth weights produce Ŵ(1/6+iv₁)\overline{Ŵ(1/6+iv₂)}, because W(r) r^{-1/2+iv} r^{-1/3} integrates to Ŵ(1/6+iv). The powers collect to Q^{5/6}(Y₁Y₂)^{-1/3}.

### 4.3 Collecting (E): the main term and the ranges where extraction fails

Errors are in units of Q/Y.

| range or piece | what happens | contribution |
|---|---|---|
| q_C > P^{1/6}Z^ε | T < Z^{-6ε}, so every a ≠ 0 has H(·) ≪ Z^{-A}. The actual contribution is negligible. Truncating the K₀ series costs Σ_{q_C>P^{1/6}} (Q/Y)P^{1/6}q_C^{-2} ≪ Q/Y. | O(Q/Y) |
| lattice error in a, i.e. b \| n₁n₂ with q_b > T^{1/6} (main-term extraction fails) | O(τ min(1,T^{1/6})) per (C,n₁,n₂). Sum: (Q/Y) Σ_C q_C^{-1} min(1,(P/q_C⁶)^{1/6}) Z^ε | O((Q/Y)Z^ε) |
| n-sum error | (Q/Y)P^{1/6} Σ_C q_C^{-2}(q_C/Y)^{1/2}. Note P^{1/6}Y^{-1/2} ≤ Y^{-1/6} since P ≤ Y². | O(Q/Y) |
| ε₀ classes with E(C,ε₀,d) = 0 | Fixed-modulus non-principal characters; only the n-sum error. | O(Q/Y) |
| coprimality (n₁,n₂) = 1 | Inside the two-variable Euler product; no separate Möbius d is needed in model (M). In setting (P), use the paper's Möbius d; the tail D₀ > Z^ε costs (Q/Y)P^{1/6}Z^{-ε}. | O((Q/Y)Z^ε) |
| P small | P < 1: T < 1 for all C, so there is no main term, consistent with the large sieve. 1 ≤ P ≤ Z^{6ε}: the main term is ≤ (Q/Y)Z^ε, so D1 holds trivially. Extraction is meaningful only for P ≥ Z^δ. | — |

### 4.4 Frequencies (P) [R-]

Again sum over a first; the a-sum depends on (a,C) only through c_C(k). Then:

- **Small C.** The n₁-sum has a non-principal character of conductor dividing 72·rad C. Convexity for Hecke L-functions on Re z = 1/2+ε gives ≪ (Y/q_C)^{1/2+ε} q_C^{1/4+ε}. The 1/L(2z) factor from μ² is harmless there.
- **Trivial bound.** |c_C(k)| ≤ q_{rad e}·q_{(C,a)}, and Σ_{e|C⁵} q_{rad e} q_e^{-1/6} ≪ q_C^{5/6+ε}. This gives (Q/Y)P^{1/6} q_C^{-4/3+ε} per C.

Combining, with the split at C* = Y^{2/3}:

  𝓟 ≪ (Q/Y) P^{1/6} Z^ε [ Y^{-1/2} Σ_{q_C≤C*} q_C^{-7/12} + Σ_{q_C>C*} q_C^{-4/3} ] ≪ (Q/Y) P^{1/6} Y^{-2/9} Z^ε.

This is ≤ (Q/Y)Z^ε, because P^{1/6}Y^{-2/9} = Y^{1/9}Q^{-1/6} ≤ 1 for Y ≤ Q^{3/2}. The lattice errors of the (P) terms total (Q/Y)Z^ε, as in §4.3.

In setting (P) the n-sums are unrestricted, so the contour can be moved to Re z = ε. The saving is then larger.

### 4.5 Frequencies (N) [R]

These are exactly the paper's non-exceptional frequencies, and its four-dimensional periodic Poisson bound gives (Q/Y)(P²/Y)Z^ε. For unequal scales Y ≤ X the same argument gives

  (Q/√(YX))(1 + P_{YX}²/Y) Z^ε,   P_{YX} = YX/Q,

because both residual columns must exceed the period.

### 4.6 Remarks

- D1 is non-vacuous only when P²/Y ≪ P^{1/6}Z^{-δ}, i.e. Y ≤ Q^{11/16−δ'}. This is the paper's regime: there Y'/P^{11/6} ≫ Z^{1/12}.
- In the numerics, P²/Y ≫ 1 everywhere. So the numerics cannot test D1's error terms. They test the stronger statement of §6.2.

---

## 5. Corollary D1′: Gram bound for scale differences [R-, given D1 extended to unequal scales via §4.5]

Let λ = u(Y)/u(Y′), where u(Y) = Y^{-1/3}Ŵ(1/6+iv) (with the same W). Then

  ‖A(Y) − λA(Y′)‖²_Ω = 𝒟_Y + |λ|²𝒟_{Y′} + O((Q/Y)(1 + P²/Y + (Y/Y′)^{1/6}(1 + Y′²Y/Q²) + (Y/Y′)^{1/3}(1 + Y′³/Q²)) Z^ε).

The rank-one main terms cancel exactly. In particular, for Y ≤ Y′ ≤ Q^{2/3}:

  ‖A(Y) − λA(Y′)‖² ≪ (Q/Y)(1 + P²/Y) Z^ε.

This is the improved Gram bound, provable with the paper's tools plus §4. Through the low separation, however, A(Y) − λA(Y′) is the two-probe combination I(X,Y,Z) − λ^{1/3} I(X/λ, λY, Z) of `IDEA1_NOTES` §2. It kills the signal, so on its own it proves nothing about I.

---

## 6. Does D1 suffice for (D)?

### 6.1 No: the identification problem [R, as a statement about what D1 contains]

Put V(Y) := A(Y)/u(Y). D1 and D1′ say that, for √Q ≤ Y ≲ Q^{11/16},

  V(Y) = U + w(Y),

- with ⟨U,U⟩ = K₀Q^{5/6}, a single vector U independent of Y;
- with ‖w(Y)‖² ≍ QY^{-1/3}, and the w's at different scales nearly orthogonal;
- up to the errors of D1.

(D) asks that U = ρ to accuracy ‖U − ρ‖² ≪ QY^{-1/3}(1 + P²/Y). Here ρ = lim_{Y→∞} V(Y), which exists by Kubota.

This is a statement about V(Y) for Y → ∞, outside the range where the paper's Poisson machinery is effective. The obstruction is precise.
- **One-sided Poisson.** Pair a scale Y in range with Y₂ → ∞, and Poisson only in the large column n₂, whose period is ≪ q_C q_k ≲ R Y Y₂/Q. This kills every (N)-frequency with q_C q_d R ≤ Q/(Y Z^ε), uniformly in Y₂.
- **What remains.** The complementary range q_C q_d R > Q/Y, where n₁ is shorter than P, is bounded only trivially. Its trivial bound grows like Y₂^{1/3}·Y^{3/2}Y₂^{1/2}/Q in V-normalisation.
- **Where ρ − U lives.** This range is exactly where ‖ρ‖² − K₀Q^{5/6} and ⟨A^rest,ρ⟩ are decided.

The task's "equivalently" also needs a correction. Given D1, (D) is equivalent to

  |u|²(‖ρ‖² − K₀Q^{5/6}) + 2 Re ū⟨ρ, A^rest(Y)⟩ = O((Q/Y)(1 + P²/Y)Z^ε).

It is not equivalent to each term separately matching to relative accuracy P^{-1/6}. Under (FE), ⟨A^rest,ρ⟩ is only known to be ≪ (Q/Y)P^{1/12}/|u|, and ‖ρ‖² is only known to match to Q^{-1/36} (§6.3). (D) still holds, because it is proved directly.

### 6.2 What suffices: the Kubota functional equation [R given (FE); this is `LoA6.md`]

**(FE).** D(w;h) = Σ_s c(s) g̃(s,h) q_s^{-w} has the following properties.
- Meromorphic continuation to Re w > 1/2 − δ. Its only pole in Re w > 1/2 is simple, at w = 2/3, with residue ρ(h).
- Polynomial growth in vertical strips.
- A functional equation w ↔ 1−w. Its conductor is ≍ q_h(1+|t|)^{O(1)}. Its dual side is a finite combination of Gauss-sum series of the same shape, at the other cusps. All of this holds uniformly in h.
- Sources: Kubota; Patterson; Brubaker–Bump, Crelle 2006. Gaps [G-a] and [G-c] of `LoA6.md` remain.

**Proof of the strong (D) from (FE).** Shift A_m = Y^{-1}(2πi)^{-1}∫ Ŵ(w−1/2+iv) Y^w D(w;−m) dw to Re w = 1/2+ε. Cauchy–Schwarz in t gives

  ‖A^rest‖² ≪ Y^{-1+2ε} sup_t (1+|t|)^{-B} M(t),   M(t) = Σ_m Ω|D(1/2+ε+it;m)|².

To bound M(t), take an approximate functional equation of length X = X′ ≍ (Q(1+|t|)^c)^{1/2}. Each piece has the form Σ_{q_s≤X^{1+ε}} b_s g(s,−m) q_s^{-1/2}, with |b_s| ≤ q_s^{-1/2-ε}.

**Addition 1.** The paper's planar additive large sieve (Lemma `planar-additive-sieve`) applies directly. Use the points x/s mod O, for primary s and units x mod s. These are distinct reduced fractions with separation ≥ X^{-1}, so

  Σ_m |·|² ≤ (Q + X²) Σ_s |b_s|² φ(s)/q_s ≪ Q^{1+ε}(1+|t|)^{O(1)}.

This uses only the definition g(s,−m) = Σ_x χ_s(x) e(−mx/s). So no sextic multiplicative large sieve and no reduction for (s,m) > 1 or non-squarefree s is needed; this removes gap [G-b] of `LoA6.md`. The dual pieces are additive character sums of the same shape.

**Addition 2.** The approximate functional equation has a residue term ≍ X^{1/6}ρ(m), contributing X^{1/3}‖ρ‖² to M. Bound it by bootstrap:
- ‖ρ‖² ≤ 2|Ŵ|^{-2}(Y^{2/3}‖A(Y)‖² + Y^{2/3}‖A^rest(Y)‖²).
- At Y = Q^{1/2+δ}, the paper's Gram bound gives Y^{2/3}‖A‖² ≪ Q^{5/6+O(δ)}.
- The A^rest term is ≪ Y^{-1/3}(Q + Q^{1/6}‖ρ‖²)Q^ε, and Y^{-1/3}Q^{1/6} = Q^{-δ/3} < 1.
- Hence ‖ρ‖² ≪ Q^{5/6+ε} and M(t) ≪ Q^{1+ε}(1+|t|)^{O(1)}.

**Result.**

  ‖A(Y) − Y^{-1/3}Ŵ(1/6+iv)ρ‖²_Ω ≪ (1+|v|)^B (Q/Y) Q^ε   for all Y ≥ 1.

There is no P²/Y term.

### 6.3 Corollary: K₀ is the theta mean-square constant [R given (FE) and D1]

Combine D1 at the scale Y with (D). The rest terms satisfy ‖A^rest‖² ≪ Q^{1+ε}/Y and |u⟨ρ,A^rest⟩| ≤ |u|‖ρ‖(Q/Y)^{1/2}. Hence

  |u|² |‖ρ‖² − K₀Q^{5/6}| ≪ (Q/Y)(1 + P²/Y)Z^ε + |u| Q^{5/12}(Q/Y)^{1/2}.

Take Y = Q^{2/3}; Y = Q^{23/34} is slightly better. This gives

  Σ_m Ω(q_m/Q) |ρ(m)|² = K₀ Q^{5/6} (1 + O(Q^{-1/36+ε})).

Note Q^{-1/36} ≈ 0.8 at Q ≈ 3000. So a 20% secondary effect at numerical sizes (§7.3) is not excluded.

### 6.4 The Rankin–Selberg route (not needed; sketch) [H]

- **‖ρ‖².** Unfold the regularised inner product of |θ₆|² against the Eisenstein series on the cover of level 72·(S). The pole at s = 5/6 gives K₀'s analogue, and the next singularity gives the error.
- **⟨A(Y),ρ⟩.** A(Y) is a smoothed truncation of the Fourier coefficients of the Kubota Eisenstein series. Pairing it with θ̄₆ amounts to Voronoi for θ₆ at the cusps d/s with q_s ≍ Y:
  - The constant term of θ₆ at d/s carries the Kubota multiplier χ_s(d)^{±1} (sign as in `V_derivation.md` §4). It cancels the character χ_s(d) in g(s,·), so the d-sum collapses to φ(s). This reproduces exactly the (E) main term with C = 1. The C ≠ 1 terms come from non-primitive cusps.
  - The dual sum has length q_s²/Q ≍ P. Its large sieve over s gives the rest.
- Both pieces need the automorphy of θ₆ at all cusps. That is the same input as (FE), with more bookkeeping, so §6.2 is preferable.

### 6.5 Other definitions of the signal component

**(a) Projection onto the span of limits.** The limit lim Y′^{1/3}A(Y′) is Ŵρ_σ. Over the finitely many ray classes and twists, the span is the finite family {ρ_σ}. So this is the same as A^sig, and it needs the same (FE).

**(b) Pure scale combinations.** These are D1′. They are provable, but they cannot be used alone. In I = Q^{-1/2}⟨A(Y),B̄⟩ = Q^{-1/2}[⟨A(Y)−λA(Y′),B̄⟩ + λ⟨A(Y′),B̄⟩], the second term is the probe at the geometry (X/λ, λY). By the scale invariance of `IDEA1_NOTES` §2, its Cauchy–Schwarz bound equals the original one exactly: λ²·‖A(Y′)‖²'s excess is Q^{5/6}Y^{-2/3} again. Telescoping Y′ → ∞ would need D1′ beyond Q^{2/3}, where the (N) bound is trivial.

**(c) Finite-scale signal (the most "elementary" option) [R for the Gram part; H for the pairing].** Take X = Q^{2/3}Z^{-δ} and define

  S(Y) := (u(Y)/u(X)) A(X).

By D1 at the scale pairs (Y,Y), (Y,X), (X,X), together with §4.5 for the unequal-scale (N)-error:

  ‖A(Y) − S(Y)‖² ≪ (Q/Y)(1 + P²/Y) Z^ε   for Y ≤ X.

This is the improved Gram from the paper's tools alone. It requires only that X³ ≤ Q² and X²Y ≤ Q².

What remains is a non-CS bound for ⟨A(X),B⟩, i.e. a "(V_X)" for the finite vector V(X) instead of ρ. A plausible elementary route:
1. **Poisson in m.** The modulus is s·μ. The product of the Gauss sum g̃(s,1) with the dual Gauss sum of χ̄_s is |g̃(s,1)|² = 1, so the sextic phases cancel.
2. **Character sum in s.** What remains is Σ_s χ_k(s) χ_μ(s)² over s ≍ X, with conductor ≍ q_k q_μ, where q_k ≲ X q_μ/Q.
3. **Poisson in s.** Its dual length is q_k q_μ/X ≍ q_μ²/Q, the length predicted by (V). Its dual Gauss sums carry the cubic structure.

If this double Poisson can be carried through with the paper's tools, then both (FE) and the metaplectic automorphy in (V) are avoided. I have not done it. [G]

---

## 7. Numerics

All of this is in `D_numerics/`. Run the scripts from that directory with the session venv.

### 7.1 Data

`expA.py` is a faster re-implementation of `exp_rank1.py`, with the same conventions:
- Gauss sums by primitive-root tables (`fastgauss.py`), checked against `zomega.Prime.gauss` to 6·10⁻¹⁶;
- 12 scales Y_i = 256·2^i, i ≤ 11, so s ≤ 2²⁰ (295,067 squarefree s);
- Q = 1500, 3000, 6000, with 1248, 2526 and 5052 rows;
- the first 9 scales of Q = 3000 reproduce `A_Q3000.npy` to 1.7·10⁻¹⁷;
- runtime ≈ 4 min.

### 7.2 Checks of the derivation

- **Poisson normalisation (1.1).** Checked with a Gaussian test function, to relative error 10⁻¹¹.
- **ĝ(0) formula.** Agrees with a direct Hankel-transform quadrature to 6 digits.
- **Brute-force exceptional sum** (`brute_exc.py`). This sums C = 1 over the actual s-lists with exact symbols and the exact lattice sum in a. Compare it with the density formula of §3/§4.2, which is "exact-lattice": it keeps the exact Λ(T) but replaces the n-sums by densities (`exactlat.py`). Both are in units of ‖ρ‖², i.e. G/(u_iu_j), at Q = 3000:

| (i,j) | brute (E), C = 1 | density formula (E), all C | data ⟨V_i,V_j⟩ |
|---|---|---|---|
| (0,1) | 3.36 | 2.46 | −0.71 |
| (0,2) | 5.59 | 4.41 | 10.33 |
| (0,3) | 5.41 | 5.49 | 10.14 |
| (1,2) | 5.87 | 5.49 | 6.42 |
| (1,3) | 6.90 | 6.43 | 5.71 |
| (2,3) | 3.92 | 3.93 | 5.79 |
| (2,4) | 4.88 | 4.95 | 7.42 |
| (3,4) | 4.82 | 4.60 | 8.12 |

  - The density formula, and hence K₀, is confirmed once the n-lists are long (i,j ≥ 2: within 1–5%).
  - At small scales, the data differ from (E) by O(1) in these units. This is the non-exceptional contribution, which is large there because P²/Y ≫ P^{1/6}.

### 7.3 Off-diagonal Gram (`compare.py`)

Here ⟨V_i,V_j⟩ estimates ‖ρ‖² + ξ_i + ξ_j, where ξ_i = ⟨w_i,ρ⟩.

| Q | K₀Q^{5/6} | mean ⟨V_i,V_j⟩, all pairs | mean, top-4 block | fit data = α·(exact-lattice pred.), pairs with i ≥ 5 | corr(data, pred.) |
|---|---|---|---|---|---|
| 1500 | 2.92 | 2.47 | 2.25 | α = 0.77 | 0.61 |
| 3000 | 5.21 | 5.35 | 4.23 | α = 0.79 | 0.74 |
| 6000 | 9.28 | 9.03 | 7.29 | α = 0.77 | 0.77 |

- The off-diagonal scales like Q^{5/6}: α is the same at all three Q.
- The oscillation of the exact lattice sum predicts which (i,j) are high or low, and the data follow it (corr ≈ 0.7).
- At large scales the data are ≈ 0.8× the exceptional prediction. This is consistent with ‖ρ‖² = K₀Q^{5/6}(1 + O(Q^{-1/36})) of §6.3, and with a non-exceptional or ξ-contribution of ≈ −20%. The numerics cannot separate ‖ρ‖² from ξ, which is the same gauge freedom as in §6.1.

### 7.4 The residual test (`final_numerics.py`, `diag_exact.py`)

Write V_i = Y_i^{1/3}A_i/Ŵ = ρ + w_i. Two estimators of ‖A_i − u_iρ‖² = u_i²‖w_i‖² are used. Both are unbiased provided ⟨w_i,w_j⟩ ≈ 0 for i ≠ j. Neither needs ⟨w_i,ρ⟩ = 0.
- **(T)** u_i² · mean over ordered j ≠ k, both ≠ i, of ⟨V_i − V_j, V_i − V_k⟩. No ρ̂ is needed.
- **(R)** Put ρ̂ = mean of V_j over the top 4 scales. The estimator is u_i²(‖V_i − ρ̂‖² − (1/16)Σ_{j∈top}‖w_j‖²_est), for i outside the top 4.

Error bars are from 300 bootstrap resamples of the rows m.

**Results, Q = 3000.** Columns are ×Y/Q. "Excess if present" is the term u²K₀Q^{5/6}·Y/Q that a surviving theta excess would add.

| Y | P^{1/6} | raw G_ii | resid (T) | resid (R) | exact k = 0 diagonal | excess if present |
|---|---|---|---|---|---|---|
| 256 | 1.9 | 0.0169 | 0.0156(5) | 0.0156(5) | 0.0161 | 0.0008 |
| 2048 | 3.8 | 0.0197 | 0.0172(6) | 0.0173(5) | 0.0180 | 0.0017 |
| 16384 | 7.7 | 0.0190 | 0.0169(6) | 0.0175(6) | 0.0178 | 0.0033 |
| 65536 | 12.2 | 0.0227 | 0.0196(7) | — | 0.0179 | 0.0053 |
| 262144 | 19.3 | 0.0253 | 0.0184(7) | — | 0.0179 | 0.0084 |
| 524288 | 24.3 | 0.0249 | 0.0156(6) | — | 0.0179 | 0.0106 |

**Weighted fit** resid·Y/Q = a + b·Ŵ²Y^{1/3}/Q. Here b is the surviving excess in units of ‖ρ‖²; b/(K₀Q^{5/6}) = 1 would mean the full excess survives.

| Q | a (T) | b/(K₀Q^{5/6}) (T) | a (R) | b/(K₀Q^{5/6}) (R) | raw-diagonal fit b/(K₀Q^{5/6}) |
|---|---|---|---|---|---|
| 1500 | 0.01692(38) | −0.002 ± 0.080 | 0.01622(58) | 0.37 ± 0.23 | 0.79 |
| 3000 | 0.01728(28) | −0.030 ± 0.061 | 0.01739(44) | −0.16 ± 0.18 | 0.75 |
| 6000 | 0.01690(21) | +0.037 ± 0.053 | 0.01713(30) | −0.11 ± 0.14 | 0.79 |

- The residual divided by the exact per-scale k = 0 diagonal has mean 0.950, 0.969 and 0.961 for the three Q. This holds scale by scale, including the small-scale wiggles that the diagonal reproduces.
- So ‖A − A^sig‖² ≈ 0.96·𝒟 ≍ Q/Y, flat from Y = 256 to 5·10⁵, i.e. P²/Y up to ~10⁸. The non-exceptional frequencies contribute only ≈ −4% net.

### 7.5 Caveats [N]

- Both estimators assume that the rests at disjoint scales are nearly orthogonal. So the data establish "A(Y) = u(Y)·(one common vector) + rests of norm² ≈ 0.96·𝒟". The identification of the common vector with ρ is the limit statement of §6.1.
- The numerics are in the regime P²/Y ≫ 1, outside D1's range. They support the strong (FE) form of §6.2, not D1's error terms.
- The 20% gap between the large-scale off-diagonals and K₀Q^{5/6} is not resolved. It could be:
  - the slow Q^{-1/36} convergence of §6.3;
  - non-asymptotic lattice effects;
  - (P)/(N) contributions at finite size.

---

## 8. Gaps

- **[R-]** The n-sum mean values of §4.2 and the (P) bound of §4.4 are sketched with exponents but not written out line by line.
- **[R-]** Setting (P): the local factors for non-squarefree s and the ray-class averages. Only the structure is given; K₀ is not computed there.
- **[G]** (FE) for the sextic Kubota series with the paper's ray twists, uniform in m. This is `LoA6.md` [G-a] and [G-c]; [G-b] is removed by §6.2. It is the only analytic input needed for (D).
- **[G]** The finite-scale (V_X) of §6.5(c), which would avoid (FE) altogether. Not attempted.
