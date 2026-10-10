# Attack B: using the Gauss-sum structure of the weights to prove (V′_β)

*9 Oct 2026. This note builds on `TARGET_VLS.md`, `REPORT.md` (CORRECTION banner), `wf/V_gaps.md` (Lemma Z), `wf/AUDIT_chain.md`, `wf2/A2_theory.md` and `wf2/A2_numerics.md` (with their referee reports), `V_derivation.md`, `paper78.tex`, and Dunn–Radziwiłł, arXiv:2109.07463 (`scratchpad/dr.txt`; read for this note: §1.3–1.4, Lemma 5.1, Lemma 5.2, Lemma 5.4, Prop. 5.2).*

**Scripts in `gram_split/vls/`:**
- `check_dr_transform.py`: N1, the Dunn–Radziwiłł twist transform and the (j, j+2) root-number rule.
- `check_diag.py`: N2, sharpness of the completed-row moment for the true Patterson column weights.
- `lp_B.py`: LP with a uniform saving s₀, and with the "ζ-opening" route of §4.
- `lp_mix.py`: LP with a different route in each low-side block.

**Labels:**
- **[R]** rigorous, given the cited results.
- **[R-]** routine but not written out.
- **[S]** taken from a secondary or literature statement without re-derivation.
- **[H]** heuristic.
- **[G]** gap.
- **[N]** numerical evidence only.

---

## 0. Verdict

**(V′_β) is not proved. Attack B gives no saving: s₀ = 0 at every M ∈ [0.80, 0.86].** With s₀ = 0 the LP returns **θ = 0.874957**, which is Liu's value (`lp_B.py`). The attack fails for a structural reason that can be stated precisely.

> **The Gauss-sum structure of β is already fully used by the paper's completed-row moment** (‖B⁰‖² ≪ Q^{1+ε}, i.e. reflection from Z columns to L columns). There is nothing left to extract from it.
> - Any argument that uses only β's structure plus norms of π is capped at about 0.89·Q^{11/12}. This is Lemma D (§3.2), with numerics N2.
> - So a proof of (V′_β) has to use the automorphy of π = ρζ.
> - For this π, the only automorphy available (§1.5, §5) is the sextic A₂ residue of `wf2/A2_theory.md`. That reduces (V′_β) to the beyond-the-edge large sieve (V′_LS). Attack B adds no tool for it.

**Results of this note:**

| | statement | status |
|---|---|---|
| **B1** | **Parity lemma.** The quadratic kernel ψ_μ(m) = χ_m(μ)³ can be absorbed only by Gauss sums of odd index. The column side has only γ₂ (cubic), so nothing is absorbed there. Twisted multiplicativity produces only even kernel exponents 2j, so it can neither create nor remove the kernel. The only non-trivial rewrite is γ₁(m)ψ_μ(m) = g̃₆(m, μ³), which is the A₂ object. In particular the task's item (ii), "move conj γ₂(m) onto the μ side", always leaves a kernel χ_m(μ)^{odd}. | [R] |
| **B2** | **Dunn–Radziwiłł Lemma 5.2 is the paper's reflection.** For a twist ψ = χ_r^j of the cubic theta, the DR Kloosterman transform equals χ_r(λ⁸)^j τ_j τ_{j+2} χ_r(u)^{−(j+2)}. This reproduces the paper's ω_{p,j} and all three cases of B_p, which is an independent derivation of Lemma Z. Hence "Type I via DR Voronoi in μ, uniform in the modulus m" is un-reflection (back to A ≍ Z). Type I with a fixed divisor d only makes the dual longer. | [R] algebra; [N] to 1.8·10⁻¹⁵ |
| **B3** | **Diagonal obstruction (Lemma D).** Every μ-side manipulation followed by Cauchy–Schwarz or Hölder in m gives at least ‖π‖·‖B⁰‖. This covers Type II splits μ = μ₁μ₂ with Heath-Brown's cubic large sieve, and DR's corrected cubic sieve. For the true Patterson weights, ‖B⁰‖² = (0.92–0.96)·Q‖β‖², B⁰_m has Gaussian moments, and random-phase β gives the same values. Type II on μ yields only the pointwise bound \|B⁰_m\| ≪ L^{1/3+ε}, which loses L^{1/3}. | [R] logic; [N] values |
| **B4** | **Type I/II in m is unavailable.** ρ has no (twisted) multiplicativity at squarefree indices (Kazhdan–Patterson non-uniqueness). The Möbius factor in ζ_m = μᾱ·conj γ₂ is the zero-detecting information, so using its cancellation is circular. | [R-]/[S] |
| **B5** | **"Opening" ζ_m.** Patterson's identity τ₃ = γ₂ gives an exact two-sided approximate functional equation of γ₂(m) in cubic characters of modulus ≍ N(m)^{1/2}. This makes π·ψ_μ Voronoi-able, but multiplies the column modulus by Q^{1/2} per opened Gauss sum. ζ_m ∝ ᾱ²γ₂(m)² needs two openings: effective modulus R = LQ > Q, so there is no saving for any M. Even an (illegitimate) single opening saves only for M < 3/5, and its LP value is 0.9605. γ₁(m) cannot be opened, since τ₆ ≠ γ₁. | [R-]; LP |
| **B6** | **Rigidity of the probe.** The detector needs the Jacobi-sum collapse γ_k(A)γ_j(A) = μ·Hecke·ray. Here k is the column Gauss index and j the row-twist exponent, and the collapse forces j + k ≡ 3 (mod 6). The reflection root number then has Gauss index 2j + k ≡ j + 3. It is a ray function only if k ∈ {0, 3}, and no squarefree-supported theta column has such k. So with Patterson's column (k = 2, j = 1), ζ ∝ γ₄ ∝ conj γ₂ is forced. No re-orientation of twist exponents removes ζ. This is F1 derived from first principles. | [R] given the paper's reflection; [S] for "no other theta column" |
| **B7** | **Bootstrap remark.** One might try to bound T through T = 𝒮 + T′ using the paper's global supremum β* (zero-freeness to the right of β*). That is vacuous, because 𝒮 *is* the detected signal. (V′) must be "blind". | [R] as logic |
| **LP** | s₀ = 0 gives 0.874957. The low-side block d = ℓ (M′ = M − 2ℓ ≈ 0.5, where the A₂ dual would even be empty) is **not binding**: an infinite saving there leaves θ = 0.874957 (`lp_mix.py`). Any useful saving must hold in the d = 0 block, at M′ = M ≈ 0.83. | [R] as computation |

**Net status.** (V′_β) ⇔ (beyond-the-edge large sieve / A₂ input), exactly as in `wf2/`. Item (i) of the task fails because of B2 and B3. Item (ii) fails because of B1 and B6. Item (iii) fails because of B5 and §4.4. The honest status stays at Liu's 0.874957.

---

## 1. The Gauss-sum content of the summand

### 1.1 Normalisations

The setting follows `TARGET_VLS.md` §1:
- Q = Z^M, L = Q²/Z = Q^{2−1/M}. At M = 0.826 this is L = Q^{0.789} and Z = Q^{1.211}.
- T = Σ_μ β_μ F(μ), with F(μ) = Σ_m π(m) ψ_μ(m) W(q_m/Q).
- π = ρζ, and ‖π‖² ≍ Q^{5/6}.
- Σ_{μ≍L} |β_μ|² ≍ 1 per dyadic block.
- The Cauchy–Schwarz scale is ‖π‖·‖B⁰‖ ≍ Q^{11/12}.

All identities below hold for squarefree, good, primary, coprime m and μ, up to the finite ray and unit factors that the paper carries (ψ₁, κ₀, R, G, ϑ, χ_n(λ)^{−2}, the λ⁴ shift).

### 1.2 Dictionary [R]

**(a) Row scalar.** Lemma Z gives ζ_m = κ₀ψ₁(m)·conj(γ₂(m))⁴. Since γ₂³ = μα and |γ₂| = 1, four equivalent forms follow:

  conj(γ₂)⁴ = μ(m)ᾱ(m)·conj(γ₂(m)) = ᾱ(m)²·γ₂(m)² = χ_m(4)·conj(χ_m(−1))·ᾱ(m)²·γ₁(m)γ₃(m).

The last form uses Lemma Z step 3, γ₁γ₃ = conj(χ(4))χ(−1)γ₂².

**(b) Column.** From paper78 eq. `reflection-infinity-coefficient` together with `reflection-cubic-sextic-conversion`:

  d(λ^{−3}nb³) = 3^{5/2}|b|·g̃₃(n),  g̃₃(n) = χ_n(λ)^{−2}γ₂(n).

So β_μ ∝ γ₂(n)α(μ)|b|/√q_μ. The column carries Gauss index **k = 2**, in the language of §5.

**(c) Kernel.** ψ_μ(m) = (λ⁴μ/m)₂ = χ_m(μ)³, with denominator m. By reciprocity it equals χ_μ(m)³·(ray), with denominator μ.

**(d) Twisted multiplicativity** (Lemma Z step 2):

  γ_j(ab) = (χ_a(b)χ_b(a))^j γ_j(a)γ_j(b) = R(a,b)^j χ_a(b)^{2j} γ_j(a)γ_j(b).

**(e) Numerator absorption.** g̃_j(c,h) := q_c^{−1/2} Σ_x χ_c(x)^j e(hx/c) = χ_c(h)^{−j}γ_j(c) for (h,c) = 1.

**(f) Jacobi factorisation.** At primes, g(χ^a)g(χ^b) = J(χ^a,χ^b)g(χ^{a+b}) for a + b ≢ 0 (mod 6). The sextic and cubic Jacobi sums are Hecke characters (Weil). So γ_aγ_b = (Hecke)·γ_{a+b}, and on squarefree arguments the sign is the Möbius factor, as in γ₂³ = μα. If a + b ≡ 0, then γ_aγ_{−a} = χ(−1).

### 1.3 Lemma P (parity) [R]

**(i)** γ_j(c)·χ_c(x)^k can be written as a single g̃_j(c,h) iff gcd(j,6) | k.
- For the quadratic kernel k = 3 this needs j odd.
- The column carries only γ₂ (j = 2). So the kernel cannot be absorbed on the μ side.
- On the m side, ζ_m contains γ₁ and γ₃, with exactly two absorptions:
  - γ₁(m)χ_m(μ)³ = g̃₆(m, μ³), the A₂ twisting parameter of `wf2/A2_theory.md` §1.3;
  - γ₃(m)χ_m(μ)³ = g̃₃^{(2)}(m, μ), a quadratic Gauss sum with numerator, which is a ray rewrite only.

**(ii) Merge.** For any j, (d) gives

  γ_j(m)γ_j(μ)·χ_m(μ)³ = R^{−j}·γ_j(mμ)·χ_m(μ)^{3−2j}.

The exponent 3 − 2j is odd. So **no merging of an m-Gauss sum with a μ-Gauss sum into a Gauss sum of the product removes the kernel.** It always leaves χ_m(μ)^{±1} or χ_m(μ)³.

*Proof.* (i) follows from (e): one needs h with χ(h)^{−j} = χ(x)^k for all characters of order 6, i.e. −je ≡ k (mod 6) solvable. (ii) is (d) rearranged. ∎

**Example: the best version of the task's item (ii).** Use ζ ∝ ᾱ²γ₂(m)². Then

  β_μ ζ_m ψ_μ(m) ∝ γ₂(μ)γ₂(m)²χ_m(μ)³ = conj(R)²·γ₂(mμ)·γ₂(m)·conj(χ_m(μ)).

This is a "Type II sum of cubic Gauss sums of the product" with a residual **sextic** kernel. It has no cubic-large-sieve structure and no Voronoi structure in either variable.

### 1.4 Consequence: every one-variable rewrite [R]

By Lemma P, the complete list of rewrites of the summand in which one variable carries an automorphic sequence (a theta or Kubota coefficient) together with a character kernel is:
1. **Unreflected.** ρ(m)·[γ₂(A)ᾱ(A)χ_A(m)], with A ≍ Z: (sextic theta) × (cubic theta) and sextic kernel. This is the paper's starting point.
2. **Reflected.** π(m)·[β_μ χ_m(μ)³], with μ ≍ L. The m-sequence π is not a theta sequence.
3. **A₂.** ρ(m)ᾱ²(m)·g̃₆(m, μ³)·β_μ. The m-sequence is the sextic A₂ residue at twisting index μ³, with conductor N(μ)³ [H] (`wf2/A2_theory.md` §2 and its referee item 1).

Forms 1 and 2 are exchanged by the μ ↔ A Voronoi (§2). Form 3 is the only one in which the long variable is automorphic and the kernel has been eliminated.

### 1.5 Which Gauss sums in m are "automorphic in m" [S]

- γ₂(m) is the m-th Fourier coefficient of Patterson's cubic theta: τ₃(m) ∝ γ₂(m) at squarefree m (paper78 eq. `reflection-infinity-coefficient`; Bröker–Hoffstein §1 "τ₃(d) = G₁(1,d)").
- γ₁(m) is **not** a coefficient of the sextic theta. τ₆(m) = q^{1/12}ρ(m) is undetermined at squarefree m (Kazhdan–Patterson non-uniqueness).
- So the only Gauss factor of ζ_m that can be traded for an automorphic expansion in an auxiliary variable is γ₂(m). This is used in §4.

---

## 2. Item (i), Type I: Dunn–Radziwiłł Lemma 5.2 is the paper's reflection

### 2.1 The transform [R]

**DR Lemma 5.2** (dr.txt, l. 1814). Let ψ be periodic mod r, with ψ̂ supported on residues coprime to r. Then the ψ-twist F₁(w;ψ)_S of the cubic theta function equals (3N(r))^{−1} times a combination of the theta functions at the three cusps. Those are evaluated at the inverted point w′ (scale |r|^{−2}) and twisted by

  ψ̃(u) = Σ_x ψ(x) S_{1/3}(x,u;r),   S_{1/3}(x,u;r) = (λ/r)₃ Σ_{(λ⁴d)(λ⁴a) ≡ 1} (a/r)₃ ě((xd + ua)/r).

**Lemma R [R].** Let r be squarefree and good, and ψ = χ_r^j. Then

  ψ̃(u) ∝ χ_r(λ⁸)^j · τ(χ_r^j) · τ(χ_r^{j+2}) · χ_r(u)^{−(j+2)},

with the two degenerate cases:
- j = 4: the a-sum is a Ramanujan sum;
- j = 0: ψ̂(d) = −1 for d ≠ 0.

*Proof.*
1. ψ̂(d) = χ_r(d)^{−j}τ(χ_r^j).
2. Substitute d = (λ⁸a)^{−1}.
3. The cubic symbol is (a/r)₃ = χ_r(a)². The a-sum is therefore τ(χ_r^{j+2})χ_r(u)^{−(j+2)}. ∎

**Check N1** (`check_dr_transform.py`). The identity holds for all j = 0,…,5 on 12 split primes, with maximum relative error 1.8·10⁻¹⁵. Run in the paper's own e(·) convention.

### 2.2 Consequences [R]

**(a) The paper's reflection is recovered term by term.**
- The row factor τ_jτ_{j+2} is ω_{p,j} = τ⁻_{p,j}τ⁺_{p,j+2}(…) in eq. `reflection-row-phase`.
- The dual column factor χ^{−(j+2)} is B_p(x) = χ_p(x)^{−j_p−2}.
- The two degenerate cases reproduce B_p for j_p = 4 (Ramanujan) and for active j_p = 0.

For the base probe the twist of the A-variable is χ_A(m) = R·χ_m(A), so j = 1. The lemma then gives:
- the row factor γ₁(m)γ₃(m), which is **Lemma Z**, now derived from DR's formula rather than the paper's local Fourier computation;
- the column twist χ^{−3}, the quadratic character.

**(b) The map j ↦ −(j+2) is an involution on twist exponents:** {1,3}, {0,4}, {2}, {5}.
- The reflection swaps the sextic twist on A ≍ Z (j = 1) with the quadratic twist on μ ≍ L (j = 3).
- The root numbers are inverse: γ₁γ₃·γ₃γ₅ = χ(−1)² = 1.

**(c) Type I sums are useless here.** The task's "Type I sums of cubic Gauss sums twisted by quadratic characters mod m" are exactly Σ_μ β_μ χ_m(μ)³W(q_μ/L), i.e. j = 3.
- DR Voronoi (Lemma 5.2 → Prop. 5.2-type functional equation, uniform in m) maps them to the A-form: length N(m)²/L = Z, sextic twist, root number γ₃γ₅ = conj(ζ-type factor). That is **un-reflection**.
- A Type I piece with a fixed divisor d | μ has, by (d), a twist by (d/μ′)₃^{±}χ_m(μ′)³ of modulus dm. Its length L/N(d) is below √(conductor) = N(dm) ≥ Q. So Voronoi lengthens it to N(dm)²N(d)/L ≥ Z. [R- given DR Prop. 5.2-type uniformity]
- Since the μ-sum is always shorter than the square root of its conductor (L < Q), no Type I sum in μ can be evaluated or shortened.

---

## 3. Item (i), Type II, and the diagonal obstruction

### 3.1 Type II on μ [R]

Take squarefree μ = μ₁μ₂ with μᵢ ≍ Lᵢ and L₁L₂ = L. Twisted multiplicativity and cubic reciprocity give

  β_{μ₁μ₂} = β_{μ₁}β_{μ₂}·conj((μ₁/μ₂)₃),

a cubic kernel between μ₁ and μ₂. (A Heath-Brown-identity-type weighting of the factorisations is routine and costs Z^ε.)

For fixed m, Heath-Brown's cubic large sieve applied to x_{μ₁} = β_{μ₁}ψ_{μ₁}(m) and y_{μ₂} = β_{μ₂}ψ_{μ₂}(m) gives

  |B⁰_m| ≪ (L₁ + L₂ + (L₁L₂)^{2/3})^{1/2}·Z^ε.

The right side is at least L^{1/3}Z^ε for every split, with equality at the balanced split L₁ = L₂ = L^{1/2}. That is a pointwise bound, against true size 1 and trivial size L^{1/2}. Summing against π gives |T| ≪ Q^{11/12}L^{1/3}, worse than Cauchy–Schwarz.

**DR's corrected cubic sieve (their (1.13), Props. 9.1–9.2).**
- It improves only the (AB)^{2/3} term of the cubic sieve, by subtracting the cubic-theta main term. Its A and B terms, the "number of moduli" and the "length", stay.
- In the Type II form above, DR's absorption (b/a)₃γ₂(b) = g̃₃(b, ·) simply re-merges μ₁μ₂ into μ, by Lemma P(ii).
- Our obstruction is the **length term in m** (§3.2), which no μ-side estimate touches.

### 3.2 Lemma D (diagonal obstruction) [R]

Let B⁰ = (Σ_μ β_μψ_μ(m))_m. Every argument of the form

> "rewrite the μ-side in any way, then apply Cauchy–Schwarz, or Hölder with any exponents 1/p + 1/p′ = 1, in m"

gives a bound at least ‖π‖_p·‖B⁰‖_{p′}. Suppose π is flat (|π(m)| ≍ Q^{−1/12}, which is the numerical unitarity of τ [N]). Then the power-mean inequality gives:
- ‖π‖_p‖B⁰‖_{p′} ≥ Q^{11/12}·(E|B⁰|²)^{1/2} when p′ ≥ 2;
- ‖π‖_p‖B⁰‖_{p′} ≥ Q^{−1/12}‖B⁰‖₁ ≈ Q^{11/12}·E|B⁰| when p′ ≤ 2.

More strongly, sup{|⟨π,B⁰⟩| : |π(m)| ≤ Q^{−1/12}} = Q^{−1/12}‖B⁰‖₁. So **no bound that uses only the size of π and the structure of β can beat Q^{11/12} by a power**, unless ‖B⁰‖₁ = o(Q^{1−η}).

*Proof.* The Cauchy–Schwarz/Hölder output is a norm of the fixed vector B⁰, however B⁰ is computed. Rewriting the μ-side can only overestimate that norm. The supremum is attained at π = Q^{−1/12}·sgn(conj B⁰). ∎

**Check N2** (`check_diag.py`). Setup:
- true weights β_μ = γ₂(μ)α(μ)N(μ)^{−1/2}·bump on squarefree primary μ ∈ (L, 2L];
- rows: primary m ∈ (Q, 2Q] with (N m, 6) = 1.

The statistic is R = ‖B⁰‖²_w/(‖w‖₁·‖β‖²); the diagonal alone gives R ≈ 1.

| Q | L | true β | conj-orientation β | random phases (3 draws) | coherent \|β\| |
|---|---|---|---|---|---|
| 2·10⁴ | Q^{0.60} = 381 | 0.963 | 0.926 | 0.919–0.958 | 1.16 |
| 2·10⁴ | Q^{0.79} = 2499 | 0.960 | 0.939 | 0.922–0.938 | 2.06 |
| 2·10⁴ | Q^{0.95} = 12189 | 0.960 | 0.949 | 0.933–0.958 | 6.90 |
| 6·10⁴ | Q^{0.79} = 5953 | 0.922 | 0.961 | 0.932–0.951 | 2.25 |

At Q = 2·10⁴ and 6·10⁴ (ℓ = 0.79) the shape of B⁰ is also measured:
- E|B⁰|/(E|B⁰|²)^{1/2} = 0.884 and 0.885 (Gaussian value 0.886);
- E|B⁰|⁴/(E|B⁰|²)² = 2.06 and 2.04 (Gaussian value 2).

**Reading [N].**
- The cubic-Gauss-sum structure of β is invisible in the m-norms of B⁰. They are indistinguishable from random phases, and the statistic does detect structure (the coherent control).
- The paper's completed moment ‖B⁰‖² ≪ Q^{1+ε} is sharp. The bound in Lemma D is ≈ 0.89·Q^{11/12}.
- So the structure of β has been used up by the reflection Z → L.

**Corollary [R + N].** A proof of (V′_β) must use structural information about π beyond its size. The only such information known is the automorphy of ρ together with the explicit ζ, i.e. the A₂ residue (§1.4, form 3), or the Voronoi of ρ in the unreflected form (§4.4(a)).

### 3.3 Type I/II in m [R-]/[S]

A combinatorial decomposition of the m-sum needs multiplicativity in m. There are two candidates.
- **ρ.** ρ(m₁m₂) has no known (twisted) factorisation at squarefree indices. This is Kazhdan–Patterson non-uniqueness of the Whittaker model for n = 6 [S]. The Bröker–Hoffstein data show "patternless" prime values [S], and `wf2/A2_numerics.md` finds no Hecke-like structure [N]. So there are no Type II sums in m.
- **ζ_m.** In the form μ(m)ᾱ(m)conj(γ₂(m)), the factor μ(m) invites a Vaughan identity. Its cancellation is Möbius cancellation twisted by finite-order Hecke data: exactly the information the method detects (`V_gaps.md` §1.6; §6 below). Using it is circular.
- In the form ᾱ²γ₂² the Möbius factor disappears, so it is not intrinsic to π. It reappears whenever ζ is split into a cubic Gauss sum times a Hecke character.

---

## 4. Item (ii)/(iii): opening ζ_m, and the other arrangements

### 4.1 Lemma E: exact expansion of γ₂(m) in cubic characters [R- given Patterson/Kubota]

Let D₃(w; m) = Σ_{t primary, good} g̃₃(t, m)N(t)^{−w} be the cubic Kubota series with numerator m and modulus t. Then:
- it has its only pole in Re w > 1/2 at w = 5/6, with residue ∝ τ₃(m)N(m)^{−1/6} ∝ γ₂(m)N(m)^{−1/6} (Patterson; DR Prop. 5.2 for the twisted version);
- its (vector-valued, over cusps) functional equation w ↔ 1 − w has factor N(m)^{1/2−w}.

Shifting contours past the pole and applying the functional equation gives, for U = N(m)^{1/2},

  γ₂(m) = c₃·N(m)^{1/6}U^{−5/6}·[ Σ_{t≍U} γ₂(t)·conj((m/t)₃)·W(N t/U) − (the same at the dual cusps, t′ ≍ N(m)/U) ] (up to ray factors),

since g̃₃(t, m) = conj((m/t)₃)γ₂(t) for (t, m) = 1. The coefficient vector has ℓ²-mass ≍ 1: N(m)^{1/3}U^{−5/3}·U = 1.

**In words.** γ₂(m) is an exact combination of ≍ N(m)^{1/2} cubic characters of m, of modulus ≍ N(m)^{1/2}. This is the cubic analogue of the Y-expansion of ρ (`D_derivation.md`), with the threshold P = U²/N(m) = 1.

### 4.2 The cost: conductor inflation [R-]

Insert Lemma E into the row sequence.

**Optimistic single opening.** Pretend ζ_m ∝ γ₂(m)·Hecke without the Möbius factor; this is false, see §1.2(a). Then

  π(m)ψ_μ(m) = ρ(m)·Hecke(m)·Σ_t a(t)·conj((m/t)₃)·χ_m(μ)³.

1. The m-sum is a θ₆ ⊗ Hecke twist by a character of modulus r = tμ, of norm R := L·Q^{1/2}. Its twist exponents are 3 at μ (cubic dual) and 4 at t (quadratic dual). Both are incoherent, so there is no main term.
2. Voronoi as in `V_derivation.md` §5 gives a dual of length R²/Q with character kernel ξ_r(n).
3. The large sieve over the moduli r (coefficient mass Σ|β_μ|²|a(t)|² ≍ 1) gives, exactly as in `V_derivation.md` §0(d) with L replaced by R,

  |T| ≪ Q^{11/12}·[R/Q + R²/Q² + R²Q^{−5/3}]^{1/2}·Z^ε.

This saves only if R < Q and R² < Q^{5/3}, i.e. L < Q^{1/3}, i.e. **M < 3/5**. Even the first term alone needs M < 2/3.

**LP.** The `zexp` mode of `lp_B.py`, with −2s = max(3M′/2 − 1, 3M′ − 2, 10M′/3 − 2), gives **θ = 0.960487**. By distributivity, θ_true = min(0.874957, 0.9605) = Liu.

**Honest version.** ζ_m ∝ ᾱ²γ₂(m)² needs **two** openings. That gives moduli t₁t₂μ of norm R = L·Q > Q, and the bound is worse than Cauchy–Schwarz for every M. The alternatives fail too:
- opening once in the form μᾱ·conj γ₂ leaves ρ(m)μ(m), which has no Voronoi formula;
- γ₁ cannot be opened at all (§1.5).

**Conclusion.** Opening the Gauss sum trades the "rough" ζ for conductor Q^{1/2} per factor. "Fewer columns than rows" (R < Q) is destroyed in every admissible version. [R- given Lemma E]

### 4.3 Comparison with the A₂ route

The same counting in A₂ language: the effective column conductor is N(μ)^{3/2} (conductor L³), and the saving condition L^{3/2} < Q is M < 3/4 (`A2_theory.md` §2.6). The referee's stronger coupling obstruction moves this to k < 1.53. The hierarchy of "effective column modulus" is:

| model | effective modulus R | saving needs |
|---|---|---|
| ζ-free (V) | L | M < 6/7 (from L²/Q^{5/3} < 1) |
| A₂ residue | L^{3/2} (conductor L³), and coupled | M < 3/4, coupled: none |
| ζ opened once (illegitimate) | LQ^{1/2} | M < 3/5 |
| ζ opened twice (legitimate) | LQ | never |

### 4.4 Other arrangements, item (iii)

**(a) DR-type dispersion in the unreflected form [R-].** In Σ_A c_Aχ_A(m) with c_A ∝ γ₂(A), absorb the cubic part of the kernel, γ₂(A)χ_A(m) = g̃₃(A, m)·(m/A)₂. Then expand |·|² over m and apply Poisson modulo A₁A₂, so that the dual length is Z²/Q > Z. The complete sums produce τ(χ_{A₁})τ(conj χ_{A₂}), and γ₂γ₁ = μαG collapses the Gauss sums:

  ‖B‖² = Q‖c‖²·(1 + o(1)) + (Q/Z)·Σ_{0<h≲Z²/Q} |Σ_A μ(A)ν(A)conj(χ_A(h))·Z^{−1/2}W|² + … .

The off-diagonal is the Möbius-twisted family (the signal family), and the number of moduli Z²/Q exceeds the length Z. So this is neither "fewer columns" nor non-circular.

**(b) Poisson in μ, Heath-Brown's quadratic-sieve duality [R].** Σ_μ|F(μ)|² = L·Σ_{mm′=□}π(m)conj(π(m′)) + (L/Q)·Σ_{0<h≤Z} F₁(h)F₂(h), with F₁ = Σ π ε ψ_h and F₂ = Σ conj(π) ε ψ_h. Cauchy–Schwarz in h returns exactly (L + Q)‖π‖². Beating it needs cancellation in the *bilinear* (non-positive) h-sum, which Poisson in h turns back into the original problem. The duality L ↔ Z is an identity of the sharp large sieve and gives nothing for free.

**(c) The A₂ "m-Kubota" form (μ ≍ L, s ≍ L³/Q) [R given KP].** Here there are fewer columns than rows (L < L³/Q). But the dual coefficient is ρ(sμ³)γ₁(s) = N(μ)^{−1/6}G₄(1,μ)conj(χ_μ(s))⁴ρ(sμ)γ₁(s), and ρ(sμ) couples s and μ through a theta coefficient at the product. It is not a character kernel (`A2_theory.md` §2.3, referee item 3). No large sieve applies.

**(d) Multivariable convexity [H].** In the (t, A, m, s) multiple Dirichlet series of `A2_theory.md` §4, Phragmén–Lindelöf in tubes interpolates *affinely in the exponents*. A linear functional on a polytope is minimised at a vertex, so the convexity bound for T is the best trivial bound among the images, i.e. vertices of the hull of the images of the absolute-convergence region. Every image identified so far has a trivial bound worse than Q^{11/12}. Cauchy–Schwarz is a positivity bound and is not of this type. A complete enumeration of the 7 objects is still [G] (as in `wf2`), but it can only reproduce some image's trivial bound.

---

## 5. Rigidity of the probe: why ζ cannot be removed (B6)

Consider a probe of the paper's shape:
- sextic Kubota rows Σ_s g̃₆(s,−m)…;
- a column Σ_A c_A χ_A(m)^j with Gauss-sum column weights c_A ∝ γ_k(A)·(Hecke);
- the zero detector obtained by Poisson in m (paper78 §6).

**Lemma H [R, given the paper's reflection formula and §1.2(f)].**

**1. Decoupling.** The complete sum modulo sA is

  Σ_{x mod sA} conj(χ_s(−x))χ_A(x)^j e(kx/sA) = R^j·χ_s(A)^{j−1}·χ_s(k)χ_A(k)^{−j}·τ(conj χ_s)τ(χ_A^j).

The s- and A-sums decouple (high side = product of a character sum in s and a sum in A) **iff j ≡ 1**. For the conjugate orientation of the whole probe, j ≡ 5.

**2. Collapse.** The A-weight after Poisson is c_Aτ(χ_A^j)/√q_A ∝ γ_kγ_j. By §1.2(f), γ_kγ_j is unimodular and of the form μ·Hecke·ray (needed for Σμη = 1/L) iff k + j ≡ 3. So the pair is {1, 2} or {4, 5}, and the detector forces (j, k) = (1, 2) up to global conjugation.

**3. Root number.** Lemma R generalises: a column of Gauss index k has multiplier χ^k (k = 2 is Kubota's cubic symbol (c/d)₃ = χ²), and its twist transform produces the root number γ_jγ_{j+k} = (Hecke)·γ_{2j+k}. For the Patterson column this is γ_jγ_{j+2}. Under k + j ≡ 3 the Gauss index is 2j + k ≡ j + 3. The root number is ray (|·| = 1 and constant on classes) iff 2j + k ∈ {0, 3} (mod 6), i.e. iff j ∈ {0, 3}, i.e. iff k ∈ {3, 0}. For k = 2 fixed, without the detector constraint, ray means 2j + 2 ≡ 0, i.e. j ∈ {2, 5}; these are the self-dual twists of §2.2(b).

**4. No such column.** k = 3 means quadratic Gauss sums, i.e. the n = 2 theta, whose coefficients live on squares. k = 0 means c_A = μ(A), which is not automorphic. Neither is a column of the required kind [S: for n ≠ 3 no theta function has Gauss-sum coefficients at squarefree indices; n = 6 is the non-unique KP case].

**5. Conclusion.** The admissible probe has j = 1, k = 2, and root number γ₁γ₃ = (Hecke)·γ₄ ∝ μᾱ·conj γ₂ (Lemma Z).

So **the same Jacobi-sum identity γ₁γ₂ = μαG that manufactures the detector forces the cubic Gauss sum into ζ.** Re-orienting the twist or the column does not help. Twists with j ∈ {2, 5} (self-dual under j ↦ −(j+2), ray root number) kill the detector. This makes the referee's "structural point" in `V_gaps.md` §1.6 precise.

---

## 6. The bootstrap remark (B7) [R as logic]

The paper already works with the global supremum β* of real parts of zeros (§`sec:detector`, "assume β* > 51/100"). One might try the following:
- write T = 𝒮 + T′ (`V_gaps.md` §1.6);
- bound 𝒮 ≪ Q^{5/6}Z^{β*−5/6+ε} via zero-freeness in Re s > β*, and bound T′ similarly;
- at M = 0.826 and β* = 0.875 this gives 𝒮 ≪ Q^{11/12}Z^{−0.027}.

**This is vacuous.** 𝒮 is the principal high-side signal: the signal coefficient is K̂(5/6), the s = 5/6 residue of the theta part (REPORT §0; `V_gaps.md` referee item 1). The contradiction "signal > low-side bound" cannot be reached if the low-side bound contains an upper bound for the signal that is valid whenever β₀ ≤ β*. In the identity H = I the term 𝒮 cancels exactly.

Hence (V′) must be **blind**: valid without zero information, like Cauchy–Schwarz and like (V′_LS). This also rules out exploiting the Möbius factor in ζ_m = μᾱ·conj γ₂ (§3.3).

---

## 7. LP checks

All runs use the model of `wf/lp_const_saving.py`, with the P²/Y piece dropped (LoA6) and the theta piece g = b/6 − 2s(M′).

| run | θ | geometry (l_x, l_y, ℓ; M) |
|---|---|---|
| s₀ = 0 (this note's result) | **0.874957** | (0.3549, 0.4783, 0.1668; 0.8332) |
| s₀ = 0.002 (reference) | 0.873332 | (0.3558, 0.4789, 0.1653; 0.8347) |
| Cauchy–Schwarz ('none') | 0.874957 | same as s₀ = 0 |
| ζ opened once (illegitimate, §4.2) | 0.960487 | (0.2956, 0.3485, 0.0004; 0.6441) |
| C–S at d = 0, **zero** theta piece at d = ℓ (`lp_mix.py`) | **0.874957** | unchanged |

The last row matters for the orchestrator. The low side has two blocks, d = 0 (M′ = M ≈ 0.833) and d = ℓ (M′ = M − 2ℓ ≈ 0.50).
- In the d = ℓ block the marked reflected length is L = Z^{2M′−1} ≈ Z⁰. There even the A₂ functional equation would make the dual empty (L³ < Q′), so a huge saving is plausible there [H].
- But the block is **not binding**: removing its theta piece entirely leaves θ = 0.874957.
- So a saving must hold in the unmarked block, at **M′ = M ∈ [0.80, 0.86]**.
- The per-block mixing that the `wf2` LP did not try therefore cannot rescue the A₂ dual route.

---

## 8. What is left

1. **The reduced problem is unchanged by Attack B.** (V′_β) needs a blind bound for ⟨π, B⁰⟩ that uses π's automorphy. By §1.4 and §5 that automorphy is the sextic A₂ residue. By Lemma D, the β-side cannot contribute beyond ‖B⁰‖ ≍ Q^{1/2}. So (V′_β) is as hard as the beyond-the-edge statement:
   - family μ ≍ L = Q^{0.79};
   - conductor ~ L³ (or with the k = 3 caveat [H]);
   - coefficients coupled through ρ(sμ³);
   - a power saving over (L + Q)‖π‖² needed in the m-length term.
2. **Positive facts recorded here** (usable by any later attempt):
   - the (j, j+2) rule (Lemma R) and its independent derivation of Lemma Z;
   - the parity lemma (Lemma P);
   - the exact expansion of γ₂(m) in cubic characters (Lemma E), with its precise conductor cost;
   - sharpness of the completed-row moment for the true column, with Gaussian statistics of B⁰ [N];
   - the non-binding d = ℓ block.
3. **Gaps that remain from `wf2`:**
   - the A₂ functional equation for n = 6 over Q(ω), including the vector-valued equation over indices (1, μ³d^f) at p | μ [G];
   - the conductor exponent k = 3 [H];
   - no large sieve for coupled theta coefficients [G].

**Status: failed for Attack B; (V′_β) is reduced to (V′_LS)/A₂, with the reasons above. θ = 0.874957 (Liu).**

---

## Referee report

*Independent skeptical referee, 9 Oct 2026. Referee scripts are in `scratchpad/refB/`: `lp_mix2.py`, `dict_check.py` and `diag_run.py`.*

### Verdict

**The note's main conclusion stands.** Attack B gives **no power saving** for (V′_β) at any M ∈ [0.80, 0.86], and the LP value stays at Liu's 0.874957. I found no hidden saving and no error that would turn the result around. The note claims no positive result, so there is nothing positive to refute (refuted = false).

Several side claims are stated too strongly or carry the wrong label. They are listed under "Objections" below. None of them changes the verdict.

### Re-runs (all reproduce)

| check | result |
|---|---|
| `lp_B.py` | const s₀ = 0: 0.874957; s₀ = 0.002: 0.873332; zexp: 0.960487; none: 0.874957 ✓ |
| `lp_mix.py cs zero` | 0.874957 for both ✓ (the d = ℓ block is not binding) |
| **new** `refB/lp_mix2.py`: per-block route mixing for the ζ-opening route | (zexp at d=0, zero at d=ℓ) = **0.960452** = (zexp, zexp); (cs, zero) = 0.874957 |
| `check_dr_transform.py` | maximum error 1.84·10⁻¹⁵ ✓, on split primes only |
| `check_diag.py`, Q = 2·10⁴, ℓ = 0.79, new seed | R_true = 0.960, R_conj = 0.939, R_coherent = 2.064, random phases 0.900 / 0.940 / 0.906; E\|B⁰\|/rms = 0.884, kurtosis 2.062 ✓ |
| **new** `refB/dict_check.py`, 768 split primes, N ≤ 6000 | γ₂³ = μα (3·10⁻¹⁵); γ₁γ₂ = μα·χ̄(4)γ₃ (2·10⁻¹⁵); γ̄₂⁴ = ᾱ²·Ḡ·γ₁ (5·10⁻¹⁵); γ₃ constant on classes mod 72 (spread 1.5·10⁻¹⁵) |

**Notes on the re-runs.**
- **LP mixing.** "θ_true = min(θ_CS, θ_zexp) by distributivity" needed a check, because a per-block minimum of two routes is not convex. The new run settles it. The LP is monotone in each block's theta piece, so every mix of (CS or zexp) at d = 0 with any route at d = ℓ is at least min(0.874957, 0.960452) = 0.874957. Minor: the note's 0.960487 sits on the S-grid boundary (S = 0.0005); the optimum as S → 0 is 0.960452. This changes nothing.
- **Dictionary.** The new run confirms §1.2. It confirms that the A₂ form ζ_m ∝ (ray)·ᾱ²(m)·γ₁(m) holds exactly (at primes; composites carry ray signs only). It also confirms that quadratic Gauss sums γ₃ are ray functions, which B6 step 3 uses.

### Bookkeeping checked (correct)

- **B5 thresholds.** R = LQ^{1/2} < Q ⇔ M < 2/3, and R² < Q^{5/3} ⇔ M < 3/5. Two openings give R = LQ, which never saves.
- **The §4.3 table.** (V): L² < Q^{5/3} ⇔ M < 6/7. A₂: L³/Q < Q ⇔ M < 3/4.
- **The twist bookkeeping in §4.2** agrees with the sextic rule of `V_derivation.md` §0(b), χ^j ↦ χ^{1−j}:
  - j = 4 at t (the conj cubic symbol) gives a quadratic dual;
  - j = 3 at μ gives a cubic dual;
  - neither twist is coherent (j = 1).
- **§4.4(b).** Since LZ = Q², the duality returns exactly (L + Q)‖π‖².
- **Lemma H step 1, re-derived by CRT.** The cross factor is χ_s(−1)^{-1}R^jχ_s(A)^{j−1}, so the s- and A-sums decouple iff j ≡ 1. ✓
- **Lemma P (i) and (ii).** ✓
- **§2.2(b).** The root-number product is γ₁γ₃·γ₃γ₅ = χ(−1)² = 1. ✓
- **Column index k = 2.** paper78 eq. `reflection-infinity-coefficient` and eq. `reflection-cubic-sextic-conversion` give d = 3^{5/2}|b|·χ_n(λ)^{-2}γ₂(n). ✓

### Attempts to rescue a saving inside Attack B (all fail, consistent with the note)

1. **Type II in μ = μ₁μ₂ with Cauchy–Schwarz over μ₂ and (m, μ₁) inside.** The kernel χ_m(μ₂)³κ(μ₁, μ₂) has moduli (m, μ₁), at least Q of them. This is worse than Cauchy–Schwarz.
2. **Absorb the kernel into the quadratic Gauss sum.** Write γ₃(m)χ_m(μ)³ = g̃(χ_m³, μ) and open it additively. The μ-sum becomes an additive twist Σ_μ β_μ e(μx/m) of the cubic theta. DR Voronoi sends it to length N(m)²/L = Z, which is un-reflection (B2).
3. **Open ρ through its Y-expansion** g̃(s, −m) and use χ_s(m)^{-1}γ₁(m)χ_m(μ)³ = g̃₆(m, sμ³)·(ray). The m-sum is then a Kubota series in m whose residue is ρ(sμ³). That is exactly the coupled dual coefficient of `wf2`, so no character kernel appears.

### Objections (labels to change, claims to weaken)

1. **"(V′_β) ⇔ (V′_LS)/A₂"** (§0 "Net status", §8.1, and "reduces to" in the summary) **is an overclaim.** Only (V′_LS) ⇒ (V′_β) is proved, by Cauchy–Schwarz in μ. What the note shows is that the tools of Attack B extract nothing beyond what (V′_LS) would give. That should be labelled **[H]**, not stated as an equivalence. (V′_β) remains formally weaker than (V′_LS).
2. **Lemma D is a tautology.** It says Cauchy–Schwarz/Hölder is sharp against a worst-case π of the given size. It is correct, and as a statement about arguments that see only |π| it is fine.
   - The gloss "the Gauss-sum structure of β is used up by the reflection" is **[H]**. Lemma D says nothing about arguments that use β *jointly* with π's structure. For example, in the A₂ form the μ-sum Σ_μ β_μ g̃₆(m, μ³) is itself a structured object.
   - The Hölder/sup version also needs π to be flat. That is [N] only, as the note says. If π were very spiky, Hölder with ‖π‖₁·‖B⁰‖_∞ (with ‖B⁰‖_∞ ≪ L^{1/3}) could in principle win. The unitarity numerics make this implausible.
3. **The circularity rationale in B4 and B7 conflates two Möbius functions.**
   - The μ(m) in ζ_m = μᾱγ̄₂ sits on the **row** variable. The detector's Möbius is μ(A) on the **column** (`V_gaps.md` §1.6).
   - Cancellation in Σ_m μ(m)·[ρ(m)γ̄₂(m)ᾱ(m)ψ_μ(m)] is not the same as zero-free information for the detected Hecke L-functions, because the bracket is not multiplicative. It is "circular" only in the heuristic sense of V_gaps' "structural point [H]".
   - B4's conclusion stands for the other reason the note gives: ρ has no twisted multiplicativity, so Vaughan/Heath-Brown identities give no usable Type I/II sums. Also, the μ(m) factor depends on the chosen representation; it is absent in the ᾱ²γ₂² and ray·ᾱ²γ₁ forms.
   - Relabel the circularity sentence **[H]**. The last sentence of §6 ("This also rules out exploiting the Möbius factor in ζ_m") does not follow from B7. B7 itself, that bounding 𝒮 by zero information is vacuous, is fine.
4. **Lemma E's identity is incomplete as written.**
   - Shifting the contour to Re w < 0 also crosses w = 1/6. That is the image of the 5/6 pole under w ↔ 1 − w, and it carries a dual-cusp theta coefficient. At U = N(m)^{1/2} its residue has the same size, N(m)^{1/4}, as the main residue.
   - The λ-power index shift must also be tracked: Patterson's ∞-support is λ^{-3}nb³.
   - Both points change only the shape and constants of the expansion. One can take U = N(m)^{1/2+δ} to make the 1/6 term lower-order, at a slightly larger conductor cost. The conclusion of B5 (a column-modulus cost of about Q^{1/2} per opening) is unaffected.
   - Keep the label [R-] and add these two items to the gap list.
5. **B6 (rigidity) depends on the probe template.** The template is: sextic Kubota rows, a Gauss-index-k column, the twist χ_A^j, and a decoupled Poisson detector that needs a μ(A) collapse.
   - "For every admissible probe orientation" should read "for every probe of this template".
   - Probes with coupled (non-decoupled) high sides, or with columns that are not Gauss-index columns (e.g. Eisenstein or Kubota-series columns), are not covered.
   - So B6 is not a proof that no probe can remove ζ.
6. **B2's "independent derivation of Lemma Z" is slightly overstated.**
   - N1 checks only the local prime-level identity, and only at split primes; inert primes are untested.
   - The global assembly (the σ_p and ε_p factors, and the CRT cross terms) is taken from `V_gaps.md`. The script does not re-check it.
   - The [R] label for the algebra is fine.
7. **μ₁₂ ⊄ Q(ω).**
   - The cubic inputs need only μ₆ ⊂ K, which holds. These are Patterson for Lemma E, DR Lemma 5.2, and Heath-Brown's cubic large sieve.
   - The only results relying on n = 6 metaplectic theory are the Kazhdan–Patterson citations [S] and the inherited A₂ functional equation [G].
   - Any future A₂ functional equation here must go through Kubota-type arguments, as in `wf/FE_LoA6_proof.md`, not through BBF/BBFH, which assume μ_{2n} ⊂ F. The note should say this next to its [G] item.
8. **N2 numerics are small.**
   - The test has 712 columns and 6036 rows at Q = 2·10⁴.
   - The random-phase spread across seeds is 0.90–0.96. So "true β is indistinguishable from random phases" is the right reading, and it is all the test can show.
   - Only squarefree μ are used, and the λ⁴ shift, ϑ, and B_p at non-squarefree rows are dropped. All of this is harmless for an [N] illustration of ‖B⁰‖² ≍ Q‖β‖².

### Strongest surviving claim

The Gauss-sum structure of the weights does not produce a saving with any of the routes (i)–(iii):
- Type I/II in μ;
- DR Voronoi, which is the paper's reflection, Lemma R;
- Hasse–Davenport and twisted-multiplicativity rewrites, which only reach the A₂ form ρ·ray·ᾱ²·γ₁;
- opening γ₂(m), which costs about Q^{1/2} of column modulus per opening.

So s₀ = 0 and θ = 0.874957. A proof of (V′) must use the joint structure of π = ρζ against the quadratic kernel, beyond its size. The only identified candidate is the A₂ residue, which remains open.
