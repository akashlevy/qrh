# Attack C: moments via multiple Dirichlet series and spectral theory

*9 Oct 2026. Target: (V′), (V′_LS) and (V′_β) of `TARGET_VLS.md`. Builds on:*
- *`REPORT.md` (CORRECTION banner);*
- *`wf/V_gaps.md` (Lemma Z, §1.6 the S + T′ split);*
- *`wf/FE_LoA6_proof.md` (Thm A);*
- *`wf2/A2_theory.md` and its referee report;*
- *`wf2/A2_numerics.md`;*
- *`IDEA1_NOTES.md` §2;*
- *`LIT_ideas234.md`;*
- *`paper78.tex`: eq. `signal` l.965, `low-separated` l.3470, `probe-poisson-identity` l.3831, Prop. `continuation-criterion` l.401.*

Scripts are in `gram_split/vls/`, with logs next to them:

| script | what it computes |
|---|---|
| `tits_hull.py` | the reflection group of the second-moment series, its Coxeter type, and the Bochner/convexity hull bound (§2.4) |
| `tits_membership.py` | which Lindelöf-on-average points the hull contains (§2.4) |
| `groupoid_moments.py` (+`.log`) | Sawin–Whitehead / Heckenberger groupoid finiteness for every moment series of this attack (§2.3, §4.1) |
| `check_lemmaC0.py` (+`.log`) | the Möbius–Kubota form of the unreflected row (Lemma C0) |
| `lp_signal_circularity.py` (+`.log`) | the LP fed with a signal-vertex bound (§4.3) |

Labels:
- **[R]** rigorous, given the cited results;
- **[R-]** routine but not written out in full;
- **[H]** heuristic;
- **[G]** gap;
- **[N]** numerical.

Citation labels: **[V]** checked in the source; **[S]** secondary (via `LIT_ideas234.md` or `A2_theory.md`).

---

## 0. Verdict

**Status: failed, with the reasons identified. Attack C does not prove (V′), (V′_LS) or (V′_β), and gives no θ below 0.874957.** The three parts fail for three different structural reasons. The reasons are new and fairly robust.

### (i) The second moment as an MDS is of hyperbolic type

Σ_μ|F′(μ)|² is a double theta residue of a 5-variable series (§2.1).
- Its Weyl groupoid is infinite [R as combinatorics].
- In the reduced 3-variable model, the functional equations generate the Coxeter triangle group:
  - (2,6,6), which is hyperbolic, when the twist conductor is N(μ)³;
  - (2,4,4), which is affine, when the conductor is N(μ)²;
  - (2,3,3) = A₃, which is finite, when it is N(μ).

  [R as group theory.]
- With the A₂ conductor exponent k = 3 [H], Bochner plus Phragmén–Lindelöf convexity gives nothing beyond the trivial bound at L = Q^λ, λ ∈ [0.67, 1) [R as computation]:
  - the hull of the images of the convergence region yields only Σ_μ|F′|² ≪ Q^{5/6}·Q·L;
  - this does not even reach the large sieve Q^{5/6}(L + Q);
  - the hull contains no point of Lindelöf-on-average type; on the diagonal slices tested (a = b, Re u = 1.01, 1.5, 2.5) it coincides with the region given by pointwise convexity of each twist.
- Only rigorous inputs [R]: rank-one continuation (Kubota 1969 + Selberg/EGM, `FE_LoA6_proof.md` Thm A) and Poisson in μ. These give even less.
- So no upper bound ≪ Q^{1−η}‖π‖² follows from any MDS argument available. With k = 2 (affine) it would follow in principle (§2.4).

### (ii) Spectral and bilinear tools do not apply

- **Louvel's Kuznetsov.** It is on the wrong cover (cubic), and Kuznetsov needs smooth weights in the modulus.
  - Lemma C0 [R, checked to 1.8·10⁻¹⁵] shows the unreflected row is a Möbius-weighted sum of conjugate sextic Gauss sums.
  - So the modulus weight in T is μ(A)·(θ₆-pairing), not smooth.
  - The relevant Kloosterman sums are degenerate (n = 0) or GL₃^{(6)} ones, for which no Kuznetsov formula exists.
- **Duke–Iwaniec / Dunn bilinear bounds.**
  - They do not apply to the primal form: the kernel is a character, not a coefficient at a product.
  - On the A₂ dual, even ignoring the coupling kernel, they never beat the trivial dual bound for λ > 0.4. A saving needs a factor Q/L² = Q^{−0.58} [R as algebra].
- **Type I/II (Vaughan)** on the Möbius form of T has Type-I main terms equal to the zero-detecting signal (§3.4).

### (iii) The first moment is finite type, but it carries the signal, so evaluating it is circular

- **Finite type.** With its natural column weights (cubic theta coefficients), T lifts to a 4-variable series whose groupoid is finite: 7 objects, 336 bases. With unweighted columns it is infinite [R as combinatorics].
- **But T contains the full zero-detecting signal.**
  - The Mellin transform of T has the polar line z = 5/6, carrying H(u)/L(u + 5/6, χ) [R-] (§4.2, Lemma C1).
  - The coefficient is exactly 1 (S = f), by a Y-scaling argument [R-] (§4.2, Lemma C2).
- **Consequences.**
  - Any evaluation of T by continuation (MDS, Voronoi chains, Type I/II with main terms) reproduces f.
  - Bounding f with an a priori zero-free region θ_k (the paper's 7/8, or Liu) returns θ = θ_k in the LP, exactly. Checked: θ_k = 0.874957 → 0.874957 and θ_k = 0.86 → 0.86, whereas a signal-blind constant s₀ = 0.027 would give 0.8586.
  - Bounds on T′ = T − S carry no zero information at all, because J − f = Y^{−1/3}Ŵ T′ + ⟨A^rest, B⟩ is the probe identity.
- **Hence (V′) must be proved by a signal-blind inequality** (large sieve / mean value type). The route in `TARGET_VLS.md` §3, "S handled together with the high side", is vacuous [R as logic, given Lemma C2].

### Net

Attack C reduces (V′) to a **signal-blind beyond-the-edge mean value**, (V′_LS) or a refinement of Cauchy–Schwarz. For that statement:
- every MDS formulation is hyperbolic (k = 3);
- no spectral large sieve exists for the relevant family (GL₃^{(6)} Whittaker coefficients at the sparse index set (1, μ³)).

The decisive analytic parameter is the effective conductor exponent k of the quadratic twist of π:
- k ≤ 2 would make the second-moment MDS affine or finite;
- the best available heuristic is k = 3 (`A2_theory.md` §2 and referee item 1).

The honest status remains Liu's 0.874957.

---

## 1. Set-up and three structural facts

### 1.1 Normalisations

- τ(m) := q_m^{1/12} ρ(m) is unitary on average: Σ_{q_m≤X}|τ|² ≍ X (`D_derivation.md` §6.3; G3 of `V_gaps.md`).
- On squarefree rows π(m) = ρ(m)ζ_m with |ζ_m| = |κ₀| constant (Lemma Z). Hence ‖π‖²_Ω ≍ Q^{5/6}.
- F′(μ) := Σ_m Ω(q_m/Q) π(m) ψ_μ(m) and T = Σ_μ β_μ F′(μ), with Σ_{q_μ≍L}|β_μ|² ≍ 1 per dyad.
- Write λ := log L/log Q = 2 − 1/M:

| M | 0.80 | 0.826 | 0.86 |
|---|---|---|---|
| λ | 0.750 | 0.789 | 0.837 |

### 1.2 How much (V′_LS) actually asks for [R]

(V′_LS) asks, **on average over μ**, for |F′(μ)|² ≤ Q^{1−λ−η}‖π‖². At M = 0.826 that is Q^{0.21−η}‖π‖².
- This is much weaker than Lindelöf on average, which would give O(‖π‖²) on average.
- It is nevertheless "beyond the large sieve": for general sequences the bound (L+Q)‖π‖² is attained by π ∝ ψ̄_{μ₀}.
- A pointwise version, sup_μ|F′(μ)|² ≤ Q^{1−λ−η}‖π‖², would also suffice. It requires |F′(μ)| ≤ Q^{5/12+(1−λ−η)/2}, a pointwise saving of Q^{(λ+η)/2} ≈ Q^{0.39} over the trivial bound Q^{11/12}. This is a subconvexity-strength statement for a twist of heuristic conductor L³ = Q^{2.37} [H]. Out of reach.

### 1.3 Poisson duality in μ does not reduce the problem [R-]

Let Φ ≥ 0 be supported in [1,2].
- Poisson in μ ∈ O (a rank-2 lattice), for the quadratic character χ = ψ_{m₁}ψ̄_{m₂} of modulus q, N(q) ≍ Q², gives

  Σ_μ Φ(q_μ/L) χ(μ) = c_O (L/N(q)) Σ_{k∈𝔡⁻¹} Φ̌(k√(L/N(q))) g_χ(k),

  where Φ̌ is the planar Fourier transform of x ↦ Φ(|x|²).
- For squarefree coprime m₁, m₂, quadratic reciprocity over Q(ω) gives g_χ(k) = (k/m₁m₂)₂ ε(m₁)ε(m₂)·(ray)·N(q)^{1/2}, where ε is the normalised quadratic Gauss sum, a ray class function.
- Hence

  Σ_μ Φ|F′(μ)|² = (diagonal, k = 0) ≍ L‖π‖²  +  (L/Q) Σ_{k≠0, N(k)≲Q²/L} Φ̌_k·|F″(k)|²   (after Mellin separation of Φ̌),

  with F″(k) = Σ_m π(m)ε(m)(k/m)₂ Ω.
- The dual family has size K = Q²/L > Q, so it is "inside the edge". But the dual weights have **exact mean zero**: by Poisson again, Σ_k Φ̌(kt) = c t^{−2} Σ_ν Φ(N(ν)/t²) = c t^{−2}Φ(0) = 0, because t² = L/N(q) < 1/2.
- So the identity is self-dual and sign-indefinite. (V′_LS) is equivalent to a cancellation, by a factor Q^{−η}, of the dual off-diagonal against the dual diagonal. Poisson takes that cancellation back to the primal off-diagonal. No positivity is gained.

### 1.4 Lemma C0: the unreflected row is a Möbius-twisted conjugate sextic Kubota sum [R]

**Statement.** For squarefree primary A and every m,

  c_A χ_A(m) = μ(A)·q_A^{−1/2} Σ_{y mod A} χ̄_A(y) e(−my/A),   c_A = γ₂(A)ᾱ(A)Ḡ(A),

with the ray factor ν carried along unchanged. Hence

  B_m = Σ_{A sf} μ(A) ν(A) V_G(q_A/Z) q_A^{−1} ğ(−m; A),   ğ(h; A) := Σ_y χ̄_A(y) e(hy/A).

**Proof.**
1. Eq. `signal` gives γ₁γ₂ = μαG. With |G| = 1, this gives γ₂ᾱḠ = μ/γ₁ = μ·γ̄₁.
2. γ̄₁(A) = q_A^{−1/2} Σ_x χ̄_A(x) e(−x/A).
3. Multiply by χ_A(m) and substitute x = my.
4. If p | (m, A), both sides vanish: the left because χ_A(m) = 0, the right because its p-factor is Σ_y χ̄_p(y) = 0. ∎

**Check [N].** `check_lemmaC0.py`: 644 pairs (A, m), with A a product of 1–3 split primes and m including multiples of a prime of A, by brute force over O/A. Maximum error 1.8·10⁻¹⁵.

**What it shows.**
- The Davenport–Hasse alignment built into the probe's columns turns B into a Möbius-weighted Kubota sum, whose numerator m is the free row variable.
- So T = Σ_A μ(A)ν(A)q_A^{−1}V_G Λ(A), with Λ(A) := Σ_m ρ(m) ğ(−m;A) K(q_m/Q) ξ(m). This is a Möbius sum over moduli of degenerate sextic Kloosterman sums S₆(−m, 0; A) = ğ(−m; A), paired with θ₆.
- This form explains §3 and §4: the Möbius weight is where the zeros enter.

---

## 2. C(i): the second moment as a multiple Dirichlet series

### 2.1 The series [R on squarefree coprime terms]

**Five-variable series.** By `A2_theory.md` §1.3, γ₁(m)(μ/m)₂ = g̃(m, μ³) for (m, μ) = 1. Hence π(m)ψ_μ(m) = ρ(m)g̃(m, μ³)ᾱ(m)²·(ray). Define

  Z₅(w₁,z₁,u,z₂,w₂) = Σ_{s₁,m₁,μ,m₂,s₂} c(s₁)g̃(s₁,−m₁) g̃(m₁,μ³)ᾱ²ψ(m₁) · conj[c(s₂)g̃(s₂,−m₂) g̃(m₂,μ³)ᾱ²ψ(m₂)] · N(s₁)^{−w₁}N(m₁)^{−z₁}N(μ)^{−u}N(m₂)^{−z₂}N(s₂)^{−w₂}.

**Residue.** Take the residues at w₁ = w₂ = 2/3 (FE_LoA6 Thm A for each inner Kubota series):

  Z₂(z₁,z₂,u) := Res Res Z₅ = Σ_μ N(μ)^{−u} R(z₁;μ) conj R(z̄₂;μ),   R(z;μ) = Σ_m π(m)ψ_μ(m)N(m)^{−z}.

**Mellin identity** [R, on lines of absolute convergence]:

  Σ_μ Φ(q_μ/L)|F′(μ)|² = (2πi)^{−3} ∫∫∫ Ω̂(z₁) Ω̂(z₂) Φ̂(u) Q^{z₁+z₂} L^u Z₂(z₁,z₂,u) dz₁ dz₂ du.

Non-squarefree p-parts differ from Chinta–Gunnells' stable p-parts (`A2_theory.md` referee item 2). These are [G], irrelevant for exponents.

### 2.2 What continuation and polar structure is provable

All of these are valid over Q(ω) without μ₁₂ ⊂ K.

| object | provable statement | label |
|---|---|---|
| D_c(w;m), rank one | Meromorphic in C, polynomial growth in vertical strips. The only pole in Re w > 1/2 is w = 2/3. Proof: Kubota 1969 (needs only μ₆ ⊂ K, totally imaginary) and Selberg/Elstrodt–Grunewald–Mennicke for Eisenstein series on Γ(N) ⊂ SL₂(O) with the finite-order Kubota character. Eckhardt–Patterson / Bröker–Hoffstein give the same functional equations for n = 6 over Q(ζ₆). | [R] (`FE_LoA6_proof.md` Thm A) |
| R(z;μ), the A₂ theta residue twisted by μ³ | Meromorphic continuation with polynomial growth. Route (a): Chinta–Gunnells' Bochner argument from rank-one equations plus local p-parts at good primes, as Kubota class vectors. Route (b): the 6-fold cover of SL₃(O) splits over a congruence subgroup Γ(N), because the congruence kernel is μ_K (Bass–Milnor–Serre; Moore). So the Borel Eisenstein series, and its residue E_P(θ₆ ⊗ ᾱ², s) on the maximal parabolic, are Eisenstein series for the arithmetic group Γ(N) ⊂ SL₃(C) with a finite-order character κ₃; Langlands' theory gives their continuation and W(A₂) functional equations. | framework [R-]. Identification of the Whittaker coefficient at (1, μ³) with R(z;μ) (Bruhat-cell computation as in BBFH, angle twist ᾱ² = non-spherical K-type) [G] |
| Functional equation of R(z;μ) | z ↔ 5/6 − z, with factor N(μ)^{3(1/2−z)}. The dual coefficients ρ(sμ³)γ₁(s)… couple s and μ. | conductor exponent k = 3 is [H] (`A2_theory.md` referee item 1) |
| Poles of R(z;μ) | None in Re z > 5/12 for non-square μ. The constant terms involve L(·, ᾱ¹²), which is entire. | [R-]/[G] (`A2_theory.md` §1.4 and referee) |
| Z₂ in u | Poisson in μ: Hecke L(u, χ_{m₁m₂}), with conductor N(m₁m₂). Maps Z₂ to a series of the same shape with ε-factors (§1.3). | [R] (bookkeeping [R-]) |
| Diagonal pole u = 1 | Residue is the Rankin–Selberg series Σ_m |π(m)|² N(m)^{−z₁−z₂}·(local factors). Its leading behaviour is the proved mean square. Continuation of θ₆ × θ̄₆ Rankin–Selberg beyond the edge needs Zagier regularisation (DFDH do this for n = 3, 4). | [R] at the edge; [G] beyond |
| A Rankin–Selberg integral producing Σ_μ R R̄ N(μ)^{−u} directly (picking out the cube indices (1, μ³)) | none known | [G] |

So the **rigorous** symmetries of Z₂ are only Poisson in μ, plus absolute convergence. The others are heuristic or not yet identified.

### 2.3 The groupoid and the reflection group

**Sawin–Whitehead groupoid** (`groupoid_moments.log`; dictionary of `A2_theory.md` §4.1, n = 6). Nodes s₁ [4], m₁ [4], μ [0], m₂ [2], s₂ [2]; edges s₁–m₁ 5, m₁–μ 3, μ–m₂ 3, m₂–s₂ 1.
- The root coefficients exceed 60 after 73 bases, so the groupoid is **infinite**.
- The reduced (m₁, μ, m₂) model is infinite too.
- Control: the ζ-free model (m-nodes plain [0]) is finite, A₃-like with 24 bases. This is the second moment of quadratic twists of a Dirichlet-type sequence, i.e. exactly why the ζ-free (V) of `V_derivation.md` was tractable.

[R as combinatorics for the blow-up itself. Reading it as "infinite" uses the fact that finite arithmetic root systems have bounded root coefficients; the cut-off 60 is a practical threshold, far above the coefficients met in the finite cases here (all ≤ 6) [S, Heckenberger's classification not re-checked]. Independently, the reduced reflection group below is *proved* infinite (hyperbolic) for k = 3.]

**Reduced reflection group** (`tits_hull.py`). Use centred unitary coordinates a = Re z₁′ − 1/2, b = Re z₂′ − 1/2, c = Re u − 1/2. Let k be the conductor exponent of R(z;μ) in N(μ). The reflections are:
- s_c: (a,b,c) ↦ (a+c, b+c, −c) (Poisson in μ, [R]);
- s_a: (a,b,c) ↦ (−a, b, c+ka) and s_b: (a,b,c) ↦ (a, −b, c+kb) (functional equations of R in z₁, z₂, [H]).

The orders of the products are |s_a s_b| = 2 and |s_a s_c| = |s_b s_c| = m_k:

| k | m_k | Coxeter type | invariant form | provable continuation |
|---|---|---|---|---|
| 1 | 3 | (2,3,3) = A₃, finite (order 24) | definite | all of C³ |
| 2 | 4 | (2,4,4) = C̃₂, affine | degenerate; null root δ = (−1,−1,2) | the Tits cone, a half-space |
| **3** | **6** | **(2,6,6), hyperbolic** | signature (2,1) | the Tits cone, a light-cone-type region |
| 4 | ∞ | (2,∞,∞), hyperbolic | (2,1) | the same |

[R as group theory; the invariant form is computed in `tits_hull.py`.]

### 2.4 Bochner + Phragmén–Lindelöf: what bound comes out [R as computation, H inputs]

**The bound.** At a point (a,b,c) of the domain of holomorphy and polynomial growth, the Mellin identity gives

  Σ_μ|F′|² ≪ Q^{5/6}·Q^{a+b}·L^{1/2+c},

so the exponent is f := a + b + λ(c + 1/2), plus polar residues.
- Absolute convergence C₀ = {a, b, c > 1/2} gives f > 1 + λ (trivial).
- The large sieve is f = 1.
- (V′_LS) is f < 1.

**The hull.** Bochner's tube theorem plus convexity gives holomorphy and polynomial growth on the convex hull of ∪_g g(C₀). The infimum of the linear f over the hull is the minimum over the image cones g(C₀) of their apex values. If some cone has a recession direction with f_lin < 0, the infimum is −∞ modulo polar terms.

**Enumeration** (5634 group elements, word length ≤ 16):

| k | λ | min apex f | negative recession direction? | reading |
|---|---|---|---|---|
| 3 | 0.75 / 0.789 / 0.837 | 1 + λ (attained at C₀ itself) | no | **trivial bound only, worse than the large sieve** |
| 3 | 0.5 ≤ λ < 2/3 | 4λ − 1 (word ab) | no | no saving over the large sieve unless λ < 1/2 |
| 3 | 0.45 | 0.70 | no | saving only for λ < 1/2, i.e. M < 2/3 (cf. `A2_theory.md` §2.6) |
| 2 | 0.80 / 0.85 / 0.88 | < 0 (f_lin(δ) = 2λ − 2 < 0) | yes | the affine group *would* give (V′_LS), modulo polar residues and the coupled-dual convergence [G] |
| 1 | any | −∞ | yes | finite type: full asymptotic |

**Lindelöf on average.** `tits_membership.py` asks which points (a, a, c) lie in the hull. By Landau's theorem applied to the positive-coefficient series Σ_μ |R(σ+it;μ)|² N(μ)^{−u}, Lindelöf on average on Re z′ = σ is equivalent to holomorphy at (σ−1/2, σ−1/2, 1/2+).

For k = 3 the hull meets c = 0.51 only for a ≥ 0.497, and c = 1 only for a ≥ 1/3. This matches the line 3a + c ≥ 2 at all three slices tested (c = 0.51, 1, 2), i.e. **pointwise convexity of each R(z;μ) plus absolute convergence in μ**: averaging over μ gains nothing. Lindelöf on average (a = 0 at c = 1/2+) is reached only for k = 1.

**One-variable version (Hardy–Gabriel convexity in σ) [R as algebra].**
- Take Σ_μ |R̃(σ+it;μ)|² ≪ L at σ = 1 (large sieve), and ≪ L^{1+3}K at σ = 0 (functional equation with K the sublattice mean-square constant; trivial K ≤ L).
- Log-convexity gives Σ|F′|² ≪ Q^{2σ−1/6} L^{4−3σ} K^{1−σ}.
- This decreases in σ when λ > 2/(3 + log K/log L). The optimum is σ = 1, which is the trivial Q^{11/6}L.

**Conclusion.** With k = 3, MDS continuation plus convexity cannot reproduce even the large sieve. The large sieve's (L + Q) is a duality (ℓ²) statement, not a contour-shift statement.

### 2.5 Smoothings and auxiliary families do not help [R as algebra]

| attempt | result |
|---|---|
| Choice of weights Φ, Ω | Changes only the Mellin weights; the obstruction is the region of continuation (§2.4). |
| Enlarging the μ-family to size L′ ≥ Q (positivity) | (L′ + Q)‖π‖² ≥ Q‖π‖² |
| Amplifier A(μ) = Σ_{ν∼N} ψ̄_{μ₀}(ν)ψ_μ(ν) | |F′(μ₀)|² ≤ (L/N + Q)‖π‖², so no gain at any N |
| Hölder / fourth moment | Needs Σ|F′|⁴ ≪ Q^{11/3}L^{−1}. The large sieve for π∗π gives ≍ Q^{11/3}, which is worse. |
| Extra archimedean average over \|v\| ≤ Z^ε | A family of size Z^ε: negligible |
| A Rankin–Selberg integral selecting the indices (1, μ³) | No candidate [G]. The indices are cubes, a sparse set of Whittaker indices of a GL₃^{(6)} Eisenstein series in its spectral variable. |

### 2.6 Verdict C(i)

**No upper bound Σ_μ|F′|² ≪ Q^{1−η}‖π‖² follows from known theory.** The rigorous inputs (Poisson in μ, rank-one Kubota theory) are far too weak. The heuristic full set has an infinite hyperbolic group (k = 3), whose Bochner hull is the pointwise-convexity region [R as computation, H inputs].

The single parameter that would change this is k:
- k = 2 (affine) would put (V′_LS) inside the Tits cone;
- k = 3 does not.

---

## 3. C(ii): spectral and bilinear tools

### 3.1 Where Kloosterman sums actually appear [R]

1. **In Σ_μ|F′|².** After Poisson in μ there are only quadratic Gauss sums g_χ(k) (§1.3). No Kloosterman sum appears, so a Bruggeman–Kuznetsov formula has nothing to act on.
2. **In T (Lemma C0).** T = Σ_A μ(A)ν(A)q_A^{−1}V_G(q_A/Z)·Σ_m ρ(m)K(q_m/Q)ξ(m)·S₆(−m, 0; A), with S₆(h,0;A) = ğ(h;A) the degenerate (n₂ = 0) Kloosterman sum for the sextic multiplier.
3. **At the A₂ level.** F′(μ) = Res_w Σ_{s,m} H((s,m); 1, μ³) is a Whittaker coefficient of a GL₃^{(6)} Borel Eisenstein series. Its Bruhat expansion involves the degenerate GL₃ Kloosterman sums of Bump–Friedberg–Goldfeld on the 6-fold cover [S].

### 3.2 Bruggeman–Kuznetsov [R]

- **Louvel, Forum Math. 2014, Thm 4.2** [S] is a sum formula for S₃(m,n;c) on the **cubic** cover with a smooth weight φ(c). It does not fit:
  - our multipliers are sextic (χ_A) and quadratic (ψ_μ); the cubic content (γ₂, the cubic theta) sits in the coefficients, not the multiplier;
  - the weight in the modulus is μ(A)·(θ₆-pairing) (item 2 above), which is not smooth;
  - the Kloosterman sums are degenerate (n = 0), whose sum formula is the Eisenstein coefficient, i.e. the Kubota series itself. That just returns D(w;m).
- **A sextic-cover SL₂(O) Kuznetsov.** Available in principle from Bruggeman–Motohashi / Lokvenec-Guleska for any unitary multiplier on a congruence subgroup [R- framework]. It fails for the same two reasons: non-smooth Möbius × θ₆ weight, and degeneracy.
- **A GL₃^{(6)} Kuznetsov formula.** None exists. Friedberg–Offen 2025 give relative-trace-formula distributions on the cubic cover of SL₃ only [S].
- **Spectral large sieve.** No spectral large sieve for metaplectic Maass forms is known; Watt's is non-metaplectic (`LIT_ideas234.md` §3) [S]. Even with one, the family needed is the cube-indexed Whittaker coefficients W_{(1,μ³)}, a sparse subsequence. The Rankin–Selberg / Lindelöf-on-average count over all indices n ≤ L³ costs L³ against the required L.

### 3.3 Duke–Iwaniec / Dunn bilinear bounds [R as algebra]

The bounds:
- DI: Σ_{a≤A, b≤B} μ²(a) α_a β_b c_g(ab) ≪ (B^{1/2} + AB^{1/4})‖α‖‖β‖ (half-integral weight);
- Dunn (cubic cusp forms): ((AB)^{1/2} + A^{3/2}B^{1/4})‖μ²α‖_∞‖β‖₂.

Both require coefficients **at products ab of a fixed form**.

- **Primal form Σ_{μ,m} β_μ π(m) ψ_μ(m).**
  - The kernel is a character, not c(μm).
  - π is not a Fourier coefficient of a fixed form on a cover of GL₂. It is a WMDS residue, a GL₃^{(6)} Whittaker coefficient in the spectral variable (`A2_theory.md` §3.3).
  - So it is not of DI shape.
- **The A₂ dual** (`A2_theory.md` §2.3, using KP: τ(sμ³) = G₄ χ̄_μ(s)⁴ τ(sμ) N(μ)^{…}). The dual is T′_dual ≈ Σ_{μ∼L, s∼S₀} x_μ y_s χ̄_μ(s)⁴ τ(sμ), with S₀ = L³/Q.
  - This is a product-coefficient form, but with an extra coupling kernel χ̄_μ(s)⁴, which DI does not allow.
  - Ignore the kernel anyway. The trivial dual bound is Q^{−1/12}L² (`A2_theory.md` §2.6), and Q^{11/12} needs a further factor Q/L² = Q^{1−2λ} = Q^{−0.58}.
  - The ratio (DI bound)/(trivial ℓ² bound) for flat x, y is:

| roles | ratio | at λ = 0.79 |
|---|---|---|
| A = #μ = L, B = #s = S₀ | L^{−1/2} + (Q/L)^{1/4} | ≥ Q^{0.05} > 1 |
| A = S₀, B = L | S₀^{−1/2} + L^{5/4} Q^{−1/2} | ≥ Q^{0.49} > 1 |

  - The second term is ≥ 1 for every λ > 0.4 (in the first row it is ≥ 1 for all λ ≤ 1). So DI-type bounds **never** beat the trivial dual bound here.
  - Two idealised cases for comparison:
    - square-root cancellation in μ alone (factor L^{−1/2}) reaches the needed Q^{1−2λ} only when λ ≤ 2/3;
    - square-root cancellation in s for each μ (pointwise Lindelöf for every twist, factor S₀^{−1/2}) would suffice for all λ ≤ 1. DI-type bounds supply that only through their first term; their second term destroys it.

### 3.4 Type I/II on the Möbius form [R + H]

By Lemma C0, T = Σ_A μ(A)·f(A) with f(A) = ν(A)q_A^{−1}V_G Λ(A).
- **Type I.** A Heath-Brown/Vaughan decomposition of μ(A) needs Type-I sums Σ_{d≤D} μ(d) Σ_e f(de). By the coherent Voronoi step, Λ has a smooth main term ∝ φ(A)q_A^{−5/6}Q^{5/6}. The Type-I main terms therefore assemble into Σ_d μ(d)·(smooth). This is the zero-sensitive signal of §4.2, not bounded without zero-free input [R].
- **Type II.** Bilinear sums Σ α_d β_e f(de) reduce, after Cauchy–Schwarz in m, to the sextic large sieve for products ğ(·,d)ğ(·,e). That sieve has the theta excess, which is sharp by DFDH-type lower bounds [H, by analogy with n = 3, 4 (DFDH [V via LIT])].

### 3.5 Verdict C(ii)

No spectral or bilinear tool applies:
- Kuznetsov: wrong cover, non-smooth Möbius × θ₆ weights, degenerate sums. GL₃^{(6)} sum formulas and metaplectic spectral large sieves do not exist.
- DI/Dunn: not applicable to the primal form, and numerically insufficient on the dual for λ > 0.4.

---

## 4. C(iii): the first moment (V′_β) as a residue of a triple series

### 4.1 The series and its groupoid [R as combinatorics] (`groupoid_moments.log`)

| series | groupoid |
|---|---|
| reflected (s, m, μ), μ = cubic theta coefficient [5] | 4 objects, 24 bases, **contains a non-admissible (Möbius) object** |
| unreflected (s, m, A), A = γ₂ᾱ [5] | the same groupoid |
| **faithful 4-node (t, A, m, s)** (t = cubic Kubota variable, A Dirichlet node before the cubic-theta residue) | **finite: 7 objects, 336 bases, no Möbius node**. Sawin–Whitehead give rationality over F_q[T] [V via `A2_theory.md` §4.2]. |
| unweighted columns β ≡ 1, (s, m, μ) with μ [0] | **infinite** (3716 bases before blow-up) |
| conjugate column orientation (s, m, A), M_AA = 1 (control) | finite, 3 objects, 48 bases, no Möbius node |

So the *natural* column weights (Patterson's cubic theta coefficients) are exactly what make the first moment finite type. In the number-field case (Q(ω), mixed cubic/sextic Kubota nodes, μ₁₂ ⊄ K), meromorphic continuation of the 4-variable series is **[G]**. The ingredients would be:
- rank-one Kubota equations for n = 3 and n = 6 [R];
- local p-parts at good primes [R-];
- vector-valued Kubota-class bookkeeping because μ₁₂ ⊄ K [G];
- the angle twist ᾱ (non-spherical K-type) [G];
- a groupoid version of Bochner [R-].

### 4.2 Why this is still useless: T carries the signal

**Lemma C1 (the zero line) [R-, via `V_gaps.md` §1.6].** Let 𝒯(z,u) be the Mellin transform of T in (Q, Z). By Lemma C0 and coherent Voronoi of θ₆ at the cusps x/A (the multiplier there is χ_A(x)^{+1}, Lemma M), 𝒯 has the polar line z = 5/6 with residue

  Σ_A μ(A)ν′(A)(φ(A)/q_A)c_{A₀} q_A^{−5/6−u} = Σ_χ c_χ H_χ(u + 5/6)/L(u + 5/6, χ),

where:
- χ runs over finitely many finite-order Hecke characters (from ν′ and the class function c_{A₀});
- H_χ(s) = Π_p (1 + χ(p)N(p)^{−s−1}/(1 − χ(p)N(p)^{−s})) is holomorphic and bounded in Re s > δ.

So 𝒯 is singular at {z = 5/6} × {u + 5/6 = ρ} for every zero ρ. The corresponding vertex term is S ≍ Q^{5/6}Z^{β−5/6}.

**Lemma C2 (S is the whole signal: c = 1) [R-].** Fix Q, Z and let Y → ∞ (so X = Q/(qY) → 0). This is legitimate because `low-separated` and `probe-poisson-identity` are identities for all (X, Y, Z).

1. **Low side.** By `low-separated`, Mellin inversion in the modulus variable s of A_{m,σ,v}(Y), and the v-integration of `V_gaps.md` §6,

   J(Y) = Q^{−1/2} Y^{−1/3} T_K(Q,Z) + b(Y),   b(Y) = Q^{−1/2}∫ dv Ŵ₀(iv) Y^{−1}∫_{(1/2+ε)} Ŵ₁(w−1/2+iv) Y^w 𝒟_v(w) dw ≪_{Q,Z} Y^{−1/2+ε}.

   Here T_K = Σ_σ Σ_m K(q_m/Q)ξ(m)ρ_σ(−m)B_{m,σ} is the theta pairing with the weight K of `V_gaps.md` §6, i.e. "Ŵ T" in the shorthand below. The bound on b holds because 𝒟_v(w) = Σ_m Ωξ(q_m/Q)^{−iv}D(w;−m)B_m is a **finite** combination of the D(w;−m). Each of these has only the pole w = 2/3 in Re w > 1/2 and polynomial growth (FE_LoA6 Thm A).
2. **High side.** In the weight 𝒲 = X^{1/2−z}Z^{x+z−1}Y^{w−1}… of `probe-poisson-identity`, the signal sits at w = 1, z = 1/6 (`IDEA1_NOTES.md` §2). So f = X^{1/3}M(1/6)Ŵ₁(1)·F_η(Z), which is exactly ∝ Y^{−1/3} at fixed (Q, Z).
3. **Comparison.** Comparing the Y^{−1/3} components of J, the signal is contained in Y^{−1/3}Ŵ T. Precisely,

   ŴT = Y^{1/3}f + T_E,

   where T_E collects the high-side terms with w + z = 7/6.

The two weight dependences agree: K̂(5/6) = Ŵ₁(1)Ŵ₀(5/6) (`IDEA1_NOTES.md` §2).

*What the argument uses.* The high side of `probe-poisson-identity` is evaluated by contour shifts in (x, w, z). It is therefore a sum of terms with definite Y-exponents (Y^{w−1}X^{1/2−z} ∝ Y^{w+z−3/2} at fixed Q), plus remainders with smaller Y-exponent. Two consequences:
- the terms with Y-exponent > −1/3 must cancel among themselves, since the low side has none;
- the Y^{−1/3} terms must equal Q^{−1/2}Y^{−1/3}T_K.

That is the [R-] content. It does not depend on the probe's geometry, because X, Y and Z are free in both identities.

### 4.3 Proposition C3: circularity [R as logic, given Lemmas C1 and C2]

**(a) Bounds through the signal vertex return the input region.** Suppose T is bounded as |T| ≤ |S| + |T′|, with S estimated through an a priori zero-free half-plane β* ≤ θ_k. That estimate is |S| ≪ Q^{5/6}Z^{θ_k−5/6+ε}, valid since 1/L ≪ |t|^ε on Re s ≥ θ_k + ε. This is a saving of

  s₀ = σ_θ − θ_k,   σ_θ = 5/6 + M/12 slot-free; in general σ_θ = base + b/12, with base = 5/6 + l_x/6 − ℓ/6 − d the LP low-side base.

In the LP of `wf/lp_const_saving.py`:
- the theta constraint becomes σ ≥ base + max(0, g, …)/2 + e/2 ≥ θ_k + e/2 ≥ θ_k;
- so θ_LP ≥ θ_k identically.

`lp_signal_circularity.log`:

| input | θ_LP |
|---|---|
| θ_k = 0.874957 | 0.874957 |
| θ_k = 0.86 (hypothetical) | 0.86 |
| a *signal-blind* constant s₀ = 0.027 (for orientation) | 0.858612 |

So the paper's own Theorem (β* ≤ 7/8 for all finite-order Hecke L over K) bounds S with a positive saving M/12 − 1/24 (= 0.027 at M = 0.826), but this buys exactly nothing.

**(b) T′ carries no zero information.** Put T′ := T − Y^{1/3}f/Ŵ = T_E/Ŵ. By Lemma C2, J − f = Y^{−1/3}Ŵ T′ + ⟨A^rest, B⟩. By the high side the left-hand side is E. So bounding T′ and the rest bounds the high-side error, which the paper already controls. It says nothing about the low-side contract |J| ≪ Z^{C(σ₀)+ω}, which still requires bounding T = Y^{1/3}f/Ŵ + T′, i.e. f itself. So no bound on T′ alone, however strong, feeds the continuation criterion (Prop. `continuation-criterion`). This holds also inside its contradiction framework, where β* may be used in the low side only with a power saving below Z^{C(β*)}: the S-vertex is exactly of size Z^{C(β*)}.

**(c) Consequence.** A usable (V′) must bound T **without evaluating its z = 5/6 polar line**: Cauchy–Schwarz, large sieve, mean value, positivity. In particular:
- finite-type MDS evaluation of the first moment (§4.1), Voronoi/Poisson chains that collect main terms, and Type I/II with main terms (§3.4) are all circular;
- `TARGET_VLS.md` §3's suggestion "S handled together with the high side" is vacuous when c = 1;
- conversely (V′) is itself a zero-free statement in disguise: (V′) + |T_E| small ⇒ β* ≤ σ_θ − s₀.

### 4.4 Two non-circular first moments, and why they do not help

- **Conjugate-column control** (M_AA = 1).
  - It is an honest finite arithmetic root system (no Möbius node). Its double residue is zero-blind and in principle computable [G in number fields].
  - But it is a *different* pairing. Without the Davenport–Hasse alignment there is no Möbius collapse, so there is no signal, and the probe needs the signal. Not usable [R].
- **Twisted first moments** M_χ = Σ_μ β_μ χ̄(μ) F′(μ), χ mod q′ with q′ > L.
  - Parseval gives Σ_μ |β_μ|²|F′(μ)|² = φ(q′)^{−1} Σ_χ |M_χ|². So (V′_LS) is equivalent to an average bound for the M_χ.
  - For χ ≠ 1 the column twist breaks the alignment, so M_χ may be zero-blind.
  - But its twisting conductor q′ ≳ L puts any MDS evaluation in the conductor-aspect regime with family ≈ conductor, which is beyond the edge again. Each M_χ is a pairing of the same "random" type as T [H]. No reduction.

### 4.5 Verdict C(iii)

- The first moment is the *right* object for MDS theory: natural weights make it finite type.
- But for that very reason it is zero-sensitive, and its evaluation is circular: it returns any input zero-free region unchanged [R as logic].
- Its continuation gives no power-saving bound usable for zero-free regions.

---

## 5. Net status

### 5.1 What Attack C establishes

| item | statement | label |
|---|---|---|
| Lemma C0 | B_m = Σ_A μ(A)ν(A)q_A^{−1}ğ(−m;A)V_G: a Möbius-twisted conjugate sextic Kubota sum | [R], [N] 1.8·10⁻¹⁵ |
| §1.3 | Poisson duality for Σ_μ|F′|² is self-dual with mean-zero dual weights; no reduction | [R-] |
| §2.3 | the second-moment reflection group is Coxeter (2, m_k, m_k): A₃ / C̃₂ / hyperbolic for k = 1 / 2 / ≥3; the SW groupoid of the true 5-node series is infinite | [R] / [R as combinatorics] |
| §2.4 | for k = 3, λ ∈ [0.67, 1): the Bochner hull gives only the trivial bound; the Lindelöf-on-average points are outside the hull; saving over the large sieve only for λ < 1/2 | [R as computation; H inputs] |
| §3.3 | DI/Dunn bilinear bounds never beat the trivial dual bound for λ > 0.4 | [R as algebra] |
| §4.1 | natural column weights make the first moment finite type (7 objects); unweighted columns do not | [R as combinatorics] |
| §4.2–4.3 | T contains the full signal (c = 1); any evaluation-based bound is circular; LP returns θ_k exactly | [R-] (Lemma C2), [R as logic and computation] |

### 5.2 What (V′) is reduced to

A **signal-blind** beyond-the-edge mean value: (V′_LS), or any refinement of Cauchy–Schwarz on ⟨ρ, B⟩ that never evaluates the z = 5/6 polar line. Any of the following would supply it:

1. **An effective conductor exponent k ≤ 2 for ψ_μ-twists of π, together with a decoupled dual.**
   - The second-moment series would then be affine (C̃₂) and (V′_LS) would lie inside the Tits cone (§2.4).
   - Against this: KP gives τ(mp³) = G₄τ(mp), with no forced vanishing (`A2_theory.md` §2.3). The referee's best heuristic is k = 3. The numerics cannot resolve k (`A2_numerics.md` §4).
2. **A spectral large sieve on the 6-fold cover of GL₃** for Whittaker coefficients at the sparse indices (1, μ³), μ ∼ L, saving over the Rankin–Selberg count. No literature, not even for GL₃ cusp forms in the non-metaplectic setting [H].
3. **A genuinely arithmetic large sieve exploiting π = ρ·γ₁·ᾱ².** None is known. The cubic analogue collapses (degree one) only because τ₃(mp²) = 0; the sextic j = 3 analogue does not (`A2_theory.md` §3.2).

### 5.3 Gaps (if one nevertheless pursued the MDS route)

| gap | content |
|---|---|
| [G1] | Identification of R(z;μ) with the (1, μ³) Whittaker coefficient of E_P(θ₆ ⊗ ᾱ², ·) on the BMS-split cover of SL₃(O); non-squarefree p-parts |
| [G2] | Number-field continuation of the 4-node (t, A, m, s) series (mixed n = 3, 6; μ₁₂ ⊄ K; angle twist) |
| [G3] | Conductor exponent k (heuristic 3) |
| [G4] | Zagier-regularised θ₆ × θ̄₆ Rankin–Selberg beyond the edge |

None of these would rescue (iii), which is circular for structural reasons. [G1], [G3] and [G4] matter for (i) only if k ≤ 2.

### 5.4 LP

Nothing in Attack C supplies a signal-blind s₀ > 0. For orientation only (`lp_signal_circularity.log`, the P²/Y piece dropped by LoA6), a signal-blind uniform saving would give:

| signal-blind s₀ | θ |
|---|---|
| 0.027 | 0.8586 |
| 0.02 | 0.8615 |
| 0.01 | 0.8668 |

By contrast, a saving of the same size obtained by bounding the signal vertex through a known zero-free region returns that region, e.g. 0.874957 (§4.3).

**Honest status: θ = 0.874957 (Liu).**

---

## 6. Files

- `vls/C_moments_spectral.md`: this note.
- `vls/tits_hull.py`: reflection group, Coxeter type, invariant form, hull bound. Arguments: λ, k, word length.
- `vls/tits_membership.py`, `vls/tits_membership.log`: hull membership of the Lindelöf-on-average points.
- `vls/groupoid_moments.py`, `vls/groupoid_moments.log`: Sawin–Whitehead groupoid finiteness for all moment series.
- `vls/check_lemmaC0.py`, `vls/check_lemmaC0.log`: Lemma C0.
- `vls/lp_signal_circularity.py`, `vls/lp_signal_circularity.log`: the LP with a signal-vertex bound (θ_LP = θ_k) against a signal-blind constant saving.

---

## Referee report

*Independent skeptical referee, 9 Oct 2026. I read this note, `TARGET_VLS.md`, `REPORT.md` (CORRECTION banner), `wf/V_gaps.md` §0, §1.6–1.7, `IDEA1_NOTES.md` §1–3, `wf2/A2_theory.md` §4, `LIT_ideas234.md`, and `paper78.tex` l.6200–6285. I re-ran every script of this note.*

### Verdict

Attack C claims no saving, and there is none: θ stays at Liu's 0.874957. Nothing positive is asserted, so there is no proof to refute (refuted = false). The negative conclusions survive. Several statements are overstated, mislabelled or not new; corrections follow.

### R1. Reproduction [N]

| script | result |
|---|---|
| `check_lemmaC0.py` | re-run: 644 pairs, max error 1.76·10⁻¹⁵. Only products of **split** primes are tested; inert-prime moduli are untested numerically (the algebraic proof covers them, given eq. `signal`). |
| `groupoid_moments.py` | re-run: all 8 lines identical to the log. |
| `tits_hull.py` | re-run at k = 3, λ = 0.45 / 0.6 / 0.79 / 0.837, and k = 2, λ = 0.79: identical (0.70 via `cab`; 4λ−1 via `ab`; 1+λ at C₀; k = 2: 107 cones with negative recession). |
| `tits_membership.py` | re-run, k = 2 and 3: identical. |
| independent BFS (referee) | 225,156 distinct elements of the k = 3 group (the note used 5,634). Minimum apex f = 1+λ at λ ∈ {0.667, 0.70, 0.79, 0.837}, and 4λ−1 at λ ∈ {0.5, 0.6}. No image generator has f_lin < 0. |
| `lp_signal_circularity.py` | as logged. **But** the `zf` mode inserts the Gram piece g = 2(θ_k − base), so the LP constraint is σ ≥ base + g/2 = θ_k *by construction*. The rows "θ_k → θ_k" are an algebraic identity, not independent evidence. The logic of Prop. C3(a) is nevertheless correct. |

### R2. Part (i), the second moment

**(a) k = 3 is not just [H] in the model the note itself uses.** Apply Heckenberger's rule (`A2_theory.md` §4.1) to the reduced object (m₁, μ, m₂), M = [[4,3,0],[3,0,3],[0,3,2]]:
- q_{m₁} = −ζ⁴ = ζ, of order 6, and q̃ = ζ³, so a_{m₁μ} = −3 and a_{μm₁} = −1;
- likewise q_{m₂} = ζ⁵ gives a_{m₂μ} = −3;
- the Cartan matrix is [[2,−3,0],[−1,2,−1],[0,−3,2]], i.e. Coxeter (2,6,6). This is exactly the k = 3 group of `tits_hull.py`.

So the groupoid line "finite=False" and the reflection-group computation are the same fact. k = 3 is forced by Lemma Z's node type: γ₁γ₃, an order-6 node with a quadratic edge. An affine (k = 2) outcome would need an order-3 (cubic-type) row node, i.e. a different row sequence.

Hence §5.2 item 1 and gap [G3] ("the decisive parameter k") should be downgraded. Within the Sawin–Whitehead dictionary, k ≤ 2 is excluded. It stays open only in the sense that the number-field identification [G1] is not written.

**(b) Rigorous form of "no negative recession direction" for k = 3.**
- In root coordinates the invariant form is B = [[1,0,−k/2],[0,1,−k/2],[−k/2,−k/2,k]].
- The target functional is f = α_a + α_b + λα_c, with B(f,f) = 2 − 6λ + 3λ².
- This is negative exactly for λ ∈ (1 − 1/√3, 1) = (0.4226, 1).
- So f lies in the imaginary cone of the compact hyperbolic (2,6,6) group, hence f ≥ 0 on the whole Tits cone, and the infimum over the hull is attained at short words.

This upgrades the enumeration-to-length-16 claim to [R as group theory] for all word lengths, given the heuristic reflections. For k = 2, B(f,f) = 2(1−λ)² > 0, so f is unbounded below on the half-space Tits cone, consistent with the note.

**(c) "The hull coincides with pointwise convexity on the slices tested" is not diagnostic.** `tits_membership.py` with k = 2 shows the same coincidence (2a + c ≥ 3/2 at c = 0.51, 1, 2), although the k = 2 hull is unbounded below in f. Only the recession / imaginary-cone analysis in (b) separates the cases.

**(d) Overstatement.** §0(i) and §2.6 say "no upper bound ≪ Q^{1−η}‖π‖² follows from any MDS argument / known theory". What is shown is narrower: it does not follow from Bochner + Phragmén–Lindelöf convexity on the (heuristic) reflection group.
- Hybrid arguments are not covered, e.g. a priori large-sieve input on vertical lines combined with the functional equations.
- They are not shown to fail, although nothing suggests they would succeed.

**(e) The k ≤ 2 escape is weaker than "modulo polar residues" suggests.**
1. By Cauchy–Schwarz, |T| ≤ (Σ|β|²)^{1/2}(Σ_μ|F′|²)^{1/2}.
2. By Lemmas C1–C2, T contains S ≍ Q^{5/6}Z^{β−5/6}.
3. So (V′_LS) with exponent η implies β* ≤ 5/6 + M/12 − Mη/2, given control of T′ at that scale. At M = 0.826 this beats Liu as soon as η > 0.066.
4. The k = 2 hull (f → −∞, so η up to 1 − λ) would give β* ≤ 1/3 + 7M/12 ≈ 0.80–0.84.

A continuation argument whose polar divisors are only explicit Hecke-L hyperplanes cannot plausibly deliver that. In any k ≤ 2 scenario the polar terms would have to carry the zeros. The "polar residues" caveat is the whole difficulty, not a technicality. The same remark applies to every route to (V′_LS): it is a zero-free-region statement (as §4.3(c) says) and must be proved signal-blind.

**(f) §2.2, route (a) (Chinta–Gunnells).** CG, like BBF, work over fields containing μ_{2n}, and here μ₁₂ ⊄ Q(ω). The Kubota-class-vector workaround is not written, so the label should be [G], not "framework [R-]". Route (b) (Bass–Milnor–Serre splitting plus Langlands) is the right framework, with the identification [G] as stated. This does not affect the negative verdict.

**(g) §1.3 bookkeeping.**
- The m₁ = m₂ terms are counted in the diagonal and again inside Σ_{k≠0}Φ̌|F″(k)|². They must be subtracted; by the mean-zero identity they contribute −Φ̌(0)-type mass.
- The conclusion is unaffected: no positivity is gained, and Poisson plus the large sieve on the dual returns (L+Q)‖π‖².

### R3. Part (ii), spectral and bilinear tools

- I re-derived the DI/Dunn ratios:
  - (A,B) = (L,S₀): L^{−1/2} + (Q/L)^{1/4} ≥ Q^{0.05};
  - (A,B) = (S₀,L): S₀^{−1/2} + L^{5/4}Q^{−1/2}, which is ≥ 1 iff λ ≥ 0.4;
  - the needed factor is Q/L² = Q^{1−2λ}.

  All correct.
- Lemma C0 and the identification S₆(h,0;A) = ğ(h;A) (a Gauss sum, i.e. a degenerate Kloosterman sum) are correct.
- Kuznetsov needs smooth modulus weights, so the μ(A)·(θ₆-pairing) weight is indeed excluded.
- Louvel Thm 4.2 and Friedberg–Offen 2025 are [V] only at abstract level in `LIT_ideas234.md`. I did not re-check them. The citations "Eckhardt–Patterson / Bröker–Hoffstein for n = 6" in §2.2 are unverified and not needed, since Kubota 1969 suffices.
- No objection to the conclusion of (ii).

### R4. Part (iii), the first moment and circularity

**Lemma C1 is well supported.** Its residue H_χ(u+5/6)/L(u+5/6,χ) has the same shape as the paper's own signal integrand Z^{C(s)}e^{(s−5/6)²}H_η(s)/L_F^S(s,η) (paper78 l.6208, l.6282), with X^{1/3}Z^{s−5/6}Z^{ℓ/6} = Z^{C(s)} (l.6275). The exponents match: Q^{1/2}Y^{1/3}·X^{1/3} ∝ Q^{5/6}, which is the z = 5/6 polar line of Lemma C1. So c = 1 is well supported.

**Lemma C2 has one unaddressed step.** The Y → ∞ argument at fixed (Q,Z) needs the high-side contour-shift expansion to stay valid as X = Q/(qY) → 0, i.e. z shifted far left with polynomial growth. The note does not address this, so [R-] is generous there. Independent support for c = 1:
- `IDEA1_NOTES.md` §2: the theta term carries the signal with coefficient K̂(5/6), and "any combination that kills the theta term kills the signal";
- `V_gaps.md` §1.6: 𝒮 is the Voronoi main term of the unreflected T.

**Novelty.** The main conclusion of (iii) is already in `IDEA1_NOTES.md` §2, in `V_gaps.md` §1.6 ("writing T = 𝒮 + T′ only re-expresses the signal") and in `REPORT.md` ("the excess cannot be cancelled, only bounded"). §0's "the reasons are new" overstates this for (iii). What is new is the explicit statement that `TARGET_VLS.md` §3's route "S handled together with the high side" is vacuous. I agree with it.

**Prop. C3(b) needs a qualification.** It says T′ is "controlled by the paper's high side". But `V_gaps.md` §1.6 finds no direct unconditional bound for T′ below the CS scale (T′ ≪ Q^{−7/12}Z^{3/2+ε}). Control of T′ = T_E/Ŵ comes only through the high-side identity and the paper's high-side error bounds, which may be used inside the contradiction hypothesis. Either way, the conclusion stands: bounds on T′ alone carry no zero information.

**Scope.** The circularity is proved for *this* probe. The general principle is that the zero-free region equals the strength of a signal-blind bound on whichever low-side term carries the signal. That leaves open probes that relocate the signal (`V_gaps.md` §1.7(b)), which Attack C does not examine.

### R5. What survives

1. **Lemma C0** [R]. B_m = Σ_A μ(A)ν(A)q_A^{−1}V_G ğ(−m;A).
2. **T contains the full signal** (Lemmas C1, C2; [R-], with the X → 0 caveat).
   - Therefore any bound on T obtained by evaluating its z = 5/6 polar line returns the input zero-free region.
   - (V′), (V′_LS) and (V′_β) must be proved by signal-blind mean-value inequalities.
   - Each of them is itself a zero-free-region statement: with T′ controlled, (V′_LS) with η > 0.066 at M = 0.826 alone beats Liu.
3. **Second-moment MDS** (k = 3, forced by the Sawin–Whitehead dictionary).
   - Bochner + convexity give only the trivial bound Q^{5/6}·Q·L for λ ∈ [2/3, 1). This is rigorous as group theory for all word lengths, by R2(b).
   - It beats the large sieve only for λ < 1/2.
4. **DI/Dunn** never beat the trivial A₂-dual bound for λ > 0.4.
5. **Honest status: θ = 0.874957 (Liu).**

### R6. Labels and wording to change

- §0(i), §2.2, §5.3 [G3]: change "k = 3 [H]" to "k = 3 [R in the Sawin–Whitehead axiomatic model; number-field identification is G1]". In §5.2 item 1, say that k ≤ 2 is excluded in that model, and that even then the residues would have to be zero-sensitive (R2(e)).
- §2.2 route (a): change [R-] to [G] (μ_{2n} hypothesis).
- §2.6: change "from known theory" to "from MDS + convexity".
- §2.4: add the imaginary-cone argument of R2(b), and drop the slice coincidence as evidence.
- Lemma C2: keep [R-], but state the X → 0 expansion step as the unchecked point.
- §4.3: note that `lp_signal_circularity.py` encodes θ_LP ≥ θ_k algebraically.
