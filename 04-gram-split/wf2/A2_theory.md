# (V′) through the sextic A₂ multiple Dirichlet series

*9 Oct 2026. Builds on `REPORT.md` (CORRECTION banner), `IDEA1_NOTES.md`, `V_derivation.md`, `D_derivation.md`, `wf/AUDIT_chain.md`, `wf/V_gaps.md` (Lemma Z), `wf/FE_LoA6_proof.md`, and Bröker–Hoffstein (`scratchpad/bh.txt`).*

Scripts in this directory:
- `lp_dual_route.py`: the LP with the theta saving that the functional-equation dual can deliver (§2.6).
- `weyl_groupoid.py`: the Heckenberger/Sawin–Whitehead Weyl groupoid of the relevant multiple Dirichlet series (§4).
- `check_twisted_meansquare.py`: mean square of quadratic twists of the true row sequence ρ·ζ at Q = 3000 (§5).

Labels:
- **[R]** rigorous, given the cited results;
- **[R-]** routine but not written out;
- **[H]** heuristic;
- **[G]** gap;
- **[N]** numerical evidence only.

Citation labels:
- **[V]** statement read in the source;
- **[S]** known only through a secondary statement.

---

## 0. Verdict

**(V′) is not reachable through the A₂ structure. The conductor is N(μ)³, not N(μ)².**

1. **Identification [R].** On squarefree rows, a(m) := ρ(m)ζ_m equals ρ(m)·γ₁(m)·ᾱ(m)²·(ray). This is the w-residue (at the theta pole w = 2/3) of the sextic twisted A₂ Weyl-group multiple Dirichlet series (WMDS), with an angle-character twist on the second variable.
   - The quadratic twist (λ⁴μ/m)₂ is not a character of conductor μ in the metaplectic sense.
   - It is the WMDS **twisting parameter m₂ = μ³**, because χ_m(μ)³ g(m,h) = g(m, hμ³). §1.3.
2. **Residual functional equation [R-].** The one-variable residue R(z; μ) = Σ_m a(m)(λ⁴μ/m)₂ N(m)^{-z} has exactly one functional equation, z ↔ 5/6 − z (centre 5/12). Its factor is N(μ)^{3(1/2−z)}.
   - The dual coefficients are ρ(sμ³)·γ₁(s)·(angle, ray).
   - Kazhdan–Patterson / Bröker–Hoffstein give τ(mp³) = G₄(m,p)τ(mp), with |G₄| = 1 and no forced vanishing, so the dual is not thinned out.
   - Hence the analytic conductor is **N(μ)³** and the dual length is **L³/Q = Z^{1.13}**, against Q = Z^{0.826}. No saving. §2.
3. **Where the "N(μ)²" heuristics fail.**
   - The second reflection σ₁ carries no μ at all.
   - The metaplectic intertwining operators see only χ_μ⁶ = 1.
   §2.5.
4. **LP [R as computation].** Using the functional-equation dual with any available bound on it gives θ = min(0.874957, 0.934) = **Liu's value**. §2.6.
5. **Cubic cross-check [R].** The A₂ residue for n = 3 is Σ τ₃(m)G₁(1,m) = Σ G₁(1,m)². This is *not* constant.
   - The constant combination τ₃·Ḡ₁ is the opposite orientation, and it does not arise.
   - In the cubic case the conductor collapses to N(μ) because τ₃(mp²) = 0.
   - The sextic analogue τ₆(mp³) ≠ 0 does not collapse.
   - Bröker–Hoffstein's data give a(π)² = (patternless root of unity)·ᾱ(π)²·γ̄₂(π)·(algebraic junk)·(ray) [H]. There is no simpler object (no Hecke character, no GL₂ form, no square-root identity at squarefree indices). §3.
6. **Groupoid check [R as combinatorics].** In the Sawin–Whitehead / Heckenberger formalism:
   - the reflected and unreflected pairings lie in one Weyl groupoid;
   - its reflection at the column node *reproduces Lemma Z* (the m-node acquires the Gauss sum γ₁γ₃);
   - the Kubota reflection at m has n_{mμ} = 3, an independent confirmation of the N(μ)³ conductor;
   - the Dirichlet reflection at m (the paper's Poisson identity) turns the column node into a Möbius node. That is the Davenport–Hasse collapse c_A/τ(χ̄_A) = μ(A), i.e. the signal.

   The faithful four-variable object (cubic Kubota t, column A, row m, sextic Kubota s) is an honest rank-4 arithmetic root system (7 objects). None of its functional-equation images that I could identify makes the bilinear form short in both variables. §4.
7. **Numerics [N].** At Q = 3000 the true sequence a(m) behaves like a random sequence against quadratic characters:
   - mean |F′(μ)|²/‖a‖² = 0.99, over N(μ) ≤ 600;
   - the coherent control gives 4–13.

   So (V′) is very plausibly *true*, even with the ideal large-sieve saving. §5.
8. **What would have to be proved.** A large sieve with few moduli for this specific sequence, (V′_LS) of §6:

   Σ_{N(μ)≍L, μ sqfree} |Σ_{N(m)≍Q} a(m)(λ⁴μ/m)₂ W(N(m)/Q)|² ≪ Q^{1−2s₀}·Σ|a|², in the range Q^{1/2} ≤ L ≤ Q.

   This is a "beyond-the-edge" mean value for quadratic twists of a degree-3-conductor object. No current technique reaches it, even for GL(3) cusp forms [H]. (V′) is **plausible but not plausibly provable with present technology.**

---

## 1. The double series and its identification

### 1.1 The row sequence (recap)

From Lemma Z (`wf/V_gaps.md` §1.2, [R]): for a squarefree good row m,

  ζ_m = κ₀ ψ₁(m) γ̄₂(m)⁴ = κ₀ ψ₁′(m) · ᾱ(m)² · γ₁(m) γ₃(m).

Here:
- γ_j(c) = q_c^{-1/2} Σ_{v mod c} χ_c(v)^j e(v/c);
- γ₃ is a quadratic Gauss sum, hence a ray function up to quadratic-reciprocity signs;
- ψ₁′ is a ray character;
- ᾱ² is the Hecke Grössencharacter m ↦ m̄/m (on primary generators).

The pairing needs

  F′(μ) = Σ_m ρ(m) ζ_m (λ⁴μ/m)₂ W(N(m)/Q),   T′ = Σ_μ β_μ F′(μ),

with μ squarefree, N(μ) ≍ L = Q²/Z, and β_μ ∝ (cubic theta coefficient at λ⁴μ)·α(μ)·N(μ)^{-1/2}.

### 1.2 The double series [R for squarefree coprime terms]

  Z(w,z) = Σ_{s,m} g̃(s,−m) γ₁(m) ψ(m) χ(m) N(s)^{-w} N(m)^{-z},   ψ = ᾱ²·(ray)·γ₃,   χ = (λ⁴μ/·)₂.

For s, m coprime and squarefree, g̃(s,−m) = χ̄_s(−m) g̃(s,1), so the coefficient is

  g̃(s)g̃(m)·(m/s)₆^{-1}·(ray).

By reciprocity (m/s) = (s,m)_S (s/m), with the Hilbert symbol at S a class function. This is exactly the coprime coefficient H(c₁,c₂) = g(1,c₁)g(1,c₂)(c₁/c₂)^{-1} of the twisted A₂ WMDS with c₁ = s and c₂ = m. I checked this against the twisted-multiplicativity rule and the p-part of Chinta–Gunnells, *WMDS of type A₂* (arXiv math/0703040, eq. (2)–(3)) [V].

**Residue in w.** Summing over s first gives D(w;−m), the full Kubota series. Its residue at w = 2/3 is ρ(−m) (`FE_LoA6_proof.md`, Thm A [R-]). So

  Res_{w=2/3} Z(w,z) = Σ_m ρ(−m) γ₁(m) ψ(m) χ(m) N(m)^{-z} = R(z;μ),

which is the Dirichlet series of F′(μ).

**p-parts [R].**
- Our series has (k,l) ∈ {(0,0), (1,0), (0,1), (2,1)}.
- The stable A₂ WMDS (CG, eq. (3)) also has (1,2) and (2,2).
- They differ only at m with p² | m. Such rows are treated by the frozen-k block structure of `V_gaps.md` §4 and do not affect exponents [R-].

### 1.3 How the quadratic twist enters [R]

χ^6 = 1, so (μ/m)₂ = χ_m(μ)³ = χ̄_m(μ³). Since g(m,h)χ̄_m(a) = g(m,ah) for (a,m) = 1,

  **γ₁(m)·(μ/m)₂ = g̃(m, μ³)**  ((m,μ) = 1).

In WMDS language, H(c; 1, μ³) = (μ³/c₂)^{-1}H(c) = (μ/c₂)₂ H(c) on c coprime to μ. **The quadratic twist is the twisting parameter (m₁,m₂) = (1, μ³).**
- Terms with p | (m, μ) are where our series and the WMDS(1,μ³) p-parts differ. They are finite Dirichlet-polynomial corrections at p | μ [R-].
- For a general twist exponent j (twist by χ_m(μ)^j), the parameter is μ^{6−j}.
- The quadratic case j = 3 is the only one that is the same in both orientation conventions.

This is the metaplectic phenomenon behind the result. A character of order dividing n acting on the modulus variable of an n-th order Gauss sum is not a ramified twist. It is a *translate*: a Whittaker coefficient at another index. For the n-fold cover the principal series is determined by the character on n-th powers, and χ_μ⁶ = 1. So the inducing data stay unramified, and the whole effect of μ is the index μ³.

### 1.4 The ray and angle factors

- **γ₃, ψ₁′, κ₀, S-adic Hilbert symbols.** These are finite-order class functions: Ψ-functions in the BBF/CG sense, or congruence-class vectors in Kubota's sense [R-]. In the Sawin–Whitehead dictionary of §4, γ₃ is a node of exponent 0 ("ray") [R].
- **ᾱ² on c₂.** This is a Grössencharacter of infinite order. It corresponds to a Borel Eisenstein series on the 6-fold cover of GL₃ induced from a torus character with nontrivial archimedean angle in the c₂-direction, i.e. inducing data (1, 1, ψ). The Weyl group permutes these data.
  - The theta hyperplane H₁ (α₁, untwisted) stays polar.
  - H₂ and H₁₂ are **not** polar: the constant terms involve L(·, ᾱ¹²), which is entire.
  - Consequence [R-]: R(z;μ) has no poles in Re z > 1/2 for any μ, including μ = 1 and square μ. So F′(μ) has no main term at all. This corrects the ζ-free picture of `V_derivation.md` §6 and agrees with §5 below.
  - Under σ₂ the angle character moves to the c₁-direction. It affects only Γ-factors, not the μ-conductor.

### 1.5 Hypotheses of the literature, and what replaces μ₂ₙ ⊂ F

**What I checked.**
- Chinta–Gunnells, *Constructing WMDS* (JAMS 2010; arXiv 0803.0691, §1 and §2) [V]: "Let F be an algebraic number field containing the 2n-th roots of unity". Their §2 says it follows BBF08 §§2.3–2.7 "essentially verbatim", so the BBF twisted series assume the same [V, via CG].
- CG's A₂ paper (math/0703040) [V]: works over F_q(t) with q ≡ 1 mod 4n, and extends to "global field[s] containing the 2n-th roots of unity".
- BB06 and DIPP: also assume μ₂ₙ ⊂ F (`FE_LoA6_proof.md` §2 [V]).
- BBFH III (Annals 166, 2007): the abstract says "F containing the n-th roots of unity" [V, abstract only]. I did not read the body; its local and global setup follows BBF [S].

For n = 6 and F = Q(ω), μ₁₂ ⊄ F, so none applies verbatim.

**Where μ₂ₙ is used.** Only in global bookkeeping:
- (−1, x)_S = 1;
- unit-independence of H·Ψ;
- the isotropy of Ω = O_S^× F_S^{×n};
- Gauss-sum identities such as g(1,ϖ;ε^t)g(1,ϖ;ε^{−t}) = |ϖ|, which become ε^t((−1/ϖ))|ϖ| in general.

The extra signs (−1/ϖ)₆ are ray characters mod 4. Bröker–Hoffstein §5 show for exactly n = 6 over Q(ζ₆) that a careless choice of coset representatives (V′ = {x ≡ 1 mod 3}) spoils clean identities that hold for a good choice V. So the issue is real, but it is bookkeeping.

**Replacements, in increasing generality.**

(a) **Rank one: [R-], already done.** Kubota 1969 needs only μ_n ⊂ F, F totally imaginary [V]. Together with Selberg/EGM theory on Γ(N₀) it gives the vector-valued functional equation of D(w;m) (`FE_LoA6_proof.md` Thm A). Independently, the Eckhardt–Patterson S-integer functional equation is applied by Bröker–Hoffstein to n = 6, K = Q(ζ₆) (BH Thm 2.4, (3.1)) [V in BH].

(b) **Rank two via rank one (Bochner): [R-], not written.** CG's proof of Thm 6.1 uses only:
- the rank-one functional equation of each inner sum (their Prop. 5.7 / Thm 5.8);
- local p-part identities at good primes;
- Bochner's tube theorem.

The local identities at p ∤ 6 are local (CG §3–4; also Sawin–Whitehead's local parts, "the same shape regardless of the choice of global field") up to the ray signs above. Replace the Ψ-space by Kubota's congruence-class vectors mod N₀ (as in `FE_LoA6_proof.md` §1.2) and the rank-one input by (a). Then the CG argument gives continuation and the W(A₂) functional equations of the twisted series, as a vector over classes.

(c) **Eisenstein series: [R-]/[G], write-up needed.**
- The 6-fold metaplectic cover of SL₃ over a totally imaginary field splits over a congruence subgroup (Bass–Milnor–Serre / Moore).
- So the Borel Eisenstein series is an Eisenstein series for the arithmetic group ker κ₃ with a finite-order character. Langlands' theory gives continuation and functional equations; Mœglin–Waldspurger also treat covers [S].
- BBFH III identify its Whittaker coefficients with the A₂ series for r = 2 [S]. Their computation is a Bruhat-cell computation that transposes to Γ(N₀).

(d) **Angle twist ᾱ².** This needs the rank-one Kubota series twisted by Grössencharacters, i.e. Eisenstein series with nontrivial archimedean type. Patterson (PLMS 1987) and Heath-Brown–Patterson use exactly such series for Gauss sums at primes [S].

(e) **Sawin–Whitehead (arXiv 2507.08662).** Axiomatic multiple Dirichlet series with mixed exponents, functional equations of Kubota and Dirichlet type, and Weyl groupoids [V]. Proved over F_q[T] only; used in §4 for the combinatorics.

(f) "Friedberg's thesis-type constructions": I could not identify a construction under this name that adds anything to (a)–(e). The Friedberg–Hoffstein–Lieman and Brubaker-thesis series are n = 3 over Q(ω), where μ₂ₙ = μ₆ ⊂ F, so they do not face the problem.

**Status of item 1.** The match is exact on squarefree coprime terms [R]. The analytic theory (continuation and the A₂ functional equations for n = 6 over Q(ω)) is [R-]: no conceptual obstacle, not in the literature.

---

## 2. Polar divisor and the residual functional equation

### 2.1 Coordinates and poles

In CG normalisation (Re s_i > 3/2, centre 1) put s₁ = w + 1/2 and s₂ = z + 1/2, with u = w − 1/2 and v = z − 1/2. The simple reflections are:
- σ₁: (s₁,s₂) ↦ (2−s₁, s₁+s₂−1), i.e. (u,v) ↦ (−u, u+v);
- σ₂: (s₁,s₂) ↦ (s₁+s₂−1, 2−s₂), i.e. (u,v) ↦ (u+v, −v).

The functional equation is (CG Thm 6.1 (6.2) [V])

  Z*(s; m) = |m_i|^{1−s_i} Z*(σ_i s; m, σ_iΨ).

The polar set lies in the W-translates of the root hyperplanes n⟨α,·⟩ = 1:
- H₁ = {u = 1/6}: the theta pole w = 2/3;
- H₂ = {v = 1/6};
- H₁₂ = {u + v = 1/6}.

**The theta residue [R given rank one].** It is the residue on H₁. On H₁ the s-sum is the full Kubota series with numerator m, whose only pole in Re w > 1/2 is w = 2/3.

### 2.2 The functional equation that survives on H₁ [R- given §1.5(b)]

No reflection fixes H₁. The residue is transported as follows:
- σ₂ takes the point (1/6, v) of H₁ to (1/6+v, −v), which lies on H₁₂;
- σ₁ takes (1/6+v, −v) to (−1/6−v, 1/6), which lies on H₂;
- the diagram automorphism c₁ ↔ c₂ (transpose-inverse on GL₃) identifies residues on H₂ with residues on H₁, with ψ ↦ ψ̄ (§1.4).

The result is a single involution of the residual line, **z ↔ 5/6 − z, with centre z = 5/12**. The factors come only from the |m_i|^{1−s_i}:
- σ₂, applied at s₂ = z + 1/2, gives |μ³|^{1−s₂} = **N(μ)^{3(1/2−z)}**;
- σ₁, applied at s₁ = z + 2/3, gives |m₁|^{...} = 1.

Schematically, as a vector over the finitely many cusps and classes,

  R(z;μ) = N(μ)^{3(1/2−z)} G(z) Σ_s ρ^{(•)}(sμ³) γ₁(s) ψ̃(s) N(s)^{−(5/6−z)},   (2.1)

where:
- G is a fixed Γ- and scattering factor, independent of μ up to unimodular class phases;
- ρ^{(•)} are theta coefficients at the dual cusps;
- ψ̃ is the transported angle and ray data.

The residue on H₂ is Res_{s₂} Σ_{c₂} g(c₁μ³, c₂)… = ρ(c₁μ³). This is where ρ(sμ³) comes from.

**Size check [R as algebra].** Take |ρ(n)|² ≈ N(n)^{-1/6} on average. Then the dual side of (2.1) with weight W at length Q gives F′(μ) ≈ Q^{1/2} Σ_{N(s) ≲ L³/Q} ρ(sμ³)γ₁(s)N(s)^{-1/3}. Its random-size is Q^{5/12}, matching the direct side. The dual is consistent.

### 2.3 The dual coefficients do not thin out [R given KP]

Kazhdan–Patterson Hecke relations, in Bröker–Hoffstein's normalisation (bh.txt §1 [V]; τ = N^{1/12}ρ, which is the paper's normalisation, `V_gaps.md` referee item 5). For (m,p) = 1:
- τ(mp^j) = G_{j+1}(m,p) τ(mp^{4−j}) for 0 ≤ j ≤ 4;
- τ(mp⁵) = 0;
- τ(mp⁶) = N(p)^{1/2}τ(m).

For j = 3, **τ(mp³) = G₄(m,p)τ(mp)**. With G₄(s,μ) = χ̄_μ(s)⁴·G₄(1,μ), this gives

  ρ(sμ³) = N(μ)^{-1/6} G₄(1,μ) χ̄_μ(s)⁴ ρ(sμ)   ((s,μ) = 1).

G₂G₄ ≡ 1, so the relations force no vanishing. So the dual of (2.1) is a sum over all s coprime to μ, of length N(s) ≲ N(μ)³/Q, with coefficients of full size. Its terms are ρ(sμ)·(cubic character of s mod μ)·γ₁(s): theta coefficients along the sublattice μO.

**Conclusion [R- given §1.5(b) and KP]: the analytic conductor in μ is N(μ)³.**

### 2.4 The table over twist exponents

The dual coefficients are ρ(sμ^{6−j}). Orientation conventions can swap j ↔ 6−j except for j = 3.

| twist χ_m(μ)^j | WMDS parameter | KP at the dual | effective conductor |
|---|---|---|---|
| j = 1 (sextic) | μ⁵ | τ(mp⁵) = 0, so the dual is supported on μ \| s; then τ(s′μ⁶) = N(μ)^{1/2}τ(s′) | N(μ)⁴ |
| j = 2 (cubic) | μ⁴ | τ(mp⁴) = G₁τ(m) | N(μ)⁴ |
| **j = 3 (quadratic)** | **μ³** | **τ(mp³) = G₄τ(mp) ≠ 0** | **N(μ)³** |
| j = 4 (cubic) | μ² | τ(mp²) = G₃τ(mp²): a support condition only | N(μ)² |
| j = 5 (sextic) | μ | τ(mp) | N(μ) |

The pairing's twist is j = 3, fixed by the reflection [R].

### 2.5 Why the two "N(μ)²" heuristics are wrong

1. **"Two reflections ⇒ N(μ)²"** (`V_gaps.md` referee 4.2). The two reflections do not contribute symmetrically.
   - σ₁ is the s-Kubota functional equation, with numerator m₁·c₂ = c₂. Its conductor N(c₂) is absorbed in the variable shift s₂ ↦ s₁ + s₂ − 1. It carries no μ.
   - All of μ enters through σ₂, as |m₂| = N(μ)³.
2. **"GL₂×GL₁ Rankin–Selberg ⇒ N(μ)²"** (the E_P(θ⊗χ_μ) picture). On the cover the Gindikin–Karpelevich factors are L(n·⟨α,s⟩, (χ_i/χ_j)ⁿ). Since χ_μ⁶ = 1, the intertwining operators see no ramification. χ_μ only changes the Whittaker index, to (1, μ³).

So the correct picture: R(z;μ) is the first-index-(1, μ³) Whittaker coefficient of the GL₃^{(6)} Eisenstein series E_P(θ₆ ⊗ ᾱ², s), viewed as a function of s. This is the referee's "GL(3)^{(6)} induced from θ₆" picture, now with the conductor derived rather than guessed [R-].

### 2.6 Consequence for (V′)

**Lengths.**
- Dual length: N(μ)³/Q ≍ L³/Q = Z^{3(2M−1)−M}. At M = 0.826 this is Z^{1.13}, against Q = Z^{0.826}.
- The dual is shorter than Q only for L < Q^{2/3}, i.e. M < 3/4.

**Even there it does not help.** The dual coefficients ρ(sμ³) couple s and μ through a theta coefficient at the product. There is no character kernel for a large sieve.

**Bound from the dual.**
- The best bound I can get from the dual is Cauchy–Schwarz in s plus a sublattice mean square Σ_{N(s)≤X}|τ(sμ)|² ≪ K·X.
- The uniform K = 1 is not known. The trivial bound is K = L.
- This gives T′ ≪ Q^{-1/12} L² K^{1/2}, against Cauchy–Schwarz Q^{11/12}.
- The ratio exponent is 3M − 2 (K = 1) or 4M − 5/2 (K = L).
- So there is a saving only for M < 2/3 or M < 5/8.

**LP (`lp_dual_route.py`, theta piece g = max(0, b/6 − 2s) with −2s = 6M′−4 or 8M′−5, P²/Y dropped):**

| model | θ | geometry |
|---|---|---|
| Cauchy–Schwarz (paper) | 0.874957 | (0.355, 0.478, ℓ = 0.167) |
| dual route, K = 1 | 0.934113 | (0.315, 0.366, 0) |
| dual route, K = L | 0.960923 | (0.295, 0.348, 0) |

By distributivity θ_true = min of the rows = **0.874957**. The functional-equation route gives exactly Liu [R as computation, given the stated bounds].

---

## 3. The cubic cross-check and Bröker–Hoffstein

### 3.1 The cubic A₂ residue is not constant [R]

Twisted A₂ for n = 3 has coefficient H(c₁,c₂) = g(c₂,c₁)g(1,c₂) on coprime squarefree terms. Its c₁-residue is τ₃(c₂) with numerator c₂. By Patterson (BH: τ₃(c) = G₁(1,c), same symbol [V]):

  Res = Σ_c τ₃(c) G₁(1,c) N(c)^{-z} = Σ_c G₁(1,c)² N(c)^{-z} = Σ_c Ḡ₁(1,c)·μ(c)α(c)·(ray) N(c)^{-z},

using G₁³ = μα·(ray) on squarefree c. This is a Möbius- and Hecke-twisted conjugate cubic Kubota series, not a constant.

The identically constant combination τ₃(d)·Ḡ₁(1,d) ≡ 1 is the *opposite* orientation.

**Which orientation does the sextic row have? [R]** By Lemma Z and the signal identity (γ₁ = Gγ₂² on squarefree m),

  a(m) = ρ(m)·γ₁(m)·ᾱ²·(ray)  and  a(m)/(ρ(m)γ̄₁(m)) = γ₂(m)·μ(m)ᾱ(m)·(ray),

which oscillates. So the sextic row has the A₂ orientation, i.e. the non-constant cubic analogue. The "constant" cubic analogy does not transfer.

### 3.2 Why the cubic conductor is small and the sextic one is not [R given KP]

For n = 3 a twist of order 3 is a parameter μ or μ².
- For μ²: KP gives τ₃(mp²) = 0. The dual is supported on μ | s, and τ₃(s′μ³) = N(μ)^{1/2}τ₃(s′). The effective conductor drops from N(μ)² to **N(μ)**.
- For μ: τ₃(sμ) = G₂τ₃(s), again conductor N(μ).

So the cubic residue is "degree one", consistent with τ₃ being a Gauss sum. For n = 6 and j = 3 the corresponding relation τ₆(mp³) = G₄τ₆(mp) relates two *undetermined* coefficients and forces nothing. The degree stays 3.

### 3.3 Is ρ(m)γ₁(m) a simpler object? No evidence for it

What is proved (KP, via BH §1 [V]):
- τ₆(p⁵) = 0;
- τ₆(p⁴) = G₁(1,p);
- τ₆(p) = G₂(1,p)τ₆(p³);
- τ₆(mp²) = G₃(m,p)τ₆(mp²).

The squarefree coefficients τ₆(m) are undetermined (non-uniqueness of Whittaker models).

Conjectures:
- **CFH**, quoted in BH: Σ τ₆(m²)N^{-u} = (Σ τ₃(m)N^{-u})(Σ G₁^{(3)}(1,d)N^{-u}). This is a "square-root" identity, but only at square indices. It is irrelevant to squarefree rows.
- **BH Conj. 5.7**, numerical: (τ₆(π)/τ₆(1))² = ζ₆^{k(π)}·g₃(1,ε²,π)/|π|^{1/2}·3^{l(π)}·(1∓3ζ₆)^{m}(…)^{n}. The elements in the last factors have norm 7 or 13. BH "have not found a pattern in the exponents k(π), nor a condition when the elements of norm 7, 13 appear".

**Transport to the paper's convention [H].** One symbol defines both the Kubota series and G₂. BH's class set V differs from the paper's primary classes, which changes ρ by S-adic root-of-unity weights. With that caveat,

  a(π)² = ζ₆^{k(π)}·(junk)·ᾱ(π)²·γ̄₂(π)·(ray).

So at primes a(π) behaves like a patternless square root of a Hecke-twisted *conjugate cubic Gauss sum*.

**Not a GL₂ object.** a(m) is not a Fourier coefficient of a GL₂ form: it is not Hecke-multiplicative at squarefree arguments. It is a Whittaker coefficient of a GL₃^{(6)} Eisenstein series taken as a function of the spectral parameter (§2.5). I see no simpler object [H].

---

## 4. The Weyl-groupoid view (Sawin–Whitehead, arXiv 2507.08662 [V])

### 4.1 Dictionary [V]

For n even and a character χ of order n, an axiomatic multiple Dirichlet series with symmetric matrix M has heuristic coefficient Π(f_i′/f_i)^{M_ii} Π_{i<j}(f_i/f_j)^{M_ij}. So:
- node i carries the Gauss sum with character ξχ^{M_ii} = χ^{M_ii+n/2};
- M_ii = 0 is a Dirichlet-type node (quadratic Gauss sum, i.e. "ray");
- **M_ii = n/2 is a Möbius node**: the Gauss sum with trivial character, Ramanujan c_f(1) = μ(f). Compare Pellet's formula.

Heckenberger's bicharacter has q_ii = −ζ^{M_ii} and q̃_ij = ζ^{M_ij} (§5.1 of the paper). The reflection uses m_ij = min{m : q_ii^{m+1} = 1 or q_ii^m q̃_ij = 1}. Admissibility requires q_ii ≠ 1, so Möbius nodes are not admissible.

Kubota-type functional equations exist at i when all q̃_ij are powers of q_ii. Then the numerator is F = Π f_j^{n_ij} with n_ij(M_ii + n/2) ≡ −M_ij (their Thm 1.3).

**Our nodes (n = 6).**
- s: g̃(s) with χ¹, so M = 4.
- m (reflected): γ₁γ₃, so M = (1−3)+(3−3) ≡ 4. Products of Gauss sums at one node add exponents; this is the function-field shadow of Davenport–Hasse / Stickelberger.
- m (unreflected): no Gauss sum, so M = 0.
- Column A (unreflected): c_A ∝ γ₂, so M = 5.
- Column μ (reflected): cubic theta ∝ γ₂ at squarefree μ, so M = 5.

**Our edges.**
- s–m: (m/s)^{-1}, so 5.
- m–μ: quadratic, so 3.
- m–A: χ_A(m)^{+1}, so 1. The coherent j = 1 orientation of `V_derivation.md` §3 is the one with opposite signs on m–A and m–s.

For the faithful four-variable object, add the cubic Kubota node t (M = 5) with edge t–A = 4. Then Res_t gives τ₃(A) ∝ γ₂(A), and A becomes a Dirichlet node (M = 0).

### 4.2 Results (`weyl_groupoid.py`) [R as finite combinatorics; orientation conventions checked as above]

**Sanity checks.**
- A₂ sextic WMDS: 6 chambers, 1 object.
- A₃: 24 chambers.
- Friedberg–Hoffstein–Lieman M = (0 1; 1 0): 3 objects, matching Sawin–Whitehead's Figure 1.

**The reflected pairing R = (s,m,μ), M = [[4,5,0],[5,4,3],[0,3,5]], and the unreflected U = (s,m,A), M = [[4,5,0],[5,0,1],[0,1,5]], lie in the same groupoid**, which has 4 objects. The reflection at the column node maps U to R with basis (e_s, e_m + 2e_A, −e_A), and it changes M_mm: 0 → 4. In words, the row node acquires exactly the Gauss-sum character of γ₁γ₃. **This is Lemma Z (F1), recovered from the Heckenberger formalism** without the paper's local Fourier computation.

**The Kubota reflection at m in R has n_{m s} = 1 and n_{m μ} = 3.** The numerator is s·μ³ and the conductor N(s)N(μ)³: an independent derivation of §1.3 and §2.

**The Dirichlet reflection at m in U** (Poisson in m: the paper's Poisson identity and high side) gives [[0,1,0],[1,0,5],[0,5,3]]. **The column node becomes a Möbius node** (M = 3). This is precisely the Davenport–Hasse collapse c_A/τ(χ̄_A) = μ(A)·(ray), with exponents −1 − 2 ≡ 3. The zero-detecting signal sits on a non-admissible node, where no functional equation exists.

**A control.** With the conjugate column orientation (M_AA = 1) the three-node series is an honest arithmetic root system: 3 objects, no Möbius node. So the degeneracy of the three-node model is caused exactly by the Davenport–Hasse alignment that makes the signal.

**The faithful four-node object** (t,A,m,s), M = [[5,4,0,0],[4,0,1,0],[0,1,0,5],[0,0,5,4]], is an honest rank-4 arithmetic root system: 336 chambers, 7 objects, every node in every object of Kubota or Dirichlet type. So Sawin–Whitehead's rationality theorem applies over F_q[T] [R in the function-field model]. The Möbius node of the three-node model appears only after the cubic-theta residue: Patterson's γ₂ times the sextic Gauss sum produced by the m-reflection.

### 4.3 What the groupoid says about (V′) [H]

The full object is "of finite type", so T′ is a double residue (t at the cubic theta, s at the sextic theta) of a series with a finite groupoid of functional equations. In the number-field case this is [G]: Sawin–Whitehead prove it over F_q[T] only, and the reflection at the column node is a theta–Voronoi step that is available here from the paper.

The images of T′ that I can identify all fail to shorten the form in *both* variables:

| form | lengths |
|---|---|
| unreflected | (A ≍ Z, m ≍ Q) |
| reflected | (μ ≍ L, m ≍ Q) |
| m-Kubota | (μ ≍ L, s ≍ L³/Q) |
| m-Poisson | the signal, plus a dual of length Z²/Q in Ramanujan sums (`V_gaps.md` §1.6) |

In each form at least one variable is below the square root of its conductor:
- Q < (L³)^{1/2};
- L < (Q²)^{1/2}.

Multivariable convexity (Phragmén–Lindelöf in tubes) interpolates trivial bounds on these images. Every trivial bound is worse than Cauchy–Schwarz, by L^{1/2} or Z^{1/2}. So the multiple-Dirichlet-series route gives nothing better than CS.

A systematic search of all 7 objects and their residue loci for a shortening image is not done [G]. I expect it to fail for the same reason, but have not checked.

**Structural caution [H].** In every functional-equation form, T′ contains the zero-detecting signal (§4.2). A multiple-Dirichlet-series asymptotic would carry it as a 1/L-type normalising factor. Any bound obtained by isolating that term makes the low-side/high-side comparison circular. A usable (V′) must therefore be a *mean-value* bound for T′ as a whole, like the paper's Cauchy–Schwarz, not an asymptotic.

---

## 5. Numerics: the true row sequence against quadratic characters [N]

`check_twisted_meansquare.py`:
- Q = 3000, 2526 rows, 2388 of them squarefree.
- ρ is estimated by the cross-scale estimator at scales Y ≥ 4096, as in `exp_sextic_jtable.py`.
- γ₂(m) is computed by twisted multiplicativity with unit correction; ||γ₂|−1| ≤ 1.1·10⁻¹⁵.
- μ runs over the 169 squarefree primary μ with N(μ) ≤ 600.

The statistic is R(μ) = |Σ Ω a ψ_μ|² / Σ Ω²|a|². A random sequence has E R = 1. Saturating the large sieve would require R ≈ Q/#μ ≈ 18 on many μ.

| sequence, twist | N(μ) ∈ [1,60) | [60,200) | [200,600] |
|---|---|---|---|
| **true a = ρ·γ̄₂⁴, quadratic (j=3)** | 1.04 (0.67) | 1.38 (0.50) | 0.83 (0.19) |
| ζ-free ρ, quadratic (j=3) | 2.34 (1.03) | 0.53 (0.23) | 0.59 (0.18) |
| true a, sextic (j=1) | 1.92 (0.83) | 0.64 (0.32) | 0.97 (0.21) |
| control: ρ, coherent sextic (j=1) | **13.3 (5.0)** | **6.3 (1.0)** | **4.0 (0.4)** |

**Overall.** For the true sequence, Σ_μ R/#μ = 0.99.

**Tails.** The per-μ tails (max 8.7, 9% above 5) are comparable to the known-incoherent ζ-free baseline (max 10.7, 4% above 5). Both reflect estimator noise of about 1–3 units per μ (`V_derivation.md` §8.3).

**Reading.**
1. F′(μ) has no main term for any μ, as §1.4 predicts.
2. Σ_{μ≲L}|F′(μ)|² ≈ #μ·‖a‖². That is the "ideal" large sieve, not the (Q+L) worst case.

If this persists asymptotically, (V′) holds with the ideal large-sieve saving s₀ = (1−M′)/2 (≈ 0.087 at M = 0.826). The ideal-large-sieve rows of `REPORT.md` / `IDEA1_NOTES.md` then give θ ≈ 0.860–0.863, depending on the variant. Q = 3000 cannot test exponents. But nothing suggests (V′) is false.

**Independent cross-check.** The parallel report `A2_numerics.md` uses the exact ζ_m (checked against the product formula to 10⁻¹⁴), 12-scale data, Q = 1500/3000/6000, and class-by-class tests. It reaches the same conclusions:
- the quadratic twists of ρ·ζ are at the random level for N(ν) ∈ [7, 4Q];
- there is no coherent or Hecke-like structure, which also excludes the cubic-type degeneration;
- the conductor exponent cannot be resolved numerically at Q ≤ 6000.

The exponent N(μ)³ therefore rests on §2, not on numerics.

---

## 6. Conclusion: is (V′) provable, and what exactly is needed

**Sufficient statement (V′_LS) [G].** For squarefree μ with N(μ) ≍ L and Q^{1/2} ≤ L ≤ Q, uniformly over the ray and angle data of §1.4 and the frozen-k blocks of `V_gaps.md` §4:

  Σ_{N(μ)≍L} |Σ_m ρ(m)ζ_m (λ⁴μ/m)₂ W(N(m)/Q)|² ≪ Q^{1−2s₁+ε} Σ_m |ρ(m)|² W,  some s₁ > 0.

By Cauchy–Schwarz in μ this gives (V′) with s₀ = s₁ in Q-units, i.e. s₀ = M·s₁ in Z-units, for *every* column vector β. The LP needs only s₀ > 0 at L = Q^{2−1/M} ≈ Q^{0.79}.

**Why it is beyond present technology [H, with the reasons proved above].**
1. *The functional equation is useless in this range* [R-]. The conductor is N(μ)³ (§2). The dual is longer than Q for L > Q^{2/3}. For Q^{1/2} < L < Q^{2/3} the dual has no bilinear character structure (§2.6). The LP confirms there is no gain.
2. *No Hecke multiplicativity.* a(m) is a metaplectic Whittaker coefficient in the spectral variable. GL₃ tools (Voronoi with Kloosterman sums, the δ-method, GL₃ Kuznetsov) presuppose Hecke eigenvalues. Even for a GL₃ cusp form F, a power saving in Σ_{d≍L}|Σ_{n≍N}λ_F(n)χ_d(n)|² over (L+N)Σ|λ|² for N^{1/3} < L < N is not known to me. This is the "beyond the edge" regime (family L, conductor L³), like the energy wall W(w) of `REPORT.md` §2.
3. *No simpler object* (§3). The cubic collapse to degree one has no sextic analogue at j = 3, and the BH data show patternless phases at primes.
4. *Multiple-Dirichlet-series methods do not beat Cauchy–Schwarz here* (§4.3). Every functional-equation image of T′ that I identified is short in at least one variable relative to its conductor, and the signal sits on the Möbius node.

**Possible directions (long shots, none established) [H].**
- (a) A large sieve over the μ-family that exploits the structure of the *columns*. β_μ are cubic theta coefficients, so T′ is a double residue of the rank-4 arithmetic root system of §4.2. A spectral theory for that object on the 6-fold cover of GL₃ does not exist.
- (b) Change the probe so that the theta component of A is ρ·ζ̄ (`V_gaps.md` §1.7(b)). This changes the high side.
- (c) Accept the F₁ structure and seek a different split of the Gram matrix, in which the theta part is paired in the *unreflected* form. But there the m-twist is coherent and the pairing *is* the signal (`V_gaps.md` §1.6).

**Bottom line.**
- The A₂ multiple Dirichlet series determines the analytic structure of z ↦ Σ_m ρ(m)ζ_m χ(m)N(m)^{-z} completely: one functional equation, conductor N(μ)³, dual length L³/Q.
- That rules out the Voronoi route, and with it every route through functional equations plus generic large sieves.
- (V′) is numerically plausible (random-like twisted mean squares), but proving it needs a new mean-value theorem of beyond-the-edge type.
- At the slot-free geometry the honest status remains **Liu's 0.874957**.

---

## Referee report

*Independent skeptical review, 9 Oct 2026. I read this file, `REPORT.md` (CORRECTION banner), `IDEA1_NOTES.md`, `wf/V_gaps.md` (Lemma Z and its referee report), `A2_numerics.md` §0–3, `bh.txt` §1, and Chinta–Gunnells, "Constructing WMDS" (arXiv 0803.0691, text extracted). I re-ran `lp_dual_route.py` and wrote two small checks in the scratchpad: `n2check.py` (n = 2 Gauss-sum p-parts) and `n6check.py` (sextic Gauss sums g(up³, p^k) at split p = 7, 13).*

### Verdict: not refuted

The main conclusion survives:
- the A₂ functional-equation route gives no power saving for (V′) at the slot-free geometry (θ stays at Liu's 0.874957);
- proving (V′) needs a new mean-value input of beyond-the-edge type.

The conclusion is in fact **more robust than the file argues**: it does not depend on the conductor exponent being exactly 3 (item 3). But one deduction is logically invalid (item 1), and one identification is mislabelled (item 2). Neither changes the verdict.

### What I verified

| check | result |
|---|---|
| CG (6.2): Z*(s;m,Ψ) = \|m_i\|^{1−s_i} Z*(σ_i s; m, σ_iΨ) | [V], arXiv 0803.0691 Thm 6.1 |
| CG hypothesis "F contains the 2n-th roots of unity" | [V], CG §1 and §2 (stated twice). BBFH III body not checked (only a search summary saying "n-th roots"); I agree with the file's [S]. |
| KP relations in BH §1: τ_n(mp^j) = G_{j+1}(m,p) τ_n(mp^{n−2−j}) for 0 ≤ j ≤ n−2; τ_n(mp^{n−1}) = 0; τ_n(mp^n) = N(p)^{1/2} τ_n(m) | [V], bh.txt |
| G₂(m,p) G₄(m,p) = 1 | [R]: G_j(m,p) = (m/p)^{−j} G_j(1,p) and G₄(1,p) = conj G₂(−1,p) |
| γ₁(m)(μ/m)₂ = g̃(m, μ³) for (m,μ) = 1 | [R], from g(am,c) = (a/c)^{−1} g(m,c) |
| Composition σ₁σ₂(s₁,s₂) = (3−s₁−s₂, s₁), total factor \|m₂\|^{1−s₂} \|m₁\|^{2−s₁−s₂} | [R] |
| On H₁ (s₁ = 7/6): z ↔ 5/6 − z with factor N(μ)^{3(1/2−z)} | [R] |
| The H₂ residue of the *same* series is Σ_{c₁} g(1,c₁) ρ(μ³c₁) | [R]. The diagram-automorphism step of §2.2 is therefore unnecessary, though harmless. |
| Orientation in §3.1 | [R]: BH give τ₃(d) = G₁(1,d) with the same g_n(m,d) convention as A_m(s), so ρ·γ₁ is the τ₃·G₁ = G₁² analogue, not the constant τ₃·Ḡ₁. Lemma Z gives a/(ργ̄₁) = (ray)·μᾱγ₂, as stated. |
| Bound T′ ≪ Q^{−1/12} L² K^{1/2} and ratio exponents 3M − 2, 4M − 5/2 | [R as algebra] |
| `lp_dual_route.py` | reproduces 0.874957 / 0.934113 / 0.960923 at the stated geometries |

### 1. "No forced vanishing ⇒ conductor N(μ)³" is not a valid inference (§2.3). Downgrade [R-] to [H].

Non-vanishing of the coprime dual coefficients does not exclude a conductor drop. A counterexample is rank one with n = 2 (`n2check.py`):
- D(s; μ³) has functional-equation factor |μ|^{3(1−s)}.
- Its coprime dual coefficients g(μ³, c) = (μ/c) g(1,c) do not vanish.
- Yet the coprime series has conductor |μ|.
- The excess |μ|² is absorbed by the p-part polynomial. Numerically, g(up³, p^k) = 1, 0, φ(p²), 0, −p³ for k = 0,…,4. After dividing by the ζ(2s−1) local factor this gives P*(s) = 1 + p^{2−2s}, which satisfies P*(s) = p^{2(1−s)} P*(2−s).

So the file's criterion detects only one of two reduction mechanisms.

**Supporting heuristic for k = 3 (mine).** The second mechanism, "the same coprime coefficients occur at a smaller index", gives no reduction here.
- On coprime parts, ρ(m) g(μ³,m) also equals:
  - the H₁ residue at index (μ⁴, μ²), via τ(mp⁴) = G₁ τ(m), with factor exponent 2 + 4 = 6;
  - the H₂ residue at (μ³, 1), with exponent 3.
- The minimum is 3. This heuristic reproduces both the n = 2 drop and the cubic drop.
- The j = 1 row shows that vanishing can push below the minimal index (4 < 5). For j = 3 neither mechanism applies.

So k = 3 is the best available heuristic, but it is [H], not [R- given KP].

### 2. §1.3, "finite Dirichlet-polynomial corrections at p | μ", is inaccurate [R-] → [G]

At p | μ, the naive index-(1, μ³) coefficients ρ(m) g(μ³, m) are nonzero only for v_p(m) ∈ {0, 4}. I checked |g(up³, p^k)| = 1, 0, 0, 0, p^{3.5}, 0 at p = 7 and 13 (`n6check.py`).

For m = m′p⁴, twisted multiplicativity gives χ_p(m′)², and ρ(m′p⁴) ∝ G₁(m′,p) gives χ_p(m′)^{−1}. So relative to the coprime series the m′-sum carries an extra **odd** power of χ_p(m′). It is a *different twist* (index (1, μ³p⁵) on coprime parts), not a polynomial multiple of the coprime series.

Consequences:
- F′(μ), which is automatically the coprime sum, is not (polynomial) × (WMDS residue).
- Transferring the functional equation to F′(μ) needs a vector-valued functional equation over the indices (1, μ³d^f), d | μ. That functional equation is not written down.

Inclusion–exclusion through such pieces lengthens the dual (N(μ)³N(d)^{≥1}/Q in the naive count) rather than shortening it. So this does not threaten the conclusion. It does mean the statement "the conductor of the residual series of F′(μ) is N(μ)³" is heuristic.

### 3. The conductor is not decisive, so the conclusion is robust [R as algebra]

Suppose the dual has length N(μ)^k/Q, with coefficients of the size the functional equation dictates (theta coefficients at an index ~ sμ^k, sublattice mean square K). The computation of §2.6 then gives:
- F′(μ) ≲ N(μ)^{k/2} Q^{−1/12} K^{1/2};
- T′ ≲ L^{(k+1)/2} Q^{−1/12} K^{1/2}, against Cauchy–Schwarz Q^{11/12}.

A saving needs L^{(k+1)/2} < Q, i.e.

  **k < 2M/(2M−1) − 1 ≈ 1.53 at M = 0.826** (K = 1), or k < 0.53 with K = L.

So **even a dual of length L²/Q (k = 2) gives no saving** unless the dual is a *bilinear form with a character kernel*, to which Heath-Brown or de Faveri apply. That was the situation in the ζ-free model, where the Voronoi dual was x_μ y_n (n/μ)₃ with y_n independent of μ.

For the true row sequence, the dual coefficients ρ(sμ^k) couple s and μ for every k ∈ {2, 3, 4}: KP leaves τ(mp), τ(mp²) and τ(mp³) undetermined. Only a degenerate "degree-one" structure (k ≤ 1) would help, and `A2_numerics.md` §3.1 excludes the cubic-type degeneration.

So the brief's premise ("dual length L²/Q ⇒ saving") holds only with bilinear structure. The file says this in §2.6 ("no character kernel") but should make it the headline: the verdict follows from the coupling and the convexity arithmetic, whether k is 2 or 3.

### 4. Overclaims and minor points (no effect on the verdict)

- **§4.2, "n_{mμ} = 3, an independent confirmation of N(μ)³".** This is not independent. It is the same Kubota-numerator computation (numerator μ³), so it confirms the functional-equation factor, not the effective conductor (item 1).
- **§2.2, size check.** It passes for *every* exponent k: the exponents give 5q/12 identically. It checks consistency of the normalisation, not that k = 3.
- **§1.4.** The ᾱ² Grössencharacter on c₂ (K-type Eisenstein series) is outside every cited framework, including CG's finite-order Ψ-spaces. The claim "R(z;μ) entire, no main term for any μ" is [R-] at best. The constant-term computation, L(·, ᾱ¹²) entire, is correct.
- **Numerics (§5, `check_twisted_meansquare.py`).**
  - The twist used is (m/μ)₆^j, with reciprocity and ray factors omitted. The docstring says this "does not change these statistics". That is true for random-level means, but it could hide coherence aligned with a ray factor (`A2_numerics.md`'s class-by-class tests cover this).
  - μ = 1 is excluded.
  - Per-μ noise is 1–3 units, so the [1, 60) bin (1.04 ± 0.67) cannot distinguish "structurally tiny" from "random". "No main term for any μ" rests on §1.4, not on the numerics.
  - The cross-scale estimator is unbiased only up to the O((Q/Y)(1 + P²/Y)) off-diagonal Gram residual of Lemma D1. That is immaterial at the stated precision.
  - As the file says, the conductor exponent cannot be resolved at Q ≤ 6000.
- **§3.3, transport of BH Conj. 5.7.** This is orientation-dependent. With the other orientation, a(π)² is (patternless root of unity) × (Hecke) × junk. Either way it gives no usable structure, and the file correctly labels it [H].
- **§1.5(c).** The Bass–Milnor–Serre splitting over a congruence subgroup for SL₃ over a totally imaginary field is correct as a reference point.

### Recommendation

Accept "no_saving" for the A₂ / functional-equation route, with these edits:
1. Relabel the conductor N(μ)³ as [H] (well supported), not [R- given KP]. Add the n = 2 caveat.
2. Replace "finite Dirichlet-polynomial corrections" in §1.3 by "other twisted coprime series (vector-valued functional equation needed) [G]".
3. State item 3 prominently: no saving unless k ≲ 1.5, or unless the dual has a μ-independent character kernel; the coupling ρ(sμ^k) prevents the latter.

The open problem is correctly identified as (V′_LS), a beyond-the-edge large sieve for a(m) = ρ(m)ζ_m against quadratic characters. The honest status at the slot-free geometry remains Liu's 0.874957.
