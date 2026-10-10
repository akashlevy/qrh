# (FE) and (LoA6) for the sextic Kubota series, and the clean Gram bound (D)

*9 Oct 2026. Companion to `gram_split/LoA6.md` and `D_derivation.md` §6.2. Both are superseded here where they disagree.*

**Status labels**
- **[R]** rigorous, given the cited theorem.
- **[R-]** routine but not written out line by line; the mechanism is given.
- **[H]** heuristic.
- **[G]** gap.

**Citation labels**
- **[V]** I read the statement in the source: TeX source, OCR text or postscript.
- **[S]** I know the statement only through another paper's restatement.

---

## 0. Verdict

**(FE), (LoA6) and Corollary D hold, with no heuristic step left.** What remains is three write-out tasks of known type (§7), none conceptual.

| claim | status |
|---|---|
| Kubota character on Γ(N₀) and its explicit form | [R] (Kubota 1966/1969 [V], local Hensel facts) |
| Continuation, vector FE, holomorphy on Re w = 1/2 | [R] given Selberg/EGM theory for cofinite Kleinian groups; κ-twist transfer [R-] |
| Structure of the scattering matrix (Lemma C) | [R-]. Kubota did the cusp pair (∞,∞) [V]; general pairs are the same computation |
| Only pole in Re w > 1/2 is the simple pole at w = 2/3 | [R] given Lemma C |
| FE of D_c, conductor ≍ N(m)(1+\|t\|)^{10}, uniform in m | [R] given the vector FE |
| Identification of D_c (ray twist) with Eisenstein coefficients | [R-] (conventions only) |
| Polynomial growth in vertical strips, uniform in m | [R] given finite order; finite order is [R-] |
| LoA6: Σ_{N(m)≤Q} \|D_c(1/2+ε+it;m)\|² ≪ Q^{1+ε}(1+\|t\|)^R, at every cusp | [R] given the rows above |
| Corollary D: ‖A − A^sig‖² ≪ (1+\|v\|)^R (Q/Y) Z^ε, no P_a^{1/6} or P_a²/Y term | [R] given LoA6 |

**Four corrections to the earlier files** (details in §6):
1. Brubaker–Bump 2006 Theorem 1 does **not** apply to n = 6 over Q(ω). It assumes −1 is an n-th power, i.e. μ₁₂ ⊂ F [V]. The correct base is Kubota's 1969 theory, which needs only μ_n ⊂ F with F totally imaginary [V], together with Selberg's general theory.
2. The constant-term normaliser in `LoA6.md` is miswritten. Moreover, at a general cusp pair the scattering entries involve Hecke L-ratios L(6s−6,ψ)/L(6s−5,ψ), not one ζ_K-ratio. The pole conclusion is unchanged.
3. The multiplicative sextic large sieve (BGL / de Faveri) is not needed. The paper's planar additive large sieve suffices, which confirms `D_derivation.md` §6.2, Addition 1.
4. The bootstrap for the theta residue (`D_derivation.md` §6.2, Addition 2) is unnecessary. A smoothing weight G_s that vanishes at the theta pole and at the ζ-pole removes both residue terms.

---

## 1. Objects and normalisations

### 1.1 Paper conventions (paper78 §"Arithmetic and analytic preliminaries", lines 557–720)

- K = Q(ω), O = Z[ω], λ = √−3, q_a = N(a) = |a|².
- e(z) = exp(2πi Tr(z/λ)) = ě(z/λ), where ě(z) := exp(2πi Tr z) = exp(2πi(z+z̄)).
- Primary means ≡ 1 mod 3.
- χ_p(a) = (a/p)₆ ≡ a^{(q_p−1)/6} mod p, extended to χ_c = Π_p χ_p^{v_p(c)}, with zero on non-units. A valuation divisible by 6 gives 1_{(a,c)=1}, not 1.
- Gauss sum: g(s,h) = Σ_{d mod s} χ_s(d) e(hd/s), and g̃ = g/√N(s).
- S is the paper's fixed excluded set. It contains the primes over 6, the primes of norm ≤ P₀, and the defining-modulus supports of η, of T̂, and of the fixed characters.

### 1.2 The series

The paper's coefficient is

  c(s) = 1_{s∈σ} χ_s(b_*) / (τ ξ(s)),  s primary, (s,S) = 1   (paper lines 3460–3466, 3606–3610).

- 1_{s∈σ} is a ray class of the fixed ray group T, whose modulus 𝔪_T is supported on S.
- s ↦ χ_s(b_*) is a fixed-numerator ray character with conductor supported on S (paper Lemma `fixed-numerator-ray`).
- ξ is a character mod b_*.

Hence c is periodic modulo a fixed 𝔣 supported on S.

**Choice of level.** Fix N₀ ∈ O with:
- 72 · 𝔪_T · b_* | N₀;
- every prime of S divides N₀;
- N₀ supported on S.

Let R(N₀) be the set of residue classes mod N₀ that are primary and prime to N₀. Then

  D_c(w;m) := Σ_{s primary, (s,S)=1} c(s) g̃(s,−m) N(s)^{−w} = Σ_{s₀∈R(N₀)} c(s₀) D_{s₀}(w;m),
  D_{s₀}(w;m) := Σ_{s ≡ s₀ (N₀)} g̃(s,−m) N(s)^{−w}.   (1.1)

Each s ≡ s₀ (N₀) is automatically primary and prime to S. This absorbs the ray-class twist completely: c is just a weight on finitely many classes.

### 1.3 Eisenstein series (Kubota 1969 §5 [V]; EGM normalisation)

- Γ := Γ(N₀) ⊂ SL₂(O) acts on H³ = {(z,v)}. It is normal in SL₂(O), and every cusp has width lattice Λ = N₀O. Its ě-dual is Λ* = (λN₀)^{−1}O.
- κ(γ) := (c/d)₆ for γ = (a b; c d) with c ≠ 0, and κ(γ) := 1 for c = 0.
- A cusp 𝔞 is *essential* if κ is trivial on its stabiliser Γ_𝔞. Fix representatives η_i of the essential cusps and scaling matrices σ_i ∈ SL₂(O).

For a unitary character ψ ∈ {κ, κ̄} define, for Re s > 2,

  E_i(w,s,ψ) = Σ_{γ∈Γ_i\Γ} ψ̄(γ) v(σ_i^{−1}γw)^s.

Its expansion at any cusp k (Kubota (36)–(42) [V]; DFDH §3.2 [V]) is

  E_i(σ_k w,s,ψ) = δ_{ik} v^s + v^{2−s} K(0,s) φ_{ik}(0,s,ψ) + Σ_{δ∈Λ*+h_k, δ≠0} v^{2−s} K(δv,s) φ_{ik}(δ,s,ψ) ě(δz),

with:
- K(u,s) = (2π)^s |u|^{s−1} K_{s−1}(4π|u|)/Γ(s) and K(0,s) = π/(s−1);
- φ_{ik}(δ,s,ψ) = V^{−1} Σ_{γ=(∗ ∗; c d) ∈ Γ_∞\σ_i^{−1}Γσ_k/Γ_∞, c≠0} ψ̄(σ_iγσ_k^{−1}) ě(δd/c) N(c)^{−s}.   (1.2)

The constant term is present only for essential k, where h_k = 0. The double cosets are in bijection with the pairs (c, d mod cΛ), c ≠ 0 (Kubota (36) [V]).

**Dictionary.** s = w + 1/2. The theta pole s = 1 + 1/6 = 7/6 corresponds to w = 2/3, and the unitary axis Re s = 1 to Re w = 1/2.

---

## 2. Literature: what is used, and verification status

| source | what it gives | status |
|---|---|---|
| **Kubota**, *On automorphic functions and the reciprocity law in a number field*, Kyoto Lectures in Math. 2 (1969). OCR text read from KURENAI, hdl 2433/84907 | §1: "F a totally imaginary number field … which contains the n-th roots of unity for a fixed n ≥ 2". No μ_{2n} hypothesis. §5: Γ = congruence group mod N with n² \| N; character χ (Prop. 2); essential cusps; E_i; Fourier expansion (38)–(42) with Dirichlet series in χ(c,d) exp(2πi tr(md/c)) N(c)^{−s}; FE Φ(s,χ)Φ(2−s,χ) = I and E = Φ(s)E(2−s) "from the general theory". §6: constant term at ∞ ∝ L(ns−n,λ)/L(ns−(n−1),λ) over class-group characters; pole at s = (n+1)/n; "the singularities of E(u,s,χ) in the half-plane Re s > 1 are singularities of the constant term"; Theorem 6: the residue Θ_χ ≠ 0 | [V] |
| **Kubota**, *Ein arithmetischer Satz über eine Matrizengruppe*, J. reine angew. Math. 222 (1966) 55–57 | κ = (c/d)_n is a character of Γ(n²) when μ_n ⊂ F | [S] (title and pages verified; content via Kubota 1969 Prop. 2 and many restatements) |
| **Elstrodt–Grunewald–Mennicke**, *Groups acting on hyperbolic space*, Springer 1998: §3.4 Thm 4.1; §6.1 Thm 1.11; §6.3 Lemma 3.5; and Cohen–Sarnak, *Notes on the trace formula*, Ch. 6 | For cofinite Γ ⊂ PSL₂(C) with trivial character: continuation; FE; Φ(s)Φ(2−s) = I; poles in Re s > 1 finitely many, real, simple, in (1,2], and those of E_i coincide with those of φ_{ii}(0,s); E and Φ holomorphic on Re s = 1, Φ unitary there; Φ bounded on Re s ≥ 1, \|Im s\| ≥ ε | [S] (restated with theorem numbers in DFDH §3.2.1, which I read) |
| **de Faveri–Dunn–Hoffstein**, arXiv:2607.07911v2, §§3–8 | The same classical framework with a unitary character over Q(i) (quartic); the FE of Fourier coefficients (their eq. `Eisenstein_functional_equation`); LoA for quartic Kubota series via Rankin–Selberg (Prop. "Second moment of ψ": Σ_{N(a)≤A} \|ψ(a,s;u)\|² ≪ \|s\|^{2000} A^{1+ε} for Re s ≥ 1+ε). They write "We expect similar arguments to work for all n ≥ 5", and note that a BGL-type large sieve "could also be applied … as in [DdFDS, Lemma 6.6]" | [V] |
| **David–de Faveri–Dunn–Stucky**, arXiv:2410.03048v2, Lemma 6.3 (Patterson's FE (6.11)), Lemma 6.6 | Cubic template. F(r,s) = N(r)^{1−s}[A(s)F_∞(r,2−s) + B(s)F(r,2−s)], with A, B explicit Dirichlet-polynomial ratios; the dual side has numerators r·νλ^b. LoA: Σ_{N(m)≤M} \|ψ(hm,σ+it)\|² ≪ N(h)^{1/2} M^{1+ε}(1+\|t\|)² for 1 < σ ≤ 5/4, via an AFE, the cubic large sieve and a sup-in-t bootstrap | [V] |
| **Brubaker–Bump**, *On Kubota's Dirichlet series*, J. reine angew. Math. 598 (2006) 159–184. Postscript (Wayback copy of sporadic.stanford.edu/bump/kubota.ps), text extracted | Theorem 1 (S-integer FE with a finite scattering matrix P_η(s)). Hypothesis, §1: F contains μ_n, "and moreover" −1 is an n-th power; §2: "Moreover, we will also impose the condition that −1 is an n-th power. If n is even, this means that μ_{2n} ⊂ F." **Fails for n = 6, F = Q(ω).** | [V] (garbled OCR, but these sentences are legible) |
| **Diaconu–Ion–Paşol–Popa**, arXiv:2607.27131, Thm 3.2 (= [BB06, Thm 1]) | Modern restatement: D̂^S(s,𝔞,ψ) has at most poles at s = 1/2 ± 1/r; FE with ∏_{v∈S_f}(1 − q_v^{rs−r/2−1})^{−1} Σ_η P_η(s)…; Γ-factor (2π)^{−(r−1)(s−1/2)} Γ(rs−r/2)/Γ(s−1/2); "(s−1/2−1/r)(s−1/2+1/r) D̃^S is entire of order 1". Also needs μ_{2r} ⊂ F | [V] |
| **Eckhardt–Patterson**, Proc. LMS (3) 64 (1992) 225–264, §5 | General-n FE in S-integer form. Bröker–Hoffstein apply it to K = Q(ζ_n), n > 2, specifically n = 6: their Thm 2.4 (continuation, at most one pole at s_B = 1/n in Re s_B ≥ 0) and FE (3.1), with Γ(ns)/Γ(s), G_f(s) = ∏_{v∈S_f}(1−\|π_v\|^{1−ns})^{−1}, T_{ij}(r,s) ∈ Q̄(\|π_v\|^s) | [S] (via Bröker–Hoffstein, Math. Comp. 85 (2016), arXiv:1312.0568, §§2–3, which I read) |
| **Müller**, *The trace class conjecture…*, Ann. of Math. 130 (1989) §4 | Finite order of rank-one (cuspidal-type) Eisenstein series | [S] (via Bernstein–Lapid, arXiv:1911.02342, intro) |
| paper78, Lemma `planar-additive-sieve` (line 3498) | Σ_{q_m≤Q} \|Σ a_j e(m z_j)\|² ≪ (Q + δ^{−2}) Σ\|a_j\|² for δ-separated z_j ∈ C/O | [V] (proof in the paper) |

**Consequence.** The basis for n = 6 over Q(ω) is Kubota 1969 plus the general theory, i.e. the classical Γ(N₀) route below. The S-integer theorems (BB06, DIPP) do not literally apply. EP92 (via BH) gives an independent statement whose shape matches ours (§3.9).

---

## 3. Theorem A: continuation, poles, functional equation, growth

**Theorem A.** Let N₀ be as in §1.2, i an essential cusp, k any cusp, and δ ∈ Λ* + h_k with δ ≠ 0.

1. φ_{ik}(δ,s,ψ) continues meromorphically to C.
2. In Re s > 1 its only possible pole is a simple pole at s = 7/6.
3. It is holomorphic on Re s = 1.
4. It satisfies the vector FE (3.4) with a scattering matrix Φ(s) independent of δ.
5. It has polynomial growth (3.8) in vertical strips, uniformly in δ.

Consequently D_c(w;m) has:
- meromorphic continuation to C;
- in Re w > 1/2, only the simple pole at w = 2/3, with residue ρ_c(m);
- holomorphy on Re w = 1/2;
- no poles in Re w > 1/2 − c₀/log(|Im w|+3);
- an FE w ↔ 1−w of conductor ≍ N(m)(1+|Im w|)^{10};
- polynomial growth on Re w ≥ 1/2, uniformly in m.

### 3.1 The Kubota character on Γ(N₀) [R]

**Lemma A1.** For γ = (a b; c d) ∈ Γ(N₀) with c ≠ 0, write c^S for the part of (c) prime to S. Then

  κ(γ) = (c/d)₆ = χ_{c^S}(d)   (up to the global orientation convention of the Hilbert symbol, fixed as in `V_derivation.md` T1).

*Proof.*
1. (c/d)₆ = Π_{v|d} (c,d)_v. By the product formula this equals Π_{v∤d} (c,d)_v^{−1}.
2. Archimedean places are complex, so their symbols are trivial.
3. At v ∈ S we have d ∈ 1 + N₀O_v ⊂ (K_v^×)⁶, so (c,d)_v = 1. This uses:
   - 1 + 8O₂ ⊂ squares, and 1 + 2O₂ ⊂ cubes since 3 ∈ O₂^×;
   - 1 + λO_λ ⊂ squares, and 1 + λ⁴O_λ ⊂ cubes, because (U^{(2)})³ = U^{(4)} as e/(p−1) = 1 for Q₃(√−3);
   - 1 + 𝔭O_𝔭 ⊂ sixth powers for 𝔭 ∤ 6, by Hensel.
4. At v ∉ S with v | c, the tame symbol is (c,d)_v = (d/𝔭_v)₆^{±v_v(c)}. The product of these is χ_{c^S}(d)^{±1}. ∎

κ is a homomorphism on Γ(N₀) because 36 | N₀ (Kubota 1966; Kubota 1969 Prop. 2 and §5 [V]).

### 3.2 General theory, via the kernel Γ₁ = ker κ [R given EGM; transfer R-]

Γ₁ := ker κ is normal of index 6 in Γ, hence cofinite. For an essential cusp η_i we have Γ_i ⊂ Γ₁, and Γ_i\Γ = ⊔_{δ∈Γ₁\Γ} Γ_i\Γ₁δ. Therefore

  E_i(w,s,ψ) = Σ_{δ∈Γ₁\Γ} ψ̄(δ) E^{Γ₁}_{δ^{−1}η_i}(w,s),   (3.1)

a fixed linear combination of trivial-character Eisenstein series of the cofinite Kleinian group Γ₁. The EGM/Selberg theory for Γ₁, restricted to the ψ-isotypic part (which the scattering operator preserves), gives:

- (a) meromorphic continuation of E_i(·,s,ψ), and hence of every φ_{ik}(δ,s,ψ), to C;
- (b) **E(w,s,ψ) = Φ(s)E(w,2−s,ψ)** and Φ(s)Φ(2−s) = I, where Φ(s) = (K(0,s)φ_{ij}(0,s,ψ))_{i,j essential};
- (c) E and Φ are holomorphic on Re s = 1, and Φ is unitary there;
- (d) in Re s > 1, the poles of E_i are real, simple, lie in (1,2], and coincide with the poles of φ_{ii}(0,s,ψ).

On (d): at a pole s₀ ∈ (1,2] the residue matrix of Φ is positive semidefinite, being a Gram matrix of the residual forms' constant terms. So a zero diagonal entry kills its row.

Kubota 1969 §5 [V] states (a) and (b) in exactly this setting, citing "the general theory". DFDH §3.2 [V] states (a)–(d) for a unitary character with the same proof skeleton. The only unwritten step is the routine isotypic transfer of the EGM statements to (3.1). [R-]

### 3.3 Lemma C: structure of the constant terms [R-]

Let Ψ be the finite group of ray class characters mod N₀ of the form θ⁶, with θ a character of (O/N₀)^× (θ⁶ is trivial on μ₆). Write L(u,ψ) for the Hecke L-function with the Euler factors at S removed.

**Lemma C.** For essential i, j,

  φ_{ij}(0,s,ψ) = Σ_{ψ'∈Ψ} P_{ij,ψ'}(s) · L(6s−6,ψ') / L(6s−5,ψ'),   (3.2)

where each P_{ij,ψ'} is a rational function of {q_v^{−s}}_{v∈S} such that:
- its denominator is a product of factors (1 − β q_v^{P(1−s)}) with |β| = 1 and P ∈ N;
- it is holomorphic in Re s > 1, given there by an absolutely convergent series;
- it is holomorphic and bounded uniformly in t on every line Re s = σ ≠ 1.

*Proof.*
1. **Multiplier.** For γ ∈ σ_i^{−1}Γσ_j with bottom row (c,d), the bottom row is ≡ (c_{ij}, d_{ij}) mod N₀. The generalisation of Lemma M (`V_derivation.md` §4) gives

     κ(σ_iγσ_j^{−1}) = χ_{c^S}(d)^{±1} · h_{ij}(c_S; d mod π^K; c^S mod N₀).

   Here c = c_S c^S, and h is unimodular. It depends on d only modulo a bounded S-power, because the Hilbert symbols at S see d only modulo (K_v^×)⁶ ⊃ 1 + π_v^{k_v}O_v. It depends on c^S only through ray classes mod N₀ (fixed-numerator symbols χ_{c^S}(c_S) and Hilbert symbols).
2. **The d-sum.** The double-coset parameter is d mod cN₀ in a fixed class mod N₀. By CRT over c_S N₀ and c^S, the sum over d splits. Σ_{d mod c^S} χ_{c^S}(d)^{±1} equals φ(c^S) when c^S is a sixth power (as an ideal) and 0 otherwise.
3. **The c^S-sum.** For c^S = (c₀)⁶, all sextic cross symbols are 1. The remaining condition is that the element c lie in its class mod N₀. That condition is a ray-class condition on c₀ through c₀ ↦ c₀⁶ mod N₀ (well defined, since η⁶ = 1 for units η). Detecting it with characters θ of (O/N₀)^× gives

     Σ_{c₀} θ⁶(c₀) φ(c₀⁶) N(c₀)^{−6s} = L(6s−6,θ⁶)/L(6s−5,θ⁶),

   using Σ φ(𝔞) N𝔞^{−u} = L(u−1)/L(u) and φ(c₀⁶) = N(c₀)⁵ φ(c₀).
4. **The c_S-sum.** c_S = u Π_{v∈S} π_v^{a_v}. For a_v ≥ a₀ the d-count is β · q_v^{a_v} with β depending only on the classes of a_v modulo a fixed period. Cross-place couplings, through π_w^{a_w} mod π_v^k and Hilbert symbols, are periodic in the a_w. So each generating function Σ_a (periodic)·q_v^{a(1−s)} is rational with denominator (1 − β q_v^{P(1−s)}). ∎

**Sanity check.** At (i,j) = (∞,∞), Kubota 1969 §6 [V] finds exactly L(ns−n,λ)/L(ns−(n−1),λ) over class-group characters λ. Q(ω) has class number 1, and with c ≡ 0 mod N₀ no class condition on c^S arises, so only ψ' = 1 survives there.

**Why Hecke L-functions appear.** The class conditions at general cusp pairs (c ≡ c_{ij} mod N₀ with c_{ij} a unit) produce the L(·,ψ') with ψ' ≠ 1. This is the sextic analogue of Huxley's L(2s−1,ψ)/L(2s,ψ) for Γ(N) ⊂ SL₂(Z). In the S-integer formulation (BB06, EP92) these terms are absorbed into the finitely many classes η ∈ K_S^×/(K_S^×)⁶O_S^×.

### 3.4 Poles: only w = 2/3 in Re w > 1/2, and the boundary Re w = 1/2 [R given 3.2, Lemma C]

- **Re s > 1.**
  - L(6s−6,ψ') is entire for ψ' ≠ 1, and has only the simple pole s = 7/6 for ψ' = 1.
  - L(6s−5,ψ') ≠ 0 for Re(6s−5) ≥ 1.
  - P_{ij,ψ'} is holomorphic there.

  So φ_{ii}(0,s) has in Re s > 1 at most the simple pole s = 7/6. By 3.2(d) the same holds for every E_i, and hence for every φ_{ik}(δ,s), since each Fourier coefficient is an integral of E over a compact set.
- **Re s = 1.** Holomorphic, by 3.2(c). In particular s = 1 (w = 1/2) is not a pole.
- **Left of the axis.** By (3.5) below, φ = Z/𝓛. In 5/6 < Re s < 1 the only possible poles are at s = (5+ρ)/6, i.e. w = (2+ρ)/6, for zeros ρ of the finitely many L(·,ψ'). The classical zero-free region excludes Re s > 1 − c₀/log(|t|+3).

Hence (task item 1): D_c(w;m) has, in Re w > 1/2, only the simple pole at w = 2/3, and is holomorphic on Re w = 1/2. No GRH is used.

### 3.5 FE of the Fourier coefficients [R given 3.2(b)]

Take the δ-th coefficient at cusp k of E_i(σ_k w,s) = Σ_j Φ_{ij}(s) E_j(σ_k w,2−s). Since K_{s−1} = K_{1−s} and the v-powers match,

  **φ_{ik}(δ,s) = (2π)^{2−2s} N(δ)^{1−s} (Γ(s)/Γ(2−s)) Σ_j Φ_{ij}(s) φ_{jk}(δ,2−s),  Φ_{ij}(s) = π(s−1)^{−1} φ_{ij}(0,s).**   (3.3)/(3.4)

Three properties:
- **Uniform in δ.** Φ does not depend on δ at all. The only δ-dependence is N(δ)^{1−s}, so the conductor is N(δ) ≍ N(m).
- **Dual side.** The same frequency δ, the same expansion cusp k, and the Eisenstein series at every essential cusp j. In Gauss-sum language, these are the Kubota series with numerator m twisted by fixed S-classes: the "Gauss sums at the other cusps" of `LoA6.md`.
- **Archimedean part.** With the ζ/L-functional equations it is ≍ (1+|t|)^{10(1−σ)} on Re s = σ (checked numerically, §3.9). So the analytic conductor is N(m)(1+|t|)^{10}.

### 3.6 The completed coefficient Z and its FE matrix M [R given Lemma C]

Put 𝓛(s) := Π_{ψ'∈Ψ} L(6s−5,ψ') = Σ_𝔢 λ(𝔢) N(𝔢)^{5−6s}, where |λ(𝔢)| ≤ τ_{|Ψ|}(𝔢), and

  Z_{ik}(δ,s) := 𝓛(s) φ_{ik}(δ,s).   (3.5)

Then

  Z_{ik}(δ,s) = (2π)^{2−2s} N(δ)^{1−s} (Γ(s)/Γ(2−s)) Σ_j M_{ij}(s) Z_{jk}(δ,2−s),   M_{ij}(s) := Φ_{ij}(s) 𝓛(s)/𝓛(2−s).   (3.6)

Written out,

  M_{ij}(s) = π(s−1)^{−1} Σ_{ψ'} P_{ij,ψ'}(s) L(6s−6,ψ') Π_{ψ''≠ψ'} L(6s−5,ψ'') / Π_{ψ'''} L(7−6s,ψ''').

**Properties.**
- **(M1) Holomorphy.** M is holomorphic in Re s < 1 and at s = 1.
  - In Re s < 1: 𝓛(s) cancels every 1/L(6s−5,·); L(7−6s,·) ≠ 0 when Re(7−6s) > 1; and P_{ij,ψ'} has poles only on Re s = 1.
  - At s = 1: Φ is holomorphic there by 3.2(c), and 𝓛(s)/𝓛(2−s) is a ratio of simple poles.
- **(M2) Trivial zeros.** Every term of M carries L(6s−6,ψ'). For K = Q(ω), L(−k,ψ') = 0 for every k ≥ 1 and every finite-order ψ', because Γ_C(u) has poles at u ∈ −N₀ while Λ(u,ψ') is finite there. (Numerically, L(−1,χ_{−3}) = 0 to 30 digits, and ζ_K(−1) = 0.)

  Hence M vanishes at s = 5/6, where the dual coefficient Z_{jk}(δ,2−s) has its theta pole, and at s = 0, −1, …, where Γ(s) has poles.
- **(M3) Bounds.** On each line Re s = σ < 1, M ≪_σ (1+|t|)^{K(σ)}. This uses Hecke FE/convexity for L(6s−6,·) and L(6s−5,·), boundedness of 1/L(7−6s,·) for Re(7−6s) > 1, and Lemma C for P.
- **(Z1) Holomorphy of Z.** Z_{ik}(δ,·) is holomorphic in −1/4 ≤ Re s ≤ 3 except at s = 1 and s = 7/6.
  - For Re s ≥ 1, use 3.4 and holomorphy of 𝓛 away from s = 1.
  - For Re s < 1, use (3.6), (M1), (M2), and the fact that the dual coefficient Z(δ,2−s) is holomorphic in 1 < Re(2−s) ≤ 9/4 except at 2−s = 7/6.
- **(Z2) Absolute convergence.** For Re s > 2, Z_{jk}(δ,s) = V^{−1} Σ_𝔢 Σ_{(c,d)} λ(𝔢) N(𝔢)⁵ ε_{jk}(c,d) ě(δd/c) (N(𝔢)⁶N(c))^{−s} converges absolutely as a triple sum over (𝔢, c, d). Here ε_{jk} = ψ̄(σ_jγσ_k^{−1}), with |ε| = 1.

### 3.7 Identification with D_c: the ray-class twist [R-, conventions only]

Take ψ = κ̄, i = ∞ and k = k(s₀), the cusp 1/s₀ with σ_k = (1 0; s₀ 1).
- This cusp is essential: on its stabiliser κ(1 − s₀μ, μ; −s₀²μ, 1 + s₀μ) = 1.
- The double cosets Γ_∞\Γσ_k/Γ_∞ are exactly the coprime (c, d mod cN₀) with c ≡ s₀ and d ≡ 1 (N₀); every such pair occurs (construction of Lemma M).
- By Lemma M, κ(γσ_k^{−1}) = χ_c(d).
- Lift d by CRT: d̃/(cN₀) ≡ x·N̄₀/c + c̄/N₀ mod O.

This gives, for δ = −m/(λN₀) ∈ Λ*,

  Σ_{d̃} χ_c(d̃) ě(δd̃/c) = e(−m s̄₀/N₀) · χ_c(N₀) · g(c,−m).

The map c ↦ χ_c(N₀) is constant on classes mod N₀: its conductor divides 72·rad N₀ | N₀, by reciprocity and the mod-72 periodicity of χ_·(λ), χ_·(2), χ_·(ζ₆) (`D_derivation.md` §2, checked to norm 3·10⁴). Hence

  **D_{s₀}(w;m) = V · e(m s̄₀/N₀) · χ_{s₀}(N₀)^{−1} · φ^{(κ̄)}_{∞,k(s₀)}(−m/(λN₀), w + 1/2),**   (3.7)

up to the orientation conventions of `V_derivation.md` §3.1/§4. The phase depends only on m mod N₀ and has modulus one.

**Uniformity in m.** (3.7) is an identity for every m ≠ 0. The Fourier coefficient is the complete sum over all c ≡ s₀, including s that share primes with m and s that are not squarefree. No reduction (`LoA6.md` gap G-b) is ever needed.

**Theta coefficients.** ρ_c(m) := Res_{w=2/3} D_c(w;m) is the corresponding combination of the Fourier coefficients of Kubota's theta function θ = Res_{s=7/6} E_∞(·,s,κ̄) at the cusps 1/s₀. This agrees with `V_derivation.md` §4.

### 3.8 Growth in vertical strips, uniform in m

**Step 1: finite order [R-].** For each fixed δ, (s−1)(s−7/6) Z_{ik}(δ,s) has finite order on −1/4 ≤ Re s ≤ 9/4. Two routes:
- **(a) Maass–Selberg.** Bound ‖Λ^T E_i(·,s)‖² for 1 ≤ σ ≤ 2 using the Maass–Selberg relations (EGM Ch. 6) and the explicit, polynomially bounded Φ of Lemma C, with its σ-derivative near σ = 1. Extract φ_{ik}(δ,s) by Parseval over a height window v ∈ [v₀, 2v₀]. In the oscillatory range 4π|δ|v₀ ≪ |t| the average of |K_{s−1}|² is ≫ e^{−π|t|}|t|^{−1}, which cancels 1/Γ(s); this gives polynomial growth on 1 ≤ σ ≤ 2. Extend to the left by (3.6) and (M3).
- **(b) Citation.** Müller 1989 §4 proves finite order for such Eisenstein series [S]. EP92 / Eckhardt's thesis (Thm 2.4, quoted by Bröker–Hoffstein) proves exponential decay of the completed series in strips [S].

**Step 2: Phragmén–Lindelöf [R given Step 1].**
- On Re s = 9/4: |Z| ≪ N(δ)^ε, by (Z2) and the Gauss-sum bound |Σ_d| ≪ N(c)^{1/2} N((c,δ))^{1/2}, or trivially, since Re s > 2.
- On Re s = −1/4: by (3.6) and (M3), |Z| ≪ N(δ)^{5/4+ε}(1+|t|)^K.
- Phragmén–Lindelöf in between (any finite-order a priori bound suffices) gives

  |(s−1)(s−7/6) Z_{ik}(δ,s)| ≪ N(δ)^{5/4+ε}(1+|t|)^{K'} on the strip.   (3.8)

On Re s ≥ 1 we have |1/((s−1)𝓛(s))| ≪ log^{|Ψ|}(|t|+3). So D_c(w;m) ≪ N(m)^{A}(1+|Im w|)^{A} on Re w ≥ 1/2, outside a fixed disc around w = 2/3, uniformly in m.

### 3.9 Cross-checks

- **Gamma factor.** The factor in (3.3) composed with ζ_K's FE is (2π)^{2−2s}Γ(s)/Γ(2−s) · π/(s−1) · Γ_C(7−6s)/Γ_C(6s−6). With s_B = s−1 this equals exactly **6 ×** Bröker–Hoffstein's Γ(s_B)Γ(−6s_B)/(Γ(6s_B)Γ(−s_B)) (checked with mpmath at t = 7.3 and 55.1: ratio 6.000000000000).
  - Its growth exponent tends to 10(1−σ). Measured log|G|/log t: 7.36 (σ = 0.3), 5.29 (σ = 0.5), 1.15 (σ = 0.9) at t = 10⁴, against predictions 7, 5, 1, with the expected O(1/log t) drift.
  - The N(r)^{1−s} conductor and the location of the poles also match EP92/BH [S].
- **Cubic analogue.** Patterson's F(r,s) = N(r)^{1−s}[A(s)F_∞ + B(s)F] (DdFDS (6.11) [V]) has the same structure:
  - a Dirichlet-polynomial matrix;
  - the dual numerators r·νλ^b are S-class shifts of r;
  - Gamma exponent 4(1−σ) = 2(n−1)(1−σ) with n = 3;
  - normaliser ζ_λ(3s−2), the n = 3 case of L(ns−n+1). For n = 6 this is L(6s−5) ✓.

---

## 4. Theorem B (LoA6)

**Theorem B.** Fix ε ∈ (0, 1/12). There are K = K(ε, N₀) and C_ε such that, for all Q ≥ 1, t ∈ R, every essential i and every cusp k,

  Σ_{δ∈Λ*+h_k, 0<N(δ)≤Q} |φ_{ik}(δ, 1+ε+it)|² ≤ C_ε Q^{1+ε}(1+|t|)^K.

By (1.1) and (3.7), the case i = ∞, k = k(s₀) (and Cauchy–Schwarz over the finitely many s₀) gives

  **Σ_{0<N(m)≤Q} |D_c(1/2+ε+it; m)|² ≪_ε Q^{1+ε}(1+|t|)^K,**

and the same at every other cusp. This is LoA6, and also the input "LoA6 at 𝔞" used for G3 in `IDEA1_NOTES.md` §7.

### 4.1 Approximate functional equation [R given Theorem A]

Fix s = 1+ε+it and X ≥ 1. Define

  G_s(u) := e^{u²} (u+s−1)(u+s−7/6) / ((s−1)(s−7/6)).

Then G_s(0) = 1, and G_s vanishes at u = 1−s (the pole of 𝓛) and at u = 7/6−s (the theta pole). Since |s−1| ≥ ε and |s−7/6| ≥ 1/6−ε,

  |G_s(u)| ≤ e^{(Re u)² − (Im u)²}(1+|u|)²/ε²,

uniformly in t.

Shift (2πi)^{−1}∫ Z_{ik}(δ,s+u) X^u G_s(u) du/u from Re u = 2 to Re u = −5/4−ε.
- By (Z1), the only poles of the integrand in between are at u = 0 (residue Z(δ,s)), at 1−s and at 7/6−s; the last two are cancelled by G_s.
- The horizontal segments vanish by (3.8) and the Gaussian.
- On Re u = 2 expand by (Z2).
- On Re u = −5/4−ε apply (3.6). Then Re(2−s−u) = 9/4 > 2, so the dual series also converges absolutely.

The result is

  Z_{ik}(δ,s) = 𝒮₁(δ) + 𝒮₂(δ),
  𝒮₁(δ) = V^{−1} Σ_𝔢 Σ_{(c,d)∈𝒞_{ik}} λ(𝔢) N(𝔢)⁵ ε_{ik}(c,d) ě(δd/c) P^{−s} V_s(P/X),
  𝒮₂(δ) = N(δ)^{1−s} V^{−1} Σ_j Σ_𝔢 Σ_{(c,d)∈𝒞_{jk}} λ(𝔢) N(𝔢)⁵ ε_{jk}(c,d) ě(δd/c) P^{−(2−s)} Ṽ_{ij}(XP/N(δ)),

where P := N(𝔢)⁶N(c), and

  V_s(y) = (2πi)^{−1} ∫_{(2)} y^{−u} G_s(u) du/u,   Ṽ_{ij}(y) = (2πi)^{−1} ∫_{(σ_u)} y^u H_{ij}(u) du,
  H_{ij}(u) = (2π)^{2−2(s+u)} (Γ(s+u)/Γ(2−s−u)) M_{ij}(s+u) G_s(u)/u.

By (M1)–(M3), H_{ij} is holomorphic in Re u < −ε and ≪ (1+|t|+|Im u|)^{K(σ_u)} e^{σ_u² − (Im u)²}. So σ_u may be any number < −ε; the series 𝒮₂ is the same for every such choice. Also V_s(y) = 1 + O(y^B) as y → 0 and V_s(y) ≪ y^{−2} as y → ∞, uniformly in s.

### 4.2 The large sieve [R]

**Lemma LS.** This is the paper's Lemma `planar-additive-sieve` after the linear change of variables Λ ↦ O (frequencies δ' = λN₀δ ∈ O, points x/N₀ ∈ C/O). Let x_j ∈ C/Λ be Δ-separated (|x_j − x_l − μ| ≥ Δ for j ≠ l and all μ ∈ Λ). Then

  Σ_{δ∈Λ*+h, N(δ)≤Q} |Σ_j a_j ě(δx_j)|² ≪_{N₀} (Q + Δ^{−2}) Σ|a_j|².

The shift h is absorbed into the a_j.

**Separation.** For pairs (c,d) ≠ (c',d') with N(c), N(c') ≤ C and μ ∈ Λ ⊂ O,

  |d/c − d'/c' − μ| = |dc' − d'c − μcc'| / |cc'| ≥ 1/C

whenever the numerator is non-zero. If it is zero, then c | c' and c' | c, so c' = ηc with η a unit and d' ≡ ηd mod cΛ. Hence each fraction class carries at most |O^×| = 6 double cosets. Group them, at a cost of a factor 6.

**Mass.** For each c there are N(c) values of d mod cN₀ in the required class mod N₀.

This is all the large sieve needs. No Gauss-sum evaluation, no multiplicativity, no squarefree or coprimality reduction, and no multiplicative sextic sieve are used.

### 4.3 The direct piece [R]

Split P ∼ C dyadically. Apply Cauchy–Schwarz over 𝔢 with N(𝔢)⁶ ≤ 2C, using weights N(𝔢)^{∓(1+ε)}. For fixed 𝔢 apply Lemma LS with c ∼ C/N(𝔢)⁶ and Δ = N(𝔢)⁶/C. The coefficient mass is ≪ (C/N(𝔢)⁶)^{−2ε} N(𝔢)^{−2−12ε}. Summing over 𝔢 and over the dyadic C, with the cutoff V_s(P/X) ≪ min(1,(X/C)²),

  Σ_{N(δ)≤Q} |𝒮₁(δ)|² ≪_ε (Q + X²) Q^{ε}.

### 4.4 The dual piece [R]

Split N(δ) ∼ Q₁ and P ∼ C dyadically, and insert the Mellin form of Ṽ on a line σ_u, chosen per block. By Cauchy–Schwarz in u (weight |H|), Cauchy–Schwarz over 𝔢, and Lemma LS, the block (Q₁, C) contributes

  ≪ (∫|H|)² · Q₁^{−2ε} (Q₁ + C²) C^{3ε} (XC/Q₁)^{2σ_u}.

- For C ≤ Q₁Q^ε/X, take σ_u = −2ε. The factor (XC/Q₁)^{−4ε} costs at most Q^{O(ε)}, so the block is ≪ Q^{O(ε)}(Q₁ + Q₁²/X²).
- For C > Q₁Q^ε/X, take σ_u = −2. The geometric decay (XC/Q₁)^{−4} sums to ≪ Q₁ + (Q₁/X)².

So Σ_{N(δ)≤Q} |𝒮₂(δ)|² ≪ Q^{O(ε)}(Q + Q²/X²)(1+|t|)^{K}. The t-power comes from ∫|H| on the lines Re(s+u) ∈ {1−ε, −1+ε}.

### 4.5 Conclusion [R]

Take X = Q^{1/2}. Both pieces are ≪ Q^{1+O(ε)}(1+|t|)^K. On Re s = 1+ε, |1/𝓛(s)| ≪_ε 1, so the same bound holds for φ_{ik}. Rename ε. ∎

**Remarks.**
1. **Verification of `D_derivation.md` §6.2, Addition 1.** In the classical formulation both pieces are additive sums over well-spaced fractions d/c with unimodular coefficients, and the paper's planar sieve gives (Q + X²)Σ|b|². The mass Σ|b|² ≍ Σ_{N(c)≤X} N(c)^{1−2Re s} is bounded at Re s = 1+ε (and ≍ log X at Re s = 1). The BGL bound, de Faveri's sextic bound (arXiv:2610.04045, not verified by me), and the `LoA6.md` §3 reductions are all unnecessary.
2. **Addition 2 (bootstrap) is unnecessary.** The theta contribution ≍ X^{1/6}ρ(δ) lives inside 𝒮₁, and the large sieve bounds it there. This is consistent: X^{1/3}Σ|ρ|² ≍ X^{1/3}Q^{5/6} ≤ Q + X².

   The AFE sketched in Addition 2 also omitted the pole of the ζ-normaliser at s+u = 1 (the analogue of 𝓛 here). G_s removes it. The alternative (keep that residue and bound Σ|φ(δ,1)|²) would need a centre-value mean square that is not available.
3. **Comparison with DdFDS 6.6.** They shift only to the reflected line and close with a sup-over-t bootstrap, using Patterson's explicit A(s), B(s). Here the dual series is taken absolutely convergent at Re = 9/4, which needs (M3) on Re s = −1/4. That is supplied by Lemma C, so no bootstrap is needed.
4. **No GRH.** Zeros of the L(·,ψ') enter only through 1/𝓛 on Re s ≥ 1+ε, where it is bounded.

---

## 5. Corollary D: the clean Gram bound [R given Theorems A, B]

**Setting.** For the paper's additive polynomial

  A_{m,σ,v}(Y) = Y^{−1} Σ_s c(s) W₁(q_s/Y) (q_s/Y)^{−1/2+iv} g̃(s,−m)   (paper line 3460),

Mellin inversion of r ↦ W₁(r) r^{−1/2+iv}, whose transform is Ŵ₁(z−1/2+iv), gives

  A_m(Y) = Y^{−1} (2πi)^{−1} ∫_{(3/2)} Ŵ₁(z−1/2+iv) Y^z D_c(z;m) dz.

The line is in absolute convergence, by (Z2) and the identification (3.7).

**Contour shift.** Shift to Re z = 1/2+ε.
- The only pole crossed is z = 2/3 (Theorem A). Its residue is A^sig_m := Y^{−1/3} Ŵ₁(1/6+iv) ρ_c(m).
- The horizontal segments vanish, by (3.8) and the rapid decay of Ŵ₁ in strips.

**Bound.** Cauchy–Schwarz in Im z gives

  |A_m − A^sig_m|² ≤ Y^{−1+2ε} ‖Ŵ₁(ε+i·)‖_{L¹} ∫ |Ŵ₁(ε+i(τ+v))| |D_c(1/2+ε+iτ; m)|² dτ.

Sum over 0 < N(m) ≤ C_Ω Q, the row ball of the smooth majorant Ω, and apply Theorem B:

  **Σ_m Ω(q_m/Q) |A_{m,σ,v}(Y) − A^sig_m|² ≪_ε (1+|v|)^K (Q/Y) Q^ε Y^{2ε},  uniformly in Y ≥ 1.**

**Consequences.**
- Since Q, Y ≤ Z^{O(1)}, the right side is ≪ (1+|v|)^K (Q/Y) Z^{ε'}. Both the P_a^{1/6} and P_a²/Y terms of Prop. `probe-gram` (eq. `gram-quantitative`) are absent for A^rest := A − A^sig.
- The factor (1+|v|)^K is integrated against Ŵ₀(iv) in `eq:low-separated` and is harmless; compare the (1+|ν|)^{J_Gr} of Prop. `probe-gram`.
- With this, claim (i) of the task is [R] modulo the [R-] items of §7.
- **Numerical support** [N]: `D_numerics/` §7.4 shows the residual ‖A − Y^{−1/3}Ŵρ‖² flat at about 0.96 × the k = 0 diagonal for P_a²/Y up to ~10⁸ (`D_derivation.md` §7.4).

---

## 6. Corrections to the earlier working files

1. **`LoA6.md` §1, normaliser.** It says "ζ_K(6w−2)/ζ_K(6w−1) … only pole … at 6w − 2 = 2". This is wrong as written.
   - At the cusp pair (∞,∞) the constant term is ∝ ζ_K^S(6w−3)/ζ_K^S(6w−2), i.e. ζ(6s−6)/ζ(6s−5); the pole is at 6w−3 = 1, i.e. w = 2/3.
   - At general cusp pairs it is a combination of L(6w−3,ψ')/L(6w−2,ψ') (Lemma C).
   - The non-constant coefficients are normalised by Π_{ψ'} L(6w−2,ψ'), not by one ζ.
   - The pole conclusion is unchanged.
2. **`LoA6.md` §1 and `D_derivation.md` §6.2, source.** Both cite Brubaker–Bump 2006 for the n = 6 FE, but BB06's Theorem 1 assumes μ_{2n} ⊂ F (−1 an n-th power), which fails for Q(ω), n = 6. DIPP's restatement has the same hypothesis. Replace with Kubota 1969 §§5–6 [V] plus EGM/Selberg (this file, §3). EP92 §5 as applied by Bröker–Hoffstein [S] is an independent S-integer statement.
3. **`LoA6.md` §2, conductor.** "Conductor ≍ N(m)(1+|t|)^{O(1)}" is confirmed. The archimedean part is (1+|t|)^{10}. The dual side is the Eisenstein coefficient family at all essential cusps, with the same frequency.
4. **`LoA6.md` §3, large sieve.** The sextic multiplicative large sieve and the (s,m) ≠ 1 / non-squarefree reductions are unnecessary; use Lemma LS. Gap G-b is gone. Gap G-c (ray class) is handled by (1.1) and (3.7). Gap G-a is reduced to Lemma C.
5. **`D_derivation.md` §6.2, Addition 2.** Replace it by the choice of G_s in §4.1. Addition 2's own AFE also had an unhandled pole at s+u = 1, the ζ-normaliser pole.

---

## 7. What is left to write out (no heuristics; all [R-])

1. **Lemma C at every pair of essential cusps of Γ(N₀).** This needs the multiplier κ(σ_iγσ_j^{−1}) = χ_{c^S}(d)^{±1}·h_{ij} at general cusps, a generalisation of Lemma M, and the bookkeeping of the S-adic generating functions. Kubota did the pair (∞,∞) [V]. What is used downstream:
   - (i) holomorphy of P_{ij,ψ'} in Re s > 1;
   - (ii) poles of P only on Re s = 1;
   - (iii) the Hecke-L structure (3.2).

   Even (iii) is needed only to obtain (M2)–(M3). An alternative proof of (M3) on the line Re s = −1/4: use the S-integer FE of EP92 (via BH, (3.1)), which has T_{ij} ∈ Q̄(q_v^{−s}) and G_f(s) = ∏(1−q_v^{1−ns})^{−1}. That route is citation-only [S].
2. **Isotypic transfer of the EGM statements to (3.1).** This includes the fact that the poles of E_i in (1,2] are poles of φ_{ii}(0,s,κ).
3. **Finite order (3.8, Step 1).** Write out the Maass–Selberg argument of route (a), or verify Müller 1989 §4 directly. Note: only finite order per δ is needed; uniformity in δ then follows from Phragmén–Lindelöf.
4. **The constant in (3.7).** Fix the orientation conventions (κ vs κ̄, e vs ě, the sign of m, the phase e(m s̄₀/N₀)). None of this affects any estimate.

**Out of scope.** (V), with its gaps G3 (now reduced to LoA6 at 𝔞 by `IDEA1_NOTES.md` §7 plus this file) and G4, and the Part II bookkeeping (P). The conditional theorem still depends on those.

---

## 8. References

**Verified ([V])**
- T. Kubota, *On automorphic functions and the reciprocity law in a number field*, Lectures in Math. Kyoto Univ. 2, Kinokuniya, 1969. hdl:2433/84907 (OCR text read; §1 hypothesis, §5 (34)–(43), §6, Theorem 6).
- B. Brubaker, D. Bump, *On Kubota's Dirichlet series*, J. reine angew. Math. 598 (2006) 159–184. Postscript of 8 Mar 2006 read: §1 hypotheses, Theorem 1.
- A. Diaconu, B. Ion, V. Paşol, A. Popa, arXiv:2607.27131, Thm 3.2 and remarks (restatement of BB06; μ_{2r} hypothesis).
- C. David, A. de Faveri, A. Dunn, J. Stucky, arXiv:2410.03048v2, §6: Lemmas 6.1, 6.3 (FE (6.11), convexity (6.13)), 6.5, 6.6.
- A. de Faveri, A. Dunn, J. Hoffstein, arXiv:2607.07911v2: §2.2, §§3, 5 (Eisenstein framework, FE of the Fourier coefficients), Prop. "Second moment of ψ".
- R. Bröker, J. Hoffstein, *Fourier coefficients of sextic theta series*, Math. Comp. 85 (2016), arXiv:1312.0568: §2 (Thm 2.4) and §3 (FE (3.1)).
- J. Bernstein, E. Lapid, arXiv:1911.02342 (introduction: finite order via Müller 1989 §4; Moeglin–Waldspurger cover covering groups).
- paper78: Lemma `planar-additive-sieve`, Prop. `probe-gram`, eq. `low-separated`, coefficient c_σ (lines 3460–3500, 3600–3670, 8340–8560).

**Secondary ([S])**
- T. Kubota, J. reine angew. Math. 222 (1966) 55–57.
- J. Elstrodt, F. Grunewald, J. Mennicke, *Groups acting on hyperbolic space*, Springer 1998 (theorem numbers as cited in DFDH).
- P. Cohen, P. Sarnak, *Notes on the trace formula*, Ch. 6.
- S. J. Patterson, *A cubic analogue of the theta series*, J. reine angew. Math. 296 (1977) 125–161.
- C. Eckhardt, S. J. Patterson, Proc. LMS (3) 64 (1992) 225–264, §5.
- W. Müller, Ann. of Math. 130 (1989) 473–529, §4.
- C. Mœglin, J.-L. Waldspurger, *Spectral decomposition and Eisenstein series*, CUP 1995.

---

## Referee report

*Independent skeptical referee, 9 Oct 2026. Labels as above: [R] rigorous, [R-] routine-not-written, [H] heuristic, [G] gap.*

### Verdict: not refuted (sound with gaps)

I found no mathematical error that changes the conclusions of this file. (FE), (LoA6) and Corollary D hold at the stated level: [R] given three unwritten standard lemmas, which the file labels honestly as [R-]. The file is, however, silent on a chain-level fatal issue (§R5): Corollary D on its own gives **no** improvement past Liu.

### R1. What I checked and confirm

**Sources.**
- **Brubaker–Bump 2006 hypothesis.** Read in `scratchpad/lit/kubota_bb.txt`, §1: "we will assume that F contains μ_n … and moreover that −1 is an n-th power in F. If n is even this implies that F actually contains μ_{2n}." So the correction in §6.2 is right: BB06 Thm 1 does not apply to n = 6 over Q(ω).
- **Kubota 1969.** In `kubota1969.txt`, §1 assumes F totally imaginary and containing μ_n, and nothing more. In §6, the (∞,∞) constant term is L(ns−n,λ)/L(ns−(n−1),λ), with the "general theory" statement that singularities in Re s > 1 come from the constant term. Both are confirmed.
- **Bröker–Hoffstein.** `bh.txt` §§2–3 confirms Thm 2.4 (at most the pole s = 1/n in Re s ≥ 0) and FE (3.1), with G_f(s) and T_ij ∈ Q(|π_v|^s). This is for n = 6 over Q(ζ₆) = Q(ω), so an S-integer FE for exactly this case exists in the literature (EP92 §5, via BH), as stated.

**Lemma A1 (multiplier on Γ(N₀)).** I re-derived it via the product formula. The local facts check out:
- 1+8O₂ ⊂ squares and 1+2O₂ ⊂ cubes;
- (U^{(2)})³ = U^{(4)} in Q₃(√−3), since e/(p−1) = 1.

**Identification (3.7).** I recomputed it independently.
- The cusp 1/s₀ is essential: the conjugated translation (1−s₀μ, μ; −s₀²μ, 1+s₀μ) has κ = 1, by the product formula.
- The double cosets are exactly the coprime (c′, d′ mod c′N₀) with c′ ≡ s₀ and d′ ≡ 1 (N₀). Every such pair is realised: choose b′ ≡ 0 (N₀) in the solution coset of a′d′ − b′c′ = 1.
- κ(γ) = ((c′−d′s₀)/d′)₆ = (c′/d′)₆ = χ_{c′}(d′)^{±1}.
- The d-sum is e(−m s̄₀/N₀) χ_{c′}(N₀) g(c′,−m).
- The identity holds for every m and every c′ ≡ s₀, with no squarefree or coprimality reduction needed. **Confirmed.**

**FE of the Fourier coefficients (3.3).** The v-powers match and K_{s−1} = K_{1−s}, so φ_ik(δ,s) = (2π)^{2−2s} N(δ)^{1−s} Γ(s)/Γ(2−s) Σ_j Φ_ij(s) φ_jk(δ,2−s). Φ is independent of δ. **Confirmed.**

**(M1)–(M3) and (Z1), given Lemma C.**
- Trivial zeros L(−k,ψ′) = 0 hold for finite-order Hecke characters of an imaginary quadratic field, because the gamma factor is Γ_C.
- At s = 5/6 they kill the dual theta pole. At s = 0, −1, … they kill the poles of Γ(s).
- P has poles only on Re s = 1.
- 𝓛(s)/𝓛(2−s) is finite at s = 1.

**Confirmed.**

**AFE (§4.1).** G_s vanishes at 1−s and at 7/6−s, so the only residue is at u = 0. The dual series on Re(2−s−u) = 9/4 converges absolutely. H_ij is holomorphic in Re u < −ε. Ṽ_ij is independent of σ_u. **Confirmed.**

**Large sieve (§§4.2–4.4).** I redid the block computation. With fixed 𝔢 and N(c) ∼ C/N(𝔢)⁶:
- the spacing is Δ = N(𝔢)⁶/C;
- the mass is ≍ N(𝔢)^{−2} C^{2ε+2σ_u};
- the δ-factor is Q₁^{−2ε}(X/Q₁)^{2σ_u};
- so the block bound is Q₁^{−2ε} C^{2ε} (Q₁ + C²)(XC/Q₁)^{2σ_u}, as stated.

With the two choices σ_u = −2ε and σ_u = −2, and X = Q^{1/2}, both pieces are ≪ Q^{1+O(ε)}(1+|t|)^K. The additive sieve suffices here exactly as in DdFDS 6.6. This is LoA "at the edge": family size equals conductor. **Confirmed.**

**Corollary D.** The Mellin inversion is correct: the residue is Y^{−1/3}Ŵ₁(1/6+iv)ρ_c(m), and the Cauchy–Schwarz and v-integration are right. **Confirmed.**

### R2. Minor errors (non-load-bearing; no conclusion changes)

1. **The Γ-factor cross-check constant is 3, not 6.** I recomputed (2π)^{2−2s}Γ(s)/Γ(2−s)·π/(s−1)·Γ_C(7−6s)/Γ_C(6s−6), with Γ_C(u) = 2(2π)^{−u}Γ(u) as written in §3.9. Divided by BH's Γ(s_B)Γ(−6s_B)/(Γ(6s_B)Γ(−s_B)), it gives 3·(2π)^{10 s_B} exactly (mpmath, σ ∈ {1, 0.7, 0.3}, t ∈ {7.3, 55.1}).
   - The factor (2π)^{10 s_B} is the archimedean part of BH's y₁^{2s}.
   - So "ratio 6.000000000000" is not reproducible from the formula as written. Most likely a different K(0,s) normalisation was used.
   - The growth exponent 10(1−σ) is correct.
2. **Choice of N₀ (§1.2).** Paper Lemma `fixed-numerator-ray` (line 647) explicitly asserts *no* bound on the conductor exponents at S of s ↦ χ_s(b_*). The file's "72·𝔪_T·b_* | N₀" therefore does not by itself guarantee that c is periodic mod N₀. The constancy of χ_c(N₀) on classes mod N₀ (§3.7) rests on a numerical mod-72 check.
   - Fix: require 𝔣_S | N₀, where 𝔣_S is the common modulus of that lemma for all numerators u·Π π_𝔭^{v_𝔭}. This is not circular, since the family depends only on v_𝔭 mod 6.
   - This is trivial, but the numerical check should not be cited as the proof.
3. **Starting line in Corollary D.** The line Re z = 3/2 is the boundary of (Z2), which needs Re s > 2. Instead, use absolute convergence for fixed m in Re z > 1: g(p^k, m) = 0 for p ∤ m and k ≥ 2. Or start at Re z = 2. Trivial.

### R3. Labels: where the rigor actually sits

The table's [R] entries for "only pole at w = 2/3", "holomorphy on Re w = 1/2", (M1)–(M3) and hence Theorem B are all conditional on **Lemma C at every pair of essential cusps**, which is [R-]. That makes them effectively [R-]. The file says this in §7, but the verdict sentence "all go through" should carry that qualifier.

Lemma C is the most substantive unwritten item. It is not only bookkeeping: it supplies holomorphy of M in Re s < 1 and its bound on Re s = −1/4. Without it, Φ(s) = Φ(2−s)^{−1} could a priori have poles at zeros of det Φ(2−s), and the dual-at-9/4 device fails.

Two remarks may help:
- **For the pole claim, only the diagonal entries (i,i) are needed** (by §3.2(d)). Since Γ(N₀) ⊲ SL₂(O), the entry φ_ii(0,s,κ) equals the (∞,∞) entry for the conjugated character κ^{g_i}. Here κ^{g_i}/κ is a character of Γ(N₀)/Γ(N₀^k) (a local-splitting discrepancy at S). So Kubota's own (∞,∞) computation carries over with an S-adic weight [R-].
- **The off-diagonal entries matter only for (M2)/(M3).** The EP92/BH S-integer FE is a citation route for those [S].

The other two [R-] items are standard: the isotypic transfer (or Selberg's theory with a unitary multiplier, which applies directly) and finite order (Müller 1989, or Maass–Selberg).

### R4. Corrections section (§6): confirmed

- `LoA6.md` puts the pole of ζ_K(6w−2) at w = 1/2, not 2/3. The correct normaliser is ζ_K^S(6w−3)/ζ_K^S(6w−2).
- `D_derivation.md` Addition 2 is incomplete. If one works with the uncompleted D, shifting left of Re s = 1 crosses the zeros of 𝓛. If one completes, the s+u = 1 pole appears. Either way, completing by 𝓛 together with the G_s weight is the right fix.

### R5. Missed chain-level fatal issue (outside this file's scope, but it bears on "past Liu")

The file ends with "the conditional theorem at 0.8665 still depends on [(V), (P)]". That understates the situation. I independently re-ran two checks from `AUDIT_chain.md`.

1. **`audit_lp_no_theta.py` (re-run).** With Corollary D alone, i.e. the P²/Y term dropped and the theta part still bounded by Cauchy–Schwarz:
   - θ = 0.894737 at zero slots;
   - θ = 0.875092 at the best slot supply, which is the paper's value.

   So **(FE)/(LoA6)/Corollary D give no improvement by themselves.** ‖A^sig‖² = Q^{5/6}Y^{−2/3} is exactly the P_a^{1/6} term, and P_a²/Y is not binding. All of the gain to 0.8662/0.8665 comes from step (ii).
2. **`check_row_phase.py` (re-run; 6028 split primes, N(p) ≤ 60000).**
   - In all 432 classes mod 72, w(p) := ᾱ(p)²τ⁻_{p,1}τ⁺_{p,3} times γ₂(p)⁴ is constant within the class.
   - w(p) itself takes 18–20 distinct arguments per class, with median within-class spread 1.99.
   - This confirms F1. The reflected row's scalar ζ (paper eq. `reflection-row-phase`, line 1832, which contains ᾱ(c)² and ∏_{active} ω_{p,j_p}) is an m-dependent cubic-Gauss-sum/Grössencharacter phase, not a ray phase.
   - The pairing that actually occurs is therefore Σ_m ρ(m)ζ_m(λ⁴μ/m)₂…, to which (V) for ρ(m)(μ/m)₂ does not apply.

**Conclusion for the headline.** This file correctly establishes ingredient (i) up to [R-] items. But (i) is not the binding ingredient. The claimed improvement past 0.874957 is unsupported until (ii) is restated as a hypothesis (V′) for ρ(m)ζ_m and proved [G]. No route is known.

### R6. Summary of referee findings

| item | referee status |
|---|---|
| Theorem A (continuation, poles, FE, growth) | [R-] overall (Lemma C, transfer, finite order); no error found |
| Theorem B (LoA6 at every cusp) | [R] given Theorem A; the computation was re-checked |
| Corollary D | [R] given Theorem B; two trivial fixes (R2.2, R2.3) |
| Literature corrections (BB06, normaliser) | confirmed from source text |
| Γ cross-check "ratio 6" | wrong constant (3), harmless |
| Bearing on 0.8665 | none on its own; the chain is blocked by F1 at step (ii) |
