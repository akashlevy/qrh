# xmark — a cross-mark saving for lem:marked? (10 Oct 2026)

Agent `xmark`. P78 = `openai/.../September-30-2026/build/paper.tex`; line numbers refer to it.
Labels: **[P]** proved (given P78's lemmas), **[B]** proved as exponent bookkeeping, **[C]** conditional,
**[H]** heuristic, **[N]** numerical (exact Z[ω] arithmetic), **[M]** model.

## 0. Verdict

* **No power saving is proved.** There is no rigorous cross-mark saving of any positive power anywhere in
  r ∈ [0.66, 0.81], z ∈ [0.13, 0.23]. Exponent payoff: 0.
* Every listed tool was tried on the exact sum (§1). Each one stops at a named step (§3), and they share one
  core obstruction (§4).
* New material:
  - **[P] Prop X1 (pair large sieve).** For z ≤ 1/4 and arbitrary bounded slot coefficients,
    Σ_u|M_uQ_u|² ≪ U^ε(U + U^{1/2}(Σ_u|M_u|⁴)^{1/2}).
    - So for z ≤ 1/4 the cross-mark problem reduces to a mark-free Möbius fourth moment.
    - [B] It is useless as a route. Any count-based bound on the fourth moment gives back only the slot-only witness
      Q², which never binds.
    - [B] The fourth-moment input is also far stronger than needed: by itself it would give counts 1 − 2rδ.
    - [N] That input is true numerically: the moment ratio is ≈ 2.1, i.e. Gaussian.
  - **[P] Exact original-side form.** On the original side, lem:marked(r,z) is equivalent to a bound for a twisted
    second moment of the inverse witness, taken over pair characters ξ_{pp′}. In this form P78 already supplies a
    cross-mark saving of U^{1−r−z}.
  - **[P/B] No-go statements** for each suggested tool (§3), with the exact failing inequality.
  - **[M] Payoff of a partial saving δ** (`xm_model.py`):
    - for small δ the worst exponent improves by ≈ −0.07δ;
    - 41/47 needs δ* = 0.159 uniformly, with slot powers and the cubic branch improved too;
    - in relative form this is φ = 0.656 of the mark length, i.e. mk = 0.744.
  - **[N] Original-side numerics** (`orig_side.py`): the cross terms are random-size. The only coherent component is
    the predicted mixed diagonal. There is no correlation with (p/p′)₃ or with γ₂(p)·conj γ₂(p′).

## 1. The exact sum, the exact range, and what is already proved

**Parameters.** Row scale U = Z^m with m = 1. Inverse witness M_u = M_u(Z^r; W). Slots Q_u of total length z.
The binding inverse splits at x = 1/2 are r ≈ 0.66–0.81. They need z up to mk(1 − r), with mk = 0.744 and slot powers,
or z = 1 − r. That gives z ≈ 0.13–0.23 and r + z < 1 < r + 2z.

### (a) Canonical form (Lemma D, `agents/marked-cap/NOTES.md` §3)

At the initialization (P78 l.11606–12340) take the principal pattern: C = d′ = s = 1, G = P₁ = V = 0. Then

- N = D′ = r + z (column length);
- M = 2D′ − m (row length);
- F = D′ (target exponent).

Define

    E = Σ_{0<q_k≤Z^M} |Z^{−N/2} Σ_n A_*(n) χ_n(k) 𝔡(n) W(q_n/Z^N)|²,   A_* = ᾱγ₂ν_*,  𝔡(n) = Σ_τ a_τ 1_{τ|n},
    X_{kτ} = A_*(τ) χ_τ(k) Θ_k(τ),   Θ_k(τ) = Z^{−N/2} Σ_{n′} A_*(n′) χ_{n′}(kτ⁴) W(q_τ q_{n′}/Z^N),
    E = E_diag + E_off,   E_off = Σ_{τ≠τ′} a_τ ā_{τ′} G_{ττ′},   G_{ττ′} = Σ_k X_{kτ} conj X_{kτ′}.

The other init patterns have the same shape, with P₁, V > 0 and a puncture.

- **Target:** |E_off| ≪ Z^{F+ε}.
- **[P] Lemma D:** E_diag ≪ Z^{F+ε} whenever r + z ≤ 1 − c.
- **[P] Trivial bound (Cauchy–Schwarz):** |E_off| ≤ #τ·Z^{F+ε} = Z^{F+z+ε}.
- **[B] Bound from P78's argument** (row branch Z^{M+z} plus the cubic branch, with the invariant violated by e):
  |E_off| ≪ Z^{F+e(r,z)+ε}, where e = max(0, r + 2z − 1, (2r + 8z − 3)/4).
- **Saving needed beyond P78:** Z^{e}. At (0.69, 0.23) this is Z^{0.15}. Along z = 0.744(1−r) it is e = 0.488(1 − r),
  i.e. 0.166 / 0.151 / 0.122 / 0.093 at r = 0.66 / 0.69 / 0.75 / 0.81.

### (b) Original-side form [P]

This is new and elementary. Write τ for tuples (one prime per slot), χ_τ(u) = Π_i(u/p_i)₆, and
ξ_{ττ′} = (χ_τ χ̄_{τ′})^{ε_χ}. Then exactly

    Σ_u Φ(u/U)|M_uQ_u|² = Z^{−z} Σ_τ |a_τ w_τ|² Σ_u Φ|M_u|² 1_{(u,τ)=1}  +  C(r,z),
    C(r,z) = Z^{−z} Σ_{τ≠τ′} a_τ ā_{τ′} w_τ w_{τ′} ν(τ)ν̄(τ′) T(τ,τ′),
    T(τ,τ′) = Σ_u Φ(u/U) |M_u(Z^r)|² ξ_{ττ′}(u).

- The first term is ≪ U^{1+ε} by lem:marked at z = 0, since r ≤ 1 − c.
- Hence **lem:marked(r,z) ⟺ |C(r,z)| ≪ U^{1+ε}.**
- So the cross-mark problem is a **twisted second moment of the inverse witness by the pair characters ξ_{ττ′}**.
  These have conductor ≤ Z^{2z}, far below the row range U.

How much is already known on this side:

- Trivially |C| ≤ U^{1+z}.
- [B] P78's argument together with Lemma D gives |C| ≪ U^{1+e}. So **P78 already supplies a cross-mark saving of
  U^{z−e} = U^{1−r−z}** (rigorously [P] inside P78's range r + 2z ≤ 1).
- The open part is the remaining U^{e}.

## 2. What a partial saving buys [M] (`xm_model.py` → `xm_model.out`, `xm_thresh.py` → `xm_thresh.out`)

Setup: weighted-numerator high side at ℓ = 0.1773, B = 41/47, fine δ-grid, b ∈ [0, 0.1].
Model of a saving δ on E_off (U-units): the excess must satisfy e(r,z) ≤ δ. Cubic-branch options:

- keep: the branch is unchanged;
- relax: (2r + 8z − 3)/4 ≤ δ;
- drop: no cubic branch.

| hypothesis | worst E |
|---|---|
| paper (δ = 0) | +0.01124 |
| δ = 0.005 / 0.01 / 0.02 / 0.03 / 0.05 (cap kept, physical supply; identical for the other options at δ ≤ 0.05) | +0.01090 / +0.01036 / +0.00982 / +0.00906 / +0.00770 |
| δ = 0.10: keep / relax | +0.00442 / +0.00408 |
| δ = 0.13: keep / relax, physical supply / relax, unlimited supply | +0.00442 / +0.00236 / +0.00198 |
| δ = 0.15: relax, unlimited supply | +0.00057 |
| **δ* = 0.159: relax, unlimited supply (slot powers)** | **0** |
| δ = 0.2: relax, physical supply | +0.00187 (saturates at the DIAG-marked value) |
| relative saving φ, z ≤ (1−r)/(2−φ), unlimited supply, no cubic branch: φ = 0.5 / 0.6 / **0.656** / 0.7 / 1 | +0.0038 / +0.0014 / **−0.00002** / −0.0010 / −0.0078 |

Reading the table:

- **Small δ.** Any positive power pays about −0.07·δ in the worst exponent (slope −0.068 to −0.088 over δ ∈ [0.005, 0.1]). So a "tiny" saving would be a genuine,
  if small, gain.
- **With the cubic branch kept**, the gain saturates at +0.0044 for δ ≥ 0.1.
- **With physical supply**, it saturates at +0.0019.
- **Two caveats for any Hecke-structured tool:**
  - Unlimited supply comes from polarised slot powers (amplitude Prop 3.3). Their colour-mask coefficients are not
    ray characters. A tool that uses the Hecke structure of the slots (§3.6) must therefore be proved directly for the
    e_k(slot) marks.
  - The cubic branch, (NL)^{2/3}, is also entirely a cross-term phenomenon: Lemma D's diagonal needs only 3M ≤ 4F.

## 3. Route by route: what each tool gives, and the exact failing step

### 3.1 Averaging over the mark pair as a variable: the pair large sieve [P], useless [B]

**Prop X1 [P].** Let 0 ≤ r ≤ 1 − c and z ≤ 1/4. Let the slots be a bounded number of disjoint prime windows with
*arbitrary* bounded coefficients (product form is not needed). Then

    Σ_u Φ|M_u(Z^r)Q_u|² ≪_ε Z^ε [ U + U^{1/2} (Σ_u Φ|M_u(Z^r)|⁴)^{1/2} ].

*Proof.* Use §1(b).

1. *Group the cross terms.* Write τ = σα and τ′ = σα′, where σ is the common part (on the slot set J where the
   tuples agree) and α ≠ α′ differ in every remaining slot.
2. *The characters are primitive and distinct.* ξ_{ττ′} = 1_{(u,σ)=1}·(χ_α χ̄_{α′})^{ε}. Each (·/p)₆ is a nontrivial
   character of the field 𝒪/p, so χ_α χ̄_{α′} is primitive mod αα′. The pair (α,α′) is recovered from the character,
   because (·/p)₆ ≠ conj(·/p)₆.
3. *Apply the multiplicative large sieve for 𝒪 = ℤ[ω].* Fix σ and put β_u = Φ|M_u|²1_{(u,σ)=1}. Huxley's sieve
   bounds Σ over primitive χ with N(cond) ≤ Q₀ of |Σ_u β_u χ(u)|² by ≪ (U + Q₀²)Σ|β_u|².
4. *Conductor range.* Here Q₀ ≍ Z^{2(z−z_J)}, so Q₀² ≤ Z^{4z} ≤ U.
5. *Cauchy–Schwarz over the pairs.* Σ_{α,α′}|a_α a_{α′} w w|² ≤ (Σ|a_α w_α|²)² ≪ Z^{2(z−z_J)}. Summing over σ, with
   weight Σ_σ|a_σ w_σ|² ≪ Z^{z_J}, gives |C| ≪ 2^{|I|}·Z^ε(U·ΣΦ|M|⁴)^{1/2}.
6. *Diagonal.* The diagonal term is ≪ U^{1+ε}, as in §1(b). ∎

**Corollary [P] (dichotomy).** Suppose z ≤ 1/4 and lem:marked fails at (r,z) by U^η. Then Σ_u|M_u(Z^r)|⁴ ≫ U^{1+2η−ε}.
So any counterexample must make the inverse witness have a large fourth moment.

**Why this does not help [B]:**

1. **Bounding the fourth moment through counts gives nothing.** On the high side, the only information about |M_u| is
   the bin structure. Restrict β to a bin row set C, where |M_u|² ≈ U^{δr} and |Q_u|² ≥ U^{2qz}. Prop X1 then gives
   #C·U^{δr+2qz} ≲ #C·U^{δr} + (U·#C)^{1/2}U^{δr}, so **#C ≤ U^{1−4qz}**.
   - That is exactly the slot-only witness Q² (sextic sieve, 2z ≤ 1/2), the log-linear interpolation.
   - It is never binding: 1 − 4qz ≥ 1 − δ/2 > R_need.
   - Using the full Σ|M|⁴ instead, bounded by the counts R(δ′), cannot do better, because the bin itself contributes
     R(δ) + 2δr. Even the paper's own counts give only S₄ = max_{δ′}(R(δ′,0) + 2δ′r) = 1.28 / 1.32 / 1.41 / 1.50 at
     r = 0.66 / 0.69 / 0.75 / 0.81 (`m4count.out`). So Σ|MQ|² ≲ U^{1.14 / 1.16 / 1.21 / 1.25}, never better than P78's
     row branch U^{r+2z} at z ≤ 0.22.
2. **The untruncated input is far stronger than needed.** Σ_u|M_u(Z^r)|⁴ ≪ U^{1+ε} at 2r ≈ 1.4 lies beyond the
   u-family's diagonal: the sextic sieve has the column term U^{2r}, and the trivial bound is U^{1+r}. If proved, it
   would give the M²-witness count 1 − 2rδ directly (0.38 < 0.571 at δ = 0.47). That alone would close 41/47 by a wide
   margin. It is a Lindelöf-on-average fourth moment of 1/L over a thin family.
3. **The pair large sieve is still best possible as a step.** z ≤ 1/4 makes {ξ_{pp′}} a subfamily of a full
   character family below its edge (Z^{4z} ≤ U). So the loss sits entirely in the fourth moment, not in the sieve.

[N] The input is true numerically. Over the u-ball with R² = 50–90·U, Σ|M|⁴/(U m₂²) = 2.06 / 2.08 / 2.10
(`orig_side.out`). This is the Gaussian value 2.

### 3.2 The cubic residue symbol (P₁/P₂)₃, Heath-Brown's cubic sieve, de Faveri's n-th order sieve [P]

**Where the cubic symbol occurs (exact).**

- **Tuple cocycle.** With s = p₁ and n = p₂ in eq. initial-common-extraction: A_*(p₁p₂) = A_*(p₁)A_*(p₂)·conj(p₁/p₂)₃.
  The marks also enter Θ only through χ_{n′}(p)⁴ = conj(p/n′)₃. So
  X_{kp} = A_*(p) Σ_{n′} c_k(n′) (kn′⁴/p)₆, with c_k(n′) = Z^{−N/2}A_*(n′)χ_{n′}(k)W.
- **Mixed diagonal.** The cross terms with n₁p = n₂p′ contain no cubic symbol. In the original-side form their size is
  ≍ U/P per pair. They are harmless, since they sum to ≤ U, and they are the only coherent part [N, §5].

**(i) Sieve in p row by row is trivial.** For fixed k, E = Σ_k|⟨a, Y_k⟩_p|² with Y_k(p) = X_{kp}.

- Any bound that separates p from k by Cauchy–Schwarz at fixed k gives Σ_k‖a‖²‖Y_k‖² = L·E_diag, the trivial bound.
- A large sieve in p only bounds Σ_p|Y_k(p)|². It therefore cannot improve on this.
- Here the conjugate variable n′ has length N − z = r ≥ 2z. So Heath-Brown's (N′ + L + (N′L)^{2/3}) and de Faveri's
  bounds are in their diagonal regime anyway: they reproduce the p-diagonal.

**(ii) The bilinear cubic form in (p₁,p₂).** For fixed (k,n″) the two-slot sum is a bilinear form with kernel
conj(p₁/p₂)₃. The cubic sieve bounds it pointwise by L₀^{5/3}, against the trivial L₀².

- That pointwise bound exceeds the form's mean-square size L₀. The sextic sieve in h = kn″⁴ is diagonal, since
  H ≫ L₀⁴.
- Inserting a pointwise bound inside |Σ_{n″}|² also loses the n″-orthogonality, a factor Z^{N−2z₀}.
- So the (P₁/P₂)₃ saving acts on a quantity that is not the bottleneck.

**(iii) A cubic sieve over the pair product P₁P₂.** The cross term would have to factor as
Σ α_{P₁}β_{P₂}(P₁/P₂)₃ × (something independent of the pair). It does not. Exact identity (h2 §2.1b):
G_{P₁P₂} is a twisted correlation of the P₁- and P₂-progression pieces of one theta sum.

[N] The original-side cross terms show no correlation with (p/p′)₃: |corr| = 0.008 / 0.021 / 0.047 / 0.006, against a
null of 0.039 / 0.039 / 0.039 / 0.017.

### 3.3 Poisson in the slot primes [P]

1. Primes must first be split into type I/II pieces (Vaughan, Heath-Brown, or the cubic-first-moment identity
   Λ = (2m − m*m*1)*log). Every smooth piece has length E ≤ P = Z^z.
2. As a function of the slot variable, the summand is a character of conductor N(h):
   - h = u on the original side (N(h) ≤ U);
   - h ≍ k·n′ on the canonical side, with N(h) ≈ Z^{M+N−z} ≈ Z^{1.51} at the binding point.
3. Poisson in e has dual length N(h)/E > U/P ≫ E. It is useful only for E > N(h)^{1/2}, i.e. z > 1/2.
4. Averaging over h does not help: Poisson in h is the u ↔ k move, which is an involution with no contraction
   (h2 §1.3).

### 3.4 Bilinear structure in (τ, τ′) [P]

- **Rank one.** E_off = Σ_k|Σ_τ a_τX_{kτ}|² − E_diag is a rank-one functional of the pair matrix aa*.
- **Every inequality on the pair matrix is trivial.** Schur, Cauchy–Schwarz in τ or τ′, and the insertion of τ′ as a
  canonical label all return L·E_diag. The label insertion fails because the partner sum Y_k^{(τ′)} is independent of
  τ′ up to one term, so Σ_{τ′}|Y^{(τ′)}|² ≈ L·E.
- **Polarisation needs 2-prime marks, which cost the same.** ⟨X_τ, X_{τ′}⟩ is a combination of 2-prime marked
  energies. Their conductor cost depends on the prime size z, not on the number of primes (mark-design), so each costs
  Z^{M+z}.
- **No off-diagonal Cauchy–Schwarz.** Beyond the wall the row Gram has kernel dimension ≥ N − K, and there the
  off-diagonal form is indefinite (h2 §1.2).
- **Adversarial Gram (mark-design L5 transposed) [P].** G_adv = (K/L)(I + L·ûû*), with û the coherent direction:
  - it is consistent with Lemma D (diagonal ≈ K) and with P78 (a*Ga ≤ KL for |a| ≤ 1, any sub-support);
  - it is consistent with the rowwise cubic and sextic sieves;
  - it is consistent with Prop X1 under the proved bound Σ|M|⁴ ≤ U^{1+r}. On the original side G_adv is
    T(p,p′) ≈ U^{r+z}û_pū_{p′} (so C ≈ U^{r+2z}, P78-saturating), and Σ_{p≠p′}|T|² ≈ L²U^{2(r+z)} ≤ U^{2+r}
    ⟺ r + 4z ≤ 2, which holds in the whole range;
  - yet it gives E_off = KL.

  Hence **a proof must use a phase-sensitive input at the coherent vector**. The same holds for every random-sign
  version a_τ → ε_τ a_τ: all proved inputs are invariant under it.
- **Interpretation.** In the original form, G_adv means |M_u|² correlates coherently with |Q_u|². That is exactly
  "both witnesses saturated on the same rows", which is the high-side adversary itself. So a cross-mark bound *is* a
  count bound; nothing weaker suffices.

### 3.5 Twisted multiplicativity [P]

- **What it does.** It turns the mark into a cubic twist: A_*(τn′) = A_*(τ)A_*(n′)χ_{n′}(τ)⁴, so the column becomes the
  theta sum at the Kubota index kτ⁴.
- **Why no gain.** The twist has conductor τ, which enters the reflection squared. That is the cost (mark-design L2:
  amplification and conductor are one symbol). It produces the cocycle of §3.2, which carries no usable cancellation.
- **Two-sided / label versions.** The partner's label differs (τ′ ≠ τ), so no canonical (same-label) family arises.

### 3.6 Heath-Brown–Patterson, Patterson's bias, the cubic first moment (OpenAI, 25 Sep 2026)

**Where Gauss sums enter [P].**

- For ray-character slots, b_p := a_pA_*(p) = (Grössencharacter ᾱ)·γ₂(p)·(finite characters).
- In X_{kp} the phase appears as γ₂(p)χ_p(h) = γ₂(p;h)(h/p)₂, with h = kn′.

**Absorption identity [P].** γ₂(p;h) = γ₂(p)·conj(h/p)₃, so the Gauss-sum phase depends on p alone. In any large sieve
over h it is a coefficient: Σ_h|Σ_p b_pγ₂(p;h)ψ|² = Σ_h|Σ_p (b_pγ₂(p))conj(h/p)₃ψ|². Hence:

- Gauss-sum equidistribution (HBP), Patterson's bias, and the angular first moment act only on sums over p with
  **fixed** h.
- Seeing the cubic theta coefficient at h from primes p ~ P needs P > N(h)^{1/2} (Voronoi length). Here
  P/N(h)^{1/2} = Z^{0.22−0.755} on the canonical side and U^{0.22−0.5} on the original side. Both are hopeless.
- The OpenAI theorem is for one fixed angular mode, with no uniformity in growing parameters. It cannot be used at
  shifts N(h) ≫ P².

**The cubic-first-moment mechanism needs three ingredients:**

- (a) coefficient independence: genuine convolutions of prime sequences;
- (b) a dominant *theta* term in the sieve at the configuration;
- (c) BHM near-orthogonality strong enough to exclude theta profiles.

**Canonical side.** Both (a) and (b) fail:

- (a): A_* is only twisted-multiplicative, A_*(n₁n₂) = A_*(n₁)A_*(n₂)χ_{n₂}(n₁)⁴. The column sum is not a product of
  independent polynomials, so moments of different lengths are unavailable.
- (b): at the lem:marked edge the dominant term is the quadratic sieve's column term T_d = 2M − N + 2z. It is attained
  by every coefficient vector in the row span, so there is no exceptional coefficient class to exclude.

**Original side.** Here (a) holds, since μ·(slot) is a convolution. And (b) holds: lem:marked beyond P78 is "no
sextic theta term" U^{5/6}C^{1/3} for C = Z^{r+z}, which exceeds U once C > U^{1/2}. But (c) fails [B]:

- A sextic theta profile has V² = C^{4/3} on R rows.
- BHM, with rms Gram entries √C, only gives R·C^{4/3} ≤ C² + R·C^{3/2}. That always holds, because 4/3 < 3/2.
- In general, theta profiles are excluded only when 2/n > 1/2, i.e. n = 3. This agrees with 09 diag.

**Aside [C/H]: the cubic cap (old-eq:4.11).** It is the cubic sieve's theta term in the mark variable, with
coefficients ã(P) ∝ ᾱ(P)²τ⁺_{P,2}·(finite). A uniform-in-shift version (n ≤ L^{2−ε}) of the angular first moment would
remove it, because the angular twist kills the Kubota pole (h2 §4.5 measured it as nearly absent). Its payoff alone
is 0, since the cap binds only for r < 1/2. It matters only alongside a δ ≳ 0.1 cross saving (§2: +0.0044 → +0.0041
at δ = 0.1).

### 3.7 Burgess / van der Corput in the row variable [B]

- Each inner sum Σ_{q_k≤K}χ_{pn₁′}χ̄_{p′n₂′}(k) is a character sum of length K = q^{(2N−1)/(2N)} ≈ q^{0.45}, inside the
  Burgess range. It saves a fixed small power q^{−η_B}.
- But using it inside Σ_{n₁′,n₂′} requires absolute values. That loses the n′-orthogonality factor Z^{N−z} = Z^{r}.
- Net: |G_{pp′}| ≤ Z^{M+N−2z−2Nη_B}, against the trivial Z^{M−z}. This is a loss of Z^{r−2Nη_B} > Z^{0.6}.

### 3.8 Zero-free half-plane / Halász–Montgomery over the pair family [B]

- **Sparse rows.** On a bin set C, the large sieve over ξ_{pp′} can be replaced by Halász–Montgomery. That needs
  prime sums Σ_p χ_p(u/u′), which are only ≪ P^{β*} with β* = 0.875.
- **No gain.** The count gains only if 2q > 0.875·(slot exponent). Since 2q ≤ δ ≤ 0.75 < 0.875, there is no gain.
- **Under GRH.** P^{1/2} would suffice, but then everything is trivial anyway (10 NOTES §3).

## 4. The common failing step

In every route the mark variable p ends up in one of three places:

- **(A) Inside a square whose outer variable is beyond its edge.** This is the quadratic sieve over k after
  reflection: the dual length T_d = 2M − N + 2z exceeds F exactly when r + 2z > 1. The column term T_d‖c‖² is attained
  by every c in the row span V = span{χ_k}. Excluding it means showing that the marked theta vector
  c = Σ_P a(P)β·χ_P^{−2} is far from V. That is an uncertainty principle ‖Π_VΠ_X‖² ≲ K/N_d, i.e. (H2-ms).
- **(B) Separated from the row by Cauchy–Schwarz.** This returns L·E_diag (§3.2(i), §3.4).
- **(C) Moved into a moment of the inverse witness beyond its diagonal.** This needs Σ|M_u(Z^r)|⁴ ≪ U at 2r > 1
  (§3.1), which is stronger than the conclusion.

The extra structures (Gauss-sum phases, the cocycle (P₁/P₂)₃, primality, product form, Hecke coefficients) either:

- are coefficients absorbed by every sieve (§3.6);
- live on fixed shifts far shorter than their conductors (§3.6);
- sit in the mixed diagonal, which is already harmless (§3.2); or
- act pointwise on a quantity that is already at its mean-square size (§3.2(ii)).

None of them acts on (A).

[H] The truth is random-matrix-like (§5; 09 h2, mark-design, edge-theta). So a proof would have to supply a
freeness/uncertainty statement between the quadratic row space and the cubic-marked theta space. No such tool exists
for thin (one-character-per-modulus) families.

## 5. Numerics [N] (`orig_side.py` → `orig_side.out`; `mixdiag.py` → `mixdiag.out`)

**Setup.**

- Exact sextic symbols in ℤ[ω]; all rows u with N(u) ≤ U and weight Φ = (1 − N/U)².
- M_u over squarefree primary n ∈ (R, 2R] with R = U^{0.7}, so R²/U = 52–91.
- Split slot primes in (P, 2P].
- Random-phase control: the same |M|² weights with random sextic values.

| U | slots | m₄/m₂² | rms O, real / random / predicted | coherent Σ_{p≠p′}T/Σ_pT(p,p): real / random / mixed-diagonal prediction | corr with (p/p′)₃ | corr with γ₂γ̄₂ | null |
|---|---|---|---|---|---|---|---|
| 2·10⁴ | 26 | 2.055 | .0083 / .0071 / .0071 | .104 / −.010 / .110 | .008 | .056 | .039 |
| 4·10⁴ | 26 | 2.078 | .0053 / .0051 / .0051 | .055 / .005 / .066 | .021 | .113 | .039 |
| 8·10⁴ | 26 | 2.103 | .0039 / .0037 / .0036 | .053 / .003 / .059 | .047 | .041 | .039 |
| 8·10⁴ (P = 400) | 60 | 2.103 | .0039 / .0036 / .0036 | .074 / −.001 / .099 | .006 | .015 | .017 |

**Reading the table.**

- **Random size.** rms of the normalized cross terms ≈ rows^{−1/2}, within 3–17% of random. The coherent excess is the
  predicted mixed diagonal n₁p = n₂p′.
- **No structured phase correlation persists.** The Gauss-phase correlation 0.113 at U = 4·10⁴ (2.9σ) does not
  replicate at 8·10⁴ (1σ, and 0.9σ at P = 400).
- **The Möbius fourth moment is Gaussian** far beyond the diagonal.

## 6. Conditional statements (for the coordinator)

- **[C] (H2-mk-δ).** For all init patterns at r + z ≤ 1 − c: |E_off| ≪ Z^{F+e(r,z)−δ}. Then:
  - lem:marked holds on e(r,z) ≤ δ (Lemma D plus §1);
  - the payoff is in §2: ≈ −0.07δ for small δ; 41/47 needs δ ≥ 0.159 for e_k(slot) / masked coefficients, together
    with the cubic branch.
- **[C] Möbius fourth moment.** Σ_u|M_u(Z^r)|⁴ ≪ U^{1+ε} for r ≤ 0.81 implies lem:marked for all z ≤ 1/4 (Prop X1).
  It also directly implies counts 1 − 2rδ, which would close 41/47 with a large margin. Truth: Gaussian [N]. Status:
  beyond the edge (Lindelöf-on-average fourth moment for a thin family).

## 7. Files

| file | content |
|---|---|
| `xm_model.py` → `xm_model.out` | count LP and high side with cross-saving δ / φ, cubic-branch options, supply options |
| `xm_thresh.py` → `xm_thresh.out` | δ* = 0.159 (relax, unlimited); small-δ slope |
| `orig_side.py` → `orig_side.out` | original-side exact numerics: T(p,p′), m₄/m₂², phase correlations |
| `mixdiag.py` → `mixdiag.out` | mixed-diagonal prediction for the coherent cross ratio |
| `m4count.py` → `m4count.out` | count-based Möbius fourth-moment exponent S₄(r) and the resulting Prop X1 bound |
| `gmodel.py`, `amp_model.py`, `vm.py` | copies of the coordinator, amplitude and marked-cap models |
