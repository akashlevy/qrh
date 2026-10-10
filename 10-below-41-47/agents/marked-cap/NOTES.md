# marked-cap — can the actual slots raise lem:marked's prime capacity? (10 Oct 2026)

Agent `marked-cap`. P78 = `openai/.../September-30-2026/build/paper.tex` (line numbers refer to it).
Labels: **[P]** proved (given P78's lemmas), **[C]** conditional, **[H]** heuristic, **[M]** model/numerics.

## 0. Verdict

* **No capacity gain is proved.** Exponent payoff 0. No exponent that would be cited moves.
* **New structural result [P] (Lemma D, §3).** In the DIAG range r + z < m, the *tuple-diagonal* part of the marked
  canonical energy needs no mark cost at all. Convert each slot tuple into a fourth-power **label** of
  lem:canonical-moment. The second z of lem:marked is therefore paid **only by the cross-tuple terms τ ≠ τ′**.
  The missing input is a cross-mark decorrelation of (H2) type for the lem:marked family. No (HLS) or (H1) input is needed here.
* **Requirement [M].** At the binding bins, the inverse witness has r ≈ 0.66–0.72 and wants z = L ≈ 0.22–0.23, so
  r + 2z ≈ 1.13–1.15. The cross-tuple sum must then save U^{r+2z−1} ≈ U^{0.13–0.15} over the
  Cauchy–Schwarz-trivial bound, i.e. about 60% of full square-root cancellation over pairs of marks (§4).
* **Model, fine grid [M].** The marked lever alone **cannot** reach 41/47: even z = L for every r < 1 leaves +0.0019.
  The BRIEF pair "0.8(1−r) + plain 0.4(1−2m)" also fails on a fine grid (+0.0003); the coarse grid's −0.0001 is a grid artefact.
  Two pairs that do close: plain 0.45(1−2m), or α ≤ 0.78.
* Every suggested option (absorption, Möbius/Mellin, sieve identities, Poisson in slots, primes on the original side,
  labels, explicit formula, other witnesses) is reduced in §5 to a precise P78 step. Each one either already
  happens in P78 or reproduces the same cross-mark terms.

## 1. Requirement (vectorized model `vm.py`, `fine3.py` → `fine3.out`; ℓ = 0.1773, B = 41/47)

The grid is δ ∈ [0.2, 0.75] in steps of 0.01, x ∈ {0, .2, .3, .4, .45, .5}, t has 101 points, r has 201 points, b ∈ [0, 0.1] in steps of 0.005.

| input | worst E |
|---|---|
| paper | +0.01124 (δ = 0.47, x = 1/2, R = 0.590) |
| marked saturated, z_M = L for all r < 1 (≡ DIAG-marked, ≡ mk = 1) | +0.00187 |
| r + 2z ≤ 1.13 / 1.045 (wall-type offset in lem:marked only) | +0.00231 / +0.00794 |
| mk 0.8(1−r) + cP 0.4 (BRIEF) | **+0.00032** (not enough) |
| mk 0.8(1−r) + cP 0.45 | −0.00018 |
| mk saturated + cP 0.5 | −0.00085 |
| mk 0.75(1−r) + α 0.75 (BRIEF) | −0.00088 |
| mk saturated + α 0.78 / 0.76 | −0.00004 / −0.00070 |
| mk saturated + cP = 1 | −0.00828 |

Anatomy (`anat2.py`): in every closing configuration the inverse witness has r ≈ 0.66–0.72. It wants
z_M = L ≈ 0.22–0.23, so r + z ≈ 0.89–0.93 < 1 (the DIAG range) and r + 2z ≈ 1.13–1.15. The plain witness has m ≈ 0.41–0.45.

**Shared wall (`joint2.py`, `joint.py`).** Suppose a lossless column-term saving Z^w were proved for lem:reflected-energy.
It feeds the probe wall (M + ℓ ≤ 1 + w). It also feeds lem:marked: the dual length there is T = 2(M+z) − N, so it gives r + 2z ≤ 1 + w/2.
41/47 then needs w ≈ 0.041 for the wall alone, and w ≈ 0.037 when lem:marked also uses it (≈ 0.033 if the transfer were 1:1).
So lem:marked is only a minor second customer of the wall estimate.

## 2. Exact structure of M_r·Q and where the cost arises (P78 steps)

Slots (l.15000–15005, 15127–15150): Q_i(u) = P_i^{−1/2}Σ_{p∈1_T} χ̄_p(u) W_i(q_p/P_i)(q_p/P_i)^{z−1}.
Relative to the row character ψ = νχ̄_•(u), the coefficient is ν̄(p)1_T(p), a fixed finite combination of θ ∈ Θ.
The inverse witness M_u has the same orientation.

1. **Initialization (l.11606–11660)** writes M·Q as a Möbius polynomial in c = n₀P₀ with μ(c) and the mark
   𝔡₀(c) = Σ_{tuples P₀|c}Π a_i(p_i)μ(p_i). Its length is D′ = r + z − 2G.
   - **The slot primes are already absorbed into μ** (μ(p) = −1 is used): eq. initial-common-extraction
     turns μ(c)γ₋₁(c) into ᾱγ₂(c) on all of c, slot primes included.
   - The residual is the divisor mark 1_{P₀|c}. That is "the support restriction k = n·p₁⋯p_K" of the brief.
2. **Init Poisson (l.11700–11760).** The conductor is z₁z₂, the coprime parts of the two columns, so a slot prime common to both copies
   sits in the gcd C and is assigned outside the square. The first Poisson charges only the length D′ = r + z.
3. **Canonical family (eq. canonical-energy).** Rows K = Z^{2D′−m−2P₁}, columns n with coefficient A_*(n)χ_n(k)χ_n(f)^4𝔡₀(n),
   label f = d′s, F = D′ − P₁.
4. **Terminal (l.9840–9925).** This is the whole-index marked reflection (cor:marked-reflection): an active mark p enters
   the conductor squared, giving T_d = 2M + 2z_a − N. The column branch (old-eq:4.10) needs T_d ≤ F, i.e. M + z₀ ≤ F.
   At the init this reads r + 2z ≤ m.
5. **[P] P78 uses no cross-mark information.** In lem:hybrid-energy at L ≤ N^{1/2} (the branch u = z_a of old-eq:4.7),
   the cubic large sieve over P returns exactly the P-diagonal: the bound N·Σ|a(P)|²/q_P is attained by the P = P′ terms.
   - So P78's terminal bound is, up to Z^ε, the sum of single-mark bounds.
   - It is unchanged if a(P) is replaced by a(P)ε_P with independent random signs.
   - The coefficients enter only through |a(P)| ≤ 1.

## 3. Lemma D — the tuple-diagonal is free in the DIAG range [P]

**Setting.** Fix one outer pattern of the initialization: C, d′, s, the ray summand, the dyads and the Fourier mode.
Let 𝒯 be the set of surviving slot tuples τ = (p_i)_{i∈I₀} (one prime per surviving slot, disjoint windows) with
coefficients a_τ = Π a_i(p_i)μ(p_i), |a_τ| ≤ 1. Write the canonical energy at the init as

E = Z^{−V₀}Σ_{f₀} w₀(f₀) Σ_k |Σ_τ a_τ X_{k,τ}(f₀)|²,
X_{k,τ}(f₀) = Z^{−N/2} Σ_{n: τ|n} A_*(n) χ_n(k) χ_n(f₀)^4 ρ(n) W(q_n/Z^N).

Then E = E_diag + E_off, with E_diag = Z^{−V₀}Σ_{f₀}w₀Σ_τ|a_τ|²Σ_k|X_{k,τ}(f₀)|² and
E_off = Z^{−V₀}Σ_{f₀}w₀Σ_{τ≠τ′}a_τ ā_{τ′}Σ_k X_{k,τ}X̄_{k,τ′}.

**Lemma D.** Assume the *unmarked* canonical invariant at the init parameters:
q_r ≤ Z^{F−M−c} and 3M + c ≤ 4F. At the init this is m − r − z ≥ c + (puncture), together with r + z ≤ 3m/2.
Then E_diag ≪ Z^{F+ε}.

**Proof.**
1. *Factor out the tuple.* For τ | n write n = τn′. Since τ is squarefree and prime to S, d′ and s, the CRT identity
   eq. initial-common-extraction (with d′ = h′ = 1, s = τ) gives
   A_*(τn′) = A_*(τ)A_*(n′)χ_{n′}(τ)^4, including all zeros. Also χ_{τn′}(k) = χ_τ(k)χ_{n′}(k).
2. *Drop the unimodular prefactor.* |A_*(τ)χ_τ(k)| = 1_{(τ,k)=1} ≤ 1, and dropping it is legitimate by positivity of the row sum. Hence
   Σ_k|X_{k,τ}|² ≤ Z^{−N}Σ_k|Σ_{n′}A_*(n′)χ_{n′}(k)χ_{n′}(f₀τ)^4ρ(n′)W(q_τq_{n′}/Z^N)|².
   The coprimality (n′, τ) = 1 is carried by the zero extension of χ_{n′}(τ)^4.
3. *Separate the weight.* Separate W(q_τq_{n′}/Z^N) by logarithmic Fourier inversion in the joint profile, exactly as P78
   does at l.11640. Then use Minkowski in the (f,k) Hilbert space. The factor (q_τ/Z^z)^{iξ} has modulus one, and the
   n′-profile acquires a norm twist with polynomial height cost.
4. *Rewrite as a label sum.* Put f = f₀τ (squarefree), V = V₀ + z, N′ = N − z, so that F₀ := N′ + V = F, and set
   w(f) = Σ_{f=f₀τ}w₀(f₀)|a_τ|² ≤ C d(f)^C. With these choices E_diag is bounded by
   Z^{−N+N′+V−V₀}·Z^{V₀}·[Z^{−V}Σ_f w(f)Σ_k|Z^{−N′/2}Σ_{n′}A_*(n′)ρ(n′)χ_{n′}(k)χ_{n′}(f)^4 W̃(q_{n′}/Z^{N′})|²].
   - The bracket is a canonical energy (eq. canonical-energy) with label f, an empty mark (z₀ = 0) and the same rows and puncture.
   - Other slots, if kept, form a subcollection of the original lists.
5. *Apply the canonical lemma.* lem:canonical-moment with z₀ = 0 applies under the stated invariant, so [·] ≪ Z^{F+ε}.
   The prefactor exponent −N + N′ + V − V₀ = 0 cancels, so E_diag ≪ Z^{F+ε}. ∎

**Remarks.**
1. *Sparse labels are free.* The label set (products of window primes) is sparse, but positivity makes this harmless.
2. *Why the probe cannot use this.* In the probe (lem:probe-row-norm) the target is the row term Z^{M′}, not the column term.
   There the same conversion is useless, and the 09 conclusion "(H1)-average needs (HLS)" stands.
   lem:marked is special because its target is the column length F.
3. *Partial version.* Convert only a subset S₁ of slots and keep S₂ as marks.
   Then E_diag(S₁) ≪ Z^F under M + z(S₂) ≤ F, by the canonical lemma with both a label and a mark.

**Corollary [P].** At (r, z) with r + z ≤ m − c, lem:marked follows from the cross-tuple bound
|E_off(S₁)| ≪ Z^{F+ε} for every init pattern. It is enough to take S₁ of total length z₁ = r + 2z − m and to keep the rest as marks.
P78's hypothesis r + 2z ≤ m is the case z₁ = 0.

## 4. Size of the needed cross-mark saving [P for the arithmetic, H for "truth"]

- **Trivial bounds.** Cauchy–Schwarz gives |E_off| ≤ (#τ − 1)E_diag. In the random model E_diag ≈ K = Z^{M} (the row
  term) and E_off fluctuates around 0.
- **What P78's range means.** With the truth-level diagonal, Cauchy–Schwarz gives E ≲ K·L = Z^{M+z}. This is ≤ Z^F exactly
  when M + z ≤ F, which is P78's range. So **P78's capacity is what "no cross-mark cancellation" buys.**
- **What is needed.** Going beyond requires the normalized L×L cross-mark Gram matrix Ĝ (L ≈ #primes in the converted
  windows) to satisfy a*Ĝa/‖a‖² ≤ Z^{F−M} = U^{1−r−z}, instead of the trivial L = U^{z}.
  - At the binding point (r = 0.69, z = 0.22) that is U^{0.09} versus U^{0.22}: a saving of U^{0.13}.
  - Equivalently |G_{PP′}| must be about U^{0.13} below trivial for typical pairs.
  - Pointwise square-root cancellation in G_{PP′} would give capacity z ≤ min(1/2, 1 − r).
- **Truth [H].** The 09 mark-design numerics on the probe analogue give coherent Rayleigh quotients of 1.02–1.04 against
  1.00 ± 0.02 for random vectors. So the decorrelation is true with a large margin, but it is pure random-matrix
  behaviour (09 agent h2): no main term, nothing to exploit.
- **Where the coefficient structure enters.** Only E_off can see a(p). Expanding θ(p) over its finitely many ray classes reduces to
  a ≡ 1 on one class, by Cauchy–Schwarz over O(1) classes. The phases in G_{PP′} are then the theta coefficients' own
  Gauss sums A_*(P)Ā_*(P′)χ_Pχ̄_{P′}(k), inseparable from the row factors Θ(kP⁴), Θ̄(kP′⁴).
  - The Gauss-sum counterexample a(p) = conj g̃₆(p) (09 diag §3.1) lives exactly in E_off; the diagonal is the same for all |a| = 1.
  - So for ray-character slots "DIAG ⟺ cross-mark decorrelation". Neither primality, product form nor windows enter the diagonal.

## 5. The suggested routes — the exact P78 step each one touches

1. **Absorb slot primes into μ.** P78 already does this: c = n₀P₀ carries μ(c) and A_* absorbs every prime of c.
   - [P] The residual divisor mark cannot be removed. F(c) = Σ_{p|c} w(p)W₁(q_c/(q_pD))W₂(q_p/P) is not a function of
     q_c (or of the smooth data of c). Compare two squarefree c of the same norm, one with and one without a prime in the window; one gives F = 0.
   - The only mark-free identity is the unwindowed Σ_{p|c} log q_p = log q_c (μ·log = −μ∗Λ). It needs no window on
     either factor, and the detector's dyadic windows and the physical slot windows both violate that.
2. **Möbius/Mellin decomposition of the support restriction.**
   - The generating function ∂_tΠ_{p|c}(1 + tθ(p)1_P(p)) gives coefficient −(1+tθ(p)) at window primes, i.e. marks of all orders.
   - Mellin in the window, w = s + τ, gives Σ_{p|c} log p·(P/p)^τ on a vertical line Re τ = σ. The kernel |Ŵ(σ+it)| ≈ 2^σ
     reproduces the cutoff p ≤ 2P, so the window and its conductor return. Nothing changes in eq. initial-overlap-polynomial.
3. **Sieve decompositions (Vaughan/Heath-Brown on p).** The pieces are of two kinds.
   - Möbius pieces: μ∗μ = μ·Σ_{f|c,f~F}1, i.e. divisor marks. Each costs conductor f² in cor:marked-reflection, the same per unit length.
   - Smooth pieces of length E ≤ z < 1: Poisson in e has dual length U/E > E.
   - The worst piece keeps the full mark length, as in 09 diag §3.3.
4. **Poisson in the slot variables.** This needs smooth pieces, which are useless as just shown.
5. **Keep the primes on the original side.** Then |Q|² = Σ_{p,p′}θθ̄χ̄_pχ_{p′}(u), and the off-diagonal is a twisted moment
   of |M_r|² by a character of conductor pp′. The u-Poisson modulus is n₁n₂p₁p₂, which is identical to step 2 of §2.
   The dual Gauss sums g(χ̄_p,h) are re-absorbed by μ(p) = −1, which recreates the same mark.
   Poisson in k returns to the original problem (move calculus).
6. **Marks → labels by Cauchy–Schwarz [P].** E ≤ #τ·E_diag ≪ Z^{F+z} (Lemma D), so Σ_u|MQ|² ≪ U^{1+z}.
   The count exponent becomes 1 − δr + (1 − δ)z, whose slope in z is positive, so this is worse than no slot. Blocks do not help either:
   the cost depends on the prime size, not on the number of primes (09 mark-design).
7. **Explicit formula / zeros.**
   - The pointwise bound |Q|² ≤ P^κ (β*, κ = 3/4) gives slope κ − δ > 0 per unit of slot [P]. Hybrids with the marked part keep this slope.
   - A bootstrap through bin counts is circular [P]. Bounding Σ_u|M_rQ|² bin by bin needs
     R(δ′) ≤ 1 − δ′(r+z) for all δ′, which is the conclusion.
8. **Other witnesses for x = 1/2 rows.**
   - Slot-only witnesses Q^k: capacity kz ≤ 3/8 (lem:marked) or 1/2 (sextic sieve), giving count ≥ 1 − δ/2 [P].
   - Copies u·a⁶ for M_rQ: count slope 5 − 3δ in the copy exponent [P].
   - The log-Möbius M^{log}_{r+z} would be mark-free (length ≤ 1, with no slots, covered by lem:marked), but the detector does not produce it.
     The log-variant of the Γ-integral, with C′ in place of C, gives the pair (M^{log}_r, S_m) with the same adversarial split, and
     its n-coefficient is Λ(n), so it is not a new identity [P].
9. **Disjoint dyadic windows and product form.** These are used only to assign slots to extracted factors. The wall sees
   Σz_i (09). No gain.
10. **Cap 2r + 8z ≤ 3** (the cubic branch old-eq:4.11). It does not bind at r ≈ 0.7: (3−2r)/8 ≈ 0.20 > (1−r)/2.
11. **Primal sextic sieve on the init's canonical family.** This gives the coefficient-blind range r + z ≤ 2m/3, which is better only for r < 1/3. Payoff 0.

## 6. Conditional statements

- **(H2-mk)** For every init pattern of lem:marked with r + z ≤ 1 − c: |E_off| ≪ Z^{F+ε}.
  - [C] (H2-mk) together with Lemma D implies DIAG-marked (capacity min(L, 1−r)). Worst E is then +0.00187, so 41/47 is not reached.
  - [C] Adding plain capacity 0.45(1−2m), or long-witness slope α ≤ 0.78, closes 41/47 (−0.0002 / −0.00004).
  - Only z₁ = r + 2z − 1 ≈ 0.13–0.15 of the slot must be decorrelated; the rest stays a mark.
- **Shared wall W(w).** It reaches 41/47 at w ≈ 0.037 (§1), versus 0.041 without the marked transfer.

## 7. Files

| file | content |
|---|---|
| `gmodel.py` | copy of the coordinator's model |
| `gm2.py` | adds knobs cM (mark cost per unit), wM (offset), zfun |
| `vm.py` | vectorized fine-grid version |
| `fine3.py` → `fine3.out` | the table in §1 |
| `anat.py`, `anat2.py` | binding witness splits |
| `joint.py`, `joint2.py` → `joint*.out` | shared wall (transfer w and w/2) |
