# R1 — Direct proof of the twisted-Möbius mean square below the diagonal (M_κ, κ < 1)

Working notes, 8 Oct 2026. Route key `R1_direct_meansquare`.
Scripts in this directory: `r1_bookkeeping.py` (exact exponent calculus), `r1_mds_group.py` (exact
Cartan matrices), `r1_fake_diagonal.py` + logs `r1_fd_2000.log`, `r1_fd_6000.log` (numerics in Z[ω]),
`r1_validate.py` (validation of the arithmetic), `r1_polar_check.py` (Patterson bias check).

Notation: paper2 = `preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex` (11/12);
paper78 = `.../September-30-2026/build/paper.tex` (7/8). Labels are the TeX labels.
`A_u(D) = Σ_{(n,S)=1} μ(n)ν(n)χ_n(u)W(Nn/D)`, `S(H,D) = Σ_{0<N(u)≤H} |A_u(D)|²`.

## 0. Verdict

**OBSTRUCTION.** I find no proof of M_κ for any κ < 1. I found no weaker variant (restricted or weighted row sets, or
`D^{1+ε}(H + D^{6/7})`) that is provable by the paper2 toolkit and still gives θ ≤ 6/7. The precise reasons:

1. *Exact location of the column term.* In paper2's own reduction, D² reappears at the positivity step of
   Prop. `prop:poisson-reduction`. This is eq. `(weighted)`: Lemma `lem:smooth-mean-square` is applied with
   coefficient μ(f), which takes absolute values of μ(f). The positive dual forms there are *equal to their
   diagonal*. I prove this under GRH for every κ < 1. It also holds unconditionally for κ < 1/6 (κ < 1/4 using the
   7/8 theorem), and numerically with ratio 1.000 ± 0.007. So the hypothesis `(auxiliary-target)` of that
   proposition is **false** for κ < 1, although M_κ itself is true under GRH. The positive-form route is not merely
   lossy. Its intermediate claim fails.
2. *Duality calculus (ideas (i)/(ii)).* Poisson in the row (P) and theta reflection in the column (R) are the only
   moves. They act on (δ, ε) = (target − diagonal, target − column term) by P: swap, and R: ε ↦ 2δ − ε. The
   quantity |ε − δ| = |1 − κ| is invariant. From κ < 1 the reachable configurations form a 6-node path
   A–B–C–C′–B′–A′. On it, min(δ, ε) ≤ −(1 − κ) < 0 at every node (exact, `r1_bookkeeping.py`). The only move
   that would succeed is R at level A, i.e. a column functional equation for μ, which is GRH-type information. The
   second moment that must be evaluated has an **affine** (infinite) group of functional equations: its Cartan
   matrix has determinant 0 (`r1_mds_group.py`). This is the multiple-Dirichlet-series signature of moments
   "beyond the edge".
3. *Cubic large-sieve sharpness (iii) is not the obstacle.* μ-absorption already removes the polar
   (Patterson-bias) terms that make the cubic sieve sharp (dFDH 2607.07911). That removal is what makes paper2
   work for κ > 1. The κ < 1 obstruction is the trivial column term Ξ ≥ B. No operator-norm bound can remove it.
   M_κ is a statement about the specific vector μν.
4. *R1's target is strictly stronger than the goal.* M_{6/7} implies zero-freeness in Re s > 6/7. But even
   "every Hecke L-function over K is zero-free in Re s > 6/7" gives only S ≲ H·D^{12/7}, not H·D. M_{6/7}
   additionally encodes a **beyond-density-hypothesis** zero-density estimate for the sextic family. GRH for all
   *non-principal* rows does not help R1 at all (§6(d)).
5. *The "weights favouring sixth powers" variant* is a prime-amplified mean square Σ|A_u|²|Σ_{p∼P}χ_p(u)|². It is
   paper78's Lemma `lem:marked` with one slot. With that lemma's proved walls (r + 2z ≤ m, 2r + 8z ≤ 3m) the
   extraction gives exactly 11/12, because amplification gains nothing. With the prior W(w) relaxation it gives at
   best 0.900. It gives 6/7 only if the marked moment holds out to the canonical diagonal at P = D^{5/7},
   H = D^{12/7}, i.e. (r, z) = (7/12, 5/12). That is true under GRH but far outside the proved walls; it is a
   marked-energy problem, not a column-term problem (§6(i), (C4)).

**Best θ via R1:** unconditionally 11/12 (paper2, κ = 1 + ϑ); 7/8 needs the separate Part II architecture.
**Conditionally:** M_{6/7}, even averaged over D, gives exactly θ = 6/7 (§7). This hypothesis cannot be weakened to
statements about non-principal rows.

---

## 1. What R1 must prove, and that there is no slack

(a) **Extraction (paper2 §3, eq. `prime-extract`).** Take rows u = p⁶ with Y/2 < N(p) ≤ Y and Y = H^{1/6}. Then
`|A_1 − A_{p⁶}| ≪ D/Y`. Averaging over J ≍ Y/log Y primes gives
`|A_1|² ≤ (2/J) Σ_p |A_{p⁶}|² + O(D²/Y²) ≤ (2/J) S(H,D) + O(D²/Y²)`.
With M_κ (S ≪ D^{1+κ+ε}) and H = D^κ, this is `|A_1|² ≪ D^{1+5κ/6+ε} + D^{2−κ/3}`, so
θ_prime(κ) = max((6+5κ)/12, 1 − κ/6). This equals **6/7 exactly at κ = 6/7**, which is also its minimum over κ.
Paper78 Lemma `lem:inverse-amplification` at U ≍ 1 uses the exact identity `sixth-power-column-identity`. It
gives θ_int(κ) = (6+5κ)/12, again 6/7 at κ = 6/7 (checked exactly in `r1_bookkeeping.py`, §6).

(b) **No slack.** Under GRH the truth is S(H,D) ≍ HD. The prime extraction needs H ≥ D^{6/7} (second term) and
S(H) ≤ H^{1/6}D^{12/7} (first term). With S ≍ HD both hold only at H = D^{6/7}. So 6/7 is the GRH floor of that
extraction, and R1 needs M_{6/7} on the nose.

(c) **D-averaged form suffices.** If `∫_D^{2D} S(D^{6/7}, y) dy/y ≪ D^{13/7+ε}`, the same extraction applied to
each y gives `∫_D^{2D}|A_1(y)|² dy/y ≪ D^{12/7+ε}`. Cauchy–Schwarz on dyadic blocks then makes
`∫ A_1(y) y^{−s} dy/y` absolutely convergent for Re s > 6/7. The Mellin argument of paper2 §3 then rules out zeros.
The averaging does not make the problem easier (Plancherel turns it into ∫|1/L(1/2+it)|², i.e. GRH-type, §6(f)).

(d) **Restricted or weighted row sets (variants asked for).**
- A bound `S ≪ D^{1+ε}(H + D^{6/7})` is M_{6/7}, by monotonicity of S in H.
- A Poisson-friendly row set must be lattice-like. Its sixth-power density is ≈ H^{−5/6}, the same as the ball.
- A sparse set {a⁶t : N(a) ≤ P, t ∈ T} has S = P·Σ_{t∈T}|A_t|²·(1+o(1)), since A_{a⁶t} ≈ A_t. The extraction then
  needs Σ_{t∈T}|A_t|² ≤ D^{12/7}. For T = {N(t) ≤ T₀} that is M at κ = 5/7, which is harder, not easier.
- A smooth weight concentrated near sixth powers has a Fourier transform given by degree-6 Weyl sums (§6(f)).
So no variant escapes the need for M_{6/7}.

## 2. The exact dual identity and where D² comes from (ideas (i), step 1)

Paper2's proof of Prop. `prop:poisson-reduction` proceeds as follows. Expand the square, set g = (n₁,n₂) and
n_j = g z_j. Apply Lemma `lem:poisson` with χ = χ̄_{z₁}χ_{z₂} (primitive mod z₁z₂) and exclusion g. Apply μ-absorption
`(convert2)`: μ(n)γ_{−1}(n) = χ_n(−1)G(n)^{−1}ᾱ(n)γ₂(n). Remove the coprimality (z₁,z₂) = 1 by Σ_{v|z₁,z₂}μ(v).
Change variables (b, f, k). This gives eq. `initial-column-output`:

```
S_ξ = (H/D²) Σ_{(b,f)=1} μ(f) N(b) Σ_{k≠0} Σ_{m1,m2} a(m1) conj(a(m2)) χ_{m1}(k f⁴) conj(χ_{m2}(k f⁴))
          · W0(N(bf m1)/D) conj(W0(N(bf m2)/D)) · Φ^(H N(k) / (N(f)² N(m1) N(m2)))
```

Here a = ā·γ₂·ξ with |a| = 1. Block scales are X = D/(BF), 𝓗 = C·D²/(HB²), XF = D/B.

**Fake diagonal and its exact cancellation.** The terms m₁ = m₂ = m contribute, via Lemma `lem:poisson` with the
principal character and exclusion m,
`Σ_{k≠0,(k,m)=1} Φ^(HNk/(N(f)²N(m)²)) = (N(f)²N(m)²/H)·Φ(0)·φ(m)/N(m) − Φ^(0)·[m = 1] + (Schwartz tail)`.
Put j = fm. The multiplicative identity
`Σ_{f|j} μ(f) N(f)² N(j/f) φ(j/f) = μ(j) N(j)` holds (at p | j: q(q−1) − q² = −q).
So the Φ(0) parts sum to `(Φ(0)/D²) Σ_b N(b) Σ_j μ(j) N(j) |W0(N(bj)/D)|²`, a Möbius sum over the gcd variable.
- Using only the zero-free region θ₀ of ζ_K, this is ≪ D^{ε} in the 𝓜_D normalisation, i.e. O(D^{1+ε}) for S.
- The Φ^(0) parts give O(HD) for S.
So **in the exact identity the fake diagonals cancel down to O(D^{1+ε}H)**. What remains is the signed off-diagonal
Σ_f μ(f)·Off(b,f).

**The positivity step reinstates them.** Eq. `(weighted)` bounds each (b,f)-block by the positive form
`Σ_{f,k}|S_U(b,f,k)|² = XF·𝓔_{(b)}(𝓗,X,F)` and then sums |μ(f)|. If 𝓔 ≍ 𝓗 (Prop. O1 below), block (B,F)
contributes `(HB²/D²)·F𝓗X = D/B` to 𝓜_D. Summed over B this is ≍ D, i.e. **S ≫ D² from this bound**: the
large-sieve column term. The proposition's hypothesis `(auxiliary-target)` 𝓔 ≪ D^ε XF requires 𝓗/X ≤ D^{ε},
i.e. D/H ≤ D^{ε} at B = F = 1. It is false for κ < 1. Equivalently, Prop. `prop:canonical`'s `(positive-gap)`
𝓗 ≤ Σ D^{−κ'} fails.

**Proposition O1 (the dual target is false below the diagonal).** Fix ξ, a real W ∈ C_c^∞((1,2)) with W ≢ 0, and
κ ∈ (0,1). Put 𝓗 = X^{2−κ} and 𝓔 = 𝓔(𝓗, X, 1; ξ, W) (paper2 eq. `energy`, F = 1).
- (a) Under GRH for finite-order Hecke L-functions over K, 𝓔 ≫_W 𝓗. Hence 𝓔/X ≫ X^{1−κ} → ∞.
- (b) Unconditionally, (a) holds whenever 2 − κ > 2θ₀, where θ₀ is a proven zero-free abscissa for all Hecke
  L-functions over K. That is κ < 1/6 with θ₀ = 11/12 (paper2), or κ < 1/4 with θ₀ = 7/8 (paper78).

*Proof.*
1. Take radial Φ ∈ C_c^∞([0,1)) with 0 ≤ Φ ≤ 1 and Φ ≢ 0. The k = 0 term vanishes because χ_n(0) = 0 for n ≠ 1.
   So 𝓔 ≥ 𝓐_Φ := (1/X) Σ_k Φ(Nk/𝓗)|C(k)|², where C(k) is the column sum.
2. Evaluate 𝓐_Φ exactly as in paper2 Lemma `lem:first-transfer`, with F = 1 and L = X. That computation uses only
   Lemma `lem:poisson`, `(convert1)`, `(quotient)` and CRT. Φ^ is Schwartz; compact support was only used to
   truncate the frequency range.
3. The zero-frequency term is Z = (𝓗/X)Φ^(0) Σ_C |W(N C/X)|² Π_{p|C}(1 − 1/Np) ≍ 𝓗, since Φ^(0) = ∫Φ > 0.
4. Every other term is w_{C,d} Σ_h P(h;U₁)·conj(P(h;U₂)) (after the Mellin separation of Lemma `lem:smooth`). Here
   P = P_{C,d,t} is the Möbius sum of eq. `(P)`. It is the Möbius sum of a finite-order Hecke character of
   conductor ≪ N(hCd), with length ℓ = X/(N(C)N(t)).
5. Under GRH, |P| ≪ ℓ^{1/2}(ℓN(h))^{ε}. (Contour shift to Re s = 1/2 + ε; |1/L| ≪ (Q(|t|+2))^ε there by
   IK Thm 5.19-type bounds.) The number of h is ≪ X²N(d)/(𝓗N(C)²) + O(1) (Schwartz tail), and
   w_{C,d} = 𝓗N(C)/(N(d)X²).
6. So each (C,t,d) contributes ≪ X^{1+ε}/(N(C)²N(t)). Summing gives |𝓐_Φ − Z| ≪ X^{1+ε} = o(𝓗).
7. For (b), use |P| ≪ ℓ^{θ₀+ε} instead. The sum becomes ≪ X^{2θ₀+ε}, which is o(𝓗) iff 2 − κ > 2θ₀. ∎

Since GRH also implies M_κ (trivially, |A_u| ≪ D^{1/2+ε}), O1(a) shows the following. **For κ < 1 the hypothesis of
Prop. `prop:poisson-reduction` is strictly stronger than its conclusion, and in fact false.** The loss happens at
the absolute values on μ(f) in eq. `(weighted)`.

**Numerics** (`r1_fake_diagonal.py`, exact sextic symbols and γ₂ in Z[ω]; validated in `r1_validate.py`:
0 symbol mismatches; γ₂ twisted multiplicativity matches brute force to 1e−13). Definitions:
- B_h(X) = Σ a(n)χ_n(h)W(Nn/X) with a = āγ₂, over squarefree primary n with (n,6) = 1 and X < Nn ≤ 2X.
- Q_v = Σ_{0<Nh≤𝓗} |Σ_{v|n} a(n)χ_n(h)W|², with 𝓗 = X^{2−κ}.
- Q₁ = the positive dual form, which is what paper2 must bound.
- Σ_v Q_v = the |μ|-weighted version that paper2 effectively bounds.
- Q_cop = Σ_v μ(v)Q_v = the true coprime-pair form, which has no diagonal.

X = 6000 (1695 columns, 2173 divisors v):

| κ | 𝓗 | Q₁/diag | Q₁/X² | Σ_v Q_v/X² | Q_cop/X² | A-level S(H)/diag |
|---|---|---|---|---|---|---|
| 0.5 | 464758 | 1.000 | 28.21 | 101.9 | +0.003 | 0.87 |
| 0.6 | 194721 | 0.997 | 11.79 | 42.7 | −0.030 | 0.93 |
| 0.7 | 81582 | 1.000 | 4.95 | 17.9 | +0.003 | 0.95 |
| 6/7 | 20792 | 1.006 | 1.27 | 4.57 | +0.008 | 1.00 |
| 0.95 | 9270 | 1.007 | 0.57 | 2.04 | +0.005 | 1.00 |
| 1.0 | 6000 | 1.005 | 0.37 | 1.32 | +0.002 | 1.00 |

Readings:
- Q₁ equals its diagonal ≍ 0.364·𝓗X. The κ = 1/2 row scales as X^{1/2}: Q₁/X² = 6.97, 16.26, 28.21 at
  X = 400, 2000, 6000. So the positive target X² fails by X^{1−κ}.
- The signed form Q_cop is 10⁻⁴–10⁻² of Q₁ and far below X². **The cancellation needed for M_κ is real, but it
  lives entirely in the signs μ(v).**
- The A-level ratio S(H)/diag ≈ 0.9–1.0 shows that M_κ itself looks true. (Compare `strengthening/mobius_meansquare.py`.)

## 3. Idea (i): keep Φ^ un-separated, use Φ(0) = 0, and theta-reflect (double reflection)

(i-1) **Φ(0) = 0 (dyadic rows).** Then ∫Φ^ = Φ(0) = 0, and each block's fake diagonal is O(1) instead of
N(f)²N(m)²/H. But Φ^ changes sign, so no positive arrangement exists. Writing Φ^ = Φ^₊ − Φ^₋ gives two positive
forms, each with fake diagonal ∝ ∫Φ^_± ≠ 0. Bounding them separately reinstates D². To use the cancellation one must
evaluate each positive form *asymptotically*, with error ≤ X² against a main term 𝓗X. That is a power-saving
asymptotic.

(i-2) **Reflect both column sums, then apply the quadratic large sieve to the coupled form.**
- Mellin-separate Φ^(HNk/(Nm₁Nm₂)) = ∫Φ̃(s)(HNk)^{−s}(Nm₁Nm₂)^{s} ds. This gives bilinear forms
  Σ_k (HNk)^{−s} B_k(W₀x^s)·conj(B_k(W₀x^{s̄})).
- Reflect each B_k by Prop. `lem:reflection`. The twist j = 1 becomes quadratic (j = 3), the column length X
  becomes N(k)²/X, and the coefficients become d(ℓ)α(ℓ)/√Nℓ.
- The quadratic large sieve (Lemma `lem:quadratic`, Goldmakher–Louvel Thm 1.1) gives Prop. `prop:R`'s
  `𝓗 + 𝓗²/X` per unit coefficient mass. Restoring the factors,
  `S ≪ (H/D)·X·(𝓗 + 𝓗²/X) = D² + D³/H = D^{max(2, 3−κ)}`.
- At κ = 6/7 this is D^{15/7}, worse than the trivial D². **The column term reappears as the reflected column
  term 𝓗²/X of eq. `(R)`.**
- Φ(0) = 0 kills only the ℓ₁ = ℓ₂ diagonal, which is the image of the fake diagonal. The large sieve bound is
  insensitive to this.

(i-3) **Asymptotic evaluation instead of a large sieve.**
- At level C (rows 𝓗, reflected columns M = 𝓗²/X), M_κ requires `Σ_h |T_h|² = (main term ≍ 𝓗) + O(X)`. That is a
  relative power saving `𝓗^{−(1−κ)/(2−κ)}`, which is **𝓗^{−1/8} at κ = 6/7**. The object is the second moment of
  quadratic twists T_h of the cubic theta function at unbalanced length X = 𝓗^{7/8}.
- The off-diagonal ℓ₁ ≠ ℓ₂, by quadratic Poisson in h, becomes level C′: rows J = M²/𝓗 = D^{4−3κ}, columns
  M = D^{3−2κ}, with coefficients d(ℓ)/Nℓ. C′ is in the diagonal regime, but its target is
  D^{2(1−κ)} below its diagonal.
- Reflecting C′ (j = 3 → j = 1) gives level B′ (rows D^{4−3κ}, columns D^{5−4κ}). Sextic Poisson then gives A′
  (rows D^{6−5κ}, columns D^{5−4κ}). That is a μ-mean square in paper2's regime, but with target D^{8−6κ}, which is
  3(1−κ) *below its diagonal* D^{11−9κ}.
- A′ admits only P, which leads back. The reachable graph is closed.

(i-4) **The calculus (exact; `r1_bookkeeping.py`).**
- State: δ = t − (r + c) and ε = t − 2c, where r and c are the log-rows and log-columns of a unimodular-coefficient
  form and t is the log-target. A sharp large sieve suffices iff min(δ, ε) ≥ 0. This generously grants an
  *optimal* large sieve at every level.
- P: (r, c, t) ↦ (2c − r, c, t − r + c). This acts by (δ, ε) ↦ (ε, δ).
- R: (r, c, t) ↦ (r, 2r − c, t + 2r − 2c). This acts by (δ, ε) ↦ (δ, 2δ − ε).
- Both are involutions, and |ε − δ| = |r − c| = |1 − κ| is invariant. R is defined only on Gauss/theta
  coefficients (levels B, C, C′, B′); P is defined everywhere.
- Positivity moves (row enlargement, Cauchy–Schwarz, separating weights) only lower δ and ε.

κ = 6/7 (exact):

| node | type | r | c | t | δ | ε |
|---|---|---|---|---|---|---|
| A | μ | 6/7 | 1 | 13/7 | 0 | −1/7 |
| B | Gauss | 8/7 | 1 | 2 | −1/7 | 0 |
| C | θ, quadratic | 8/7 | 9/7 | 16/7 | −1/7 | −2/7 |
| C′ | θ, quadratic | 10/7 | 9/7 | 17/7 | −2/7 | −1/7 |
| B′ | Gauss | 10/7 | 11/7 | 19/7 | −2/7 | −3/7 |
| A′ | μ | 12/7 | 11/7 | 20/7 | −3/7 | −2/7 |

- For every κ < 1 tested (1/100, 1/2, 6/7, 99/100) the maximum over the path of min(δ, ε) is exactly −(1−κ).
- For κ = 11/10 the calculus reproduces paper2: min ≥ 0 from node C on, which is where Goldmakher–Louvel applies.
- **The single move that would succeed is R at A**: (0, −ρ) ↦ (0, +ρ). That is a column functional equation for
  μν, i.e. analytic information on 1/L(s, νψ_u). It is not available.
- In matrix form PR = [[3, −2], [2, −1]] (trace 2, det 1) is parabolic, with fixed line r = c (κ = 1). Iterating
  never reaches the trivially-diagonal regime rows ≥ cols².

(i-5) **Multiple-Dirichlet-series signature (heuristic, exact algebra in `r1_mds_group.py`).**
- For Z(s,w) = Σ_h N(h)^{−w} L(s, θ⊗ψ_h), the functional equations R: (s′,w′) ↦ (−s′, w′+2s′) and
  P: (s′,w′) ↦ (s′+w′, −w′) have Cartan matrix [[2,−2],[−1,2]] (det 2, finite type C₂). RP has order 4, so the
  group is the dihedral group of order 8. **The first moment is finite type.**
- For the second moment Σ_h N(h)^{−w} L(s₁, θ⊗ψ_h) L(s₂, θ̄⊗ψ_h), the Cartan matrix is
  [[2,0,−2],[0,2,−2],[−1,−1,2]] with **det 0**. This is affine type (D₃^{(2)}, the transpose of C₂^{(1)}). On the
  diagonal slice R₁R₂·P is unipotent, the same parabolic as the scale calculus.
- Affine (infinite) groups of functional equations are exactly the "beyond the edge" moments, e.g. the fourth
  moment of quadratic Dirichlet L-functions (D₄^{(1)}, Bucur–Diaconu over F_q(t)). There MDS methods give no
  continuation to all of ℂⁿ, and power-saving asymptotics are not available.
- This is a heuristic, not a theorem. Γ-factors, non-squarefree terms and the lack of an Euler product for θ are
  ignored. It is consistent with (i-4).

**Exponent outcome for (i):** no κ < 1.
- Positivity arrangement: S ≪ D^{max(1+κ, 2)}, i.e. D² at eq. `(weighted)`.
- Reflected and coupled arrangement: S ≪ D^{max(2, 3−κ)} (D^{15/7} at κ = 6/7), from the term 𝓗²/X of eq. `(R)`.
- Asymptotic arrangement: needs a 𝓗^{−1/8} power saving in an affine-type second moment, which is open.

## 4. Idea (ii): μ(n)χ_n(u) = G(n)^{−1} ᾱ(n) γ₂(n) γ₁(n; u⁵)

The identity is verified as follows. γ₁(n; u⁵) = χ̄_n(u⁵)γ₁(n) = χ_n(u)γ₁(n), and both sides vanish if (u,n) ≠ 1.
Combine this with `(convert1)`: ᾱγ₂γ₁ = μG. So A_u = Σ_n β(n)γ₁(n; u⁵)W(Nn/D) with β = νG^{−1}āγ₂, a twisted
cubic-theta coefficient. Possible uses:
- (ii-a) Cauchy–Schwarz over n: |A_u|² ≤ ‖βW‖²·‖γ₁(·;u⁵)W‖² ≍ D². Trivial (κ_eff = 1, θ = 11/12 at best).
- (ii-b) Reflect β with γ₁(n;u⁵) as the twist. This is impossible: n ↦ γ₁(n) is not a character (it is a
  twisted-multiplicative sextic Gauss sum), and βγ₁ = μν. Any reflection of the n-sum would be a reflection of μ.
- (ii-c) View γ₁(n; r) as the n-th Dirichlet coefficient of the sextic Kubota series ψ₆(r, s) = Σ_c g₆(r,c)Nc^{−s}.
  Then A_u is a Rankin–Selberg-type pairing of the cubic theta function (3-fold cover) with a Fourier coefficient
  of the sextic Eisenstein series (6-fold cover) at index u⁵. By the identity this pairing *is* 1/L(s, νψ_u).
  There is no known structure for such mixed-cover pairings.
- (ii-d) The mean square over u: γ₁(n₁;u⁵)·conj(γ₁(n₂;u⁵)) = γ₁(n₁)conj(γ₁(n₂))·χ_{n₁}(u)χ̄_{n₂}(u). So Poisson in u
  is literally step 1 of (i).

**Outcome for (ii):** either trivial (D², which gives 11/12), or identical to (i) and blocked at the same steps.
Using u⁵ ↦ ū (inverse mod n) gives Kloosterman-type phases e(xū/n). The large sieve in u then re-indexes x = yu,
which is the same object again.

## 5. Idea (iii): Heath-Brown's cubic large sieve is sharp (de Faveri–Dunn–Hoffstein, arXiv 2607.07911)

- dFDH Thm 1.1 (unconditional): Ξ₃(A,B) ≫ (AB)^{−ε}(A + B + (AB)^{2/3}). The extremal vector is β_b = g̃₃(b), via
  the polar term P_B(a) ≈ g̃₃(a)N(a)^{−1/6}B^{5/6} (Patterson residue at 5/6). Their Conj. (1.8) for n = 6 reads
  Ξ₆(A,B) ≍ A + B + A^{5/6}B^{1/3} + A^{1/3}B^{5/6}.
- **What μ's structure does.** μ-absorption turns the sextic problem into ā·γ₂-coefficients with sextic twist. These
  reflect to *quadratic* twists, and the twisted theta Θ_k has vanishing constant term (paper2 Lemma
  `lem:theta-realization`). So no polar term appears; Prop. `prop:R` has none. Numerically (`r1_polar_check.py`):
  |Σ γ₂W| = 48, 109, 222 at X = 2000, 6000, 15000, which is Patterson bias of rough growth X^{0.76–0.8} toward
  X^{5/6}. In contrast |Σ āγ₂W| = 26, 17, 0.5, with no bias.
- This is why paper2 can beat large-sieve bounds for κ > 1. With H = D^{1+ϑ}, *if* dFDH Conj. (1.8) holds for
  n = 6 (conjectural; dFDH prove lower bounds only for n = 3, 4), then Ξ₆ ≍ H^{5/6}D^{1/3} = D^{7/6+5ϑ/6} > H for
  ϑ < 1. General-coefficient methods would then stop at S ≪ D^{13/6+5ϑ/6} > D^{2+ϑ}. Paper2's Prop. `thm:ms` is a
  statement about the specific vector μν, not a large-sieve statement.
- **Why it does not help for κ < 1.** The obstacle is the column term B, which satisfies Ξ ≥ B trivially (dFDH
  (1.1)). It is attained by β_b = conj(χ_b(u₀)) for any single row. M_κ says that μν is nearly orthogonal to the
  H-dimensional span of the row vectors (χ_n(u))_n; its projection has the "random" size (H/D)‖μν‖². That is a
  property of the vector, outside any operator-norm (large sieve) argument.
- At H = D^{6/7}: Ξ₆ ≥ D (trivially), and conjecturally Ξ₆ ≍ D^{47/42}. Target per unit mass: H = D^{6/7}.

**Outcome for (iii):** no κ < 1. The dFDH mechanism is already neutralised by μ-absorption and is not the bottleneck.

## 6. Idea (iv): other mechanisms checked

- (a) **Heath-Brown/Vaughan identity + Hölder.** Type II pieces need M, N ≤ √H with MN = D > H, which is
  impossible. Fourth moments give `√((H+M²)M²(H+N²)N²) ≥ MN·√(D²) = D²`. Column term again.
- (b) **Zero-free region + large sieve / Halász–Montgomery.** S ≤ min(H·D^{2θ₀+ε}, D^{2+ε}) with θ₀ = 7/8. This
  never reaches HD. Large-values estimates cannot separate the sixth-power rows, which are all ≈ the principal row.
- (c) **Other families.**
  - Cubic family: μγ̄₂ = āγ₂², which is not a theta coefficient, so there is no absorption.
  - Quadratic family: μ·(quadratic Gauss sum) = μ·(unit), so the dual is μ-type again. The square rows give
    √H|A₁|², which is circular.
  - Pair family χ_n(u)χ̄_n(u′): the multiplicity ≈ H of the trivial character reinstates the D² column term.
  - Sublattices, residue classes, Größencharacter (t-aspect) hybrids: the density of copies of A₁ does not
    improve, while the family size grows.
  - General principle: |A₁|² ≲ (|𝓕| + D)·D/#copies, with #copies ≤ |𝓕|^{1/6} for sextic families. So θ ≥ 11/12
    unless the column term is beaten.
- (d) **GRH for all non-principal rows does not help R1.**
  - Write S(H) = Σ_{non-principal} + Σ_{N(v)≤H^{1/6}}|A_{v⁶}|² + (the 6 unit classes). The hypothesis controls
    only the first sum.
  - The extraction needs an *upper bound on S(H)* that is independent of A₁. For H < D the only one available is
    the large-sieve D^{2+ε}. That gives H^{1/6}|A₁|² ≤ D^{2+ε}, i.e. θ ≥ 11/12 (at H = D).
  - The exact Poisson identity does not separate principal from non-principal rows. Applied twice it is the
    identity. Schematically it yields only a consistency relation of the shape
    0 = (H/D)Σ_v μ(v)MT_v + Σ_{v≠1} μ(v)(S_v − diag_v), where the S_v are v⁴-twisted copies of S at shifted
    scales. That relation carries no upper bound for S.
  - So every hypothesis implying M_κ must directly control the principal row with weight H^{1/6}.
  - (In the paper78 probe architecture GRH-for-rows does help, up to 13/15. That is a different route.)
- (e) **Asymptotic large sieve (Hilbert / Montgomery–Vaughan type).** Its error term is the operator norm of the
  off-diagonal kernel K(n₁,n₂) = Σ_h χ_{n₁}χ̄_{n₂}(h) = (𝓗/X)·(A-level operator with rows X²/𝓗). That norm is
  ≍ 𝓗 (the dual column term), as large as the main term. Applied to the vector a it is small iff M_κ. Circular.
- (f) **Averaging over D and Plancherel.** ∫_D^{2D} S dy/y ≤ D^{2σ} Σ_u ∫|Ŵ(σ+it)|²|L(σ+it,·)|^{−2} dt requires
  1/L to be analytic on Re s ≥ σ with σ = 1/2 (needed to reach HD). That is GRH. Weights concentrated on sixth
  powers have Fourier transforms Σ_v e(v⁶ξ) (degree-6 Weyl sums), which give only Weyl-type savings.
- (g) **Mollifier / amplifier inside the mean square.** A_u·L_u^{(y)} is a coefficient sum supported on n > y of
  length Dy, which is worse. Sixth powers are the only available amplifier.
- (h) **The probe architecture of paper78.** Out of scope for R1 (other routes). Note that M_{6/7} would also
  replace the paper78 raw moment `(raw-moment)` at H = D^{6/7}.
- (i) **Weights favouring sixth powers via a prime amplifier (the natural "weighted" variant; `r1_amplified.py`).**
  - Weight the rows by |M_u|², with M_u = Σ_{p∼P} χ_p(u). Then |M_{v⁶}|² ≈ P² on sixth powers versus ≈ P on
    average, so sixth powers are favoured by the factor P.
  - In the normalisation of paper78 Lemma `lem:marked` (rows Z^m, inverse length Z^r = D, one slot Z^z = P,
    m = 1), a bound S_amp := Σ_{N u≤Z}|M_u(Z^r)Q_u|² ≪ Z^{1+ε} together with the rows u = q⁶ gives
    `|A₁|² ≪ D·Z/(YP) + D²/Y²` with Y = Z^{1/6}. Hence θ = max((r + 5/6 − z)/r, (2r − 1/3)/r)/2.
  - Optimising over the admissible (r, z) (z ≤ r), exact:

    | admissible range | best θ |
    |---|---|
    | Lemma marked *as proved* (r + 2z ≤ 1, 2r + 8z ≤ 3) | **11/12**, at z = 0 (amplification useless) |
    | prior W(w) capacity z ≤ min((1+w/2−r)/2, 1−r, (3−2r)/8) | 0.9123 (w = 1/10), 0.9074 (w = 1/5), 0.9000 (w = 1/2) |
    | case-B wall dropped (r + z ≤ 1, 2r + 8z ≤ 3) | 9/10, at (5/6, 1/6) |
    | canonical diagonal only (r + z ≤ 1): "lossless" | **6/7 at (r, z) = (7/12, 5/12)**, i.e. P = D^{5/7}, H = DP = D^{12/7}; 5/6 at P = D |

  - The lossless one-slot marked moment at (7/12, 5/12) is **true under GRH**. There |M_u|² ≪ D^ε. Moreover
    Σ_u|Q_u|² = Z(1 + o(1)): P² = Z^{5/6} < Z, so Poisson in u leaves only the diagonal p₁ = p₂. It lies on the good side of the canonical diagonal
    (r + z = 1), so the below-diagonal obstruction of §3 does not apply to it.
  - But it lies far outside Lemma marked's proved range: r + 2z = 17/12 (needs ≤ 1) and 2r + 8z = 9/2 (needs ≤ 3).
  - Its obstruction is the marked-energy wall of the paper78 architecture (other routes). In the crude P/R
    calculus with marks, reflecting B_{kp⁴} costs N(p), which makes the reflected column term P·𝓗²P/D. Even with an
    optimal large sieve this needs 2r + 5z ≤ 2m, which fails at (7/12, 5/12): 13/4 > 2.
  - Like M_κ, it contains the principal row: here with weight H^{1/6}P.

## 7. Conditional results (sharply stated)

**(C1) Theorem (conditional).** Suppose that for every finite-order Hecke character ν of K, every W ∈ C_c^∞((0,∞)),
and every ε > 0,
`∫_D^{2D} Σ_{0<N(u)≤D^{6/7}} |A_u(y)|² dy/y ≪_{ν,W,ε} D^{13/7+ε}`  (D ≥ 2).
This is implied by the pointwise-in-D M_{6/7}. Then every finite-order Hecke L-function over K, hence ζ and every
Dirichlet L-function, is zero-free in Re s > 6/7.

*Proof.*
1. Fix D and take H = D^{6/7} and Y = H^{1/6} = D^{1/7}. For y ∈ [D,2D] and each of the J ≍ Y/log Y primes p with
   Y/2 < N(p) ≤ Y, (p,S) = 1, we have |A_1(y) − A_{p⁶}(y)| ≪_W y/Y (paper2 §3). The rows p⁶ are distinct with
   N(p⁶) ≤ H.
2. Hence |A₁(y)|² ≤ (2/J)S(H,y) + O(D²/Y²). Integrate dy/y over [D,2D]:
   `∫|A₁|² ≪ D^{13/7+ε}·D^{−1/7} + D^{2−2/7} = D^{12/7+ε}`.
3. Cauchy–Schwarz on each dyadic block gives ∫_D^{2D}|A₁(y)|y^{−σ}dy/y ≪ D^{6/7+ε−σ}. This is summable over
   dyadic D for σ > 6/7.
4. So 𝓜_W(s) = ∫A₁(D)D^{−s}dD/D is holomorphic on Re s > 6/7 and equals Ŵ(s)/L_K^S(s,ν) there. The choice
   W = y^{−ϱ}φ(y) of paper2 §3 yields the contradiction 0 = Ŵ(ϱ) > 0 at a putative zero ϱ.
5. The Dirichlet case follows by base change, as in paper2. ∎

Exponent check: (6+5κ)/12 = 1 − κ/6 = 6/7 at κ = 6/7 (exact). The integer extraction of paper78 needs M at all
H ≥ D'^{6/7}, D′ ≤ D (the sup over scales).

**(C2) Zero-density form.**
- Let N_ν(σ;Q) = #{u : 0 < N(u) ≤ Q, L(s,(νψ_u)*) has a zero ρ with Re ρ ≥ σ and |Im ρ| ≤ Q^{ε}}. Rows are counted
  individually, so the principal-type rows u ∈ μ₆·K^{×6} count with multiplicity ≍ Q^{1/6}.
- Here ψ_u(𝔭) = (u/𝔭)₆ is a finite-order Hecke character by Artin reciprocity for K(u^{1/6})/K, with conductor
  dividing 6^a·rad(u).
- **Claim:** if N_ν(σ;Q) ≪ Q^{1−(2σ−1)/κ+ε} for σ ∈ [1/2, θ₀], then M_κ holds.
  1. Shift the contour of A_u(D) = (2πi)^{−1}∫Ŵ(s)D^s E_u(s)/L(s,(νψ_u)*) ds to Re s = β_u + ε on |t| ≤ D^{ε'}.
     Here β_u is the largest real part of a zero in |t| ≤ D^{ε'} + 2. Use |1/L| ≪ (N(u)(|t|+2))^ε there
     (Borel–Carathéodory).
  2. Shift to Re s = θ₀ + ε on |t| > D^{ε'}, using the proven half-plane θ₀. The Ŵ-decay makes that part
     ≪ D^{1/2}.
  3. So |A_u(D)| ≪ D^{β_u+ε}N(u)^ε + D^{1/2}, and S ≪ Σ_σ D^{2σ}N_ν(σ;H) ≪ D^{1+κ+ε} at H = D^κ.
- At κ = 6/7 this is N ≪ Q^{(13−14σ)/6}, **strictly stronger than the density hypothesis Q^{2(1−σ)} for every
  σ > 1/2**.
- It contains the conclusion directly: the Q^{1/6} principal-type copies force every zero of L(s,ν) to satisfy
  Re ρ ≤ (6+5κ)/12 = 6/7.

**(C3) Equivalence (exact Poisson).** M_κ is equivalent to the following asymptotic, which I call ASM_κ:
`Σ_h Φ(Nh/𝓗)|B_h(X)|² = FD(𝓗,X) + O(X^{2+ε})` for 𝓗 = X^{2−κ}. It must hold uniformly for the auxiliary twists f,
the exclusions b, and polarised pairs of weights. Here FD is the explicit diagonal.
- The fake diagonals cancel exactly under Σμ(f) up to a Möbius sum over the gcd variable (§2).
- Via reflection this is a **power-saving (𝓗^{−(1−κ)/(2−κ)}) asymptotic for the second moment of quadratic twists
  of the cubic theta function**. Its functional-equation group is affine (§3, (i-5)).
- So it is a reformulation, not a weakening.

**(C4) Alternative (above-diagonal, marked) hypothesis.**
- Statement: for each ν, Σ_{0<N(u)≤Z} |M_u(Z^{7/12})·Z^{−5/24}Σ_{p∼Z^{5/12}} χ_p(u)|² ≪ Z^{1+ε}. This is paper78
  Lemma `lem:marked` with one slot of length z = 5/12 at r = 7/12, outside its conditions.
- Consequence: |A₁(D)|² ≪ D^{12/7+ε} by the prime extraction (exact: θ(7/12, 5/12) = 6/7), hence θ = 6/7.
- It is true under GRH and incomparable with M_{6/7}. It is not below the diagonal, but it requires removing both
  of Lemma marked's walls. Under its proved walls the same extraction gives exactly 11/12.

**Weakest hypothesis found:** the D-averaged M_{6/7} in (C1). It cannot be reduced to statements about non-principal
rows alone (§6(d)). The goal itself, zero-freeness for all Hecke L-functions in Re s > 6/7, does **not** imply it: it
gives only S ≪ H·D^{12/7+ε}. So R1 aims at a strictly stronger statement than the goal.

## 8. Where a genuinely new idea would have to enter

Only one of the following would rescue R1:
- (1) A "column functional equation" for the Möbius coefficients against sextic characters, i.e. an R-move at level
  A. By the calculus this alone would give M_κ for all κ < 1 with an optimal large sieve.
- (2) A power-saving asymptotic for an affine-type second moment, with saving 𝓗^{−1/8} at κ = 6/7, for quadratic
  twists of the cubic theta function at unbalanced length X = 𝓗^{7/8}.
- (3) A non-dual input that sees the signs μ(v) of §2, e.g. a sieve-theoretic or spectral treatment of Q_cop that
  does not pass through positivity.
- (4) A lossless one-slot marked moment at (r, z, m) = (7/12, 5/12, 1) (§6(i), (C4)). This is above the
  canonical diagonal, so the column-term obstruction does not apply. It requires removing both marked-energy
  walls of paper78 Lemma `lem:marked`, which is the subject of the Part II routes.

None is available in either paper or in the literature surveyed (R5 notes). Numerically, Q_cop is 10⁻⁴–10⁻² of the
positive form, so the needed cancellation exists. No method in hand can prove it.
