# Pushing OpenAI's quasi-Riemann bound below 7/8, and where the method stops

*Analysis of `openai/math` family 003: "The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane Re s > 7/8" (30 Sep 2026) and the companion "alternate 11/12 proof" (5 Oct 2026). Written 8 Oct 2026.*

## Summary

1. **A strictly better exponent from the paper's own lemmas.** Keep every lemma of Part II and change only the fixed geometry (`old-eq:5.1`) to

   > b = 1/8, ℓ = 1/6 + 1/10⁴, M + ℓ = 1, so l_x = 21247/60000, l_y = 28747/60000, h = 16253/20000.

   The argument then closes at

   > **θ = 7/8 − 1/40000 = 0.874975.**

   Every geometry-dependent inequality is re-checked in exact rational arithmetic in `certificate.py` (24/24 pass). That includes a new Bernstein-polynomial certificate for the nonlinear endpoint inequality, which replaces the paper's Lemma "Compensated endpoint certificate". The checker is discriminating:
   - It re-certifies the paper's published margin 49/440640 at the paper's own geometry.
   - It fails when asked for a 3·10⁻⁴ margin there.
   - It fails for ℓ = 1/6 + 1/5000.

2. **That gain is essentially all the paper's lemmas allow.** An exact LP-duality certificate (`barriers.py`, part A) shows that **no choice of geometry** gives θ < 33216105/37963129 ≈ **0.874957** from these lemmas. The optimum of this lemma set therefore lies in [0.874957, 0.874975]. 7/8 is the end of the line for this exact proof, up to about 4·10⁻⁵.

3. **Even with GRH for the auxiliary rows the method stops at 13/15.** Suppose every auxiliary "row" L-function satisfied GRH, so the zero-detector and row counting become perfect. Keep the low-side mean-value lemmas as they are. Then the architecture still cannot go below **θ = 13/15 ≈ 0.8667** (`barriers.py`, part B, exact multipliers 4/5 and 1/5).

4. **θ = 1/2 is out of reach for this method.** Getting below 13/15 needs a longer prime "amplifier" (ℓ). That forces the completed cubic-theta energy into the range where every available mean-value tool loses. The reflection has a Q²/Z dual-length term. The sextic large sieve has a (QN)^{2/3} term, which is open: removing it is the cubic large sieve conjecture. Sixth-power rows contribute a genuine main term. Approaching 1/2 needs an amplifier of unbounded length together with power-saving asymptotics for these energies. A consistency argument (§6.3) shows that is GRH-strength input. So the method cannot reach 1/2 without, in effect, assuming something as strong as what it is trying to prove.

## 1. How the 7/8 proof works

Work over F = ℚ(√−3) with sextic residue symbols χ_n(u) = (u/n)₆. The key identity is that a sextic Gauss sum times μ(n) becomes a *cubic* Gauss sum, i.e. a Fourier coefficient of Kubota–Patterson's cubic theta function.

- **The probe I_η(X, Y, Z)** is one finite character sum with two exact evaluations.
  - *Low (reflected) side:* cubic-theta reflection converts the sextic twists into quadratic ones. Heath-Brown/Goldmakher–Louvel's quadratic large sieve and a planar additive large sieve then give |I| ≲ X^{1/2}·P_a^{1/12}, with P_a = Y/X.
  - *High (Poisson) side:* Poisson summation in the row variable writes I as a sum over sixth-power-free rows u of ζ(6z)·L(w, χ_u)/L(s, η χ̄_u) times a holomorphic Euler correction.
- **The principal row u = 1** contains 1/L(s, η). Its residues at w = 1 and z = 1/6 give the "signal" Z^{C(s)}·∫H/L, with C(s) = s − 5/6 + l_x/3 + ℓ/6.
- **Non-principal rows** are moved to the line Re s = a + O(e), where a ≤ β* is the real part of the rightmost zero in a buffered box. A zero detector turns each row with a zero into a large inverse Dirichlet polynomial and a large plain one. Two moment estimates then count such rows: a marked inverse second moment, and a fourth moment with prime factors, together with sixth-power amplification.
- **Prime compensation (Part II).** The probe is modified at K disjoint prime "slots" of total log-length ℓ. This amplifies the principal row by Z^{ℓ/6} (from Σ_p p^{−5/6}). For non-principal rows the slots become twisted prime sums Σ_p χ̄_p(u)p^{−1/2}. These either cancel, or are large and earn a counting gain.
- **Continuation criterion.** If low exponent ≤ C(σ_0) and all non-principal rows are below Z^{C(β*)−σ}, then β* > σ_0 is contradicted. This happens uniformly over the whole Hecke family, which is why the family-wide supremum β* is used.

Part I is the same scheme with no slots and X = Y = Z^{1/2}, giving 11/12. The separate 11/12 paper instead extracts A₁(D) from a mean square over u ≤ D^{1+ϑ}, using the rows u = p⁶ (§6.4).

## 2. The exponent bookkeeping (transcribed from Part II)

The geometry has X = Z^{l_x}, Y = Z^{l_y}, b = l_y − l_x, M = l_x + l_y, slot length ℓ, and row range h = 1 − l_x + ℓ. With δ = 2a − 1 and prime amplitude x = q/δ ∈ [0, 1/2]:

| Piece | Exponent requirement (relative to C(σ_0)) | Source in paper |
|---|---|---|
| low side | σ_0 ≥ 5/6 + l_x/6 − ℓ/6 + ½·max(0, b/6, 2b − l_y) + ½·e | Prop. `probe-low`, `probe-gram` |
| completed-row energy | e = max(0, (1+3ℓ−2M)/4, M+ℓ−1). This is the paper's own bound, the true maximum of eq. 4.6 over block data. An earlier draft used the weaker (1+2ℓ−2M)/3, which assumes no powerful row part; LP lower bounds computed with it remain valid, since the true constraint is stronger. | Lemma `reflected-energy`, Lemma `probe-row-norm` |
| high side, d = h | (1+δ)(1−l_y)/2 − σ_0 − ℓ/2 + xδℓ + h(R_* + δ/2 − 1/6) < 0 | eq. `common-high-exponent` |
| row count | R_* = 1 − δ + (α−δ)δP_x/(2J), J = (α−δ)D_x + δP_x, D_x = 3 − 17x/9, P_x = (2 − 8x/9)(1−x), α = 5/6 | Prop. `detector-counts` |
| floor bin, intermediate rows, small rows, supply ℓ/h > 7/37 | linear checks | §§ "seven-eighths bound", "frequency ranges" |

At the paper's geometry this reproduces the published numbers exactly:

- the low side gives exactly 7/8;
- the high side's worst point is δ ≈ 79/204, x = 1/2;
- the endpoint margin is ≥ 49/440640.

Two structural facts follow:

- **The signal depends only on b.** Under M + ℓ = 1, C(s) = s − 2/3 − b/6, and the low side is **σ_low = 11/12 − ℓ/4, independent of b**. Raising ℓ lowers the target, but it lengthens the row range h = (1 + 3ℓ + b)/2, and the high side pays about 1.08 per unit of ℓ.
- **The binding constraints** are low = high at the interior point δ ≈ 0.39, x = 1/2, together with M + ℓ = 1, the energy wall in §6.2.

## 3. The improvement θ = 7/8 − 1/40000

Changes to Part II, everything else verbatim:

- **Geometry.** Replace `old-eq:5.1` by b = 1/8, ℓ = 5003/30000, h = 16253/20000, l_x = 21247/60000, l_y = 28747/60000, M = 24997/30000 (still M + ℓ = 1).
  - The Mellin exponent is unchanged: C_II(s) = s − 11/16.
  - The low bound becomes Z^{3/16 − 1/40000} = Z^{C(σ_0)}.
- **Target and bootstrap.** Set σ_0 = 34999/40000 and Δ = β* − σ_0. Part I still gives β* ≤ 11/12, so δ ≤ κ ≤ 5/6, the range certified.
- **Plain fourth moment.** Use κ = max(3/4, 2β* − 1) in Lemma `plain`; that lemma requires κ ∈ [3/4, 1]. The capacity loss is then ≤ max(0, β* − 7/8)/4 ≤ Δ/4, absorbed exactly as in the paper because 1 − h/4 > 0.
- **Euler region two and the outer-row lemma.** Extend region two from Re s ≥ 7/8 down to Re s ≥ σ_0. The four principal local error exponents stay ≤ −σ_0, and the good-prime defects stay below −1. The checks are in the script. Lemma `outer-row-tails` uses D₁(3/8) only through β* + 1/2 > 1, which still holds.
- **New endpoint certificate.** −E(h) ≥ 9·10⁻⁵ uniformly on δ ∈ [0, 5/6], x ∈ [0, 1/2], certified exactly by Bernstein expansion (35 boxes).
- **Other margins:**
  - floor bin: −0.0057;
  - intermediate rows: −0.0033;
  - small rows: −0.086;
  - supply: ℓ/h − 7/37 = 0.016;
  - room above h: 5ℓ − h = 0.021;
  - principal margins: l_y/20 and h/600.

As a corollary, every consequence in the paper improves by the same amount, for example π(x; q, a) − Li(x)/φ(q) ≪ x^{0.874975} log x.

**Caveat.** This is a modification at the level of the paper's own exponent bookkeeping. It inherits the correctness of the paper's (unrefereed) lemmas, which are stated for general geometry. It is not Lean-checked: 19 files of the 487k-line formalization hard-code 17/48, 23/48 and 13/16.

## 4. Why the paper's lemmas stop at about 0.87496

`barriers.py` (A) takes a nonnegative rational combination of:

- the low-side constraint (weight 0.812);
- its Gram/energy epigraphs (0.406, 0.344, 0.0625);
- the high-side constraint at δ = 0.387 and 0.39 (weights 0.089, 0.099).

The combination eliminates every parameter and leaves σ_0 ≥ 33216105/37963129. Since it uses only some of the high-side constraints, it is a valid lower bound for every geometry. So there is no hidden better parameter choice. The authors' 7/8 sits within 4.3·10⁻⁵ of the true optimum of their own machinery.

## 5. Which ingredient would have to improve

The best σ_0 after re-optimizing the geometry in each case (`sensitivity.py`):

| Change to one ingredient | best θ |
|---|---|
| none (paper) | 0.87496 |
| plain-moment capacity 2/9 → 1/3 (Lemma `plain`) | 0.87471 |
| amplification slope α: 5/6 → 2/3 (Lemma `inverse-amplification`) | 0.87407 |
| perfect row counting, R = 1 − δ (density hypothesis for every row family) | 0.86979 |
| Gram bound without the P_a^{1/6} sixth-power term | 0.86398 |
| GRH for all rows, detector floor → 1/2 | **13/15 = 0.86667** (exact, §6.1) |

So improving the zero-detector and counting side alone is worth at most about 0.008. The real bottleneck is the low side: the completed-theta energy and the Gram term.

## 6. Barriers to θ = 1/2

### 6.1 Even GRH for all rows only reaches 13/15

With GRH for every row, rows contribute trivially (count U, |L(1/2)| ≈ U^{o(1)}). The high side then reads σ_0 > (1 − l_y)/2 − ℓ/2 + 5h/6.

Combine it with the low side using the multipliers 4/5 (low) and 1/5 (high), plus the Gram and energy epigraphs (1/5, 1/5, 1/3, 1/15). Everything cancels except **σ_0 ≥ 13/15**. Equality holds at l_x = l_y = 2/5, ℓ = 1/5, which is the vertex where the energy wall M + ℓ = 1 binds. The problem is convex, so this is the global optimum.

### 6.2 The energy wall (M + ℓ ≤ 1)

The marked completed-row energy is "diagonal" (≈ Q) only when the row length plus the active mark length is at most the completed length plus the mark length. Past that point:

- **Reflection** adds a dual-length term, a loss of M + ℓ − 1.
- **The sextic large sieve** (Blomer–Goldmakher–Louvel) has a (QN)^{2/3} term.
- **Sixth-power rows** m = k⁶j add N·Q^{1/6} no matter which large sieve is used.

Removing (QN)^{2/3} is the open cubic/sextic large sieve conjecture. Even that would not remove the sixth-power-row term.

### 6.3 Why no mean-value input can carry it to 1/2

Below 13/15 the bookkeeping demands ℓ → ∞, i.e. an amplifier of unbounded length, while still keeping the low bound under the signal.

But the probe *is* "principal row + other rows". Take a target η with zeros approaching the family supremum β* ≥ 1/2; such zeros exist unconditionally, since L-functions have zeros on the critical line. For that target the principal row is Ω(Z^{C(β*)−ε}), and the high-side machinery keeps the other rows smaller. Any valid low bound is therefore at least the principal size. That means the energies have main terms governed by the same zeros. Pushing toward 1/2 requires power-saving asymptotics for amplified (mollifier-length) energies far past the diagonal. That is the familiar mollifier-length barrier, and it is equivalent in strength to the conclusion. With the exponent model, the "optimal large sieve" idealization produces nonsense bounds (σ_0 < 0), which confirms that those main terms are real.

The architecture's other floor comes from its sextic shape: the signal carries Z^{σ − 5/6} from the ζ(6z) pole, against a Cauchy–Schwarz low bound of X^{1/2}. Without prime amplification, θ ≥ 5/6 + l_x/6, so the 5/6 floor is structural to the sextic probe.

### 6.4 The 11/12 companion has its own wall

There θ = (6 + 5κ)/12, where the mean square Σ_{N(u)≤D^κ} |A_u(D)|² ≪ D^{1+κ+ε} is known only for κ > 1, beyond the large-sieve diagonal. Getting 7/8 that way needs κ = 9/10, and getting 1/2 needs κ → 0: square-root cancellation over families smaller than the sum length, again GRH-strength.

## 7. Bottom line

- A strictly smaller exponent than 7/8 follows from OpenAI's own lemmas: θ = 7/8 − 1/40000, exactly certified here, and that is essentially the best those lemmas can do.
- Any substantial improvement must come from new mean-value input on the low side, not from the zero-detector or counting side.
- Even perfect row control caps the architecture at 13/15.
- The method cannot reach θ = 1/2 without GRH-strength mean-value estimates.
- Your friend's unpublished lower θ can't be checked from here. If it is substantially below about 0.8749, it must change one of the ingredients in §5, most plausibly the Gram/energy side.

## Files

| File | What it does |
|---|---|
| `certificate.py` | Exact (stdlib-only) check of all Part II inequalities at θ = 7/8 − 1/40000, including the Bernstein endpoint certificate. |
| `barriers.py` | Exact LP-dual certificates: θ ≥ 0.8749570 for the paper's lemmas, θ ≥ 13/15 under GRH for rows. |
| `sensitivity.py` | Re-optimizes the geometry under each modified ingredient (numpy/scipy). |
