# new-architectures: alternative probes and architectures (10 Oct 2026)

Agent: `new-architectures`. Direction: architectures that change one of the P78 scalings
(copies, slot amplification, row suppression, low side, energy wall). Status labels:
**PROVED** (a finite algebraic identity, standard), **MODEL** (exponent bookkeeping in the
P78 Part II style, the LP being exact for the stated constraints), **HEUR** (an analytic input
assumed by analogy, not proved), and **[S]** (a literature fact taken from a secondary source,
since arXiv and numdam are blocked here).

## 0. Verdict

None of the four directions gives an architecture below 0.875 with paper-type lemmas. Even
under GRH for rows, none of them goes below 13/15. Each one fails for a precise reason:

| idea | best exponent (GRH for rows, optimistic lemmas) | obstruction |
|---|---|---|
| (1) higher-rank theta duals (KP unique Whittaker) | ≥ 0.90 in the P78 shape (quartic ≥ 7/8); θ ≥ 1 in the P2 shape | A degree-r dual has reflected dual length (QP)^r/(ZP). The wall becomes (r−1)(M+ℓ) ≤ 1. |
| (2) trilinear probes | 13/15 (A·A′·B is the same as P78); 11/12 (B·B′); worse (A·B·B′) | Gauss sums are twisted-multiplicative. An extra Möbius row has device ratio ≈ 11–14 > 4. |
| (3) denser copies (c > 1/6) | 1 − 4c/5 would be great (c = 1/3 gives 11/15), but no realization exists | c = 1/N. For N < 6 the absorbed dual has degree ≥ 3 or none exists. |
| (4) other placements of the marks | 13/15 at best (quadratic m-twist ≡ P78); otherwise 8/9 or 14/15 | Placement dichotomy: the B side pays the wall, the A side pays Gram (exact cancellation). |

Four structural facts are new here. They are useful as a map for any future attempt.

1. **13/15 is pinned twice** (MODEL, `template_lp.py`).
   - The GRH floor is the common vertex of (wall ∩ case A) on the low side and of (wall ∩ rows).
   - Improving the row side alone does not move it. Row suppression ρ ≥ 2/3, or row-range growth κ_U ≤ 0.9, still gives exactly 13/15.
   - Removing case A alone also does not move it.
   - Closed forms, with copy exponent c (and a = c, Gram slope c) and wall M+ℓ ≤ W (LP-checked on a grid):
     - for W ≥ 1 (the wall ∩ case-A vertex): **θ_GRH = 1 − c(7+W)/10**;
     - for W ≤ 1 (the wall ∩ rows vertex): θ_GRH = 1 − 2c(1+W)/5.
     - At W = 1 both give 1 − 4c/5: c = 1/6 gives 13/15, c = 1/5 gives 21/25, c = 1/4 gives 4/5, c = 1/3 gives 11/15.
     - At W = 1+w the first is 13/15 − w/60, which matches 02-strengthening.
2. **Device criterion** (MODEL, `device.py`).
   - Consider any new structure ("device") of log-length λ that multiplies signal/low by Z^{g_L λ} and rows/signal by Z^{g_R λ}, and touches neither the energy wall nor case A.
   - It lowers the GRH floor **iff g_R/g_L < 4**. At ratio 3 the floor is 41/48; at ratio 2 it is 5/6.
   - The P78 marks have ratio 2 (g_L = 1/6, g_R = 5/6 − 1/2 = 1/3). Their gain is capped only because they consume the wall.
   - A second theta row has ratio (4/3 − β)/(β − 5/6), which is 11 at β = 7/8.
3. **Absorption classification** (PROVED for the identities; MODEL/HEUR for the conductors). See §2.
   - P2's analogue for order-N twists with a degree-r dual gives θ = 1/2 + (1 − 1/N)(1 − 1/r).
   - Only (N, r) = (6, 2) is realizable, giving 11/12.
4. **Amplifier placement dichotomy** (MODEL). See §5. The energy wall QP ≤ Z does not depend on where the marks are placed.

## 1. The template (`template_lp.py`, validated)

Exponents are in units of log Z. The geometry is X = Z^{l_x}, Y = Z^{l_y}, P = Z^ℓ and Q = XY, with M = l_x + l_y and b = l_y − l_x. The template has five parameters:
- c, the copy exponent (the principal frequency has copies ua^N with c = 1/N);
- a, the slot amplification;
- ρ, the relative slot suppression;
- κ_U, the row-range growth;
- r, the degree of the reflected dual.

**Low side:**
σ ≥ (1−c) + c·l_x − aℓ − d + ½·max(0, c·b, 2b − l_y + d) + ½·max(0, caseA, (r−1)(M′+ℓ′) − 1, (M′+4ℓ′−2)/3), for d ∈ {0, ℓ}.

**High side under GRH:**
σ ≥ (1 − l_y)/2 − ρℓ + (1−c)(1 − l_x + κ_U ℓ).

The signal is X^{1/2}·F^{c}·Z^{β−1}, with F = ZP/X the frequency range. So the amplification is a = c·κ_U: the marks add copies because they extend F.

The base parameters reproduce `r4/model.py` with `grh=True`: 13/15 at (l_x, l_y, ℓ) = (2/5, 2/5, 1/5). The case W(0.1) gives 0.865 = 13/15 − 1/600. The sensitivity table is in `scenarios_output.txt`.

Two further checks:
- Removing the reflected wall and using the coefficient-blind BGL energy, ‖B‖²/Q ≪ 1 + ZP/Q + (ZP)^{2/3}Q^{−1/3}, gives θ = 1 for every c. Theta reflection is therefore indispensable on the low side.
- With paper counts, case A does not bind: the r4 model gives 0.8749568 without case A and 0.8749570 with it.

## 2. Idea 1: higher-rank theta functions

**μ-absorption is a statement about Gauss sums, not about the group** (PROVED, classical).

Setup: let π be a primary prime and χ = χ_π a residue character of order N ≥ 3, with γ = q^{−1/2}g.

The identity: from J(χ, χ^{−2}) = g(χ)g(χ^{−2})/g(χ̄) and the Stickelberger/Iwasawa normalization J ≡ −1 mod (1−ζ_N)², the quotient μ(π)√q/J(χ,χ^{−2}) is a Hecke character value. Hence

  **μ(n)·γ(χ̄_n) = Hecke(n)·γ(χ_n)·γ(χ_n^{−2})**  (the cross terms match for squarefree n).

The Poisson dual of an order-N twist family therefore always absorbs μ into a product of two Gauss sums. The cases are:
- **N = 6.** Davenport–Hasse at p = 2 gives γ(χ)γ(χ⁴) ∝ γ(χ²)γ(χ³). Since γ(χ³) is a ray factor, the result is a single cubic Gauss sum (P78/P2).
- **N = 4.** Here γ(χ²) is a ray factor, giving a single quartic Gauss sum of χ. This is the refuted quartic transplant.
- **N = 3.** The result is γ(χ)². Inverting the P2 identity, it is ∝ γ₆(φ) with φ² = χ̄, a single sextic Gauss sum.
- **N = 10 ↔ 5, 14 ↔ 7, …** These come from Davenport–Hasse for 2m with m odd (`04-gram-split/LIT_ideas234.md` §5).

**Consequence: GL_r^{(n)} thetas offer no new identities, only new realizations of known duals.**
- **Order-3 Gauss sums** are realized on GL₂^{(3)} (degree 2), and possibly on GL₃^{(3)} first-row coefficients (degree 3).
- **Order-4 Gauss sums** are realized on GL₃^{(4)} (degree 3), or as GL₂^{(4)} Kubota–Eisenstein coefficients.
  - For the Eisenstein realization, the twist is the shift ν = h³ with conductor N(h)³. This was verified in `03-six-sevenths/quartic-transplant`.
- **Order-5 Gauss sums** are realized on GL₄^{(5)} or GL₅^{(5)}.
- **Order-6 Gauss sums** are realized on GL₅^{(6)} or GL₆^{(6)}, or as GL₂^{(6)} Eisenstein coefficients.
  - For the Eisenstein realization, the twist is a shift, γ(φ)φ(h)² = γ(φ, ν) with ν ≡ h^{−2} ≡ h⁴ mod sixth powers. The conductor is N(h)⁴ (HEUR, by the quartic analogue).
- [S] Kazhdan–Patterson: the theta on GL_r^{(n)} has a unique Whittaker model iff n ∈ {r, r+1}.
  - That its first-row coefficients at squarefree indices are order-n Gauss sums is recalled, not verified. The conclusions below do not depend on it.

**Exponents** (MODEL, optimistic: a diagonal large sieve for every reflected dual, no case A, side-condition margin 0). The high side is identical to P78 with c = 1/N. Only the low-side wall changes, because a row of length ZP twisted with modulus QP has reflected dual length (QP)^r/(ZP), and that dual length must be at most Q.

| family → dual | c | r | θ_GRH (P78 shape) | θ (P2 shape) = 1/2 + (1−1/N)(1−1/r) |
|---|---|---|---|---|
| sextic → GL₂^{(3)} theta (paper) | 1/6 | 2 | 13/15 | 11/12 |
| sextic → GL₃^{(3)} first row | 1/6 | 3 | 11/12 | 19/18 |
| quartic → Kubota N(h)³ | 1/4 | 3 | 7/8 (true value 1: the quartic dual twist has no good sieve) | 1 |
| cubic → sextic Kubota N(h)⁴ | 1/3 | 4 | 9/10 | 1 |
| cubic → GL₅^{(6)} / GL₆^{(6)} | 1/3 | 5 / 6 | 17/18 / 0.976 | 31/30 / 19/18 |
| order 10 → GL₄^{(5)} / GL₅^{(5)} | 1/10 | 4 / 5 | 1.04 / 1.10 | 47/40 / 61/50 |
| *hypothetical* degree-2 dual, cubic family | 1/3 | 2 | 11/15 | 5/6 |
| *hypothetical* GL₂ quintic theta, order 10 | 1/10 | 2 | 23/25 | 19/20 |

**Obstruction.**
- The copy gain of a lower-order family is always paid for by a dual of degree ≥ 3, or ≥ 4 for the cubic family. In the P2 shape the extraction then needs H < D^{…} while the transfer needs H ≥ D^{2−2/r}. These are exactly complementary for (4,3) and (3,4).
- Order 10 is worse on both counts: fewer copies (c = 1/10) and degree 4–5.

## 3. Idea 2: trilinear probes

- **A·A′·B (two Gauss polynomials, one theta row)** (PROVED reduction).
  - For coprime s, s′: g(χ_s,m)g(χ_{s′},m) = conj(χ_s(s′)χ_{s′}(s))·g(χ_{ss′},m).
  - So A·A′ is a single Gauss-sum polynomial of length YY′ with convolution coefficients: it is P78 with Y ↦ YY′. The exponents are identical.
  - The numerator becomes L(w,χ_u)², which is only worse off GRH.
- **B·B̄′ (two theta rows, no A)** (MODEL, analytic).
  - The signal is Q^{1/2−c}(ZZ′)^{β−1+c}, the low side is Q^{1/2}, and the rows are (ZZ′/Q)^{1/2}.
  - With t = M/(1+z′) ≤ 1/2 we get β > max(1−c+ct, 3/2−c−(1−c)t), so the best is 1 − c/2 = 11/12 at t = 1/2.
  - Marks only worsen it.
- **A·B·B′.**
  - The extra Möbius row of length Z′ multiplies signal/low by Z′^{β−5/6} and rows/signal by Z′^{4/3−β}.
  - That is a device ratio of 11 at β = 7/8 and 14 at β = 13/15. The criterion needs < 4.
  - It also needs its own energy condition Q ≤ Z′.
- **Other Cauchy–Schwarz groupings.**
  - Cauchy–Schwarz in A or in s gives the same X^{1/2}: it is coefficient-blind on one side, or it is the operator norm of the same Gram matrix.
  - Moving part of the s-length to the theta side of Cauchy–Schwarz trades Gram (gain y₁/6) for the wall (loss ≥ y₁/6, and the s-part is sextic, so BGL). The net is ≤ 0.
- **Where the losses are.** In the random model |I| ≈ Y^{−1/2}. Cauchy–Schwarz loses Q^{1/2} and the absolute row sum loses F^{1/2}. At the GRH vertex, low ≈ rows ≈ signal, so a gain must come from the high side and the low side together.

## 4. Idea 3: denser copies

- **Copies are tied to the twist order.** Copies of the principal frequency are exactly the u with χ_n(u) = 1 for all n, so c = 1/N with N the twist order. By §2, N < 6 forces a dual of degree ≥ 3 (N = 3, 4) or gives no absorption at all (N = 2, where μγ₂ = μ·ray).
- **Subfamilies do not help.**
  - u ∈ squares makes the twists cubic; u ∈ cubes makes them quadratic. Both inherit those obstructions.
  - The Poisson dual over u = b² is μ·conj(γ₂), which has no structure. Over u = b³ it is μ·ray, the move calculus again.
- **Larger fields do not help.** Over K′ ⊇ K with rows u ∈ O_K, the identity (u/𝔫)_{K′} = (u/N𝔫)_K keeps the ratio of copies to family size at H^{1/6}.
- **Hybrids lose.** A sextic theta row with a cubic Gauss side creates "cube rows" u = b³ whose numerator has a pole and so gets no Y-suppression. That costs F^{1/3} unsuppressed rows.
- **The value of c, had it been realizable.** With a degree-2 dual and a quadratic dual twist:
  - θ_GRH = 1 − 4c/5, which is 0.84 for c = 1/5 and 4/5 for c = 1/4;
  - in the P2 shape, 1/2 + (1−c)/2, which is 7/8 for c = 1/4 and 5/6 for c = 1/3.
- **Why it is not realizable.** It would need a degree-2 automorphic object whose squarefree coefficients are Gauss sums of order 4 or 6. Kazhdan–Patterson non-uniqueness on GL₂^{(n)} for n ≥ 4, and the Bröker–Hoffstein data (τ₆(p) is not a Gauss sum), rule that out.

## 5. Idea 4: other placements of the marks (amplifier Σ_p w_p ψ_p)

Rows can only see an amplifier ψ_p(u) if p enters the Poisson modulus in m. That forces F ↦ FP, i.e. κ_U = 1, and hence a = c. The κ-marks invariance in ../../NOTES.md is the special case of rescaling Z.

| placement | effect (MODEL) | θ_GRH |
|---|---|---|
| A-marks 1_{p\|A}, Z→Zp (P78) | a = 1/6, ρ = 1/2, wall Q²P/Z ≤ Q | 13/15 |
| restriction p\|m | χ̄_s(p)χ_A(p) sends the s- and A-characters to different frequencies (u/p vs up). The target row 1/L(x,η) occurs only at H = p⁵a⁶, where the numerator is L(w, χ_p⁴·), with no pole. **The signal loses its w = 1 residue (≈ Y^{−1/2}).** | worse |
| p\|s, Y→Yp | Identical to the unmarked probe at (X, YP, Z): a reparametrization of l_y | 13/15 (no new freedom) |
| p³\|A (cube part) | Rough: the local factor is p^{1−3x−3z} = p^{1/2−3β} at the signal point, with no μ(p) coherence (the case-A mechanism) | no amplification |
| twist m by χ_p(m), weights conj γ₆(p), grouped with B | Same amplifier, range and suppression as A-marks. But after reflecting B the columns are pk with a **sextic** p-part, so BGL's (QN)^{2/3} gives the wall Q^{3/2}P ≤ Z. | 8/9 (ℓ* = 0) |
| twist m by (m/p)₂, weights conj γ₃(p), grouped with B | The reflected columns pk are purely quadratic, so the wall is exactly QP ≤ Z, with no cube-part marks (no case A). Square rows u = b² are coherent but thin (F^{1/2} ≤ FP^{−1/2}). | 13/15; 0.8749568 with paper counts (≡ P78) |
| m-twist grouped with A | A·M is a Gauss polynomial of modulus sp on fixed Q, i.e. (X, Y) → (X/P, YP). The Gram excess (YP²/X)^{1/6} costs exactly P^{1/6}. If YP² ≤ X instead, the optimum is 8/9 for every ℓ. (`aside_marks.py`) | 8/9 |
| m-twist, sup-norm on B (Z ≥ Q²) | The wall is P-free, but M ≤ 1/2 | 14/15 |

- **Dichotomy.** Put the amplifier on the theta side of Cauchy–Schwarz and the reflected dual length becomes Q²P/Z, giving the wall QP ≤ Z wherever the mark sits. Put it on the Gauss side and the Gram shadow cancels the copy gain exactly. **So the energy wall cannot be removed by moving the marks.**
- **On the earlier χ_p(m) check.** Without the conj γ₆(p) weights the signal is incoherent; that is the earlier finding. With the weights the amplification P^{1/6} comes back, but the energy is worse (sextic). The quadratic version is exponent-equivalent to P78.
  - It is structurally cleaner: no cube-part marks, and quadratic slot polynomials. That might help lemma-level constants, but per r4, c_P = 1/3 gives only 0.87471.

## 6. What a genuinely new architecture would have to supply

Each of the following is necessary at the exponent level. Under GRH for rows each would beat 13/15; with paper counts each must be paired with the row machinery.
- **(a)** A degree-2 μ-absorbing dual with a quadratic dual twist for a twist order N < 6. This is impossible by §2 and §4, as far as Kazhdan–Patterson/Patterson are known.
- **(b)** A wall-free and case-A-free device with g_R/g_L < 4. Every candidate examined has either ratio ≥ 11, or ratio 2 but consumes the wall.
- **(c)** Energy beyond the wall, W(w): the known lemma-level lever, outside this direction.

## 7. Files (this directory)

| file | content |
|---|---|
| `template_lp.py` | parametrized GRH template LP (c, a, ρ, κ_U, Gram, degree r, case A, wall); validated against r4 |
| `scenarios.py`, `scenarios_output.txt` | every number quoted above (`python3 scenarios.py`) |
| `device.py` | generic wall-free device; threshold g_R/g_L = 4 |
| `aside_marks.py` | amplifier grouped with the Gauss polynomial (Gram cancellation) |
