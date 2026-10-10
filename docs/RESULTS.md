# Results: proved, certified, conditional, partial, retracted

Every result from the project, sorted by strength of evidence. θ always means a zero-free half-plane Re s > θ.

The L-function family is OpenAI's:
- every finite-order Hecke L-function over Q(√−3) = Q(ζ₃);
- hence, by base change L_K(s, χ∘N) = L(s,χ)L(s,χχ₋₃), every Dirichlet L-function and ζ.

Notation follows OpenAI's 7/8 paper (P78) and its 11/12 companion (P2).

## A. Machine-checked in Lean (standard axioms only)

| # | Statement | Evidence | Where |
|---|---|---|---|
| A1 | **β\* ≤ 10499/12000 = 0.8749166…** ζ, all Dirichlet and all Q(ζ₃) Hecke L-functions are nonzero on Re s > 10499/12000. This is 4.04·10⁻⁵ below the public cubic-root endpoint. | Isolated rebuild; Lake build; Comparator ×3 with negative control; con-ron + NanoDa + Lean three-kernel PASS (author-side, unsandboxed); submitted as [tracker PR #6](https://github.com/GTimG/QuasiRiemannTracker/pull/6) | [08](../08-weighted-numerator-10499-12000) |
| A2 | **ζ(s) ≠ 0 for Re s > 23/24.** A standalone, tree-shaken proof of 292 modules and 121,523 code lines, with a paper-numbered spine and a matching 46-page paper. | Comparator; import audit | [05](../05-lean-23-24) |
| A3 | **L(s,χ) ≠ 0 for Re s > 23/24** for every Dirichlet χ (ζ as a corollary), on the pruned openai/math fork | Comparator | [05](../05-lean-23-24), [akashlevy/math#1](https://github.com/akashlevy/math/pull/1) |
| A4 | OpenAI's 7/8 theorems (ζ, Dirichlet, Hecke) still verify after in-place pruning of `DirichletL` (−25.7% code lines) | Comparator, all 5 QRH challenges | [05](../05-lean-23-24) |
| A5 | The energy induction of P78 extends to κ ≥ 13/18. This is **not** a zero-free theorem. | `#print axioms`; standard three | [06](../06-kappa-bootstrap) |
| A6 | `WeightedQRH.Parameters.physical_low_barrier`: the 13/15 low-majorant identity, outside the proof path | Lean | [08](../08-weighted-numerator-10499-12000) |
| A7 | Formal mean square for every σ > 1/20 (`outside_mean_square_of_canonical_density_at_radius`). This is the prerequisite of the 15/16 lead. | Lean, already in openai/math | [05](../05-lean-23-24) |

## B. Paper-level proofs with exact rational certificates (no Lean)

These reuse P78's analytic lemmas at new parameters. Those lemmas are unrefereed and not Lean-checked at the new geometry.

| # | Statement | Certificate | Where |
|---|---|---|---|
| B1 | **θ = 34999/40000 = 7/8 − 1/40000.** Geometry b = 1/8, ℓ = 1/6 + 10⁻⁴. Found independently of and on the same day as Liu. | `certificate.py`: 24/24 checks, Bernstein endpoint certificate | [01](../01-seven-eighths-lemma-set) |
| B2 | κ-bootstrap: 17499141/20000000, then **43747851/50000000**. Not novel: stage 2 equals Guo's published rational, and the mechanism is publicly known. | `bootstrap_certificate.py` | [06](../06-kappa-bootstrap) |

## C. Conditional results (hypothesis ⇒ conclusion, proved or exactly certified)

| # | Hypothesis ⇒ conclusion | Where |
|---|---|---|
| C1 | **Theorem C1.** D-averaged M_{6/7}, i.e. ∫Σ_{N(u)≤D^{6/7}}\|A_u\|² ≪ D^{13/7+ε}, ⇒ θ = 6/7. Full proof. The hypothesis is strictly stronger than the conclusion. | [03](../03-six-sevenths) |
| C2 | **Theorem C.** M_κ ⇒ β\* ≤ (6+5κ)/12: 6/7 at κ = 6/7, 3/4 at κ = 3/5, 1/2 as κ → 0. Uses Lemma B. | [03](../03-six-sevenths) R2 |
| C3 | BMS, a weaker averaged mean-square statement, ⇒ 6/7 | [03](../03-six-sevenths) R2 |
| C4 | Lossless one-slot marked moment at (r, z, m) = (7/12, 5/12, 1) ⇒ 6/7 (true under GRH) | [03](../03-six-sevenths) R1 |
| C5 | W(1/5) ⇒ **108/125** (2181 boxes); W(4/7) ⇒ 6/7; low side ≥ 13/15 − w/60 for every w | [02](../02-strengthening), [03](../03-six-sevenths) |
| C6 | Refined count + **W(3/10)** ⇒ 6/7 (17321 boxes, robust). Refined count + W(2/7) ⇒ 6/7, which passes only at ZETA = 10⁻⁶ and is fragile. | [03](../03-six-sevenths) R4 |
| C7 | GRH for rows + Gram bound without P_a^{1/6} ⇒ 6/7. GRH rows + W(1/7) + refined count ⇒ 6/7. GRH rows + wall w ⇒ 13/15 − w/60. | [03](../03-six-sevenths) |
| C8 | Vertex B: Gram fix + raw inverse moment for H ≥ D^{8/9} + floor 1/1000 ⇒ 6/7 + 1/1000. Vertex C: W(1/7) + refined count + κ_A = 2/3 ⇒ 6/7 + 1/1000. | [03](../03-six-sevenths) R4 |
| C9 | **Corollary D** (clean Gram): ‖A − A^sig‖² ≪ (Q/Y)Z^ε. Proved modulo three unwritten lemmas of known type (Lemma C at all cusp pairs, the isotypic transfer, finite order). On its own it gains nothing: 0.874957. | [04](../04-gram-split) |
| C10 | (V′) with saving s₀ ⇒ θ(s₀). Slotted: 0.8733 at s₀ = 0.002 … 0.8532 at s₀ = 0.04. Random rate ⇒ 0.851918. LoA_π(θ_π < M/(2M−1)) ⇒ (V′). | [04](../04-gram-split) |
| C11 | Single-ingredient payoffs (float model): Gram exponent 0 ⇒ 0.86398; perfect counts ⇒ 0.86979; c_P = 1/3 ⇒ 0.87471; α = 2/3 ⇒ 0.87407 | [01](../01-seven-eighths-lemma-set) |

## D. Partial results and new lemmas

**Lemmas about the 7/8 and 11/12 machinery**
- **Refined powerful-row count** (new, elementary). The case-A loss improves from (1+3ℓ′−2M′)/4 to (1+2ℓ′−2M′)/3. [03](../03-six-sevenths) R4.
- **Detector floor** 51/100 → 1/2 + η is legitimate. [03](../03-six-sevenths)
- **Lemma B** (bootstrapped prime extraction): |A₁ − A_{p⁶}| ≪ (D/Np)^{β+v} under H(β). [03](../03-six-sevenths) R2.
- **`lem:plain` holds for κ ∈ [13/18, 1]**; the low side at general ℓ ≤ 1/5 gives θ_low = 11/12 − ℓ/4; Euler regions extend down to 437/500. [06](../06-kappa-bootstrap)
- **Weighted original-numerator fourth moment** with selected primes, amplitude-class Hölder, and the quotient-free coefficient expansion. [08](../08-weighted-numerator-10499-12000)
- **Signed dual of P2's mean square:** exact common-divisor cancellation, the fake diagonal removed, and a large-j tail ≪ D³/(B²J^{17/8−δ}). [07](../07-method-barriers)

**Lemmas about the reflected theta structure** ([04](../04-gram-split))
- **Lemma Z:** ζ_m = κ₀ψ₁(m)·conj(γ₂(m))⁴ on squarefree m, and its generalization: Möbius absorption iff j = 1, ζ-free iff j ∈ {2, 5}.
- Lemma C0, Lemma R (the Dunn–Radziwiłł transform equals the paper's reflection), Lemma P (parity).
- Zero-slot count R₀(δ) = (15−16δ)/(15−6δ).
- Lemma D1, G3/G4 closure, Lemma H / Theorem 6.1 (probe rigidity).

**Counting and exact identities** ([03](../03-six-sevenths))
- Density exponents for distinct sextic rows: A(6/7) = 7/3, A(7/8) = 16/7, A(11/12) = 2.
- A_{u₀} = A_{u₀p⁶} − ψ(p)A_{u₀p⁶}(D/Np); the fake-diagonal Möbius identity; μ(n)χ_n(u) = G(n)⁻¹ᾱ(n)γ₂(n)γ₁(n; u⁵).

## E. Numerical evidence (does not establish provability)

- **M_κ:** per-row energy / diagonal ≈ 0.91–0.98 at κ = 6/7 (D = 3000, 10⁴). [02](../02-strengthening)
- **W(w):** case-B energies stay diagonal 8–128× past the wall. [02](../02-strengthening)
- **(V′_LS)** sits at the random rate up to Q = 9.6·10⁴, and saturation is excluded by a factor of 30–50. [04](../04-gram-split)/vls/E
- **Fake-diagonal cancellation exists in the signs:** Q_cop ≈ 10⁻⁴–10⁻² of Q₁. [03](../03-six-sevenths) R1
- **Quartic:** the dual shift is h³ (spread 10⁻¹³, against 2.00 for h); no polar terms up to norm 2·10⁶. [03](../03-six-sevenths)/quartic-transplant

## F. Retracted claims and errata

Do not cite these. They are kept for the record.

| Claim | Why retracted | Where |
|---|---|---|
| W(1/10) ⇒ 0.865; W(1/5) ⇒ 0.858 or 6/7; W(17/100), δ₀ = 1/1000 ⇒ 6/7 | The W(w) model was wrong: case-A energy, the canonical diagonal r + z < 1, and α are not relaxed. The arithmetic was sound. | [02](../02-strengthening) |
| Gram split ⇒ 0.8665 / 0.871 / 0.864 | F1: the pairing is against π = ρ·ζ_m, not ρ | [04](../04-gram-split) |
| R4: "θ = 6/7 ⟺ the argument closes with h = 6/7" | Refuted by the referee; only θ ≥ 1 − h/6 holds | [03](../03-six-sevenths) |
| Quartic transplant ≈ 0.854 (and ≈ 0.81 theta route) | It used the wrong N(h)² conductor; the true N(h)³ gives θ = 1 | [02](../02-strengthening), [03](../03-six-sevenths) |
| "Zero-slot gap" in `lem:plain` (weighted audit) | `lem:plain` item 1 has no length restriction | [08](../08-weighted-numerator-10499-12000)/audits |
| `c388e88e6` called an "original OpenAI commit" | It is a commit of the fork akashlevy/math | [08](../08-weighted-numerator-10499-12000)/history |

Errata in working files, from the referee reports:
- `04-gram-split/REPORT.md`: "t = 1 → 0.8751" should be 0.877911; "0.86624" should be 123/142.
- `04-gram-split/V_derivation.md` §0(c) is wrong: square-μ main terms are O(Z^{C(5/6)}).
- `04-gram-split/LoA6.md` has the normaliser wrong (it should be ζ_K(6w−3)/ζ_K(6w−2)). Brubaker–Bump does not apply because μ₁₂ ⊄ Q(ω).
- `04-gram-split/vls/B_bilinear.md` overclaims "(V′_β) ⇔ (V′_LS)"; only ⇐ holds.
- `02-strengthening/conditional_certificate.py`'s defaults encode the retracted W(1/10) claim and fail when run directly.
- `02-strengthening/optimize_w.py` uses the superseded model; use `audit/optimize_corrected.py`.
- `03-six-sevenths/SIX_SEVENTHS.md` §5 states the quartic closing threshold loosely as "H > D^{4/3}". The scripts give D^{15/11} (dF) and D^{11/8} (BGL), so the conclusion holds a fortiori.
