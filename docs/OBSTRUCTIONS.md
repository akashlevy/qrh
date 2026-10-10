# Obstructions and barriers

Every barrier found, grouped by what it blocks. A barrier here always bounds what a **stated package of estimates or lemmas** can prove. None of them says anything about where zeros actually are.

The status labels are:
- **exact**: a rational certificate or LP dual;
- **proved**: a written proof;
- **model / heur**: an exponent-level model or heuristic.

## 1. Ceilings of the 7/8 architecture (OpenAI's cubic-theta probe)

| Barrier | Statement | Status | Where |
|---|---|---|---|
| Lemma-set floor | With P78's own lemmas and any geometry, θ ≥ 33216105/37963129 ≈ 0.8749570. Liu's optimum (1507−2√921)/1653 ≈ 0.874957069799 lies 8·10⁻⁸ above it. | exact | [01](../01-seven-eighths-lemma-set) `barriers.py` |
| κ-feedback ceiling | Geometry + feeding proved bounds into κ cannot certify θ < θ\* ≈ 0.874957019420099, the root of 7884θ³ − 18819θ² + 14643θ − 3686 = 0. This is the public "cubic-root" endpoint. | exact, sharp for the 3-constraint relaxation | [07](../07-method-barriers) |
| **13/15 floor** | Part II gives θ ≥ 13/15 for every geometry and **every row-count input, GRH for rows included** (Theorem D). With the paper's detector floor 51/100 the floor is 167/192. The density hypothesis for rows already attains 13/15. | exact | [01](../01-seven-eighths-lemma-set), [03](../03-six-sevenths) R3 |
| 13/15 low-majorant barrier | The unchanged direct-Gram and reflected-row low majorants force B ≥ 13/15 (also in Lean as `physical_low_barrier`). This still caps the weighted-numerator route. | exact | [07](../07-method-barriers) |
| Energy-wall low side | With the paper's case A, W(w) gives low ≥ 13/15 − w/60, so 6/7 needs w ≥ 4/7. | exact | [02](../02-strengthening) |
| Row-range bound | θ ≥ 1 − h/6 always (d = 0 low constraint). The maximal h is 4/5 for the paper's lemmas and 6/7 with a Gram fix. | exact | [03](../03-six-sevenths) R4 |
| 5/6 structural floor | Without prime amplification the sextic probe has θ ≥ 5/6 + l_x/6. The residual-Cauchy–Schwarz template has σ₀ ≥ 5/6 + (l_x − ℓ)/6 + (g+e)/2, even if the theta pairing vanished. | proved | [01](../01-seven-eighths-lemma-set) §6.3, [04](../04-gram-split) |
| **No route to 1/2** | Below 13/15 the bookkeeping needs ℓ → ∞, i.e. mollifier-length energies with power-saving asymptotics, which is GRH-strength (a consistency argument). The energy wall has three sources: the reflection dual length; the BGL (QN)^{2/3} term, which is the open cubic/sextic large sieve conjecture; and sixth-power rows NQ^{1/6}. | §6.1 exact, §6.3 heuristic | [01](../01-seven-eighths-lemma-set) `QRH-analysis.md` §6 |

## 2. Structural main terms on the low side

- **The Gram P_a^{1/6} term is a genuine main term.**
  - It is the energy of one rank-one theta vector; the exceptional frequencies k = a⁶e give it.
  - It cannot be *cancelled*, only bounded, because the signal has the same Mellin coefficient K̂(5/6).
  - de Faveri–Dunn–Hoffstein show the cubic large sieve's (AB)^{2/3} term is sharp, and the conjectured n = 6 truth reproduces exactly this excess.
  - Sources: [04](../04-gram-split); [03](../03-six-sevenths) R4 §6.
- **Case A** (1+2ℓ′−2M′)/3 after the refined count: a pure-cube dual where the mark sum does not oscillate, plausibly genuine. [03](../03-six-sevenths)
- **The energy-wall excess** is the column term of the quadratic large sieve, sharp for general coefficients. [03](../03-six-sevenths)

## 3. Mean values below the large-sieve diagonal (M_κ, κ < 1)

All of these explain why 6/7 (Theorem C1) is unproved.

| Barrier | Statement | Status | Where |
|---|---|---|---|
| P2's dual target is false below the diagonal (Prop O1) | For 𝓗 = X^{2−κ}, 𝓔/X ≫ X^{1−κ}. This is GRH-conditional in general and unconditional for κ < 1/4. The loss is the absolute value on μ(f). | proved | [03](../03-six-sevenths) R1 |
| Move calculus | Under the available moves P and R, \|ε − δ\| = 1 − κ is invariant. Every reachable path has max min(δ, ε) = −(1 − κ), even with optimal large sieves. The only escape is a "column functional equation" for μ, which is GRH-type. | proved, as a model of the moves | [03](../03-six-sevenths) R1 |
| MDS signature | The first moment is finite type (C₂). The second moment is **affine** (det 0), the "beyond the edge" signature. | heuristic | [03](../03-six-sevenths) R1 |
| Coefficient-blind impossibility | Any large-values bound valid for all \|c_n\| ≤ 1 allows ≥ H^{1/6} large rows (the copies u = a⁶). This rules out the large sieve, Halász–Montgomery, Huxley and Guth–Maynard. The multiplicity-beating LV hypothesis is *equivalent* to the conclusion. | proved | [03](../03-six-sevenths) R3 |
| Target-copy detection | At u = a⁶ only inverse (Möbius) witnesses of length ≥ U^{1−o(1)} detect zeros. R3 at the target *is* M_κ. | sketch | [03](../03-six-sevenths) R3 |
| Bootstrap consistency | One target with a zero at β ≤ 13/15, plus GRH elsewhere, is consistent with every lemma conclusion of both papers. At σ = 6/7 one needs N_fam < H^{1/6}, while DH gives H^{2/7}. | exponent-level | [03](../03-six-sevenths) R2 |
| Prime-sixth-power extraction floor | max((6+5κ)/12, 1 − κ/6) ≥ 6/7, whatever the mean square. Only Lemma B (bootstrap) removes the D²/Y² term. | proved | [03](../03-six-sevenths), [07](../07-method-barriers) |
| M_{6/7} is stronger than its conclusion | It contains row u = 1, giving 13/14. It is equivalent to N(σ;Q) ≪ Q^{(13−14σ)/6}, beyond DH. | proved | [03](../03-six-sevenths) |
| Signed-dual unit block | The unit j-block needs a coefficient-sensitive O(D^{2+ε}) bound for the coprime-pair form. Positive or absolute completion loses D/H. | identified gap | [07](../07-method-barriers) |

## 4. Gram/theta-split obstructions (removing P_a^{1/6})

All in [04](../04-gram-split).

**The split itself**
1. **F1, the reflection scalar.** The pairing is against π = ρ·ζ_m, not ρ. This was fatal to 0.8665 / 0.871 / 0.864.
2. **Norm-only obstruction.** ‖ρ‖² ≤ Q^{5/6+ε} and ‖B‖² ≤ Q^{1+ε} give exactly Q^{11/12}, and this is sharp for norm-only arguments.
3. **Probe rigidity.** Zero detection (j = 1) and a ζ-free row (j ∈ {2, 5}) are mutually exclusive in the template; a scan of 50 designs confirms it.
4. **Signal circularity.** T contains the whole zero-detecting signal (c = 1). So (V′) must be signal-blind, which makes it a Möbius-randomness (zero-free-region) statement.

**Missing tools**

5. **Conductor coupling.** The A₂ dual has conductor N(μ)³ and coupled coefficients ρ(sμ^k), so no large sieve applies. A saving would need k < 1.53.
6. **Hyperbolic second moment.** The reflection group is Coxeter (2,6,6) for k = 3, so Bochner + Phragmén–Lindelöf gives only the trivial bound.
7. **β-structure used up.** Arguments using only |π| and β are capped near 0.89·Q^{11/12}.
8. **Beyond the edge.** (V′_LS) needs a large sieve with family L ≍ Q^{0.79} against conductor about L³, with no automorphic kernel. This is not known even for GL₃ cusp forms.
9. **Spectral tools.** There is no GL₃^{(6)} Kuznetsov and no metaplectic spectral large sieve. Louvel's Kuznetsov lives on the cubic cover.
10. **μ₁₂ ⊄ Q(ω).** Brubaker–Bump, BBF/BBFH, Chinta–Gunnells and DIPP all assume μ_{2n} ⊂ F, so every n = 6 multiple-Dirichlet-series statement is heuristic here.
11. **No Type I/II structure** in m (Kazhdan–Patterson non-uniqueness).

## 5. Other low-side no-gos

From `independent_low.txt`, in [07](../07-method-barriers):
- An auxiliary residue character with ξ⁶ ≠ 1 removes the z = 1/6 pole.
- Varying-Q kernel cancellation kills the signal.
- Sums of separately bounded probes cannot beat min_j(L_j − c_j).
- A direct reciprocal-L per-column strategy needs λ < 2(1 − B) ≈ 0.25, but has λ ≈ 4/5.

## 6. Other fields

- **Quartic transplant (Q(i)) is refuted.**
  - The row twist enters as h³, so the reflection conductor is N(h)³ (Kubota–Diaconu functional equation).
  - The dual twist stays quartic: quartic Gauss sums are Eisenstein, not theta, coefficients.
  - There is no polar term.
  - The transplanted Part II gives **θ = 1**. In P2's architecture, closing the argument needs H ≥ D^{15/11}, while the extraction needs H < D^{4/3}.
  - Sources: [03](../03-six-sevenths) §5 and `quartic-transplant/`.
- **The sextic case is unique.** Davenport–Hasse linearization of Möbius × Gauss sum works only for odd m in order 2m. The quintic theta is unknown. Larger fields containing μ₆ give the same exponents. [04](../04-gram-split) `LIT_ideas234.md`.

## 7. Formalization obstructions

All in [05](../05-lean-23-24).

| Obstruction | Detail |
|---|---|
| Comparator proof-term identity | Restated constants are compared including their proof terms. Deleting declarations renames shared `_proof_N` auxiliaries, so modules restated by a challenge must be exempt from pruning. |
| Hard-coded descent constants | 1/20, 1/40, 1/1000 and 3000 appear in 66 declarations, which blocks a parametric 11/12. |
| Hard-coded 7/8 geometry | 19 Lean files hard-code 17/48, 23/48 and 13/16, so the 34999/40000 re-optimization is not mechanically checkable without edits. |
| Stale build objects | Pruning leaves 253 stale `.olean` files that can silently satisfy imports. |
| Automation golf fails | One-tactic replacement saves 0.4%. Merging α-equivalent duplicates saves about 157 lines. |
| Paper-form descent rewrite | At most about 15% saving for weeks of work; declined. |
| macOS verification | Landrun and NanoDa sandboxing are unavailable, so the runs were unsandboxed. A Linux sandboxed rerun is needed for adversarial-grade trust. |
