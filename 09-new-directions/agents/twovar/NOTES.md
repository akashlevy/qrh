# twovar: splitting the cubic triple (1,1,1) over several Poisson variables (10 Oct 2026)

Labels: **PROVED** (finite algebra), **R-** (rigorous modulo standard theory not written out),
**MODEL** (exponent bookkeeping in P78 Part II shape; the LP is exact for the stated pieces),
**HEUR** (assumed by analogy). Files are listed in sec 8.

## 0. Verdict

**Dead.** Every split with k ≥ 2 Poisson variables has a GRH-for-rows floor above 1 in P78 shape.
Counts can only raise θ, because the floor bin of the counts model is the GRH constraint. So θ with counts is also > 1.

| k | triple split | B-row (degree, c_B) | low-side base 1/2 + k/6 | GRH floor, unmarked, perfect energy = 1/2 + 5k/18 | **GRH floor (model)** | θ with counts |
|---|---|---|---|---|---|---|
| 1 | (P, θ⊗θ) | γ₁², deg 4, **c_B = 4** | 2/3 | 7/9 | **10/11** (13/17 if c_B were 2) | 10/11 (0.797 if c_B = 2; alt-detect) |
| 2 | (P, P, θ₃) | γ₁ theta, deg 2, c_B = 2 | 5/6 | 19/18 | **13/12** | ≥ 13/12 (vacuous) |
| 3 | (P, P, P) | λ̄η Hecke, deg 1, c_B = 1 | **1** | 4/3 | **4/3** | ≥ 4/3 (vacuous) |

The trade-off is lopsided:
- **Lowering the row degree** only relaxes the energy (the wall). At k ≥ 2 the energy is not what binds: perfect energy still gives 19/18 and 4/3.
- **Each extra Poisson variable** costs (ZP)^{1/2−c} = Z^{1/6} in signal/low (sec 2.3). This is exact normalisation bookkeeping and does not depend on the energy, Gram or row bounds.
  - With GRH rows and the Gram exceptional term, each variable adds 5/18 to the floor.

**What kills k = 2 (MODEL).** The vertex is l₁ = l₂ = 0, y₁ = y₂ = 3/4 (M = 3/2), where three constraints meet:
- the low side, = base 5/6 plus the Gram exceptional (cube-frequency) terms (Y_i/X_i)^{1/3} of both A-polynomials, giving 5/6 + 1/4;
- the generic GRH rows, 11/6 − M/2 − (l₁+l₂)/6, numbering F₁F₂ = Z²/(X₁X₂);
- the energy, which becomes diagonal at M ≥ 3/2 through the Poisson node, and costs more for M < 3/2.

Removing the Gram exceptional term would give 5/6. Reaching 7/8 needs that term cut to 13% of its size (`runs/gram_scan.out`). It cannot be cut: it is the signal's shadow, the same coherent s-sum (sec 2.4).

**What kills k = 3 (MODEL, robust).** The low side alone forces β > 1 + Σ l_i/3 ≥ 1, for any geometry, marks, energy or rows.

## 1. The exact structure (PROVED, algebra)

Let κ_A(n) = (n/A)₃ be completely multiplicative in n. Then any split of the triple (1,1,1) into k Poisson variables plus a row is a **tensor product**

  I_k = Σ_A b(A) q_A^{−1/2} Π_{i=1}^k S_i(A),   S_i(A) = Q_i^{−1/2} Σ_{m~Q_i} A^{(i)}_m κ_A(m),

where A^{(i)}_m = Y_i^{−3/2} Σ_{s~Y_i} c_i(s) g(ψ_s, −m) are cubic Gauss-sum polynomials. This holds because B_{m₁…m_k} = Σ_A b(A) q_A^{−1/2} Π κ_A(m_i).

**Absorption.** Poisson in each m_i (modulus s_iA, Gauss sum γ(ψ̄_s κ_A) = γ(ψ̄_s)γ(κ_A)·ψ̄_s(A)κ_A(s)) gives:
- the cross term = 1, by cubic reciprocity;
- g(ψ_s)γ(ψ̄_s) = N(s)^{1/2}ψ_s(−1), so the s-sum is coherent.

The row therefore needs b(A)·γ(κ_A)^k = μ(A)·Hecke(A). Using γ₁³ = μλ:
- k = 1: b = γ₁²;
- k = 2: b = γ₁, the single cubic theta;
- k = 3: b = λ̄.

Each S_i(A) is *exactly* the one-variable P78-type object. So its principal part is the one validated in the k = 1 template. **The k-variable probe is k copies of the k = 1 Poisson step, multiplied.**

**Corollary (HEUR, structural).** S₂(A) = F₂^{−1/6}γ(κ_A)·(const) + (non-principal part). Hence I₂ = F₂^{−1/6}·I₁[θ⊗θ] + (a probe with no principal frequency).
- Any analysis that keeps the signal sees either the c_B = 4 row (Cauchy–Schwarz in m₁ alone) or pays the joint Cauchy–Schwarz (sec 2).
- The second Poisson variable just manufactures the second Gauss sum of θ⊗θ inside the row.

## 2. Exponent derivation (MODEL; units log Z)

**Geometry.**
- X_i = Z^{l_i}, Y_i = Z^{y_i}, Q_i = X_iY_i, M_i = l_i + y_i, M = Σ M_i.
- Marks P = Z^ℓ on A (Z → Zp), as in P78.
- Frequency ranges F_i = ZP/X_i, with h_i = 1 + ℓ − l_i.

**2.1 One Poisson step (as in P78, re-derived).**
- Σ_{m~Q} ψ̄_s(m)κ_A(m) = Q·C^{−1/2}·γ(ψ̄_s)γ(κ_A)·Σ_{q_u ≤ F} κ̄_A(u)ψ_s(u)…, with C = Y·ZP.
- With the normalisations Q^{−1/2}Y^{−3/2} and the coherent s-sum (≈ Y at the copies), the principal part is S(A) ∋ (X/(ZP))^{1/2}·F^{c}·γ(κ_A) = F^{c−1/2}γ(κ_A).
- c = 1/3: copies u = a³ ≤ F, count F^{1/3}.
- Check at k = 1: signal = F^{c−1/2}P^{1/2}Z^{β−1/2} = X^{1/2}F^{c}Z^{β−1}. This is P78's X^{1/3}Z^{β−5/6}P^{1/6} for c = 1/6 ✓.

**2.2 Signal (k variables).**
- **S = Π_i F_i^{c−1/2} · P^{1/2} · Z^{β−1/2}.**
  - P^{1/2} is the coherent slot sum Σ_p p^{−1/2}.
  - Z^{β−1/2} is Σ_{A~Z} μλη(A) q_A^{−1/2} at a zero β.
- **Copy sharing.** Each variable has its own copies u_i = a_i³ ≤ F_i, so the copy count is Π F_i^{1/3}. Two kinds of row are not copies:
  - **Half-copy rows**: u₁u₂ a cube but u_i not cubes. The A-row is principal, so these carry the zero, but they lack the s-poles. They are smaller than the signal by Π Y_i^{−1/2}, and never bind.
  - **Partial-copy rows**: u_i a cube for i ∈ J only. These have s-poles for i ∈ J. They go into the row side (2.5).

**2.3 Low side: Cauchy–Schwarz in (m₁,…,m_k).**
- The bound is |I| ≤ Π Q_i^{−1/2} Π‖A^{(i)}‖ (Σ_m d_k(m)|B_m|²)^{1/2}.
- The divisor-weighted energy over m ~ ΠQ_i (family M, row length 1+ℓ) is ΠQ_i·E, with mean square O(1) after the marks. So

  **L = Π_i X_i^{1/2} G_i^{1/2} · E^{1/2}.**

- **Per-variable cost (MODEL, robust).** S_i/L_i = F_i^{c}(ZP)^{−1/2}G_i^{−1/2} ≤ (ZP)^{c−1/2}.
  - This uses F_i ≤ ZP. If X_i < 1, the planar Gram term ‖A‖² ≥ Y/X gives the same bound.
  - Hence β > 1/2 − ℓ/2 + Σ_i (1/2 − c)(1+ℓ) + Σ c·l_i + ½ΣG + ½E.
  - For cubic this is **5/6 − ℓ/6** at k = 2 and **1** at k = 3; the k = 3 bound is independent of ℓ.
  - So the mark amplification is a_k = 1/2 − k/6: 1/3, 1/6, 0.

**2.4 Gram norms of the cubic Gauss-sum polynomials.**
- ‖A^{(i)}‖² ≪ (Q_i/Y_i)(1 + P_{a,i}^{1/3} [+ P_{a,i}²/Y_i]), with P_a = Y/X.
- Poisson in m gives dual v ≤ P_a. At cube v the dual s-sum is |Σ_s c̃(s)|², the **same** coherent sum that produces the signal's s-pole. So the exceptional term (Q/Y)P_a^{1/3} is the signal's shadow (as in NOTES §2/§4).
- Non-cube v are negligible for smooth s-weights.
- The P_a²/Y piece is optional (gq). It never binds.

**2.5 High side (GRH for rows).**
- **Rows.** Row (u₁,…,u_k) is Π_i L(w_i, ψ_{u_i}ξ_i) / L(x, λη·κ̄_{u₁⋯u_k}) times the copy ζ-factors and the slot sum Σ_p w_p λη(p) κ̄_{u₁⋯u_k}(p) p^{z−1}.
  - Ranges: q_{u_i} ≤ F_i.
  - The distinct denominators are indexed by u₁⋯u_k mod cubes, up to Π F_i = (ZP)^k/ΠX_i. This is P78's F times (ZP/X₂)⋯: the theory agent's "extra factor ≈ Z".
- **Per-variable row sizes.**
  - Generic variable (GRH numerator, Y_i^{1/2}): X_i^{1/2}Y_i^{−1/2}(ZP)^{−1/2} per row, F_i rows, total (ZP/Q_i)^{1/2}.
  - Copy variable: F_i^{c−1/2}.
- **Constraint.** For every J ⊊ {1..k} (J = ∅ generic, J ≠ ∅ partial-copy):

  σ ≥ 1/2 − ℓ/2 + Σ_{i∉J} [ (1+ℓ−M_i)/2 + (1/2 − c)(1+ℓ−l_i) ].

  - k = 1 is exactly the template / P78 GRH row constraint ✓.
  - k = 2, J = ∅ is 11/6 − M/2 − (l₁+l₂)/6 (unmarked).
  - J = {1} is 7/6 − (2/3)l₂ − y₂/2: the one-variable cubic constraint in variable 2.
- **Marks** grow the rows by P^{k/2} and the signal by P^{1/2 − k/6}. The device ratio is 5 at k = 2, which is > 4 (new-architectures criterion), so marks hurt. LP: ℓ = 0 is optimal; ℓ = 0.1 gives 1.094 and ℓ = 0.2 gives 17/15.
- **Counts version (HEUR, not needed).** Take a zero at a = (1+δ)/2, contours x = a and w_i = 1−a, and charge each numerator its reflected cost q_{u_i}^{δ/2}. Then

  σ ≥ 1/2 + δ/2 − ℓ/2 + xδℓ + Σ_i [ −(1+δ)y_i/2 + h_i(R + δ/2 − c) ],

  with R counted in base U = Z^{h₁+⋯+h_k}. This reduces to P78 `common-high-exponent` at k = 1, and to the GRH row at δ = 0, R = 1. So θ_counts ≥ θ_GRH.

**2.6 Energy E** (MODEL, alt-detect's pluggable nodes).
- The nodes are min over {direct cubic LS, Poisson, reflect (conductor (QP)^{c_B}) + cubic LS, reflect + Poisson} of the max of their pieces, with the optional slot piece (M'+4l'−2)/3.
- k = 2: c_B = 2 (P78 j_p = 2, R-). k = 3: c_B = 1, giving e ≡ 0 (no wall), as the coordinator expected.
- The compensation term (d = ℓ) is modelled at (X_i/P, Y_i/P) with relative weight P^{−(k+1)/2} (HEUR for k ≥ 2). It never binds.

## 3. The LP (`kvar_lp.py`)

**Variables:** l_1..l_k, y_1..y_k, ℓ, σ.

**Constraints:**
- (low, d ∈ {0, ℓ}) the base of 2.3, plus ½Σ_i g_i with g_i ∈ {0, (y_i−l_i)/3, gq(2b_i − y_i + d)}, plus ½e for every piece e of the chosen energy node at (M' = M − 2kd, l' = ℓ − d);
- (rows) 2.5 for every J ⊊ [k];
- (side) l_i ≥ ℓ, y_i ≥ ℓ + margin, M − 2kℓ ≥ margin, F_i ≥ 1.

The node choice is enumerated, and scipy HiGHS solves each LP.

## 4. Validation (`validate.py`, `runs/validate.out`)

At k = 1 the LP agrees with the reference to 1e−9 in every case:
- sextic P78 floor: **13/15** ✓ (template_lp);
- sextic W(0.1): 0.865 ✓;
- cubic one-variable template, with the third piece:
  - c_B = 2: **13/17** ✓;
  - c_B = 3: 0.859615 ✓;
  - c_B = 4: 10/11 ✓;
- cubic, without the third piece: 3/4 ✓ (all against alt-detect/p78_cubic_lp).

`closed_forms.py` checks θ = 1/2 + k(1/2 − 2c/3) against the LP (unmarked, perfect energy). All agree:
- 7/9, 19/18 and 4/3 for cubic k = 1, 2, 3;
- 8/9 for the sextic copies at k = 1;
- 3/2 for sextic (1,1,4) with k = 3.

For unequal c_i the partial-copy rows bind, and the LP lies above the formula.

## 5. Results for k = 2 (`variants.py`, `runs/variants.out`; all MODEL, GRH rows)

| configuration | θ_GRH |
|---|---|
| default (third piece, marks, compensation, partial rows) | **13/12** (M = 3/2; degenerate face, e.g. l = 0, y = (3/4, 3/4)) |
| no third piece / + gq / no marks / symmetric / big box / margin 0 | 13/12 |
| energy only below the reflection wall (nodes reflect, reflect+Poisson) | 109/90 = 1.211 |
| perfect energy (no wall, no excess) | 19/18 |
| Gram exceptional ×0.5 / ×0.25 / ×0.1 / removed | 41/42 / 0.910 / 0.866 / 5/6 |
| Gram exceptional removed on variable 2 only | 49/54 = 0.907 |

## 6. Variants (task item 3)

- **(a) Y₁ ≠ Y₂.** The optimum face already contains very asymmetric points, e.g. y = (0.071, 1.357), and the floor stays at 13/12. Fixing y₂ ∈ {0.1, 0.25, 0.5}: 13/12.
- **(b) Short second variable.**
  - Q₂ tiny (l₂ = y₂ = 0.02): 1.143. The variable-2 rows (ZP/Q₂)^{1/2} are then unsuppressed.
  - l₂ = 0.75: 1.099. l₂ = 1: 47/40.
- **(c) m₂ smooth, no s₂-polynomial.**
  - If m₂ is Poisson-transformed, it is the y₂ = 0 case: **11/10**. The Gauss sum γ(κ_A) still comes from Poisson in m₂.
  - If it is *not* Poisson-transformed (Q₂ < Z, a short character sum), no third Gauss sum appears. The row then has γ₁·γ₁ = γ₁² and no μ, so there is no detection (PROVED: no μ-factor).
  - Supplying the third factor from the row instead is k = 1 θ⊗θ, with c_B = 4 → 10/11.
- **(d) One Poisson plus Gauss-sum orthogonality.**
  - Make m₂ complete (F₂ ≤ 1, Q₂ ≥ Y₂ZP). Character orthogonality then forces s₂ = A, and the third γ₁ is the polynomial's own g(ψ_A).
  - The cost is S₂/L₂ = (ZP)^{−1/2} (the F₂ = 1 end of 2.3), so β > 7/6 already on the low side. LP: **47/40** (71/60 with the third piece).
- **Shared s.** Take both Poissons in one s-sum, A_{m₁,m₂} = Σ_s c(s)g(ψ_s,−m₁)g(ψ_s,−m₂).
  - The copies are then pairs with u₁u₂ a cube, but rows and low each lose Y^{1/2} against two independent s-sums. Worse (HEUR).
- **Sextic analogues (cheap, MODEL).** The copy exponent is 1/ord(χ^e) per variable:
  - three Poisson variables with a Hecke row:
    - (1,1,4): 325/216 ≈ 1.505;
    - (2,2,2): 4/3;
    - (1,2,3): 37/30;
  - (1,2) + ray, two Poisson variables: 1.118.
  - The low bases are 4/3, 1, 1, 1. A quadratic (e = 3) Poisson variable has c = 1/2 and costs nothing on the low side, but γ₃ is a ray quantity, so it is pointless.
  - The thin two-Poisson split (1,3 | θ₃, quadratic residual) gives 1.016.
  - As the coordinator corrected, (1,1,4) with a θ row is thin; it was not run.

## 7. Bottom line

The (Poisson, Poisson, θ₃) split does achieve its aim: a single-Gauss-sum, full-index, c_B = 2 row (R-). But the second Poisson step is bought at an exact normalisation price:
- (ZP)^{−1/6} on signal/low;
- F₂ extra rows on the high side;
- a second Gram exceptional term.

Together these give 1/2 + 5k/18 even with perfect energy. The energy wall, which was the only thing c_B = 2 was meant to fix, does not bind at k ≥ 2. The cubic family has no split with floor < 7/8. Its best case is k = 1 with c_B = 4, at 10/11.

## 8. Files

| file | content |
|---|---|
| `kvar_lp.py` | k-variable P78-shape GRH LP (signal, joint C–S low with per-variable Gram, energy nodes, partial-copy rows, marks, compensation) |
| `validate.py`, `runs/validate.out` | k = 1 validation against template_lp.py / p78_cubic_lp.py |
| `variants.py`, `runs/variants.out` | all k = 1, 2, 3 runs, variants (a)–(d), sextic analogues |
| `closed_forms.py`, `runs/closed_forms.out` | θ = 1/2 + k(1/2 − 2c/3) vs LP |
| `binding_k2.py`, `runs/binding_k2.out` | constraint values at the k = 2 vertex |
| `runs/gram_scan.out` | k = 2 floor vs Gram exceptional factor (7/8 needs factor 0.13) |
