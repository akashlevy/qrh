# Attack A: functional-equation geometry for (V′)

*9 Oct 2026. Target: `TARGET_VLS.md` (V′), (V′_LS), (V′_β). Builds on `wf2/A2_theory.md` (§4 and referee report), `wf/V_gaps.md` (Lemma Z, §1.6), `V_derivation.md`, and paper78 (eq. `reflection` l.1763, `reflection-row-phase` l.1832, `signal` l.965, Poisson identity l.3831).*

**Scripts** (all in `gram_split/vls/`, run with the project venv):

| script | purpose |
|---|---|
| `groupoid_search.py`, `groupoid_bfs.py` | Heckenberger/Sawin–Whitehead Weyl groupoid; DFS and shortest-word BFS enumeration of chambers |
| `orientation_test.py`, `model_check.py` | which encodings reproduce the paper's *proven* reflection rule and the signal (§2) |
| `chamber_analysis.py` | residue-adapted chambers, lengths of all summation variables, coefficient structure, signal-line test, bounds (§3–4) |
| `u_full_table.py` | full 24-chamber table of the exact 3-node model U (§3.1) |
| `m1_deep.py` (log `m1_deep.log`) | deep search of the faithful 4-node object, up to 10 476 chambers (§3.2) |
| `design_scan.py` | 50 "modified probe" variants: zero detection vs. shortening (§5) |
| `lp_loa_threshold.py` | LP consequence of the Mellin-shift / Lindelöf-on-average reformulation (§6.3) |

**Labels:**
- **[R]** rigorous, given the cited results;
- **[R-]** routine but not written out;
- **[RC]** rigorous as a finite combinatorial or linear-algebra computation (the script is the proof);
- **[H]** heuristic;
- **[G]** gap;
- **[N]** numerical evidence.

---

## 0. Verdict

**Best achievable saving from functional-equation geometry: s₀ = 0. The LP then gives θ = 0.874957, Liu's value.** Proof status of "no saving": [RC] for the exact model U, plus the labelled [R-]/[H] items in §8. A genuine non-FE input is still required (§7).

1. **The exhaustive search finds exactly four functional-equation images of T. All four were already known.** [RC]
   - In the exact 3-node model U = (s, m, A), the variable t is eliminated by Patterson's identity: the column carries γ₂(A).
   - U has 24 chambers and 4 objects. Exactly 8 chambers are residue-adapted, and they fall into 4 representation types:
     - **R0** = unreflected;
     - **R1** = the paper's reflection;
     - **R2** = A2_theory's m-Kubota dual ρ(sμ³);
     - **R3** = the m-Voronoi, which contains the signal.
   - Past a Möbius node there is no functional equation, so this list is complete.
   - The faithful 4-node encoding (t, A, m, s) has infinite type. A search of 1 136, 3 352 and 10 476 of its chambers finds the same four types and nothing else. [N, combinatorial]
2. **None of the four beats Cauchy–Schwarz for M ∈ (2/3, 1).** Exponents are beyond Q^{5/12}, in Z-units, at M = 0.826; the CS target is M/2 = 0.413.

   | rep. | effective lengths | best available bound | vs CS |
   |---|---|---|---|
   | R0 | (Z, Q) | 0.554, sextic large sieve (de Faveri) | worse |
   | R1 | (L, Q) | 0.413, quadratic large sieve | = CS |
   | R2 | (L³/Q, L), coupled | 0.891, trivial | worse unless M < 2/3 |
   | R3 | (Z, Z/Q), coupled | 0.587, trivial | worse, and the contour must cross the coherent pole line, whose residue is the Möbius signal S |

3. **Correction to `wf2/A2_theory.md` §4.2.** [RC]
   - Its 4-node object M0 = [[5,4,0,0],[4,0,1,0],[0,1,0,5],[0,0,5,4]] ("honest, 336 chambers, 7 objects") encodes the *conjugate* cubic theta.
   - Its column reflection gives a sextic kernel and leaves the row node bare. That contradicts the paper's proven rule B_p = χ_p^{−j−2} and Lemma Z.
   - The encoding that reproduces the paper is M1 = [[1,2,0,0],[2,0,1,0],[0,1,0,5],[0,0,5,4]]. It is **not** a finite arithmetic root system: coefficients grow without bound.
   - M0 has 8 residue-adapted types. Five have an empty free variable, and two of those avoid the signal line: lengths (L, ∅) and (∅, ∅).
   - These would be *false positives* for (V′); they belong to the anti-aligned pairing.
4. **The obstruction is structural (§5).** [RC for the scan; interpretation H]
   - I scanned 50 modified-probe variants (row twist exponent j, column Gauss type, theta orientation).
   - Four are viable probes: a cubic-theta column together with a row twist that is coherent for the row theta. They form a 2×2.
     - The two Davenport–Hasse-aligned variants, the paper's and its global conjugate, are exactly the ones whose signal column is Möbius, i.e. which detect zeros. For both the geometry is rigid: 4 images, no shortening.
     - The two anti-aligned variants have short images, empty in both variables, but detect no zeros.
   - Zero detection (the Möbius node) and shortening are mutually exclusive in this formalism.
5. **Partial and smoothed applications give nothing (§6), with one exact reformulation.**
   - Applying the functional equation to dyadic pieces, splitting the reflected and unreflected ranges, or Mellin-separating a kernel leaves the critical block μ ≍ L, m ≍ Q unchanged.
   - Mellin-shifting the row variable to the centre z = 5/12 of A2_theory's residual functional equation gives:
     - **(V′) ⟸ (LoA_π):** Σ_{N(μ)≍L} |Π(5/12+it; μ)|² ≪ L^{θ_π+ε}(1+|t|)^A, with θ_π < M/(2M−1) (= 1.267 at M = 0.826), where Π(z;μ) = Σ_m π(m)(μ/m)₂ N(m)^{−z};
     - the saving is then s₀(M) = (M − (2M−1)θ_π)/2;
     - Heath-Brown's quadratic large sieve gives only θ_π = 3/2, which is sharp for general coefficients;
     - LP (§6.3, conservative marked version): θ beats Liu iff θ_π ≲ 1.31. For example θ_π = 1.2 gives 0.8688, θ_π = 1.1 gives 0.8626, and θ_π = 1 (Lindelöf on average) gives 0.8548.
     - This is a *reformulation*, not a proof: LoA_π is a beyond-the-edge second moment.

---

## 1. Set-up: the residual function, lengths, and the signal line

### 1.1 The residual function [R-]

With the Patterson identification of the column (squarefree A; non-squarefree blocks as in `V_gaps.md` §4–5),

  T = Σ_{A,m} γ₂(A)ᾱ(A) q_A^{−1/2} χ_A(m) ρ(m) V(q_A/Z) Ω(q_m/Q) ≍ ∫∫ R(u_A, u_m) Z^{u_A} Q^{u_m} (Mellin weights) du_A du_m,

  R(u_A, u_m) = Res_{u_s = 1/6} Z_U(u_s, u_m, u_A).

Here:
- Z_U is the 3-variable series with coefficient g̃(s,−m)·χ_A(m)·γ₂(A), in centred variables u = (Mellin variable) − 1/2.
- The s-residue is the sextic theta pole w = 2/3, i.e. u_s = 1/6.
- In the 4-node encoding the cubic theta pole sits at u_t = 1/3.

### 1.2 Length calculus [RC as linear algebra; H for the effective conductor]

A chamber is a basis E of the root lattice: row E_k is the k-th new simple root in old coordinates. The new centred variables are u′ = Eu. Writing Σ_i u_i ℓ_i = Σ_k u′_k ℓ′_k gives

  **ℓ′ = E^{−T} ℓ**   (ℓ = log-lengths in units of log Z; Q = Z^M).

For a simple reflection at i this is exactly "dual length = conductor/length", with conductor Π_j N_j^{m_ij}:

  ℓ′_i = Σ_j m_ij ℓ_j − ℓ_i,  ℓ′_j = ℓ_j (j ≠ i).

**Residue-adapted chambers.** These are chambers in which the residue root e_s (and e_t in 4-node models) is ± a new simple root. Then:
- the double residue is a simple residue of the new series;
- the free variables' lengths do not depend on the (infinite) residue length.

This is the precise sense in which "T is a sum over two variables" in that chamber. In a non-adapted chamber the residue length Λ appears in the free lengths (§3.1 table), so the representation is not a finite sum.

**[H].** The Heckenberger exponent m_ij = min(Kubota exponent, ord(q_ii) − 1) is used as the effective conductor exponent. The A2 referee (item 1) notes that conductor drops are possible in principle. A drop could only *shorten* a dual, so this is the one point where the negative conclusion is heuristic (§8).

### 1.3 Normalisation of bounds [R as algebra, given unitary coefficients on average]

Functional equations preserve the random size of a sum. Hence in any representation with free lengths N₁, N₂ (unitary coefficients on average):
- random size: Q^{5/12} (representation-independent);
- trivial bound: Q^{5/12}(N₁N₂)^{1/2};
- large sieve (bilinear forms only): Q^{5/12} Δ(N₁,N₂)^{1/2};
- CS target: Q^{5/12}·Q^{1/2} = Q^{11/12}.

Check: in R1 the trivial bound is Σ|β||π| = L^{1/2}Q^{11/12} = Q^{5/12}(LQ)^{1/2}. ✓

When coefficients sit on a sublattice (theta values at sμ³), the sublattice mean-square factor K multiplies the bound; K = 1 is used, and K ≤ L trivially (A2 §2.6).

### 1.4 The signal line [R for the residue, R- for the crossing criterion]

For each A, the m-sum Σ_m ρ(m)χ_A(m)N(m)^{−1/2−u_m} is coherent (j = 1). It has a pole on

  **P1 = {u_m = 1/3}**,

and its residue, summed against γ₂(A), is the Möbius signal S (Davenport–Hasse, `V_gaps.md` §1.6 [R]). The residue is zero-sensitive.

**Crossing criterion.** In a chamber with free variables u′₁, u′₂, write u_m = κ₁u′₁ + κ₂u′₂ + κ₀. If κ₁, κ₂ ≤ 0 (not both 0), the convergence cone {Re u′_k > σ_k} lies in Re u_m < 1/3. The contour must then cross P1, and the representation contains S ("P1 forced").

Any bound obtained from such a representation bounds T − S, not T. That is circular for the probe, because Q^{−1/2}Y^{−1/3}Ŵ·S is exactly the high-side signal (`IDEA1_NOTES.md` §2, `A2_theory.md` §4.3).

---

## 2. Which encoding is faithful? (orientation audit)

### 2.1 Two proven anchors

- **(F-refl)** Paper Prop. `completed-reflection` (l.1763ff) [R]: a column twist χ_p^j gives the dual column factor B_p = χ_p^{−j−2}, and the row acquires τ⁻_{p,j}τ⁺_{p,j+2}. For j = 1 this is the quadratic kernel and γ₁γ₃(m), i.e. Lemma Z.
  - In Sawin–Whitehead (S-W) exponents: new edge 3, and the row node gets (1−3)+(3−3) ≡ 4.
- **(F-sig)** `V_gaps.md` §1.6 [R]: the m-Voronoi turns γ₂(A) into μ(A), i.e. a Möbius node (exponent 3).

### 2.2 Test results (`model_check.py`, `orientation_test.py`) [RC]

| encoding | column reflection gives | m-reflection gives | faithful? |
|---|---|---|---|
| **U** = (s,m,A) [[4,5,0],[5,0,1],[0,1,5]] (A = Patterson γ₂) | [[4,5,0],[5,4,3],[0,3,5]]: edge 3, row node 4 ✓ | column → 3 (Möbius) ✓ | **yes** |
| U₁ (column γ₄, M_AA = 1) | [[4,5,0],[5,0,1],[0,1,1]]: edge 1, row bare ✗ | column → 5 ✗ | no (anti-aligned) |
| **M0** (A2_theory's 4-node) | word `At`: edge 1, row node 0 ✗ | — | **no**: behaves as U₁ |
| **M1** = [[1,2,0,0],[2,0,1,0],[0,1,0,5],[0,0,5,4]] | word `At`: edge 3, row node 4 ✓ | word `ms`: γ₁ × t-residue → Möbius ✓ | **yes** |

**The 3-node scan.** Model (t, A, m) with t-edge e_tA and column twist j; the composite `At` gives dual edge (2 − j) for e_tA = 4:

| j | 1 | 3 | 4 | 5 |
|---|---|---|---|---|
| dual edge (`At`, e_tA = 4) | 1 | 5 | 4 | 3 |
| paper's −j − 2 | 3 | 1 | 0 | 5 |

So M0 has multiplier exponent k = 2 against the paper's k = 4. Equivalently, M0 is the paper's configuration with *conjugated* column.

**Finiteness** (`weyl_groupoid.run`):

| encoding | finite? | chambers | objects | Möbius node |
|---|---|---|---|---|
| M0 | yes | 336 | 7 | none |
| **M1** | **no** (coefficients exceed 40 / 200 after 86 / 139 bases) | — | — | — |
| U | yes | 24 | 4 | yes (truncation point) |

**Consequences.**
- A2_theory §4.2's sentences "the faithful four-node object ... is an honest rank-4 arithmetic root system: 336 chambers, 7 objects" and "Sawin–Whitehead's rationality theorem applies" refer to the anti-aligned (conjugate-theta) pairing. They should be relabelled.
- The faithful object is of infinite type. Its finiteness in the 3-node form U comes only from the Möbius (non-admissible) node, where the signal sits. [RC]
- **[H]: why the S-W dictionary flips.** In the S-W dictionary, the residue of a Kubota node with Gauss character χ^c and edge −c behaves at its neighbour like a theta of the opposite orientation, relative to Bröker–Hoffstein's "τ₃(d) = G₁(1,d), same symbol". I did not resolve this at the level of the function-field dictionary. The audit uses only the two proven anchors.

---

## 3. Exhaustive search

### 3.1 The exact model U: complete list [RC]

`u_full_table.py`:
- 24 chambers, 4 objects:
  - O0 = [[4,5,0],[5,0,1],[0,1,5]];
  - O1 = [[0,1,0],[1,0,5],[0,5,3]], with a Möbius node;
  - O2 = [[4,5,0],[5,4,3],[0,3,5]];
  - O3 = [[0,5,0],[5,4,5],[0,5,3]], with a Möbius node.
- Lengths of all three variables in every chamber, with Λ the residue length (→ ∞):

| chamber | object | lengths of the chamber's nodes 0, 1, 2 (the node carrying ±e_s has the Λ-length) | adapted |
|---|---|---|---|
| id, s | O0 | (Λ, M, 1), (M−Λ, M, 1) | ✓ R0 |
| A, sA | O2 | (Λ, M, 2M−1), (M−Λ, M, 2M−1) | ✓ R1 |
| Ams, sAms | O2 | (5M−3, Λ+5M−3, 2M−1), (5M−3, 6M−3−Λ, 2M−1) | ✓ R2 |
| ms, sms | O3 | (1−M, Λ−M+1, 1), (1−M, 1−Λ, 1) | ✓ R3 |
| m (= paper's high side) | O1 | (Λ, Λ+1−M, 1) | ✗ (m-dual of length YZ/Q at finite Y = Poisson identity) |
| 15 others | O0–O3 | Λ appears in at least one free length, with coefficient ±1 or ±2 | ✗ |

**Free-variable maps** in the 4 adapted types (centred, on u_s = 1/6), with node types (S-W exponents) and residue numerators:

| type | free variables u′ | lengths | κ (u_m-coeffs) | node types | edge | s-residue at |
|---|---|---|---|---|---|---|
| R0 `id` | (u_m, u_A) | (M, 1) | (1, 0) | m: 0, A: 5 (γ₂) | 1 (sextic) | ρ(m) |
| R1 `A` | (u_m + 2u_A, −u_A) | (M, 2M−1) | (1, 2) | m: 4 (γ₁γ₃), μ: 5 (cubic theta) | 3 (quadratic) | ρ(m) |
| R2 `Ams` | (−2u_A − u_m − 1/6, 5u_A + 3u_m) | (5M−3, 2M−1) | (5, 2) | s′: 4, μ: 5 | 0 | ρ(s′·μ³) |
| R3 `ms` | (−u_m − 1/6, u_A + u_m) | (1−M, 1) | (−1, 0) | b: 0, A: **3 (Möbius)** | 0 | ρ(b·A) |

**Cross-checks against independent derivations** (this is what makes U trustworthy):
- R1 is the paper's reflection, with L = Q²/Z, the quadratic kernel, and γ₁γ₃ on the row, which is Lemma Z. [R]
- R2 reproduces A2_theory §2.3's dual coefficients ρ(sμ³) and dual length L³/Q (5M−3 = 3(2M−1) − M). This is an *independent* derivation of A2's conductor N(μ)³ [RC]. It is still [H] as an effective conductor, per the A2 referee's item 1.
- R3 reproduces `V_gaps.md` §1.6:
  - the Möbius column;
  - the dual index n = A·b of length Z·Z/Q = Z²/Q;
  - the forced crossing of the coherent pole, i.e. the signal S.

**Completeness [RC + R].**
- Every chamber reachable by admissible reflections is listed.
- The exploration stops only at Möbius nodes. A Möbius node carries Σ μ(c)χ(c)|c|^{−u} = 1/L(u,χ), which has no functional equation; its poles are the zeros of L.
- So no further functional-equation image exists in the exact model.

**Analytic realisation of the reflections.** Each reflection is realised by a proven or routine functional equation:
- s: Kubota (`FE_LoA6_proof.md` Thm A) [R-];
- m in O0: Poisson [R];
- A: the paper's reflection [R];
- m in O2: rank-one Kubota functional equation of the γ₁γ₃-node, A2 §1.5(b) [R-].

Composites need the continuation of the 3-variable series in the intermediate tubes: [R-]/[G] (A2 §1.5(b), "Bochner, not written").

### 3.2 The faithful 4-node object M1 (infinite type) [N, combinatorial]

`chamber_analysis.py` and `m1_deep.py`:

| search | chambers explored | residue-adapted types | chambers per type |
|---|---|---|---|
| maxcoef 20 | 1 136 | 4 | 4 |
| maxcoef 40 | 3 352 | 4 | 4 |
| maxcoef 80 | 10 476 | 4 | 4 |

The 4 adapted types (each with 4 sign choices of the two residue roots) are exactly R0 (`id`), R1 (`At`), R2 (`Atms`: lengths (2M−1, 5M−3), t-residue at μ, s-residue at μ³·s′) and R3 (`ms`, P1 forced). They have the same lengths and coefficient structure as in U.

Because the object is infinite, completeness is [H]. The stability over a 10× increase in search depth, and agreement with the complete U list, make a missed image unlikely.

### 3.3 A2_theory's M0 (anti-aligned): 8 adapted types, 32 chambers [RC]

| word | free lengths | P1 forced | structure |
|---|---|---|---|
| id | (1, M) | no | sextic kernel |
| At | (2M−1, M) | no | **sextic** kernel, row bare (≠ paper) |
| Atms | (2M−1, M−1 < 0) | no | coupled; one variable empty |
| AtmtAtsmtAt | (−1, M−1) | no | both empty |
| ms | (1, 1−M) | yes | coupled; column *not* Möbius |
| mAtAsmAtA | (1−2M, 1−M) | yes | coupled |
| mAtAmsmAtAm | (1−2M, −M) | yes | both empty |
| AtmtAtmsmtAtm | (−1, −M) | yes | both empty |

Had the search been run on M0 as A2_theory proposed, it would have "found" representations with T = (polar terms) + (sum of length ≤ L), i.e. apparent savings up to the ideal (1−M)/2. **These are spurious for the true pairing.** They describe Σ_A γ̄₂-type(A)χ_A(m)ρ(m), whose coherent residue is a Gauss-sum series, not μ(A), and which detects no zeros.

---

## 4. Bounds in each faithful representation, as functions of M

Exponents are beyond Q^{5/12}, in Z-units; the CS target is M/2.

| rep. | coefficients | best provable bound | saving s₀(M) = M/2 − bound | in [0.80, 0.86] | status |
|---|---|---|---|---|---|
| R0 | bilinear, sextic kernel χ_A(m); ρ coherent (main term S per A) | de Faveri Θ₆(Z,Q): ½·max(1, M, 5/6 + M/3, 1/3 + 5M/6) | M/3 − 5/12 < 0 on (1/2, 1), where the term 5/6 + M/3 dominates | −0.141 at 0.826 | [R] given de Faveri Thm 1.1 |
| R1 | bilinear, quadratic kernel, π = ρζ, β = cubic theta | Goldmakher–Louvel / Heath-Brown: ½·max(M, 2M−1) = M/2 | **0 exactly** | 0 | [R]; beating it is (V′_LS) |
| R2 | coupled: ρ(s′μ³) = N(μ)^{−1/6}G₄(1,μ)χ̄_μ(s′)⁴ρ(s′μ) (KP) | trivial: (7M−4)/2 (K = 1); (9M−5)/2 (K = L) | (4−6M)/2 > 0 iff M < 2/3 (K = 1) | < 0 | [R] algebra; length [H] (conductor) |
| R3 | coupled ρ(bA), Möbius column | trivial on regular part: (2−M)/2 | (2M−2)/2 < 0; **and contains S** | < 0 | [R-]; circular |

**Notes.**
- **R1 is the only bilinear form whose large-sieve constant equals CS.** Generic π saturates it, so any saving must come from π's arithmetic. That is (V′_LS), with the same "beyond the edge" status as in A2 §6.
- **R2 at M < 2/3 is the dual route of A2 §2.6** (`lp_dual_route.py`: θ = 0.934 > Liu).
- **Is T ≈ S?** If R3's regular part were provably ≪ Q^{11/12}Z^{−δ}, (V′) would be *equivalent* to a Möbius-sum bound, i.e. a zero-free statement. Its trivial bound is Q^{−1/12}Z ≫ Q^{11/12}, so this is **not** established. The non-oscillating A-sum (μ(A)c_A(n) = μ(A)² for (A,n) = 1) has squarefree-counting errors of size Z^{1/2}; see `V_gaps.md` §1.6. [H]

---

## 5. Why: zero detection excludes shortening (design scan)

`design_scan.py` varies:
- the s-orientation: (M_ss, M_sm) ∈ {(4,5) paper, (2,1) conjugate};
- the column type: M_AA ∈ {5, 1, 0, 2, 4};
- the row twist: χ_A(m)^j, j ∈ {1,…,5}.

That is 50 variants. For each it records the groupoid, the Möbius indicator in the signal chamber, the residue-adapted types, and the best bound among representations that do **not** force the signal line.

**Viable probes.** A probe needs two things:
- (i) a cubic-theta column, M_AA ∈ {5, 1} (Patterson). Only then does the row admit the reflection behind ‖B‖² ≪ Q.
- (ii) a row twist that is coherent for the row theta. This is j = 1 for the paper's s-orientation (V_derivation §3.2 [R]) and j = 5 for the conjugate one. Only then does the coherent pole line P1 exist, so that a signal can arise at all.

That leaves four viable variants, a clean 2×2. **[RC]:**

| viable variant | alignment | signal column after the m-Voronoi | residue-adapted types | best bound avoiding P1 | short images avoiding P1 |
|---|---|---|---|---|---|
| **paper: s (4,5), column γ₂ (5), j = 1** | aligned (Davenport–Hasse) | **Möbius (3): detects zeros** | 4 | **0.413 = CS** | **none** |
| global conjugate: s (2,1), column (1), j = 5 | aligned | Möbius: detects zeros | 4 | 0.413 = CS | none |
| s (4,5), column γ₄-type (1), j = 1 (= U₁ ≅ A2's M0) | anti-aligned | Gauss type (5): no zero detection | 8 | 0.000 (both free lengths empty) | `Ams` (M−1, 2M−1); `AmsAmA` (M−1, −1) |
| s (2,1), column (5), j = 5 | anti-aligned | Gauss type: no zero detection | 8 | 0.000 | `Ams`, `AmsAmA` |

**Non-viable variants (the other 46).** For example:
- the cubic row twist (5, 2) has an image `Ams` with lengths (L²/Q, L) and a sextic kernel, with large-sieve exponent 0.351 < 0.413. This is ζ-free-like.
- (5, 4) and (5, 5) have empty-length images.
- a few variants have a Möbius node in a free variable: (0, 2), (0, 5), (2, 4), (4, 3) and their conjugates. These have no theta column, or an incoherent twist, so no pole lies on P1. The "P1 forced" flag then carries no signal, and these are not zero-detecting probes in the paper's sense.

**Reading.**
- [RC] Among viable probes, zero detection occurs exactly for the Davenport–Hasse-aligned pair. For that pair, every functional-equation image either keeps a variable of length ≥ Q with a large-sieve constant ≥ Q, or forces the signal line.
- [RC] The anti-aligned viable probes have representations that are entirely polar terms plus an empty or short sum, but they detect no zeros. The coherent residue there is Σ_A γ₄γ₁-type, a Jacobi-sum (Hecke) × Gauss series, not μ(A).
- [H] The mechanism. The same alignment γ₁γ₂ = μαG that makes the column Möbius after the m-Voronoi (the signal) makes the column reflection deposit γ₁γ₃ on the row (Lemma Z). The row node then has exponent 4, so the quadratic twist enters as the A₂ twisting parameter μ³ (R2, conductor L³). With the anti-aligned column the row stays bare, the kernel stays sextic, and further reflections shorten. The ζ-free model of `V_derivation.md` (dual length L²/Q, cubic dual) corresponds to such a non-zero-detecting configuration.
- So the failure of (V) → (V′) is not a missed chamber. **Within functional-equation geometry, the arithmetic that lets the probe see zeros is the arithmetic that blocks every shortening.**

---

## 6. Partial and smoothed applications

### 6.1 Functional equation on part of the μ-range (dyadic blocks) [R-]

Split R1 into blocks N(μ) ≍ L′ ≤ L. The column mass of a block is ≍ min(1, (L′/L)^{1/2}), because V♯(x) ≪ x^{1/4} (paper eq. `reflection-kernel`). Per block:
- CS: Q^{11/12}(L′/L)^{1/4};
- R2 on the block: Q^{−1/12}L′²(L′/L)^{1/4}K^{1/2}, which beats the block's CS only for L′ < Q^{1/2}.

Blocks with L′ ≤ L·Z^{−4s₀} are already below the target by the V♯ decay alone. The critical blocks are L′ ≍ L, where neither applies. Hence s₀ = 0.

Using the sharper V♯(x) ≪ x^{5/6−ε} (the first pole of R(t) is at t = −5/6) changes only the non-critical blocks.

### 6.2 Mixed reflected/unreflected splits [R-]

- **A-blocks.** Each block A ≍ Z′ of the column has mass ≍ 1.
  - Unreflected: the large sieve gives ≥ Q (sextic, de Faveri).
  - Reflected: the length is Q²/Z′ ≥ Q for Z′ ≤ Q.
  - V_G(x) = (4π)^{−1/2}e^{−(log x)²/4} suppresses Z′ ≠ Z^{1+o(1)} only super-polynomially in log.
  - So every block of non-negligible weight is at the CS scale.
- **Arithmetic splits.** For example A = A₁A₂, reflecting only A₁ for fixed A₂. The twist becomes χ_m χ_{A₂}^{(·)} with conductor ≥ q_m q_{A₂}, so the reflected length grows to ≥ Q²q_{A₂}/Z₁. That is worse.

### 6.3 Mellin-shifting a single variable

**(a) Separating the coupling weight** V♯(q_μZ/q_m²) by a Mellin integral on Re w = σ multiplies the large-sieve bound by (Q²/(ZL))^σ = 1, since L = Q²/Z exactly. The lengths are invariant, so there is no gain. [R]

**(b) Shifting the row variable to the centre** [R- given A2's residual functional equation and no poles for squarefree μ ≠ 1 (A2 §1.4, R-); conductor exponent 3 is H]. With Π(z; μ) = Σ_m π(m)(μ/m)₂N(m)^{−z},

  F′(μ) = (2πi)^{−1}∫_{(5/12)} Π(z;μ)Q^zΩ̂(z)dz,  |T| ≤ ‖β‖(Σ_μ|F′(μ)|²)^{1/2} ≪ Q^{5/12}(sup_t Σ_{μ≍L}|Π(5/12+it;μ)|²)^{1/2}.

So **(LoA_π(θ_π)): Σ_{N(μ)≍L}|Π(5/12+it;μ)|² ≪ L^{θ_π+ε}(1+|t|)^A ⇒ (V′) with s₀(M) = (M − (2M−1)θ_π)/2.** This is positive iff θ_π < M/(2M−1):
- 4/3 at M = 0.80;
- 1.267 at M = 0.826;
- 1.194 at M = 0.86.

**Benchmarks:**
- the approximate functional equation (length L^{3/2}) plus Heath-Brown's quadratic large sieve gives θ_π = 3/2, sharp for general coefficients;
- Lindelöf on average is θ_π = 1, which gives the ideal saving (1−M)/2.

**LP** (`lp_loa_threshold.py`; P²/Y piece dropped as in `lp_const_saving.py`; θ_true = min(Liu, LP) by distributivity):

| θ_π | θ (unmarked lengths) | θ (marked lengths L_J = Q²Z^{ℓ′−1} + Q^{5/6} floor; conservative) | optimum geometry (marked): l_x, l_y, ℓ, M |
|---|---|---|---|
| 1.00 (LoA) | 0.851918 | 0.854755 | 0.129, 0.765, 0, 0.894 |
| 1.10 | 0.859037 | 0.862599 | 0.176, 0.672, 0, 0.848 |
| 1.15 | 0.862100 | 0.865870 | 0.196, 0.633, 0, 0.828 |
| 1.20 | 0.865152 | 0.868797 | 0.213, 0.598, 0, 0.811 |
| 1.25 | 0.869736 | 0.871433 | 0.229, 0.567, 0, 0.796 |
| 1.30 | 0.873800 | 0.873870 | 0.244, 0.538, 0, 0.782 |
| 1.35 | 0.877460 → **0.874957** | 0.877528 → **0.874957** | (Liu's geometry) |
| 1.40 | 0.881243 → 0.874957 | 0.881316 → 0.874957 | |
| 1.50 (large sieve) | 0.889585 → 0.874957 | 0.889679 → 0.874957 | |

Checks: s₀ = 0 reproduces 0.874957 at Liu's geometry (0.355, 0.478, 0.167). θ_π = 1 unmarked reproduces `IDEA1_NOTES.md`'s "ideal large sieve, P²/Y dropped" value 0.851918.

So any power saving over the large sieve's 3/2 that reaches θ_π below about **1.31** (both variants; the marked optimum is slot-free) would beat Liu. This is a beyond-the-edge second moment: family size L, conductor L³. It is the central-value form of (V′_LS), not a new route.

### 6.4 Not taking the s-residue (finite Y) [R]

With ℓ_s = y finite, the chambers of U are the paper's own geometry:
- `id`: low side, CS in m;
- `s`: Kubota dual of length Q/Y; for Y > Q it is empty, which is the theta split itself;
- `A`: the reflected row;
- `m`: Poisson in m, length YZ/Q, Möbius column; this is the paper's high side.

No other adapted chamber exists (§3.1 table). The probe already uses every functional-equation image of the 3-variable object.

### 6.5 Opening the Gauss sum β_μ = γ₂(μ)·…

This introduces an additive variable x mod μ, giving Kloosterman-type sums after Poisson in μ. That is the Kuznetsov route (Idea 2 in `LIT_ideas234.md`), not functional-equation geometry. Not pursued here.

---

## 7. What remains (the minimal non-FE input)

Any of the following would give (V′) at M ≈ 0.826. None follows from functional equations plus generic large sieves.

1. **(V′_LS)**, from `TARGET_VLS.md`: the large-sieve constant of R1, beaten for the specific π.
2. **(LoA_π(θ_π))** with θ_π < M/(2M−1) (§6.3(b)). This is the same statement at central values.
   - In the LP the break-even exponent is θ_π ≈ 1.31.
   - Large sieve: 3/2. Lindelöf on average: 1.
3. **R2 with a conductor drop to k < 2M/(2M−1) − 1 ≈ 1.53,** together with K = 1. This needs a "degree-one" collapse of ρ(sμ³)γ₁(s) that KP does not force (A2 §2.3, referee item 3). The cubic-type degeneracy is excluded numerically (`A2_numerics.md` §0.2).

In the M-range [0.80, 0.86] items 1 and 2 are equivalent up to ε. All three are "beyond the edge" (family L, conductor ≥ L³ or Q²). This is the same type of input that blocks Ideas 2 and 3.

---

## 8. Proof-status ledger

| step | label |
|---|---|
| S-W/Heckenberger groupoid computations: chambers, objects, Möbius nodes, residue-adapted chambers, ℓ′ = E^{−T}ℓ | [RC] |
| U is the faithful 3-node model (reproduces paper Prop. completed-reflection and V_gaps §1.6) | [R] for the anchors; [RC] for the match |
| M0 (A2 §4.2) is the conjugate-theta encoding; M1 is faithful and of infinite type | [RC] |
| The list R0–R3 is complete for U | [RC] (exploration closes at Möbius nodes, where no FE exists [R]) |
| The list is complete for M1 | [H] (searched to 10 476 chambers) |
| Analytic realisation of each reflection | s, A, m(O0): [R]/[R-]; m(O2): [R-] (A2 §1.5(b)); composites: [R-]/[G] (3-variable continuation not written) |
| Effective conductor = Heckenberger exponent (lengths L³/Q in R2 etc.) | [H] (A2 referee item 1); conclusion robust unless k < 1.53 (§7.3) |
| Bounds in §4 | [R] given de Faveri Thm 1.1, Goldmakher–Louvel, KP (abs τ(mp³) = abs τ(mp)) |
| P1-forced criterion ⇒ representation contains S | [R-] |
| "Zero detection excludes shortening" (§5) | [RC] for the 50-variant scan; [H] as a general principle |
| (V′) ⟸ LoA_π(θ_π), θ_π < M/(2M−1) | [R-]: needs only continuation of Π to Re z = 5/12 with polynomial growth and no poles (A2 §1.4 [R-]); the large-sieve benchmark θ_π = 3/2 uses conductor N(μ)³ [H] |
| LP thresholds | [R as computation] given the LP model |

**Net.** Attack A is complete and negative.
- The systematic search that A2_theory left open (its §4.3 [G]) is done.
- In the exact model it returns only the four known images; the faithful 4-node object has infinite type and, over 10 476 chambers searched, also returns only these four.
- No image gives (V′); s₀ = 0 and θ = 0.874957.
- The obstruction is the Davenport–Hasse alignment that gives the probe its zero detection.
- The single remaining input is a beyond-the-edge mean value: (V′_LS), or equivalently LoA_π with θ_π < M/(2M−1).

---

## Appendix: reproduction

All commands run from `gram_split/vls/` with `PY=/private/tmp/claude-501/.../scratchpad/venv/bin/python`.

| command | reproduces | approx. time |
|---|---|---|
| `$PY model_check.py` | §2.2 anchors and M0/M1 `At` | seconds |
| `$PY orientation_test.py` | 3-node scan e_tA ∈ {2, 4}, j = 0…5 | seconds |
| `$PY ../wf2/weyl_groupoid.py`, and `run(M1, 6, cap, maxcoef)` | finiteness | seconds |
| `$PY chamber_analysis.py U M0` | §3.1, §3.3 | seconds |
| `$PY m1_deep.py` | §3.2 | ~minutes |
| `$PY u_full_table.py` | full U table | seconds |
| `$PY design_scan.py` | §5 | ~1 min |
| `$PY lp_loa_threshold.py loa` | §6.3, unmarked | ~10 min |
| `$PY lp_loa_threshold.py loamarked` | §6.3, marked | ~10 min |

`groupoid_search.py` (DFS) prints all residue-adapted chambers of M0 with their basis matrices E. `chamber_analysis.analyse(model, cap, maxcoef)` returns, for every residue-adapted chamber type:
- word;
- free-variable maps;
- lengths;
- κ;
- node types;
- edge;
- residue numerators;
- bounds.

---

## Referee report

*Independent skeptical review, 9 Oct 2026. I read this file, `TARGET_VLS.md`, `wf2/A2_theory.md` (including its referee report), `wf/V_gaps.md` §1, `V_derivation.md` §3 (coherence), paper78 Prop. `completed-reflection` (l.1677–1860), and every script listed in the Appendix. I re-ran `model_check.py`, `orientation_test.py`, `u_full_table.py`, `chamber_analysis.py U M0`, `design_scan.py`, and `weyl_groupoid.run` on M0, M1 and their global conjugates (M1 up to maxcoef 1000). I read `m1_deep.log` and `lp_loa_marked.log`. New checks, saved in `vls/`:*
- *`referee_A_orient_all_j.py`: the orientation audit for every twist exponent j;*
- *`referee_A_lp_loa_sub.py` and its `.log`: the unmarked LoA LP at θ_π = 1.0, 1.30, 1.32, 1.35.*

### Verdict: not refuted

The file claims no saving, and that negative result stands:
- **In the exact 3-node model U, the functional-equation geometry gives s₀ = 0 at M ∈ (2/3, 1). So θ stays at Liu's 0.874957.**
- The A2_theory §4.2 orientation correction is right. My all-j check makes it considerably stronger than the file's j = 1 test.

Several side claims are overstated (items 1–5). None of them changes the verdict.

### What I verified

| claim | check | result |
|---|---|---|
| U: 24 chambers, 4 objects (2 with a Möbius node), 8 residue-adapted chambers, 4 types | `u_full_table.py`, `chamber_analysis.py U` | reproduced |
| Lengths: R1 (M, 2M−1), R2 (5M−3, 2M−1), R3 (1−M, 1), R0 (M, 1) | E^{−T}ℓ by hand for a simple reflection (u_i = −u′_i, u_j = u′_j + m_ij u′_i ⇒ ℓ′_i = Σ m_ij ℓ_j − ℓ_i); E for A, Ams, ms read off the table | correct. They agree with the paper (L = Q²/Z), A2 §2.6 (L³/Q = Z^{5M−3}) and V_gaps §1.6 (n = Ab, length Z²/Q) |
| Coefficient structure | s-residue numerators: R1 {m:1, μ:0} (bilinear; row node 4 = γ₁γ₃, so **Lemma Z is accounted for**); R2 {s′:1, μ:3} → ρ(s′μ³); R3 {b:1, A:1} → ρ(bA) | correct; R2 and R3 are genuinely coupled (KP leaves ρ(s′μ) undetermined) |
| Bound exponents at M = 0.826 | R0 ½·(1/3 + 5/6) = 0.554; R1 M/2 = 0.413; R2 (7M−4)/2 = 0.891; R3 (2−M)/2 = 0.587 | arithmetic correct. R0 is ≥ CS for the simpler reason that a column of length Z > Q forces a large-sieve constant ≥ Z |
| M0 finite (336 chambers, 7 objects); M1 infinite | `weyl_groupoid.run` | M1 coefficients keep growing (86 / 139 / 192 bases at maxcoef 40 / 200 / 1000). Global conjugates behave identically |
| M0's 8 types, including (2M−1, M−1) and (−1, M−1) off the signal line | `chamber_analysis.py M0` | reproduced; they coincide with the anti-aligned U₁ = design-scan row (paper, M_AA = 1, j = 1) |
| Design-scan 2×2 | `design_scan.py` | reproduced row by row |
| LP | `lp_loa_marked.log`; my unmarked rerun | unmarked θ_π = 1.00 → 0.851918, 1.30 → 0.873800 and 1.35 → 0.877460 (as in the file). New: 1.32 → 0.875256 > Liu, so break-even lies in (1.30, 1.32), consistent with "≈ 1.31" |

### The load-bearing step: orientation (M0 vs M1) — confirmed and strengthened [RC within the formalism]

The negative verdict depends on M0 *not* being the faithful 4-node object. M0 has images with lengths (L, ∅) and (∅, ∅) that avoid the signal line. Had they been images of T, they would give s₀ = (1−M)/2 ≈ 0.087 and θ ≈ 0.852–0.855. So I tested the orientation against the paper's full reflection rule for every j, not only j = 1. Prop. `completed-reflection` gives, for j ∉ {0, 4}:
- dual column factor B_p = χ_p^{−j−2}, i.e. S-W edge −j−2;
- row factor τ_j τ_{j+2}, i.e. S-W exponent 2j − 4.

`referee_A_orient_all_j.py` takes 3-node (t, A, m) encodings and looks for the chamber with e_t simple and free lengths (2ℓ_m − ℓ_A, ℓ_m):

| t-encoding (M_tt, e_tA) | j = 1 | j = 2 | j = 3 | j = 5 |
|---|---|---|---|---|
| paper predicts (edge, row) | (3, 4) | (2, 0) | (1, 2) | (5, 0) |
| **M1-type (1, 2)** | **(3, 4)** | **(2, 0)** | **(1, 2)** | **(5, 0)** |
| M0-type (5, 4) | (1, 0) | none | (5, 4) | (3, 2) |
| direct U-type (A Kubota, M_AA = 5) | (3, 4) | (2, 0) | (1, 2) | (5, 0) |
| mixed encodings (5, 2), (1, 4) | none | none | none | none |

So:
- M1 reproduces the paper's rule at four independent exponents.
- M0 reproduces it at none. Its j = 5 output (3, 2) is the global conjugate of the paper's j = 1, exactly as the file says ("paper's configuration with conjugated column").

At the level of the S-W/Heckenberger formalism, A2_theory §4.2's "honest rank-4 arithmetic root system, 336 chambers, 7 objects" therefore describes the anti-aligned pairing, and its empty-length images are not images of T.
- The tension with Bröker–Hoffstein's "same symbol" is a convention question, correctly labelled [H].
- It plausibly reflects the classical conjugation between Kubota-series residues and Gauss sums in Patterson's evaluation. I did not check this against the sources.
- It does not affect the audit, which uses only proven anchors.

### Issues (none changes the verdict)

1. **"Past a Möbius node there is no functional equation" [R] is overstated (§3.1).**
   - 1/L(u,χ) *does* satisfy a formal functional equation.
   - What fails is that the contour cannot be moved across the critical strip without picking up the zeros of L.
   - Restate it as "no reflection at a Möbius node is available without zero-dependent residues". That is [R], and it gives the same conclusion: any image beyond such a node contains signal-type, zero-sensitive terms.
2. **The θ_π = 3/2 benchmark is not established (§0 item 5, §6.3(b)).**
   - At Re z = 5/12, the dual half of the approximate functional equation has μ-coupled coefficients ρ(sμ³)γ₁(s). Heath-Brown's quadratic large sieve therefore applies only to the direct half.
   - The only *provable* LoA_π exponent is the one coming from the length-Q large sieve, θ_π = M/(2M−1). That is exactly Cauchy–Schwarz, which is consistent with the rest of the file.
   - Relabel "large sieve gives 3/2" as [H].
3. **"(V′_LS) and LoA_π are equivalent up to ε" (§7) is one-directional.**
   - LoA_π(θ_π) ⇒ (V′_LS) at the same M is shown.
   - The converse (from length-Q twisted mean squares at all scales to central values) needs an approximate functional equation whose dual half is again coupled: [H].
4. **The reduction (V′) ⟸ LoA_π is labelled [R-]; it should be [R-]/[G].**
   - It needs Π(z;μ) continued to Re z = 5/12, with no poles and polynomial growth, *for the true π*.
   - The A2 referee's item 2 (coprimality at p | μ needs a vector-valued functional equation: [G]) and item 4 (the ᾱ² angle twist lies outside every cited framework) apply here.
   - Π converges absolutely only for Re z > 11/12, so the A₂ continuation is genuinely needed.
5. **R3's "trivial" bound Q^{−1/12}Z (§4) is the random-size normalisation, not a proved bound.**
   - V_gaps §1.6's actual unconditional bound for the regular part is Q^{−7/12}Z^{3/2+ε}, which is worse.
   - No effect: both exceed CS.
6. **Design scan (§5) is [H] as a principle, as labelled.**
   - Two variants have a Möbius node in a P1-forced chamber *and* short images off the signal line, with best exponent 0.000: (paper, M_AA = 0, j = 2) and (conj, M_AA = 0, j = 4).
   - They are excluded because the twist is not coherent, so P1 carries no pole, and the column has no theta, so ‖B‖² ≪ Q is lost. That exclusion is reasonable.
   - But "zero detection" is operationalised as a formal criterion: a Möbius node in the P1-forced chamber, plus a coherent pole.
   - The design space is also narrow: 3 nodes, a bare row, two s-orientations. "Mutually exclusive" holds for this 50-element family only.
7. **M1 completeness [H] is correctly labelled.**
   - In infinite type, free lengths a + bM with large |a|, |b| could in principle nearly cancel (e.g. 6M − 5 < 0 at M = 0.826). A search to maxcoef 80 cannot exclude this.
   - Supporting heuristic, not in the file: chambers with e_t simple should be exactly the chambers of the t-residual object, which is U, and U is complete. Indeed M1's four types coincide with U's four.
8. **Cosmetic.**
   - In `m1_deep.log` type R0 (word `t`) has multiplicity 3, not 4.
   - §6.2: V_G(Z′/Z) at Z′ = Z^{1+δ} is Z^{−δ² log Z/4}, which decays faster than any power. The conclusion (only Z′ = Z^{1+o(1)} matters) is right; the wording "only super-polynomially in log" is not.
9. **Number-field transfer.**
   - Every [RC] statement is a statement about the S-W/Heckenberger model, which is a theorem over F_q[T] only.
   - For K = Q(ω), n = 6, μ₁₂ ⊄ K, so extra ray signs (−1/ϖ)₆ appear (A2 §1.5). They do not alter the Heckenberger exponents, so the negative verdict is not affected.
   - "No image beats CS" is [RC] in the model and [H] as a statement about T. The effective-conductor caveat (§1.2) is the main residual risk; it was already quantified by the A2 referee (a drop would have to reach k < 1.53 with K = 1).

Not applicable:
- No probe modification is claimed to work, so there is no high side to re-check. The shortening (anti-aligned) variants detect no zeros.
- The file contains no numerics beyond the LP, so estimator bias does not arise.

### Recommendation

Accept "s₀ = 0 from functional-equation geometry", with these relabels:
- item 1: reword the Möbius-node statement;
- item 2: θ_π = 3/2 → [H];
- item 3: "equivalent" → "⇐";
- item 4: LoA reduction → [R-]/[G].

Record the all-j orientation table as the evidence for the A2_theory §4.2 correction.

The open input remains (V′_LS), or sufficiently LoA_π with θ_π < M/(2M−1). Its LP break-even is θ_π ∈ (1.30, 1.32), reached at M ≈ 0.78, below the [0.80, 0.86] window. The only provable value is θ_π = M/(2M−1), which is exactly Cauchy–Schwarz.
