# mark-design: can the choice of slot primes, slots and weights control the energy beyond the wall?

Agent `mark-design`, 10 Oct 2026.

**Verdict: NEGATIVE for all five ideas.** No choice of prime set, slot layout, weights, two-sided placement or labels
gives a provable W(w) with w > 0. θ stays at 0.874957. Each idea ends in a proved structural lemma (§1) or a bookkeeping
no-go. The numerics show that G − D has no low-rank or block structure to exploit, even in truth (§7). The exact
requirement (§8): the marks' extra conductor per unit length would have to drop from P¹ to P^{0.687}
(dual conductor P^{1.69} instead of P²). Lemma L2 shows that the amplifying cubic symbol forces P².

Labels:
- **[P]** proved, given the P78 identities;
- **[B]** exponent bookkeeping;
- **[C]** conditional;
- **[H]** heuristic;
- **[N]** numerical.

## 0. Objects and the Gram dictionary

The binding dyad is Lemma `reflected-energy` / `hybrid-energy` (P78 l.7510, l.7747), with K = Z^M, N = K²L/Z,
L = Z^ℓ and wall excess w = M + ℓ − 1 = log_Z(N/K).

    S_{k,P} = Σ_n β(n) ψ_k(n) χ_P(n)^{-2},   G = S*S  (L'×L', PSD),   E(c) = c*Gc,   signal = ⟨c, loc⟩.

**Exact dual phase of the coherent vector [P].** Combine eq. `marked-reflection-phase` (l.7396–7405) with the
definitions at l.1818–1830 (σ_q = λ²c/q, ε_q = −((λ³c/q)σ_q)^{-1}, ω_{q,0} = −τ⁺_{q,2}χ_q(ε_q)^{-2}). The weight of a
marked tuple in the dual basis χ_P^{-2} is

    ã(P₁⋯P_r) ∝ η̄(P₁⋯P_r) · ᾱ(P₁⋯P_r)² · Π_i τ⁺_{P_i,2} · χ_{P_i}(λ³c_F)² · Π_{i≠j} χ_{P_i}(P_j)²,

and Π_{i≠j}χ_{P_i}(P_j)² = Π_{i<j}(P_i/P_j)₃² = Π_{i<j} conj (P_i/P_j)₃ by cubic reciprocity. **The multi-slot dual weight carries the cubic
theta cocycle.** Up to unimodular one-prime factors it is β(P₁⋯P_r) (twisted multiplicativity gives the same
conj (P_i/P_j)₃), not a tensor product. Numerics use the single-slot
case, c_P = ᾱ(P)²γ₂(P) with η = 1 (label `p78`).

**What is proved about G:**
- (F1) G ≼ (K+N)(N + L + (NL)^{2/3})·I. This is `hybrid-energy`: the quadratic large sieve in k, then the cubic large sieve in P.
- (F2) Each G_PP ≤ (K+N)N, and no better: the single-mark energy already has the wall (NOTES §9). Hence only tr G ≤ L'(K+N)N.
- (F3) [C] tr G ≪ L'(K + N/L)N·Z^ε, *if* the optimal large sieve **(HLS)** holds for the thin hybrid family
  {ψ_kχ_P^{-2} : k ≤ K, P ~ L} (KL characters; N ≤ KL, so it sits at or below its own edge).
  - (HLS) is coefficient-blind. It is exactly what energy-wall's O2 "break-even" silently assumed.
  - I know no proof. It does not follow by positivity from the cubic or sextic large sieves used in this project: their
    (QN)^{2/3}-type terms exceed the family size Q = KL at the 27/31 geometry, by about Z^{0.26} if taken in that form.

## 1. Five structural lemmas (used throughout)

**L1 (Cauchy–Schwarz separation) [P].** The low side is ‖A‖·‖B‖, so the wall is a property of ‖B‖ alone.
- Anything done to A leaves the wall unchanged.
- So does any unimodular factor in m, and any cancellation inside the product A_m B_m.
- Only the theta-side marks (prime set, slots, weights) enter ‖B‖.

**L2 (amplification–conductor coupling) [P].** Three identities combine at a marked prime:
- twisted multiplicativity gives β(PA′) = β(P)β(A′)χ_P(A′)⁴;
- the Poisson step in m produces g(χ_{PA′},u) = χ_P(A′)χ_{A′}(P)g(χ_P,u)g(χ_{A′},u) = 𝔯χ_P(A′)²·g g;
- χ_P(A′)⁴·χ_P(A′)² = 1.

Consequences:
- The coherent principal-row factor −η(P)𝔥(P)p^{z−1} therefore **requires the m-sum to carry χ_P(m)**. Here 𝔥 is the fixed Hecke factor of γ₁γ₂ = μαG.
- On the low side, the same χ_P⁴(A′) makes B a cubic twist of θ, with reflected length Q²L/Z. That is the wall.
- So the amplification and the conductor are one symbol. This sharpens NOTES §20 to the level of the identity.

**L3 (no coherent finite labels) [P, P78 probe class].** By `lem:local-euler`, the principal row is an Euler product that
treats all primes of A alike.
- A character Ψ inserted in A, directly or produced by the CRT from an m-twist ξ mod q₀ (g(ξχ̄_sχ_A,u) ∋ ξ(sA)χ̄_s(q₀)χ_A(q₀)), turns the target η into ηΨ. On the s-side, ζ(w) becomes a twisted L(w,·) with no pole.
- A label λ(P) on the weights is compensated only if Ψ = λ̄ on the slot primes and Ψ = 1 on every other prime. For a Hecke character this forces Ψ = 1.
- For the m-twist the only target-preserving choice is ξ = χ̄_{q₀}. Its label is ξ(P)χ_P(q₀) = 𝔯 (a fixed sign) by sextic reciprocity.
- **The only coherent labels are archimedean**, N(P)^{it}, because L(s, η|·|^{it}) has the same zeros. That is Mellin-height averaging (NOTES §9, dead).

**L4 (cubic clique bound) [P, via Heath-Brown's cubic large sieve].** Suppose a fraction ≥ 1/3+δ of ordered pairs in a
set 𝒜 of primes of norm ~L has (P₁/P₂)₃ = 1.
- Then |𝒜| ≪ L^{2/3+ε}/δ, since |Σ_{𝒜×𝒜}(P₁/P₂)₃| ≪ |𝒜|·L^{2/3+ε}.
- Across two slots the bound is A₁A₂ ≪ (L₁+L₂+(L₁L₂)^{2/3})^{1+ε}/δ.
- A slot with A = L^α primes has net amplification L^{1/6}(A/L)^{1/2} = L^{(α−2/3)/2}. So **cubic cliques amplify nothing.**

**L5 (optimal weights; the adversarial Gram matrix) [P].**
- The best signal²/energy over all designs is loc*G^{-1}loc, attained at c* = G^{-1}loc.
- Any inequality G ≼ KN·Z^ε·I + R, with R of explicit range W, gives loc*G^{-1}loc ≥ ‖Proj_{W^⊥}loc‖²/(KN·Z^ε).
- Every P78 lemma holds uniformly for |a(P)| ≤ 1. So G_adv := KN·I + (N²−KN)ûû* with û = loc/‖loc‖ satisfies (F1), (F2) and (F3).
- For G_adv, every design with |⟨c,loc⟩| ≥ α‖c‖‖loc‖ has E(c) ≥ α²(N²−KN)‖c‖².
- **Hence no design gains without a new, phase-sensitive inequality about G at the coherent vector.**

## 2. Idea 1: restricted prime sets

| restriction | result |
|---|---|
| ray classes, splitting in a *fixed* abelian extension, norms from a fixed quadratic extension | Each is a fixed ray-character restriction, which P78's lemmas allow. G is replaced by a principal submatrix; only constants change [P/B]. Numerically there is no block structure: normalised \|G_{P₁P₂}\| within/across classes of (π₀/P)₃ is 0.01586/0.01623 for θ and 0.01588/0.01612 for random-phase columns [N]. |
| growing auxiliary conductor q₀ = L^γ | The density L^{-γ} loses amplification. Model (`dens`, optimistic, high side unchanged): dens 0.9 → 0.88200, 0.8 → 0.88961, 2/3 → 0.89547 [B]. The character components of 1_class are labels, so L3 applies. |
| mutual cubic symbols (P₁/P₂)₃ = 1 within or across slots | Impossible beyond L^{2/3} (L4) [P]. The best greedy clique is 7 of 184 primes (≈ L′^{0.37}). Entries with symbol 1 versus ≠ 1 give 0.01632 vs 0.01600, the same as random columns (0.01622 vs 0.01595) [N]. |
| "χ_{P₁}χ̄_{P₂} degenerate on the dual support" | Impossible. The dual support is all squarefree n ~ N with N ≫ L², which is equidistributed mod P₁P₂. The n-diagonal of G is already mark-orthogonal. The column term comes from the k-off-diagonal, which no restriction on P touches [P]. |
| exact block structure (blocks of A primes) | Useless for A > 1. The wall depends on the size L of the primes, not on their number, so each block is beyond its own wall. Only A = 1 (full cross-orthogonality) helps, and even then one still needs Σ\|c_P\|²G_PP ≪ KN, i.e. (H1) on average ⟸ (HLS) [B/C]. |
| conjugate pairs {P, P̄} | Normalised \|G_{P,P̄}\| is 0.0219 vs 0.0161 overall. Random-phase columns show the same (0.0215 vs 0.0160), so this is a family effect, not θ: slightly worse [N]. |

## 3. Idea 2: structured weights (c ⟂ bad subspace)

- **The optimum is c* = G^{-1}loc (L5).**
  - It helps only if the column term is carried by an *explicitly identified* subspace W with loc ∉ W.
  - The proof of the column term is isotropic. In the column space the extremal set is all of V = span{ψ_k : k ≤ K} (dimension K ≫ L′). Pulled back to c-space, it singles out no subspace.
- **Explicit rank-one candidates:**
  - (i) **Single rows s_{k₀}**, the quadratic large sieve's single-character extremals. Removing them, or re-phasing rows, leaves the large-sieve bound unchanged [P]. The generic adversary has every row at |S_k|² ≈ N²/K ≤ L′N, i.e. ordinary row size, so there is nothing to remove [B]. Numerically row norms are max/mean = 1.33, and the best single-row overlap with the top eigenvector is 0.053 (controls 0.050, 0.049) [N].
  - (ii) **The Patterson-bias direction** (wall-spectral §4.3) is explicit and rank one on the principal rows. It can be nulled by design, e.g. with two ray classes of opposite Kubota residue [H]. That only lifts the 19/22 cap for large w; it is irrelevant for 27/31.
- **Rank of the bad part [C].** Under (F3) the eigenspace of eigenvalues ≥ KN·Z^{w′} has dimension ≤ L′Z^{-w′}. It is relatively low-rank but unidentified, and G_adv puts loc inside it.
- **In truth there is nothing to orthogonalise against** [N, §7]. The optimal weights gain only a factor 1.04–1.09 over the coherent choice.

## 4. Idea 3: multi-scale slots

- **[P] Slots do not factor, already at the level of weights.** The dual weight of a multi-slot mark carries Π_{i<j}(P_i/P_j)₃ (§0). Making it a tensor product needs cross-slot cliques, which L4 rules out.
- **[P] Kubota-index form.** B^{(𝒫)}_m ∝ χ_𝒫(m)·Θ(m𝒫⁴), where 𝒫 = ΠP_i and Θ is the unmarked row. So slots are index shifts of one theta row, not independent factors. In addition ψ_k(n₁n₂) and χ_P(n₁n₂) couple all the prime factors of n.
- **[B] The wall sees only the product.** The reflected length is N = K²(ΠP_i)/Z, so the wall is M + Σℓ_i ≤ 1. Each slot being below its own wall (M + ℓ_i ≤ 1) is irrelevant.
- **[B] Exchange rate.** Moving a slot into the rows by Cauchy–Schwarz breaks even only under (HLS). With the proved per-P quadratic sieve it is a strict loss of v.
  - A "better than break-even" mechanism must keep part of the slot sum's coherence inside the square, i.e. bound c*Gc at a coherent c. That is circular.
  - Hölder and slot moments are phase-blind, so O1 applies.

## 5. Idea 4: marks on the Gauss-sum side and the theta side at once

- **Same prime P on both sides (P | s and P | A).**
  - χ̄_P(m)χ_P(m) = 1 removes P from the Poisson modulus. So g(χ_P,u) never appears, and the cross symbols survive: χ_P(s′)²𝔯 from the Gauss-sum CRT and χ_P(A′)⁴ from the theta twisted multiplicativity [P].
  - The principal row therefore becomes L(w, χ_P²·)/L(x, ηχ_P⁴·): no w = 1 pole, and a twisted target.
  - The coherent term moves to rows u ∈ P⁴(F^×)⁶, a different row for each P (cf. new-architectures, the p | m case) [B].
  - Meanwhile ‖B‖ is literally P78's marked B-norm, because χ_P(m) is unimodular. So **the conductor effect does not cancel at all** (L1).
- **Distinct primes on the two sides.**
  - A-side marks are the reparametrisation (X,Y) → (X/P_A, YP_A), since Gauss sums are twisted-multiplicative, with Gram exponent (b+2ℓ_A)/6 [P/B] (new-architectures §5).
  - The combination is P78 at another geometry, which the LP already optimises.
- **No gain.**

## 6. Idea 5: exact orthogonality by averaging over labels

- **Orthogonality itself is trivial [P].** If the label group separates the slot primes (N(q₀) ≳ L), then avg_λ λ(P₁)λ̄(P₂) = δ_{P₁P₂}, so E_avg = Σ|c_P|²G_PP.
  - This still needs (H1) on average ⟸ (HLS) [C]. So exact cross-orthogonality alone is **not** W.
- **Coherence fails (L3) [P].** Each realisation fails as follows:
  - a theta-row twist changes the target;
  - an m-twist ξ produces the label ξ(P)χ_P(q₀), but also twists all of A and s, and is target-preserving only when the label is trivial;
  - a Gauss-side label (P | s) is the double mark of §5;
  - archimedean labels are the Mellin-height averaging of NOTES §9.
- **Averaging m-twists with Gauss-sum weights.** Σ_ξ conj g(ξ)ξ(m)ξ̄(P) = φ(q₀)e(−mP̄/q₀) [P]. So the label average only re-phases the family variable m, and hence the rows after reflection [B]. The quadratic large sieve is invariant under row phases [P], so nothing provable is gained.
- **Bookkeeping, even ignoring coherence [B].** A label twist of conductor q₀ ≥ L on the theta row gives the wall M + 2ℓ ≤ 1 (`wall_l = 2`): θ = 0.88169, which is worse. This matches wall-spectral §2.3.

## 7. Numerics: structure of G (`gram_design.py`)

Setup:
- exact Z[ω] arithmetic: β = γ₂ᾱ, quadratic rows ψ_k, marks χ_P^{-2};
- all windows are beyond the wall (reflection proves diagonality only for K ≤ N₀/L ≈ 13–27);
- controls are random-phase columns (same |βw|) and an iid complex Gaussian S of the same shape;
- Ĝ = D^{-1/2}GD^{-1/2}, and R(c) = c*Gc/(‖c‖²·mean D).

Results for N₀ = 4·10⁴ (11311 columns), L′ = 184 marks of norm in [1500, 3000):

| K/N₀ (rows) | spectrum of Ĝ: θ / random columns / Gaussian [MP edges] | top-1 share of ‖Ĝ−I‖²_F | stable rank (of 184) | R(p78) | R(random c): mean ± sd (max of 200) | best flat c / sup | optimal-weights gain |
|---|---|---|---|---|---|---|---|
| 1/4 (2831) | [0.51,1.70] / [0.49,1.72] / [0.57,1.57] [0.56,1.58] | 0.027 / 0.029 / 0.027 | 36.5 / 34.6 / 36.8 | 1.039 | 1.000 ± 0.023 (1.056) | 1.56 / 1.70 | 1.087 |
| 1 (11279) | [0.62,1.49] / [0.62,1.48] / [0.76,1.26] [0.76,1.27] | 0.027 / 0.026 / 0.023 | 37.0 / 38.4 / 44.2 | 1.019 | 1.001 ± 0.016 (1.057) | 1.40 / 1.49 | 1.040 |

Notes on the table:
- N₀ = 2·10⁴ with L′ = 60 agrees: at K = N₀, R(p78) = 1.010, sd 0.020, sup 1.41 and gain 1.033.
- The coherent vector's overlap with the top eigenvector is 0.001–0.012, i.e. about 1/L′.
- The θ-Gram is indistinguishable from the random-phase-column Gram. Both are a little wider than Marchenko–Pastur, which is a property of the {ψ_k, χ_P} family and not of θ.
- G − D is **not** low rank: the top eigenvalue carries ≈ 2.7% of the Frobenius mass, and the stable rank is ≈ L′/5.
- **Restricted designs** are in the §2 table: there is no block structure.
- For K′ < L′ (K = N₀/64), G is singular and the optimal-weights gain is not defined (recorded as nan).

Files: `gram_N40000_L1500.{out,json}`, `gram_N20000_L400.{out,json}`.

## 8. Model payoffs (`model_md.py` = wall-spectral `model_ws.py` + knobs `wall_l`, `dens`; `run_scan.py`, `runs/`)

| hypothesis | θ |
|---|---|
| baseline (wall_l = 1) | 0.8749568 |
| W(0.064) (check) | 0.8709554 |
| marks' extra conductor per unit length wall_l = 0.9 / 0.8 / 0.75 / 0.7 | 0.873849 / 0.872595 / 0.871904 / 0.871164 |
| **wall_l = 0.69 / 0.685 / 0.68** | 0.871009 / **0.870931** / 0.870853 (27/31 at **wall_l ≈ 0.687**) |
| wall_l = 0.6 / 0.5 / 0 | 0.869514 / 0.867593 / 0.861111 |
| labels (wall_l = 2) | 0.881693 |
| restricted density dens = 0.9 / 0.8 / 2/3 | 0.882 / 0.890 / 0.895 |

So 27/31 needs marks whose dual twist costs at most P^{1.687} per length P, instead of P². In the P78 mark class
this is impossible (L2), because the P^{1/6} gain and the P² conductor come from one cubic symbol.

## 9. Bottom line

| idea | status | obstruction |
|---|---|---|
| 1 restricted prime sets | dead | fixed data is neutral [P]; growing data are labels (L3) or lose density [B]; cliques ≤ L^{2/3} (L4) [P]; no block structure [N]; blocks of size > 1 are useless [B] |
| 2 structured weights | dead | c* = G^{-1}loc needs an explicit bad subspace; proved bounds are phase-uniform, so G_adv is consistent with them (L5) [P]; truth has none (gain 1.04) [N] |
| 3 multi-scale slots | dead | dual weights carry the theta cocycle [P]; the wall sees ΠP_i [B]; exchange rate ≥ break-even, and break-even only under (HLS) [B/C] |
| 4 both sides | dead | Cauchy–Schwarz separation (L1); a double mark kills absorption and keeps the B-norm (L2) [P]; distinct primes are a reparametrisation [B] |
| 5 labels | dead | only archimedean labels are coherent (L3) [P]; label orthogonality leaves (H1)-average ⟸ (HLS) [C]; labelled walls are worse (wall_l = 2) [B] |

Two by-products for the coordinator:
1. **(H1) can be weakened.** In wall-spectral's conditional theorem, (H1) can be replaced by its average over P. That average is implied by (HLS), a *coefficient-blind* large sieve for the thin hybrid family {ψ_kχ_P^{-2}} at its own edge. So the only theta-specific input left is (H2).
2. **O2 assumed (HLS).** Energy-wall O2's "break-even" for slot transfers assumes (HLS). With proved tools the transfer is a strict loss.

A design can only help through a new phase-sensitive inequality at the coherent vector, i.e. (H2) with saving ≥ L^{-0.93}, or
through an amplification source other than the cubic symbol.
