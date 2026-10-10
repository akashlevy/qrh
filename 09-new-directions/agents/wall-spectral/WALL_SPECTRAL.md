# wall-spectral: automorphic and spectral routes to W(w)

Agent `wall-spectral`, 10 Oct 2026.

**Verdict.** None of the four routes gives W(w) for any w > 0. Unconditionally w = 0, so θ stays at 0.874957.

For each route I give either a proved no-go, stated as exponent bookkeeping, or the precise missing input. One correction to the brief:
- Quadratic twists of the cubic theta are *not* Grössencharacter twists of sextic Eisenstein coefficients.
- They are **Möbius** twists of them (identity (∗) in §1.1, proved and checked numerically).
- So the "Eisenstein reformulation" is P78's own Poisson/Davenport–Hasse duality (NOTES §6), not a new object.

Labels:
- **[P]** proved, given P78 lemmas;
- **[B]** proved as exponent bookkeeping;
- **[C]** conditional;
- **[N]** numerical;
- **[H]** heuristic;
- **[L]** literature, from abstracts only (arXiv is blocked).

Files in this directory:

| file | content |
|---|---|
| `identity_check.py`, `.out` | identity (∗) in §1.1, 540 cases |
| `fe_group.py`, `fe_group_compare.py`, `.out` | FE groups: ours (affine), SY/Li (affine), cubic family (finite) |
| `mds_reach.py`, `.out` | invariant I, the polar orbit, contour bookkeeping |
| `opnorm.py`, `opnorm_*.out/json` | operator norm of the (k,P) matrix, i.e. W uniformly in the marks |
| `model_ws.py`, `ms_*.json`, `mh_*.json` | copy of energy-wall's `model_wall.py` with knobs `bias_relax` and `Hoff` |

## 0. Normalisation

- Write a(n) = τ(n)X(n)V(q_n/N), with |τ| ≍ 1 on squarefree n, so ‖a‖² ≍ N.
- X(n) = Σ_P c_P χ_P(n)^{−2}, with c_P = ã(P)q_P^{−1/2}.
- E = Σ_{k ≤ K} |S_k|², where S_k = Σ_P c_P S_{k,P} and S_{k,P} = Σ_n τ(n)V ψ_k(n)χ_P(n)^{−2}.
- Target: E ≪ K N Z^ε.
- Scales: K = Z^M, N = K²L/Z = K Z^w, L = Z^ℓ.
- Each twist θ⊗ψ_kχ_P^{−2} has conductor (q_k q_P)², so its approximate-functional-equation (AFE) length is KL. Since N/(KL) = Z^{M−1} < 1, **the sums S_{k,P} are shorter than the AFE length.**
- The k-family (K members) is therefore beyond the edge by the factor L. The (k,P)-family (KL members) sits exactly at the edge, but only if the P-sum is incoherent, and in E it is coherent.

## 1. Route 1: Eisenstein reformulation and Rankin–Selberg

### 1.1 The identity [P; numerically checked]

Let β(n) = γ₂(n)ᾱ(n) be P78's theta coefficient on squarefree n, with no μ. Let γ₁(n;h) = q_n^{−1/2} Σ_x χ_n(x) e(hx/n) be the sextic Gauss sum. This is the numerator of the h-th Fourier coefficient of the sextic Kubota Eisenstein series.

For squarefree n with (n, 6kP) = 1:

    (∗)  β(n) ψ_k(n) χ_P(n)^{−2} = μ(n) · G(n) 𝔯(n,k) · conj γ₁(n; k³P⁴),

where:
- G(n) = conj χ_n(4) Γ(n), a function on a fixed ray class group mod 12;
- 𝔯 is the sign bicharacter mod 4.

Proof sketch:
1. P78 eq. `signal` gives γ₁γ₂ = μαG.
2. Hence β = μ G conj γ₁. **The Grössencharacter α cancels completely.**
3. Then use γ₁(n;h) = conj χ_n(h) γ₁(n).
4. Sextic reciprocity gives ψ_k(n) = 𝔯(n,k) χ_n(k)³.
5. Cubic reciprocity gives χ_P(n)^{−2} = χ_n(P)⁴.

`identity_check.py` confirms (∗) in 540 cases with ω(n) ∈ {1,2,3} and 0 failures. Dropping μ(n) makes every case with ω(n) odd fail.

Consequences:
- The factor μ(n) = (−1)^{ω(n)} is not a Hecke character or Grössencharacter.
- So E = Σ_k |Σ_P c_P F(k³P⁴)|², where F(h) = Σ_n μ(n) G 𝔯 conj γ₁(n;h) V(q_n/N). This is the second moment of **Möbius-twisted** sextic Kubota sums at the thin index set h = k³P⁴.
- The marks are the P⁴ part of the index; on the 3-fold side they are U_P-translates, τ(n)χ_P(n)^{−2} = τ(Pn)/τ(P).
- By (∗) again, the genuine Eisenstein moment Σ_k |Σ_n γ₁(n;k³P⁴) V|² is the energy of μ·θ. That is the marked **Möbius** family, the Poisson/Davenport–Hasse dual of E (NOTES §6).
- So the 6-fold reformulation is P78's own R∘P. The automorphic object is θ⊗ψ_k on the 3-fold cover, at level k².

### 1.2 What Rankin–Selberg can and cannot see

The de Faveri–Dunn–Hoffstein / David–de Faveri–Dunn–Stucky inputs, from their abstracts [L]:
- a Lindelöf-on-average bound for the second moment of Kubota series over **all** indices h ≤ H;
- Rankin–Selberg (Zagier-regularised) lower bounds for coefficients of θ×θ.

Why these do not reach E:
- **Positivity loses too much.** Restricting the all-index bound to h = k³P⁴ loses H/#{indices} = K²L³.
- **Rankin–Selberg sees only the diagonal.** On the 3-fold cover, RS of |θ⊗ψ_k|² sees only Σ|τ|². No RS integral produces the k-off-diagonal ψ_k(n₁n₂):
  - ∫θθ̄·E^{(2)}(·,w) over the double cover vanishes, since it pairs a non-genuine function with a genuine one [P];
  - inserting the quadratic theta θ₂ yields additive convolutions (n₁ − n₂ + m² = r), not ψ_k(n₁n₂).
- **The only generating object is an affine multiple Dirichlet series (MDS).** It is Z_{P₁P₂}(s₁,s₂,w) = Σ_k q_k^{−w} D_{P₁}(s₁,k) conj D_{P₂}(s₂,k).
  - Its FE group is infinite, with quadratic growth and the degenerate invariant form (x₁+x₂+y)², so it is **affine C̃₂** [P].
  - It is the same group as the Soundararajan–Young / Li moment Σ_d L(s₁,f⊗χ_d) L(s₂,f⊗χ_d) d^{−w}, for which only log savings are known [L].
  - By contrast, the cubic-family second moment that David–de Faveri–Dunn–Stucky handle with a power saving has a **finite** FE group (order 24 or 48), per `fe_group_compare.out`.

### 1.3 The missing input, and why the MDS/RS machinery cannot supply it [P/B/C]

The group fixes I = x₁ + x₂ + y [P]. Consequences:
- Bochner/Phragmén–Lindelöf hulls of any starting region never go below its inf I.
- Absolute convergence (and convexity) gives I ≥ 3/2.
- Large-sieve or Lindelöf-on-average input gives I ≥ 1/2, with level cost L¹.
- The orbit of the diagonal polar divisor is {aI ± y = ½} and similar, accumulating at I = 0, the natural boundary.

Contour bookkeeping at the 27/31 geometry (M = 0.860, ℓ = 0.204, w = 0.064; `mds_reach.out`):

| mark pairs | what W needs |
|---|---|
| diagonal, P₁ = P₂ | Lindelöf on average at I = ½, **uniformly in the level** |
| cross, P₁ ≠ P₂ | a bound at I ≤ ½ − ℓ/M = 0.263, strictly below the large-sieve line, so unreachable |

**[C] Theorem.** Assume, for |t_i| ≤ Z^ε:
- **(H1)** Σ_k W(q_k/K) |D_P(½+it,k)|² ≪ K^{1+ε} L^ε;
- **(H2)** Σ_k W(q_k/K) D_{P₁}(½+it₁,k) conj D_{P₂}(½+it₂,k) ≪ K^{1+ε}/L for P₁ ≠ P₂.

Then E ≪ KN Z^ε for every N ≤ KL. This is W(∞), apart from the genuine bias term of §4.3.

Proof: Mellin on Re s = ½, then sum over (P₁, P₂).

Notes on the hypotheses:
- Both are true in the random model, since L ≤ K^{1/2} ⟺ ℓ ≤ M/2.
- (H1) alone is useless: Cauchy–Schwarz over P loses ℓ > w.
- (H2) with a saving L^{−θ₂} gives a uniform energy excess (1−θ₂)ℓ (`Hoff` runs, §5).
- (H1) is Lindelöf on average for a *short* family: K quadratic twists of a level-L² form. It is open, like short-t mean values of a fixed L(s,χ).
- (H2) is a cross-correlation statement in the same family, also open.

**Answer to the brief's question.** RS/MDS gives no power saving that is uniform in the length or level. The spectral input needed is not a spectral gap or a triple-product bound: there is no integral representation to put them into.

## 2. Route 2: shifted convolutions and exact orthogonality

Three no-gos.

1. **Completion lemma [P].**
   - Suppose a family F dominates the Gram matrix G = Σ_k ψ_kψ_k* on the columns and is exactly orthogonal there (Σ_F f f* = c·Id).
   - Then c ≥ ‖G‖_op ≥ ‖ψ_{k₀}‖² ≍ N, which is the column term.
   - A full character group containing every ψ_k (k ≤ K) has modulus divisible by lcm(k), which is ≥ e^{(1+o(1))K} ≫ N. Its congruence is then equality on the columns: no shifted convolution sums appear and the bound is flat.
   - BFKMM works only because their family *is* the full group mod q; there one evaluates, rather than bounds.
2. **Exact identity [P].**
   - Expanding ψ_k by Gauss sums gives E = Σ_k q_k^{−1} Σ_{t mod k} (t/k)₂ Σ_{n ≡ tm (k)} a(n)ā(m).
   - These are *multiplicatively* shifted correlations n − tm = hk, i.e. coefficients of θ(z)·conj θ(t̄z) at level ≍ q_t.
   - The t-sum weighted by (t/k)₂ re-creates ψ_k(nm) exactly. Only t = 1 (a genuine shifted convolution) is spectrally accessible, and isolating it means dropping (t/k), which is completion as in item 1.
3. **Extra averaging [B].**
   - Twisting B by ξ mod q₀ with a coherent signal multiplies the rows by q₀.
   - It multiplies the reflected dual length by q₀² (conductor q₀²). So N/K grows by q₀ and the wall moves the wrong way.
   - Dilations of m are the same twist; Mellin-height averaging is NOTES §9.

Result: w = 0.

## 3. Route 3: Kuznetsov / amplification on the cubic cover

**No-go for individual bounds [B].**
- Assume even Lindelöf for every L(½+it, θ⊗ψ_kχ_P^{−2}).
- Then |S_{k,P}|² ≪ N^{1+ε}, and with Σ_P |c_P| ≍ L^{1/2} one gets E ≪ K·L·N·Z^ε.
- This is *worse* than the large sieve's N²: KL/N = Z^{1−M} > 1.
- Individual bounds lose the coherence of the P-sum, and Λ_k = Σ_P c_P D_P is not an L-function: its level is ΠP².
- No subconvexity exponent helps. Proskurin/Louvel Kuznetsov amplification (a family of size ≍ K² per level k², with non-multiplicative coefficients, so amplifiers must use T_{p³}) yields at best individual subconvexity.

Result: w = 0.

## 4. Route 4: special structure of β·X

1. **Hecke structure [P].**
   - The marks are internal (U_P; j = 4 Ramanujan factor), so S_{k,P} = F(k³P⁴) is a "multiplicative Hankel" matrix.
   - F is not multiplicative in h, so it has no factorisation to exploit.
2. **W must hold for all slot coefficients**, which is equivalent to σ_max(S)² ≲ K′D Z^ε [N].
   - Setup: N₀ = 40000, 184 marks of norm 1500–3000. The wall would allow K ≤ N₀/P ≈ 13–27.
   - For K′ = 222 … 14110, σ_max²/(K′D) = 4.54 … 1.50, against the random-phase control 4.52 … 1.48 and the random-matrix-theory (RMT) edge (√K′ + √L′)²/K′ = 3.65 … 1.24.
   - The second singular value is equally close to the control. Max per-P energy is 1.03–1.19 times the diagonal.
   - An aligned adversary is detected (max row 52–68 times the diagonal).
   - The N₀ = 20000 run (70 marks) agrees.
   - So W looks true **uniformly in the marks**.
3. **Patterson's bias [H, consistent with de Faveri–Dunn–Hoffstein sharpness].**
   - For unweighted slots, Patterson's bias inside X(n) gives X(n) ⊃ −c·conj τ(εn)(L²/N)^{1/6}/log L, with ε a fixed unit class.
     It comes from the dual marks ã(P)χ_P(n)^{−2} = −a(P)·g̃₃(P; ε_P n), a prime-modulus Kubota sum.
   - So βX has a non-oscillating part ∝ |τ(n)|². This part is orthogonal to every non-principal ψ_k (up to Pólya–Vinogradov), so it lives only on principal rows (k_res = 1 dyads).
   - There it has size N^{5/3}L^{2/3}. That is exactly the (NL)^{2/3} term of Lemma `hybrid-energy`: an excess of (M+4ℓ−2)/3 = −0.108 at the target geometry.
   - So it is a genuine main term that never binds near 27/31.
   - Model check (`ms_*.json`):

     | wE | piece relaxed (energy-wall convention) | piece kept |
     |---|---|---|
     | 0.064 | 0.8709554 | 0.8709554 |
     | 0.1 | 0.868705 | 0.868705 |
     | 0.2 | 0.863333 | 0.863636 |
     | 0.3 | 0.861667 | 0.863636 |
     | ∞ | 0.85 | **0.863636 = 19/22** |

   - So if this piece is kept, the energy lever saturates at 19/22. That is a ceiling for any W(w) with large w, unless the (NL)^{2/3} term is shown to be an artefact on the generic dyads (q_k ≍ K, non-principal ψ_k), which I expect is true.
   - It is irrelevant for 27/31. It cannot be exploited either: it is a positive main term on a few rows.
4. **Spectral expansion of X [B].**
   - For smooth slots this is Kuznetsov with one index equal to 0. It turns βX into Rankin–Selberg coefficients θ×f̄_j, which are degree-4 objects.
   - The quadratic twists then have AFE length K², which is further beyond the edge.

## 5. Model numbers (`model_ws.py`, a copy of energy-wall's `model_wall.py`)

| hypothesis | θ |
|---|---|
| routes 1–4 unconditionally (w = 0) | **0.874957** (unchanged) |
| W(0.064), with the bias/(NL)^{2/3} piece kept | 0.870955 |
| (H1)+(H2): E ≪ KN Z^ε for all N ≤ KL, bias kept (`Hoff=0`) | **0.863636 = 19/22** |
| (H1)+(H2_θ₂): excess (1−θ₂)ℓ′, `Hoff` = 0.03 / 0.05 / **0.07** / 0.1 / ≥0.25 | 0.866648 / 0.868761 / **0.870942** / 0.874352 / no gain |

The `Hoff` rows are min(paper, offset-branch LP). The union of the two convex branches is exact here, because at d = ℓ the energy piece vanishes. So 27/31 needs:
- (H1) with level loss at most L^{0.07};
- (H2) with saving at least L^{−0.93}.

The optimum then sits far beyond the wall: (l_x, l_y, ℓ) = (0.385, 0.522, 0.288), M+ℓ = 1.195. It satisfies ℓ ≤ M/2, so (H2) is still true in the random model there.

## 6. Bottom line

| route | status | w |
|---|---|---|
| 1 Eisenstein / RS / MDS | identity (∗) [P]; affine C̃₂ with invariant I [P]; needs (H1)+(H2) [C] | 0 |
| 2 shifted convolutions | completion cap at ‖G‖ ≥ N [P]; no congruence structure [P]; averaging worsens N/K [B] | 0 |
| 3 Kuznetsov / amplification | even individual Lindelöf gives KLN > N² [B] | 0 |
| 4 β·X structure | marks internal (U_P), no factorisation; Patterson bias = (NL)^{2/3} term, harmless at 27/31 [H]; W numerically uniform in marks [N] | 0 |

The precise missing input is (H1)+(H2): Lindelöf on average, plus cross-mark decorrelation, for the K quadratic twists of the level-L² cubic theta.
- This is a *short* family: K against an AFE length of KL.
- It is the same affine class as power-saving second moments of quadratic twists of GL2 forms, uniformly in the level.
- No spectral tool known to me reaches it. RS and MDS methods cannot reach it by the I-invariance.
