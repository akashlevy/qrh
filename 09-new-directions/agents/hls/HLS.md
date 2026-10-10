# hls: the hybrid large sieve for {ψ_k χ_P^{-2}} — truth, proofs, payoff

Agent `hls`, 10 Oct 2026. Labels: **[P]** proved (standard analytic details: smooth weights, gcd/Möbius
separation, Mellin separation, all costing (KLN)^ε); **[C]** conditional; **[H]** heuristic; **[N]** numerical.

**Verdict.**
1. (HLS) is **true numerically** with a random-matrix constant [N]: λ_max/(rows+cols) = 1.6–1.9, flat in scale, and equal
   to a random product control within 1.6%.
   - There is **no theta obstruction** in the needed range [N + H].
   - Cubic and quadratic Gauss-sum vectors cannot absorb ψ_k for k ≠ 1 (verified: no Kubota pole). Their theta mass sits on the
     k = 1 rows, where Heath-Brown gives ≤ KL + N because L ≤ K.
   - The sextic Gauss-sum vectors do see the rows, at Kubota index k³P^{2 or 4}. But the 6-fold theta coefficients at cubes obey
     |τ₆(p³h)| ≈ N(p)^{−1/6}|τ₆(ph)| (measured over N(p) = 18 … 2170), and they are multiplicative at composite cubes.
     The family is therefore theta-poor, with margin K^{1/6} at the edge under the Rankin–Selberg average law.
     Without that suppression (HLS) would fail at the edge.
2. **Proved:** D(N) ≪ (KLN)^ε·min{L(K+N), K(L+N+(LN)^{2/3}), N + KL² + KL√N} [P].
   - This is optimal only in the trivial ranges N ≤ K and N ≥ K²L².
   - It beats P78's sequential bound only for N > K².
   - In the needed range K < N < KL it is exactly the sequential bound: the loss is N/K, which is the wall.
   - Obstructions [P, structural]:
     - Poisson in n and in k are involutions with fixed point N = KL;
     - the smoothed (all-n) family that every Poisson argument needs is genuinely larger than KL + N, since it is ≥ K√N from squares;
     - the off-P-diagonal blocks need coherent cancellation across blocks, which is an (H2)-type decorrelation, but coefficient-free.
3. **Payoff alone: zero.** (HLS) with Cauchy–Schwarz over the marks reproduces P78 exactly (break-even,
   0.8749568). (H2) alone is also 0.8749568. Only (HLS)+(H2) together give 19/22. No other P78 lemma
   (reflected energy, Lemma marked's 2z, DIAG) improves, because their mark sums are coherent [P, model].

## 0. Objects

- S_{k,P}(a) = Σ_{n~N} a_n ξ_{kP}(n), ξ_{kP} = ψ_kη_P, ψ_k = χ_k³ (quadratic, conductor k), η_P = χ_P^{-2} (cubic, conductor P).
- D(N) = D(N;K,L) := sup_a Σ_{k≤K,P~L}|S_{k,P}(a)|²/‖a‖²; (HLS) is D ≪ (KL+N)^{1+ε}.
- Sextic form: ξ_{kP}(n) = 𝔯(n,k)·χ_n(k³P⁴) (sextic reciprocity, 𝔯 = sign bicharacter mod 4).
- Gram forms:
  - row side G = SS^*, KL×KL, with K×K blocks G^{PP'} = Ψ diag(η_Pη̄_{P'}) Ψ^T;
  - column side S^*S = A∘B (Schur product), A = Ψ^TΨ (quadratic Gram, rank ≤ K), B = H^*H (cubic Gram, rank ≤ L).
- Needed geometry (19/22 optimum of the conditional theorem, `model_hls.py`):
  - M = 0.909, ℓ = 0.273, so K = Z^{0.909}, L = Z^{0.273}, L ≈ K^{0.30};
  - N = K·Z^{0.18} = KL·Z^{−0.09}, so N sits Z^{0.09} below the edge KL.

## 1. Is (HLS) true?

### 1.1 Gram spectrum [N] (`gram_scan.py`, `runs/scan1.{log,json}`)

Exact Z[ω] symbols (via `edge-theta/edge_lib.py`), λ = top eigenvalue of SS^* by Lanczos. Controls:
- 'prod': random signs for ψ, random cube roots for η, same zero pattern;
- RMT: (√rows + √cols)²;
- HB-quad: Heath-Brown's family alone with the same number of rows;
- seq: P78's sequential bound min(#P·λ(ΨΨ^T), #k·λ(HH^*)).

| N0 | K | L0 | cols | k's × P's = rows | λ/(rows+cols) | λ / prod-ctrl | λ / RMT | seq/λ | HB-quad λ/(rows+cols) |
|---|---|---|---|---|---|---|---|---|---|
| 5000 | 400 | 50 | 1416 | 109×10 = 1090 | 1.869 | 1.0070 | 0.939 | 5.0 | 1.835 |
| 20000 | 1000 | 100 | 5632 | 282×20 = 5640 | 1.888 | 1.0007 | 0.944 | 8.8 | 1.873 |
| 20000 | 200 | 100 | 5632 | 57×20 = 1140 | 1.636 | 0.9840 | 0.936 | 12.3 | 1.942 |
| 20000 | 200 | 300 | 5632 | 57×44 = 2508 | 1.817 | 0.9957 | 0.945 | 20.2 | 2.004 |
| 20000 | 200 | 1000 | 5632 | 57×136 = 7752 | 1.877 | 0.9939 | 0.945 | 17.9 | 1.805 |
| 20000 | 60 | 3000 | 5632 | 15×354 = 5310 | 1.879 | 0.9965 | 0.940 | 6.7 | — |
| 80000 | 2500 | 200 | 22579 | 701×32 = 22432 | 1.901 | 1.0005 | 0.950 | 14.3 | — |
| 10000 | 1250 | 40 | 2816 | 349×10 = 3490 | 1.867 | 0.9997 | 0.936 | 4.7 | — |
| 20000 | 2500 | 40 | 5632 | 701×10 = 7010 | 1.877 | 1.0033 | 0.941 | 5.0 | — |
| 40000 | 5000 | 40 | 11311 | 1415×10 = 14150 | 1.882 | 1.0071 | 0.944 | 5.1 | — |
| 80000 | 10000 | 40 | 22579 | 2831×10 = 28310 | 1.883 | 1.0069 | 0.944 | 5.1 | — |

- **The arithmetic family is indistinguishable from the random product control** (to ≤ 1.6%), and λ/(rows+cols) is flat at 1.6–1.9 over a 16-fold range of N0. The sequential bound is 5–20 times too large.
- The edge rows (rows ≈ cols) are flat in N0: 1.869 → 1.888 → 1.901.
- The needed shape N < KL (rows = 1.24·cols, L0 = 40, the last four rows, from `runs/scan2`) is also flat: 1.867 → 1.877 → 1.882 → 1.883 for N0 = 10⁴ … 8·10⁴.
- The top eigenvector is delocalised: max row share ≤ 0.7%, k = 1 share ≤ 2% (7% when K = 60).

### 1.2 Theta (Gauss-sum) vectors [N]

Normalised as ρ(a) = Σ|S a|²/Σ_{k,P}Σ_n|a_n|²|ξ_{kP}(n)|², with ρ = 1 the diagonal.

Vectors tested:
- γ_j(n), the normalised Gauss sum of χ_n^j, j = 1..5: sextic for j = 1, 5, cubic for j = 2, 4, quadratic for j = 3;
- μγ_j;
- β = γ₂ᾱ (P78's theta), μβ;
- the column extremal conj ξ_{k0P0};
- the top eigenvector.

Results:
- Every Gauss-sum vector has ρ = 0.89–1.07 in every configuration.
- The only systematic excess is γ₂/γ₄/β on the k = 1 rows: ρ(γ₂) = 1.15 with k = 1 against 0.91 without, at N0 = 2·10⁴, K = 200. This is the cubic theta term, which Heath-Brown bounds by (L + N + (LN)^{2/3}) ≤ KL + N because L ≤ K [P].
- For comparison, ρ(top) = 3.4–10 and ρ(column) = 1.75–6.1.

### 1.3 Which theta vectors can absorb the hybrid twist [P + N] (`sextic_kubota.c`, `kubota_analyze.py`, `kubota_fit.py`)

A Gauss-sum vector a_n = γ_m(n) turns S_{k,P} into a Kubota sum Σ_n γ_m(n; h) with γ_m(n;h) = conj χ_n(h)^m γ_m(n). This needs
conj χ_n(h)^m = χ_n(k³P⁴), i.e. −m·e_h ≡ (3, 4) on the exponents (mod 6).
- The cubic (m = 2, 4) and quadratic (m = 3) vectors **cannot** absorb the quadratic part ψ_k for k ≠ 1: 2e ≡ 3 (mod 6) has no solution.
  For them, ψ_k is a genuine quadratic twist of the cubic Kubota series, which has no pole.
- Only the sextic vectors absorb it: m = 1 gives index k³P², and m = 5 gives index k³P⁴. The theta residue then needs the
  6-fold theta coefficient τ₆ at a cube p³ for every p | k.

Numerical test with exact Gauss sums for all squarefree odd primary c, N(c) ≤ X:
- T_{m,j}(X) = Σ γ_m(c)·conj(π₀/c)_6^{jm}, the Kubota sum at index π₀^j;
- test primes π₀ with N(π₀) = 7 … 313 and 1009 … 2389;
- statistics: growth of mean |T|²/#c (pole ⇒ X^{2/n}), and the coherence C of the fitted pole amplitude
  between [X/64, X/8] and [X/8, X] (C ≈ 1 for a pole, ≈ 0 ± T^{−1/2} without).

**Pole detection** (`kubota_analyze.py`, `kubota_fit.py`; X = 2·10⁶, 12 test primes; `runs/kubota_*2e6.json`). C is the coherence.

| vector | index p^j: j = 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| cubic γ₂ (C) | pole (1.00) | pole (1.00) | **none** (0.02) | — | — | — |
| cubic γ₂ × quadratic twist (s odd) | none | none | none | | | |
| sextic γ₁ (C) | pole (0.99) | pole (0.85) | pole (0.87) | weak (0.65) | weak (0.50) | **none** (0.21) |
| sextic γ₅ (C) | pole (0.995) | pole (0.61) | pole (0.71) | weak (0.39) | weak (0.59) | **none** |

**Residue sizes** ρ_j(p) = |τ(p^j)/τ(1)|², from the unbiased cross-range estimator (`kubota_prime.py`, `kubota_pairs.py`).
Data:
- X = 6·10⁶ with 30 small primes N(p) ∈ [7, 313] and 95 large primes N(p) ∈ [1009, 2389] (Gauss sums reused via `-PRIMES=`);
- X = 2·10⁶ with 12 small and the same 95 large primes;
- files `runs/kubota_prime_{2e6,6e6}{,_big}.json`.

The j = 0 normalisation drifts by about 30% between X = 2·10⁶ and 6·10⁶, so ratios carry roughly that systematic error.

| X = 6·10⁶ | N(p) ≈ 18 | 76 | 162 | 259 | 1175 | 1645 | 2170 |
|---|---|---|---|---|---|---|---|
| cubic m = 2, j = 1 | 0.349 | 0.239 | 0.191 | 0.162 | 0.0971 | 0.0845 | 0.0782 |
| N(p)^{−1/3} (Patterson) | 0.38 | 0.236 | 0.183 | 0.157 | 0.095 | 0.085 | 0.077 |
| sextic m = 1, j = 1 (squarefree index) | 0.46 | 0.66 | 0.57 | 0.45 | 0.31 | 0.41 | 0.36 |
| m = 1, j = 2 | 0.48 | 0.54 | 0.22 | 0.22 | 0.19 | 0.20 | 0.076 |
| **m = 1, j = 3** | **0.185** | **0.153** | **0.027** | **0.071** | **0.050** | **0.047** | **0.024** |
| ± | 0.057 | 0.064 | 0.040 | 0.047 | 0.015 | 0.019 | 0.010 |
| **ρ₃/ρ₁** | **0.40** | **0.23** | 0.05 | **0.16** | **0.16** | **0.12** | **0.067** |
| m = 1, j = 4 | 0.19 | 0.11 | −0.015 | 0.039 | −0.032 | 0.005 | −0.015 |
| m = 1, j = 5 | 0.02 | 0.00 | −0.03 | −0.06 | −0.02 | 0.00 | −0.01 |
| m = 5, j = 3 | 0.27 ± 0.17 | 0.14 ± 0.10 | 0.33 ± 0.19 | 0.36 ± 0.22 | 0.02 ± 0.04 | 0.04 ± 0.04 | 0.01 ± 0.04 |

At X = 2·10⁶ the large primes give m = 1, j = 3: 0.004 ± 0.010, −0.009 ± 0.013, −0.018 ± 0.014, which is smaller still.

Findings:
- **The estimator is calibrated.** The cubic case reproduces Patterson's |τ₃(p)|² = N(p)^{−1/3}|τ₃(1)|² to 1–3% over N(p) = 18 … 2170, and τ₃(p²) = 0 exactly.
- **Quadratic twists kill the cubic theta.** The cubic Kubota sums twisted by s odd show no pole at any index [N]; the theory is standard [H].
- **Sextic coefficients at cubes are suppressed by a pairing law.** [N] ρ₃(p)/ρ₁(p) ≈ N(p)^{−1/3}: 0.40/0.38, 0.23/0.24, 0.16/0.16, 0.16/0.095, 0.12/0.085, 0.067/0.077, with one low outlier at N(p) ≈ 162.
  That is, |τ₆(p³h)| ≈ N(p)^{−1/6}|τ₆(ph)|.
  - Heuristic reason [H]: the Hecke (T_{p⁶}) relation pairs the exponents j ↔ 4 − j (mod 6); with the Kubota local factors |c_j| = N(p)^{−2/3−j/6}
    it gives |τ(p³h)/τ(ph)| = |c₃/c₁|^{1/2} = N(p)^{−1/6}.
  - The self-paired index p⁵ is 0, as observed.
- **The squarefree coefficients ρ₁(p) are O(1)** with weak decay (0.46–0.66 down to 0.31–0.41). The Rankin–Selberg average law Σ_{N(h)≤H}|τ₆(h)|² ≍ H^{5/6}, i.e. mean N(h)^{−1/6} (the source of de Faveri's A^{5/6}B^{1/3}), predicts N(p)^{−1/6} [H].
  - Hence ρ₃(p) ≈ N(p)^{−γ₃} with γ₃ = 1/3 + γ₁: γ₃ ≥ 1/3 from the pairing law alone, and 1/2 with the average law.
  - The fitted value over N(p) = 18 … 2170 is 0.3–0.45.
- **Composite cubes are multiplicative** [N]. At p₁³p₂³ the mean residue is 7.2·10⁻⁴ ± 4.4·10⁻⁴ against the product prediction 6.6·10⁻⁴; at p₁³p₂² it is 1.6·10⁻³ ± 0.5·10⁻³ against 1.5·10⁻³ (m = 1, 66 pairs, X = 2·10⁶).

**Consequence [H].** The only theta vectors that see the k ≠ 1 rows are γ₁ and γ₅. Their Rayleigh excess is
≈ N^{1/3}·Σ_{k≤K}Π_{p|k}ρ₃(p)·Σ_P ρ_{2 or 4}(P) ≍ N^{1/3}K^{1−γ₃}L^{1−γ₂}.
- Measured γ₂ ≈ 0.25 (index P² for γ₁); index P⁴ for γ₅ decays faster.
- (HLS) survives the theta test iff N^{1/3} ≲ K^{γ₃}L^{γ₂}.
  - **At the 19/22 point** (K = Z^{0.909}, L = Z^{0.273}, N = Z^{1.09}) this reads 0.909γ₃ + 0.273γ₂ ≥ 0.364.
    The pairing law alone (γ₃ = 1/3, γ₂ = 0.25) gives 0.371, with no margin to spare. With the average law (γ₃ = 1/2) it gives 0.52, a margin of Z^{0.16}.
  - **At the exact edge N = KL** it needs γ₃ > 1/3 (for γ₂ ≈ 1/3), i.e. a genuine decay of the squarefree coefficients beyond the pairing law. The margin is K^{γ₃−1/3} = K^{1/6} under the average law.
- **Verdict on theta.** (HLS) is theta-poor, but only by the margin K^{γ₁} that comes from the decay of the squarefree sextic theta coefficients. Without the pairing law (γ₃ = 0) it would **fail** at the edge by K^{1/3}. This is the one place where (HLS) is genuinely sensitive to arithmetic.
### 1.4 Spectral argument [P/H]

- **Lower bound.** λ_max ≥ max(KL, N)(1 − o(1)): the diagonal gives KL, and a = conj ξ_{k0P0} gives N.
- **Rank.** S^*S = A∘B has rank ≤ KL and trace KLN.
- **Reformulation.** Since A = Ψ^TΨ ≤ (K+N)Π_V, with V = span{ψ_k} of dimension K, we get
  S^*S ≤ (K+N)·Σ_P Π_{η_P V}. So
  **(HLS) ⟸ ‖Σ_P Π_{η_PV}‖ ≪ (1 + KL/N)Z^ε,** i.e. the L cubic rotations of the quadratic subspace are jointly nearly
  orthogonal. Random K-dimensional subspaces give (1 + √(KL/N))².
- **Theta mass.** A theta extremal would need a vector whose Rayleigh quotient carries Kubota residues on a positive
  proportion of rows. The family's total theta mass is Σ_{k,P}|τ₆(k³P^j)|²·N^{1/3} ≍ N^{1/3}K^{1−γ₃}L^{1−γ₂} (§1.3).
  With γ₃ = 1/2 and γ₂ = 1/3 this is ≪ KL + N throughout N ≤ K^{3/2}L, which contains the needed range [H].
- **Contrast.** Ordered by the norm of a, the full sextic family has theta mass A^{5/6}N^{1/3}, which is de Faveri's term.
  **The hybrid family is theta-poor**: its rows are cube-indexed (k³), and theta coefficients decay at cubes.

## 2. What can be proved

### 2.1 Proposition 1 [P]

For N(P) ~ L, N(k) ≤ K, N(n) ~ N, and every ε > 0:

    D(N; K, L) ≪_ε (KLN)^ε · min{ L(K+N),  K(L + N + (LN)^{2/3}),  N + KL² + KL√N }.

In particular (HLS) holds for N ≤ K (sequential) and for N ≥ K²L² (as for the general large sieve). The third bound
beats the sequential one exactly when N > K²; at the edge N ≈ KL (< K²) nothing better than L(K+N) is proved.

*Proof.* The first two bounds are the sequential ones (Heath-Brown quadratic, then cubic; or the reverse). For the
third, majorise the squarefree n-sum by the smooth sum over **all** n (positivity) and split G into K×K blocks:
- **Diagonal blocks.** Σ_{all n}W|Σ_k x_kψ_k(n)|² ≤ (N + K√N)‖x‖², because ψ_k(n₁n₂²) = ψ_k(n₁)1_{(k,n₂)=1}. Apply Heath-Brown on each square class.
- **Off-diagonal blocks.** For P ≠ P', the block form is Σ_n W(N(n)/N)·η_Pη̄_{P'}(n)·X(n)·Y(n), with X = Σx_kψ_k and Y = Σy_kψ_k.
  - Separate gcd(k,k') by Möbius, and apply Poisson modulo k₁k₁'PP'. The zero frequency vanishes because P ≠ P'.
  - The Gauss sum factors by CRT. **For fixed (P,P') every cross factor is separable**: (P/k₁)₂(k₁/P)₃², (P'/k₁)₂(k₁/P')₃,
    and their k₁' analogues. The (k₁,k₁') cross factor is a quadratic-reciprocity sign fixed by classes mod 4.
  - The Fourier weight separates by Mellin inversion.
  - The block becomes (N/KL)·Σ_{N(u) ≲ U}(η̄_Pη_{P'})(u)X'(u)Y'(u), with U = K²L²/N, unimodular re-phasing x→x', y→y',
    and X' = Σx'_kψ_k(u).
  - Cauchy–Schwarz, then Heath-Brown on each square class u = u₁u₂² of the dual variable, gives
    Σ_u|X'(u)|² ≪ (U + K√U)‖x‖². So ‖G^{PP'}‖ ≪ (N/KL)(U + K√U) = KL + N/L + K√N.
- **Schur test over blocks.** ‖G‖ ≤ max_P Σ_{P'}‖G^{PP'}‖ ≪ N + K√N + L(KL + N/L + K√N). ∎

**Remark (the square part is essential) [P].** For the smooth all-n family, c_{kP} = δ_{P,P₀}K^{−1/2} gives
|Σ c ξ(n)|² ≍ K on every square n = n₂², so D_all(N) ≫ K√N. This exceeds KL + N whenever L² < N < K², and in particular at
the edge. Any Poisson proof of (HLS) at the edge must therefore remove the non-squarefree n first, as Heath-Brown does for
the quadratic sieve. Heuristically the K√N terms of the dual blocks are cancelled by the cubic twist (η̄_Pη_{P'})(u₂²);
if so, the third bound improves to N + KL², still useless at the edge.

### 2.2 Proposition 2 (self-duality; why the needed range is untouched) [P]

- (a) **Poisson in n**, per off-diagonal block of the smooth all-n family: G^{PP'}_N = (N/KL)·D₁ G^{(P,P')}_{U} D₂ with
  U = K²L²/N, where D_i are diagonal unitaries and G^{(P,P')}_U is the block of the conjugate family at length U. So block norms satisfy β(N) = (N/KL)β(K²L²/N).
  This is an involution with fixed point N = KL.
- (b) **Poisson in k** (the cubic factor rides along, because it does not depend on k):
  G_K = KL·I + G^{off}_K with G^{off}_K ≅ (K/N)·G'^{off}_{N²/K}. This is the same family at quadratic length K* = N²/K and
  the same N and L. A bound D(N;K*) ≤ αK*L + βN returns D(N;K) ≤ KL + αNL + βK. Since the dual diagonal is K*L
  (so α ≥ 1), this never beats sequential.
- (c) **Theta reflection R** (for a = theta): length N ↔ (KL)²/N. It maps the pre-wall range N ≤ K onto N' ≥ KL², which is
  exactly where the block bound of Prop 1, with the square part removed, would be sharp. This is "R stops at the wall" in
  sieve language.
- So the needed range K < N < KL is invariant under every available duality, and its centre N = KL is the common fixed point.

### 2.3 Proposition 3 (block-wise arguments lose; the cocycle) [P + N]

- **Lower bound on blocks.** Any bound that uses only block norms is ≥ L·min_{P'≠P}‖G^{PP'}‖.
  - Heuristically ‖G^{PP'}‖ ≳ K: x = (ψ_k(n₀))/√K concentrates X at n₀ [H].
  - In the random model ‖G^{PP'}‖ ≈ 2√(KN) [H]; the measured block norms are in the table below [N].
- **Consequences.**
  - Block-Schur is at best ≈ N + 2L√(KN), i.e. 2KL^{3/2} at N = KL, a loss of L^{1/2}.
  - Cotlar–Stein is at best ≈ KL^{7/4}.
  - Even with square-root block norms (an unproved "twisted quadratic sieve"), block-Schur gives a slope-1/2 wall, which
    buys nothing in the model (§3).
- **The cocycle.** Without blocking, the full Poisson dual is (N/KL)·⟨c', (M∘G^{(U)})c'⟩ with
  M_{(kP),(k'P')} = conj(k/P')₆·(k'/P)₆.
  - M is the **partial transpose** of a rank-one matrix of sextic symbols.
  - Its diagonal reproduces exactly the KL term of (HLS).
  - (HLS) ⟺ ‖M∘G'^{(U)}‖ ≲ U. The Schur-multiplier norm of M is L, and the factorisation is explicit.
  - M acts on block (P,P') by the diagonal unitaries diag(Φ̄(·,P')) and diag(Φ(·,P)), so it preserves block norms but
    scrambles inter-block phases.
  - The dual is therefore not a square, and no Heath-Brown recursion closes.
  - (HLS) at the edge is equivalent to **coherent cancellation across P-blocks**. That is an (H2)-type decorrelation
    statement, for characters instead of theta coefficients.

**Block numerics [N]** (`block_decor.py`, `runs/block_*.json`). Columns: norm of the full G; mean and max off-diagonal block norm against 2√(K'N'); the block-Schur sum against ‖G‖; ‖G_off‖ = ‖G − blockdiag‖.

| family (k' × P', cols) | ‖G‖ | diag blocks | off blocks mean / max (2√(K'N')) | block-Schur sum | ‖G_off‖ |
|---|---|---|---|---|---|
| arithmetic 282 × 20, 5632 | 21287 | 9356 | 2413 / 2625 (2520) | 55424 (2.6‖G‖) | 15230 |
| random product 282 × 20 | 21283 | 7966 | 2413 / 2491 | 54154 | 15269 |
| arithmetic 57 × 44, 5632 | 14788 | 6785 | 1040 / 1136 (1133) | 51783 (3.5‖G‖) | 9292 |
| random product 57 × 44 | 14880 | 6461 | 1044 / 1148 | 51838 | 9336 |

- The off-diagonal blocks have exactly the random-model size 2√(KN). The twisted quadratic Gram behaves like a random matrix, so a "twisted quadratic sieve" with square-root saving is the truth.
- The block-Schur sum overestimates ‖G‖ by a factor growing like √L' (2.6 at L' = 20, 3.5 at L' = 44).
- The true norm comes from phase cancellation **across** blocks, which the arithmetic family shares with the random control.

### 2.4 Other routes checked (all ≥ sequential at the edge) [P/B]

- **Hölder between the two sieves.** Use the SVD of the dual coefficient matrix and fourth moments via Heath-Brown at
  length K² and the cubic sieve at length L². This gives D ≤ (‖s‖₁²/‖s‖₂²)·((N+K²)(N+L²+(NL²)^{2/3}))^{1/2} ≤ L·K√N
  at the edge, i.e. K^{3/2}L^{3/2} at N = KL. That is worse than the sequential KL² whenever K > L.
- **Moments of the family / Halász–Montgomery.**
  - The Gram row ℓ² norm is computable: Poisson again, and the cross phases are unimodular per row, so |·|² kills them.
    It gives Σ_{(k',P')}|G_{(kP),(k'P')}|² ≪ K²L² + KL²N.
  - Halász–Montgomery with Cauchy–Schwarz then gives D ≤ N + (KL)^{1/2}(K²L²+KL²N)^{1/2} = K^{3/2}L² at the edge. Worse.
- **"Moments of different lengths"** (OpenAI cubic first moment, `dual.tex`):
  - Its near-orthogonality lemma is exactly the Gram-row bound above.
  - The large-values contradiction needs coefficients that are bounded convolutions of prime sequences.
  - For (HLS) the coefficients are arbitrary. For the (H1) application they are μ(n)·(sextic Gauss sum) by identity (∗),
    which does not factor (energy-wall O1/O3).
- **de Faveri's optimal n-th order sieve.** In de Faveri's parametrisation the family is χ_a with a = k³P⁴, of norm
  K³L⁴ ≫ KL, so the bound is ≥ K³L⁴. It is useless unless the method can be re-proved family-by-conductor for a fixed
  exponent pattern. The cubic sub-family of the full sextic family shows that such a statement must exclude pure-cubic
  patterns.
- **Shift by a large prime p₀ (N(p₀) > K).** It gives monotonicity D(N) ≤ D(N·N(p₀)), which only transports the
  KL² + N bound downward.

## 3. Exponent model (`model_hls.py` = wall-spectral `model_ws.py` + knob `hls_a`; `run_model.sh`, `runs/model_hls.out`)

Knob `hls_a = a`: a partial hybrid sieve D ≪ KL + N·L^a in N ≤ KL.
- The diagonal-in-P energy Σ_P|c_P|²G_PP has excess M' + aℓ' − 1.
- The off-diagonal part is (H2) with saving L^{−θ₂} (Hoff = 1 − θ₂), or Cauchy–Schwarz over P (excess M' + (1+a)ℓ' − 1).
- The Patterson-bias piece (M'+4ℓ'−2)/3 is kept (bias_relax = false).
- a = 1 is P78's sequential bound; a = 0 is (HLS).

| inputs | θ |
|---|---|
| P78, none | 0.8749568 |
| (HLS) alone (a = 0, Cauchy–Schwarz over P) | **0.8749568**, exactly break-even (= energy-wall O2) |
| Prop 1 in N ≤ KL (a = 1) with Cauchy–Schwarz | 0.8816931; the paper's own bound is better, so 0.8749568 |
| (H2) alone (a = 1, Hoff = 0) | **0.8749568** |
| (HLS) + (H2) | **0.8636364 = 19/22** (geometry M = 0.909, ℓ = 0.273, N = KL·Z^{−0.09}) |
| partial sieve a + perfect (H2): a ≤ 1/3 | 19/22 (at the optimum the diagonal piece is slack once a ≤ 0.09/0.273) |
| a = 0.5 / 0.687 / 0.75 / 0.9 / 1 | 0.867593 / 27/31 / 0.871904 / 0.873849 / 0.874957 (identical to mark-design's `wall_l`) |
| (HLS) + (H2) at saving L^{−0.93} (Hoff = 0.07) | 0.870942 |
| a = 0.5 + Hoff = 0.07 | 0.873083 |
| block-wise bound with square-root block norms (an unproved "twisted quadratic sieve", §2.3) + perfect (H2) | **0.8749568**: diagonal excess max(M'−1, (M'+ℓ'−1)/2) is a slope-1/2 wall, useless (energy-wall: slope ≲ 0.10 needed) |

So the conditional theorem is an AND:
- (H2) is useless without a hybrid sieve that is **lossless** on a range w ≤ (1−a)ℓ;
- (HLS) is useless without (H2);
- even an idealised block-wise sieve buys nothing.

**Other P78 lemmas.**
- Lemma reflected-energy (l.7747) applies hybrid-energy to Σ_k|Σ_P a(𝐩)q_P^{−1/2}(…)|² with a(𝐩) arbitrary: a coherent mark sum.
- Replacing hybrid-energy by Cauchy–Schwarz over P plus (HLS) changes the exponent
  from max(H, v+ℓ_b) + max(v − z_a, 0, (2v − z_a)/3)
  to max(H + z_a, v + ℓ_b) + v.
  This is never smaller, and strictly larger whenever z_a > 0.
- Lemma marked's 2z (terminal KP ≤ N, DIAG notes §2.1) and the DIAG partial versions use only this lemma with coherent slots.
- Lemma plain's loss is a sup-norm/β* effect.
- Slot transfers are exactly break-even under (HLS) (energy-wall O2; reproduced here as hls_a = 0 without (H2) → 0.8749568).
- **So (HLS) improves no P78 lemma on its own.**

## 4. Files

| file | content |
|---|---|
| `hls_lib.py` | Family class: exact Z[ω] rows ψ_kη_P, columns, Gauss-sum vectors γ_j via CRT; fast S, S^*; Lanczos top eigenvalues; random controls |
| `test_hls_lib.py` | checks: β = edge-theta β, γ₁γ₂ = μαG, adjointness, Lanczos = brute force |
| `gram_scan.py`, `runs/scan*.{log,json}` | §1.1–1.2 |
| `sextic_kubota.c` (+ binary), `sk_read.py`, `sk_validate.py` | exact γ₁..γ₅ on all squarefree odd primary c ≤ X (validated against Python at 300 primes and all composites ≤ 4000) |
| `kubota_analyze.py`, `kubota_fit.py`, `kubota_prime.py`, `kubota_pairs.py`, `runs/kubota_*.json`, `runs/sk_*.log` | §1.3 (data binaries are in the session scratchpad, regenerable: `./sextic_kubota X out.bin threads p0...`, or `-PRIMES=prev.bin` to reuse stored Gauss sums) |
| `make_tables.py`, `run_scan2.sh` | table rows; second scan (N < KL shapes) |
| `block_decor.py`, `runs/block_*.json` | §2.3 |
| `model_hls.py` (knobs `hls_a`, `hls_f`), `run_model.sh`, `runs/model_hls{,_first,_f}.out` | §3 |
