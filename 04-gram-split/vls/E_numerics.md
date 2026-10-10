# Attack E: numerical ground truth for (V′_LS), (V′_β), Type I and Type II structure of π = ρ·ζ

*9 Oct 2026. Scripts and data: `gram_split/vls/` (data in `vls/data/`). Target and normalisations: `TARGET_VLS.md`. Builds on `wf2/A2_numerics.md` and its referee report.*

Labels:
- **[R]** exact identity, or exact implementation verified to machine precision;
- **[N]** statistical numerical observation, with its noise level stated;
- **[H]** heuristic interpretation;
- **[G]** gap or untested assumption.

Nothing in this note is a proof of (V′). Numerics of this kind can show that a statement is *false* or *plausible*. They cannot show that it is *provable* (see §8).

---

## 0. Bottom line

1. **New data, 16× larger in Q** [R].
   - ρ(m) is now estimated at Q = 6000, 12000, 24000, 48000 and 96000, up to 81,342 rows.
   - The data use dyadic scales Y = 2^8 … 2^21, so s runs up to 2^22. All 295,890 primes of norm ≤ 2^22 have exact sextic Gauss sums.
   - The new pipeline reproduces the existing `A_Q6000_K12` data to 4·10⁻¹⁷ on all 12 common scales.
   - ζ_m (Lemma Z) agrees with the paper's product formula to 1.6·10⁻¹⁴ on random rows at every Q.
   - Only scales Y ≥ 4Q are used in the main analysis (Y ≥ 5.5Q at Q = 96000). Windows Y ≥ 1.5Q, 2.5Q and 16Q serve as bias checks.
2. **(V′_LS) behaves exactly as for a random sequence: R(L,Q) ∝ L/Q. There is no saturation and no intermediate regime** [N].
   - **Ratio to the random model.** Write R_rand(L) = Σ_{N(ν)≤L} Σ_c |F′_c(ν)|² divided by its random-model expectation. Then R_rand(L) = 1.00 ± 0.05 for L ≥ Q^{0.7} at every Q.
   - **At L = Q^{0.79}** the values are 1.009 / 0.971 / 1.011 / 1.007 / 1.039, with jackknife s.e. 0.02–0.05, for Q = 6000 … 96000.
   - **Small L.** For L ∈ [Q^{0.5}, Q^{0.7}] the ratio stays within [0.90, 1.22]. That is the estimator's documented systematic scatter of ±20% (§2.3).
   - **Large-sieve ratio.** In the task normalisation, R_LS(L,Q) = S(L)/(Q Σ Ω|π|²) = κ·#ν(L)/Q with κ between 0.64 and 0.80 (mean 0.72) at every L and Q. At L = Q^{0.79} this gives 0.0330 → 0.0190 from Q = 6000 to 96000, i.e. ∝ Q^{−0.21}.
   - **Saturation (R_LS ≍ 1) is excluded** by a factor of 30 to 50 at L = Q^{0.79}, and by 4·10² to 10³ at L = Q^{0.5}.
   - **Power-law fits.** A fit R_LS ∝ (L/Q)^b Q^c over all (L, Q) gives b = 1.008 ± 0.007 and c = 0.002 ± 0.015. Independent dyadic bins give Rat(ν) ∝ (N(ν)/Q)^{−γ} with γ = 0.000 ± 0.004.
   - **Allowing for systematics,** b ∈ [0.94, 1.01]. Hence the effective (V′_LS) exponent at L = Q^{0.79} is η = 0.21·b ≈ 0.20–0.21, i.e. η ≈ 1 − log L/log Q.
3. **No large-sieve-extremal columns** [N].
   - max_ν Rat(ν) is 3.7 / 4.3 / 4.8 / 5.9 / 6.8 for Q = 6000 … 96000. Random-sign nulls on the same rows give a mean maximum of 3.7 / 3.8 / 4.4 / 4.9 / 6.2.
   - Per-ν ratios at N(ν) ≤ 512 are uncorrelated across Q (|corr| ≤ 0.13 over 146 columns).
   - Saturation would need some column with Rat ~ Q/L ≳ 30.
4. **(V′_β), the true bilinear form with Patterson's cubic theta coefficients, also sits at the random level** [N].
   - Σ_c |Σ_m Ω π B⁰_m|² divided by its random expectation is 0.75–1.38 (jackknife s.e. 0.1–0.35) in every dyadic window L ∈ [√Q, Q], including L = Q^{0.79}, for Q = 6000 … 48000.
   - Measured against Cauchy–Schwarz the ratio is ≈ 1.5–2.5/#ideals, i.e. full square-root cancellation.
   - The noise-free test "row scalar ζ against cubic-theta columns" is also at the random level up to Q = 96000. The column side has no hidden coherence with ζ.
5. **Type I sums are at the random level on both sides** [R-computation / N].
   - **μ side (exact).** S_d(R) = Σ_{d|x} β_x w(N x/L)(λx/R)₂ at L = Q^{0.79} has ⟨|S_d|²⟩_R / Σ_{d|x}|βw|² ∈ [0.71, 0.99] (mean 0.87–0.89) for every modulus d with N(d) ≤ 400 (up to 106 of them), at Q = 6000 … 48000. The max over R is about log #R.
   - **m side (noisy).** Rows restricted to d | R give a twist ratio of 0.98–1.01 for every N(d) ≤ 43 and every Q.
6. **Type II: π has no bilinear structure beyond a known class-function piece** [N].
   - **The test.** A twisted-multiplicative test (4-cycle statistic on M[a,b] = π(ab), for all six reciprocity twists j) finds exactly one component: a j = 4 twisted rank-one piece. Its size is Γ/E ≈ 0.04–0.07 (|corr| ≈ 0.2–0.25) for Q = 12000 … 48000, with no visible decay.
   - **It is ρ's class-function part.** π(j) = ρ(j − 4) identically, because ζ is exactly j = 4 twisted-multiplicative. So the piece is ρ's untwisted correlation transported by ζ.
   - **Removing it removes the signal.** At Q = 48000, removing the per-scale class means mod 36 of ρ (ρ's constant-term component, together with the estimator's deterministic small-Y component) reduces it to ≤ 0.013 in every band. It is not the conjugation symmetry V(m̄) = conj V(m): dropping conjugate prime pairs leaves it unchanged.
   - **Nothing else.** With that piece removed, |Γ/E| ≤ 0.03 at every j. So π(m₁m₂) is not ≈ f(m₁)g(m₂) under any reciprocity twist, beyond ζ × (ray class function).
7. **Still no power beyond degree 1; the conductor exponent of π's twists remains undetermined** [N/H].
   - **ρ shows a deficit; π does not.** At Q ≥ 24000 the ζ-free ρ has quadratic twists below random at small conductor: pooled Rat ≈ 0.75 for N(ν) < 512, in both windows. π shows no deficit (0.90–1.13; pooled 0.995 ± 0.027).
   - **ρ's deficit is most likely (at least partly) its degree-1 constant-term component [H; §3.1].** Multiplying by ζ turns that component into a random-twist object. At Q = 48000, removing class means mod 12 and then mod 36 moves ρ's small-ν ratios steadily toward 1 (0.74 → 0.79 → 0.87 at N(ν) ∈ [128, 512)).
   - **The decisive noise-free test.** ζ, γ₂ and γ₂ᾱ, sequences of GL(2)-metaplectic (cubic-theta) type with a Voronoi dual of length C N(ν)²/Q, show **no** small-conductor smallness up to Q = 3.84·10⁵ (N(ν)²/Q down to 3·10⁻⁴). Over the same range "1" falls to 0.05. So degree-2 duals have a large constant (C ≳ 10³–10⁴) and are invisible at accessible Q.
   - **Consequence.** The data cannot tell L²/Q from L³/Q duals for π, confirming the A2 verdict out to Q ≈ 10⁵.
8. **What it means** [H].
   - (V′_LS) is numerically true with the random exponent η ≈ 1 − log L/log Q. It holds uniformly for L ≤ Q, not only at L = Q^{0.79}.
   - If it were *provable* at that rate, the saving would be s₀ = (1 − M)/2. In the LP this is `theta_split_model` mode "random", with the P²/Y term dropped by (LoA6), and it gives **θ = 0.8519** (0.8548 with marks). Compare 0.874957 (Liu).
   - The numerics give no handle on provability. π looks like a generic sequence: random-level twists, no Type I or Type II structure to exploit, and no visible short dual.

---

## 1. Data generation (new) [R]

**Pipeline** (`vls/E_primes.py`, `vls/E_build.py`).
- **Primes.** All primary primes of O with N ≤ 2^22 (295,734 split, 156 inert) come from a vectorised lattice search rather than zomega's per-prime search. The output is identical to `zomega.primes_upto` on N ≤ 2^14 (all 1893 primes, including r and the sixth-root tables).
- **Gauss sums.** g̃(π) = g(π,1)/√Nπ by `fastgauss.gauss_split`, using one representative of each conjugate pair and g(π̄) = conj g(π). The relation is checked to 1.8·10⁻¹⁵ on 2000 primes; the run took 48 min wall time on 11 cores.
  - Agreement with the existing cache `gt_1048576.npy`: 2.3·10⁻¹⁵ on all 82,058 common primes.
- **Moduli s.** There are 1,183,821 squarefree s < 2^22, built by vectorised layering. g̃(s) uses twisted multiplicativity with pair-symbol tables, since every prime factor but the largest has N ≤ 2048.
- **Associates trick.** By Euler's criterion, (ζ₆/p)₆ = ζ₆^{(Np−1)/6}; this is verified on all primes up to 2^14. Hence
  χ_s(−ζ₆^c m₀) = ζ₆^{c t(s)} χ_s(−m₀), with t(s) = Σ_{p|s}(Np−1)/6 mod 6.
  So A_k is accumulated in six t-components on primary representatives m₀ only, and all six associates are recovered exactly. This makes the symbol table 6× smaller.
- **Validation.**
  - A direct, expA-style computation with zomega at Q = 600 and s < 2^14 agrees to 2·10⁻¹⁷.
  - At Q = 6000 the new data reproduce `D_numerics/data/A_Q6000_K12.npy` to ≤ 3.7·10⁻¹⁷ on all 12 common scales.
  - Script: `E_validate.py`.

**Data** (`vls/data/A_Q{Q}_E21.npy`, `m_Q{Q}_E21.npy`, 14 scales Y = 2^8 … 2^21):

| Q | rows (sf) | sf ideals | scales used (Y ≥ 4Q) | noise/signal power per scale (smallest → largest Y) | ‖ρ‖²_Ω |
|---|---|---|---|---|---|
| 6000 | 5052 (4728) | 788 | 2^15 … 2^21 (7) | 4.90 → 1.56 | 7.44 |
| 12000 | 10176 (9492) | 1582 | 2^16 … 2^21 (6) | 4.75 → 1.59 | 13.66 |
| 24000 | 20388 (19032) | 3172 | 2^17 … 2^21 (5) | 4.56 → 1.76 | 23.88 |
| 48000 | 40668 (37932) | 6322 | 2^18 … 2^21 (4) | 4.38 → 2.06 | 42.23 |
| 96000 | 81342 (75948) | 12658 | 2^19 … 2^21 (3) | 3.65 → 2.64 | 72.42 |

Notes on the table:
- ‖ρ‖²_Ω grows like Q^{0.82}, against the proved Q^{5/6} (D_derivation §6.3) [N, consistency check].
- **Per-row resolution.** The per-row signal is still below the noise at every scale, so individual ρ(m) cannot be resolved. Only aggregate statistics are meaningful (as in A2_numerics).
- **Limit on Q.** Y ≥ 4Q at Q = 96000 leaves only 3 scales. Going further needs s beyond 2^22, i.e. Gauss sums of norm up to 10⁷. That is a few hours of CPU, feasible but not done.

## 2. Estimators, noise and bias

### 2.1 Estimators [R as algebra]

Write V_k(m) = Y_k^{1/3} A_m(Y_k)/Ŵ = ρ(m) + w_k(m), and π_k = V_k ζ_m.

- **Twisted sums.** For each column ν and unit class c, set T_k(ν,c) = Σ_R Ω π_k(ζ₆^c R)(ν/R)₂. Then |F′_c(ν)|² is estimated by the inverse-noise-weighted distinct-scale product (`a2lib.cross_est`). This is unbiased if the w_k are mutually uncorrelated.
- **Random-model denominator.** The same estimator, row by row, gives den(ν) = Σ_m Ω²|π|²|(ν/m)₂|². That is E|F′|² for a random sequence with the same |π|.
- **Errors** come from three sources:
  1. leave-one-scale-out jackknife;
  2. random-sign nulls π → π·ε_m (6–8 draws). These keep |π| and the noise structure and destroy arithmetic structure, and E_ε[num] = den holds exactly;
  3. pairwise-scale estimates (§2.3).

### 2.2 Window dependence (the A2 referee's bias concern)

The cross-scale estimator is biased for *periodic* test functions while Y ≲ 3Q, because the deterministic k = 1 Poisson term survives.

**ρ is strongly affected; π is not** [N], per `E_pairs.py`:
- **ρ.** At Q = 96000 the untwisted per-scale unit-class sums of ρ, in units of √(random level), are:
  - −19 … −14 at Y = 2.7Q;
  - +5 … +6 at 5.5Q;
  - +19 … +26 at 11Q;
  - +27 … +31 at 22Q.

  So ρ's constant-term coherence is real, but its per-scale estimates have still not reached a plateau at Y = 22Q.
- **π.** The same sums for π are |·| ≲ 3.7 at every scale: random, with no trend.

**The (V′_LS) statistics for π are window-independent within errors:**

R_rand (± jackknife; random-sign null mean in brackets):

| Q | window (#scales) | L = Q^{0.5} | Q^{0.7} | **Q^{0.79}** | Q |
|---|---|---|---|---|---|
| 6000 | Y ≥ 1.5Q (8) | 1.13 ± 0.07 [1.05] | 1.02 ± 0.04 | **1.015 ± 0.048** | 0.992 |
| 6000 | Y ≥ 4Q (7) | 1.11 ± 0.08 [1.05] | 1.01 ± 0.03 | **1.009 ± 0.046** | 0.992 |
| 6000 | Y ≥ 16Q (5) | 1.17 ± 0.10 [1.04] | 1.02 ± 0.05 | **1.024 ± 0.061** | 0.989 |
| 6000 | K12 data, Y ≥ 1.3Q (7) | 1.24 ± 0.19 [1.03] | 1.06 ± 0.03 | **1.070 ± 0.057** | 1.000 |
| 12000 | Y ≥ 1.5Q (7) | 0.98 ± 0.11 [1.03] | 0.95 ± 0.05 | **0.979 ± 0.023** | 0.994 |
| 12000 | Y ≥ 4Q (6) | 0.94 ± 0.11 [1.03] | 0.95 ± 0.05 | **0.971 ± 0.024** | 0.995 |
| 12000 | Y ≥ 16Q (4) | 1.04 ± 0.07 [1.02] | 0.99 ± 0.07 | **0.997 ± 0.025** | 0.996 |
| 24000 | Y ≥ 4Q (5) | 0.90 ± 0.04 [1.04] | 1.02 ± 0.04 | **1.011 ± 0.035** | 1.001 |
| 24000 | Y ≥ 16Q (3) | 0.87 ± 0.14 [1.03] | 1.02 ± 0.04 | **1.022 ± 0.050** | 1.003 |
| 96000 | Y ≥ 2.5Q (4) | 1.25 ± 0.04 [0.97] | 1.07 ± 0.02 | **1.039 ± 0.011** | 1.004 |
| 96000 | Y ≥ 5.5Q (3) | 1.22 ± 0.08 [0.92] | 1.07 ± 0.04 | **1.039 ± 0.020** | 1.004 |

There is no window dependence beyond the noise. In particular the Y ≲ 3Q bias that the referee found for ρ's class sums does not appear in π's twisted statistics.

### 2.3 Systematic scatter of the estimator [N]

**The scatter.** The row-level off-diagonal scale Gram Σ_m Ω π_k conj(π_l)/Σ Ω|π|², which should be 1 for each pair k ≠ l, scatters between 0.70 and 1.35 at every Q.
- If the noises were independent across rows and scales, the scatter would be ≈ 0.02–0.07.
- The excess has two sources:
  - **ρ's deterministic class-function component D_k(m).** Since |ζ| = 1, π's row-level products inherit it.
  - **Associate and conjugation structure.** V(m̄) = conj V(m), and the six associates are DFT combinations of six components.

**Why the large-L ratio is still reliable.** The twisted numerator and the denominator see the *same* pair weights. Averaged over many ν, orthogonality of the ψ_ν makes num/den self-normalising. That is why R_rand(L) is pinned to 1.00 ± 0.01 at L ≥ Q^{0.9}.

**At small L it is not.** Pairwise-scale estimates, each normalised by its own value at L = Q^{0.8}, give:

| Q | L = Q^{0.5} | Q^{0.6} | Q^{0.7} | Q^{0.79} |
|---|---|---|---|---|
| 48000, π | 0.94 ± 0.11 | 0.97 ± 0.05 | 0.98 ± 0.04 | 0.989 ± 0.018 |
| 96000, π | 1.21 ± 0.05 | 1.08 ± 0.05 | 1.03 ± 0.03 | 1.003 ± 0.006 |

(± is the spread over pairs.)

**Is the Q = 96000 excess at small L real?** The +20% at Q^{0.5} is shared by all pairs (89 columns × 6 classes), but it does not replicate:
- the per-ν ratios at N(ν) ≤ 512 are uncorrelated with those at Q = 48000 (corr −0.02);
- the small-ν mean wanders 0.91, 1.03, 0.92, 1.13 over Q = 12k … 96k, an empirical s.d. of about 0.1.

I therefore treat ±20% as the honest uncertainty of R_rand at L ≤ Q^{0.6}, and ±5% at L ≥ Q^{0.7} [N].

## 3. (V′_LS): the large-sieve ratio R(L,Q)

**Statistic** (`E_ls.py`). The columns ν are all squarefree primary ν with (ν,6) = 1 and 7 ≤ N(ν) ≤ Q: 1695 / 3390 / 6766 / 13536 / 27095 columns. The full twist matrix (ν/R)₂ = Π_{p|R} χ_p(ν)³ is exact.

Unit classes c (f = 3) are kept separate, which absorbs κ₀ and the ray factor. Note that |Σ_c κ(c)F′_c|² ≤ 6 Σ_c|F′_c|².

Definitions:
- S(L) = Σ_{N(ν)≤L} Σ_c |F′_c(ν)|²;
- R_LS = S(L)/(Q Σ Ω|π|²), the task normalisation;
- R_rand = S(L)/Σ_{N(ν)≤L} den(ν).

**Cumulative ratios** (Y ≥ 4Q; ± is the jackknife; "null" is the random-sign mean (s.d.)):

| L | Q = 6000 R_LS / R_rand | 12000 | 24000 | 48000 | 96000 |
|---|---|---|---|---|---|
| Q^{0.50} | 0.00268 / 1.11 ± 0.08 (null 1.05) | 0.00182 / 0.94 ± 0.11 | 0.00115 / 0.90 ± 0.04 | 0.00085 / 0.92 ± 0.08 | 0.00081 / 1.22 ± 0.08 (null 0.92) |
| Q^{0.60} | 0.00664 / 1.06 ± 0.04 | 0.00437 / 0.90 ± 0.07 | 0.00357 / 1.03 ± 0.04 | 0.00263 / 0.96 ± 0.03 | 0.00225 / 1.09 ± 0.02 |
| Q^{0.70} | 0.0152 / 1.01 ± 0.03 | 0.0115 / 0.95 ± 0.05 | 0.0101 / 1.02 ± 0.04 | 0.0079 / 0.98 ± 0.03 | 0.0070 / 1.07 ± 0.04 |
| **Q^{0.79}** | **0.0330 / 1.009 ± 0.046** | **0.0278 / 0.971 ± 0.024** | **0.0247 / 1.011 ± 0.035** | **0.0214 / 1.007 ± 0.029** | **0.0190 / 1.039 ± 0.020** |
| Q^{0.90} | 0.0845 / 0.99 ± 0.02 | 0.0773 / 0.96 ± 0.01 | 0.0741 / 1.00 ± 0.01 | 0.0700 / 1.01 ± 0.01 | 0.0646 / 1.00 ± 0.02 |
| Q | 0.2034 / 0.992 ± 0.011 | 0.2047 / 0.995 ± 0.002 | 0.2024 / 1.001 ± 0.003 | 0.2051 / 1.008 ± 0.006 | 0.2043 / 1.004 ± 0.009 |

What the table shows:
- **The random-model expectation is R_LS ≈ 0.72·#ν(L)/Q.** The constant 0.72 is ΣΩ²|π|²/ΣΩ|π|² times the coprimality fraction. The measured R_LS·Q/#ν(L) is 0.64–0.80 at every L and Q (`E_summary.py`).
- **Saturation would mean R_LS ≍ 1** (the large-sieve Q-term, i.e. Cauchy–Schwarz). At L = Q^{0.79} the measured R_LS is 0.019–0.033, so saturation is excluded by a factor of 30–50.
- **Scaling in Q.** log R_LS(Q^{0.79}) has slope −0.208 in log Q (from 6000 to 48000) and −0.20 (to 96000). The random prediction is −(1 − 0.79) = −0.21.
- **Fits** (`E_summary.py`):
  - **Cumulative fit.** log R_LS = a + b log(L/Q) + c log Q gives b = 1.008 ± 0.007 and c = 0.002 ± 0.015 (45 nested points, jackknife weights).
  - **Independent dyadic bins** (60 bins over 5 Q). log Rat = −γ log(N(ν)/Q) + const gives γ = 0.000 ± 0.004 (χ²/dof = 0.91).
  - **Bins with L < Q^{0.79} only** (42 bins): γ = 0.037 ± 0.010. This is driven by the Q = 96000 small-ν bins and lies within the §2.3 systematics.
  - **Conclusion.** b = 1 − γ ∈ [0.94, 1.01]. "In between" behaviour R ~ (L/Q)^b with b ≤ 0.9 is excluded [N].
- **Extremes.** max_ν Rat is 3.7 / 4.3 / 4.8 / 5.9 / 6.8, against a null mean maximum of 3.7 / 3.8 / 4.4 / 4.9 / 6.2.
  - The fraction of columns with Rat > 5 is 10⁻⁴ for π and the same for the null.
  - There is no column with Rat > 10 at any Q.
- **Calibrations** (exact, noise-free, same rows and columns; `E_ls.py` "exact" columns):
  - **Random phases:** 1.00 in every bin.
  - **Degree 1:** the constant sequence "1" is structurally small at small conductor, and increasingly so with Q. At N(ν) ∈ [7, 64) its Rat is 0.10–0.63 at Q = 12000, 0.01–0.08 at Q = 48000, and 0.01–0.12 at Q = 96000. ᾱ² gives 0.06–0.35 at Q = 96000, N(ν) < 64.
  - **ζ alone:** about 1. (A scaling test of ζ alone is in §3.1.)
  - So the test **does** see short duals of degree-1 objects at these Q. π is not degree-1-like at any Q.

### 3.1 Small conductor: ρ against π, and the exact ζ scaling test [N/H]

Pooled over the dyadic bins with N(ν) < 512 (`E_summary.py`):

| Q | π (Y ≥ 4Q) | ρ, ζ-free (Y ≥ 4Q) | ρ (Y ≥ 16Q) |
|---|---|---|---|
| 6000 | 1.00 ± 0.05 | — | 0.90 ± 0.06 |
| 12000 | 0.90 ± 0.07 | — | 1.05 ± 0.04 |
| 24000 | 1.02 ± 0.07 | **0.75 ± 0.06** | **0.68 ± 0.07** |
| 48000 | 0.90 ± 0.07 | **0.75 ± 0.03** | — |
| 96000 | 1.13 ± 0.06 | 0.90 ± 0.13 | — |

**Observations.**
- From Q = 24000 on, ρ's quadratic twists at small conductor are 25–30% below the random level. This holds in both windows, so it is not the Y ≲ 3Q bias.
- π's twists are not: their pooled value over all Q is 0.995 ± 0.027.
- The A2 note found no difference at Q ≤ 6000. At Q ≥ 24000 the test separates ρ from π statistically. As explained below, this is **not** evidence about the conductor of π.

**Exact scaling test for cubic-Gauss-sum sequences** (`E_zeta_twist.py`, noise-free). Classes are R mod 4; Rat is the bin mean, for Q = 6000 / 24000 / 96000 / 384000:

| N(ν) | N(ν)²/Q at Q = 3.84·10⁵ | ζ = ψ₁ γ̄₂⁴ | γ̄₂ ᾱ | γ₂ | "1" (degree 1) |
|---|---|---|---|---|---|
| [7, 16) | 3·10⁻⁴ | 0.73 / 0.97 / 0.91 / 1.29 | 0.64 / 1.40 / 0.63 / 0.97 | 1.22 / 1.00 / 1.38 / 0.91 | 0.47 / 0.40 / 0.06 / **0.05** |
| [16, 32) | 1.3·10⁻³ | 0.96 / 1.16 / 0.79 / 0.92 | 0.97 / 0.74 / 1.20 / 1.26 | 0.89 / 1.12 / 0.76 / 0.86 | 0.90 / 0.45 / 0.23 / **0.05** |
| [32, 64) | 5·10⁻³ | 0.79 / 1.24 / 0.76 / 0.90 | 0.80 / 0.88 / 0.95 / 0.82 | 0.90 / 0.82 / 0.98 / 1.03 | 1.12 / 0.75 / 0.30 / **0.06** |
| [64, 128) | 0.02 | 0.95 / 1.03 / 1.10 / 1.13 | 0.91 / 1.04 / 1.09 / 1.10 | 0.84 / 0.99 / 0.82 / 0.92 | 0.91 / 0.79 / 0.61 / **0.15** |
| [128, 256) | 0.085 | 1.20 / 1.08 / 0.98 / 0.93 | 1.05 / 1.02 / 0.86 / 1.02 | 0.95 / 0.86 / 1.04 / 0.92 | 0.93 / 0.82 / 0.73 / 0.32 |
| [256, 512) | 0.34 | 0.98 / 0.99 / 1.09 / 1.01 | 0.98 / 0.93 / 0.98 / 1.07 | 1.01 / 1.05 / 1.00 / 0.95 | 1.00 / 0.94 / 0.85 / 0.48 |

What the table shows:
- The degree-1 calibration becomes structurally small as Q grows, as expected for dual length ∝ N(ν)/Q.
- The cubic-Gauss-sum sequences are GL(2)-metaplectic; θ₃ coefficients at squarefree arguments have a Voronoi formula with dual length C N(ν)²/Q. They stay at the random level everywhere, even at N(ν)²/Q = 3·10⁻⁴.
- So C ≳ 10³–10⁴, in line with the cubic numerics of V_derivation §8.2. **Degree-2 duals are invisible at every accessible Q** [N].

**Interpretation of ρ's deficit** [H; not settled].

There are three candidate explanations.

- **(a) ρ's degree-1 part.** ρ has a class-function (constant-term) component carrying ≈ 20% of ‖ρ‖² in excess of random (A2 §3.1, mod 36). Being degree-1, its quadratic twists vanish like those of "1", while it still counts in the denominator. That would give Rat_ρ → 1 − φ_h ≈ 0.75–0.8 once "1"-type twists are small, i.e. at Q ≳ 2·10⁴ by the table above. This is exactly where the deficit appears. In π = ρζ the component becomes ζ·h, whose twists are at the random level, so π shows no deficit.

  *Direct check* (`E_rho_classproj.py`, Y ≥ 4Q). Removing the per-scale class means of ρ gives:

  | Q | N(ν) | none | mod 12 | mod 36 |
  |---|---|---|---|---|
  | 24000 | [7, 128) | 0.62 | 0.55 | 0.65 |
  | 24000 | [128, 512) | 0.86 | 0.91 | 0.93 |
  | 48000 | [7, 128) | 0.68 | 0.73 | 0.78 |
  | 48000 | [128, 512) | 0.74 | 0.79 | 0.87 |
  | 48000 | [512, 1024) | 0.92 | 0.97 | 1.09 |

  π's ratios stay at 0.91–1.10 throughout. At Q = 48000, where there are about 60 rows per class mod 36, the class projection moves ρ steadily toward 1. This supports (a) at least in part. The per-scale class means are themselves contaminated by the deterministic small-Y component (§2.2), so the check cannot be sharper.
- **(b) A genuine θ₆ Voronoi effect.** This would require a much smaller constant C than in the cubic case.
- **(c) Estimator systematics specific to ρ.** ρ, unlike π, carries the deterministic periodic component. However, the deficit is the same in the Y ≥ 4Q and Y ≥ 16Q windows (0.75 against 0.68 at Q = 24000).

**Either way, it does not touch π.** π has no deficit at any Q. The exact ζ test above shows that degree-2 (cubic-theta-type) duals are invisible at these sizes. So the ρ-versus-π difference cannot be read as evidence about the conductor of π's twists.

The A2 verdict stands, now out to Q ≈ 10⁵ for π and 3.8·10⁵ for ζ: **no power beyond degree 1; L²/Q versus L³/Q is undeterminable numerically.**

## 4. (V′_β): the true bilinear form with Patterson's cubic theta coefficients

**Columns** (`E_beta.py`, from `cubic_theta.py` = paper78 eq. reflection-infinity-coefficient):
- x = n b³ primary and odd, a(x) = |b| g̃₃(n), μ = λ^{−3}x, λ⁴μ = λx;
- β_x = a(x) α(x)/√N(x);
- the twist (λx/R)₂ is exact;
- dyadic windows w_L(x) = bump(N(x)/L, ½, 1).

The paper's additive character ϑ(λ⁴μ) is emulated by 8 additive characters e(hx/9). The fixed ray factor is absorbed by the 72 class functions (unit class c) × (R mod 4).

B⁰_R = Σ_x β_x w_L(x)(λx/R)₂ is exact. Only π is noisy.

| Q | windows L (from √Q to Q, plus Q^{0.79}) | π: Σ|T|²/random, range over windows | typical jackknife s.e. | π: Σ|T|²/CS | 1/#ideals | ⟨|B⁰_R|²⟩/Σ|βw|² |
|---|---|---|---|---|---|---|
| 6000 | 93 … 6000 (8 windows) | 0.92 – 1.10 | 0.08 – 0.25 | 2.2–2.7·10⁻³ | 1.27·10⁻³ | 0.88 – 0.99 |
| 12000 | 187 … 12000 (8) | 0.77 – 1.18 | 0.12 – 0.31 | 0.9–1.4·10⁻³ | 6.3·10⁻⁴ | 0.83 – 0.98 |
| 24000 | 187 … 24000 (9) | 0.89 – 1.38 | 0.12 – 0.32 | 5.3–8.1·10⁻⁴ | 3.2·10⁻⁴ | 0.89 – 0.99 |
| 48000 | 375 … 48000 (9) | 0.76 – 1.07 | 0.20 – 0.35 | 2.2–3.1·10⁻⁴ | 1.6·10⁻⁴ | 0.93 – 0.97 |

Reading the table [N]:
- **Random level everywhere.** The bilinear form is at the random level in every window, and at L = Q^{0.79}: 1.10 / 1.18 / 0.89 / 0.92.
- **Square-root cancellation.** Measured against Cauchy–Schwarz it is ≈ 2/#ideals, i.e. |T| ≈ ‖π‖‖B⁰‖·#rows^{−1/2}. In Q-units, |T| ≈ Q^{5/12}·O(1), not Q^{11/12}.
- **ζ-free control.** The same form with ρ in place of π (the sum that (V) treats) is also at the random level.
- **Column side.** ⟨|B⁰_R|²⟩_R ≈ Σ_x|β_x w|², so the reflected column sum is itself a random-level sum for typical R. Its m-mean square is what large-sieve duality predicts (#R ≫ #x).

**Noise-free check: does the row scalar alone correlate with the columns?** (`E_zeta_beta.py`, exact, Q up to 96000.)

Σ_{cl}|Σ_R Ω ζ_R B⁰_R|²/Σ Ω²|B⁰|², over classes R mod 36:

| Q | L = Q/8 | Q/4 | Q/2 | Q |
|---|---|---|---|---|
| 6000 | 0.90 | 0.95 | 0.96 | 1.05 |
| 24000 | 1.01 | 0.99 | 0.84 | 0.86 |
| 96000 | 1.28 | 1.01 | 1.00 | 0.87 |

- The same holds for conj(γ₂), γ₂, conj(γ₂)ᾱ and random phases.
- There is **no** main term in the ζ × cubic-theta bilinear form.
- An earlier F = 1 reading of 2.6–3.9 at Q = 6000 was a single-sample fluctuation: it is 0.03–4.4 at larger Q and ≈ 1 once averaged over classes.

## 5. Type I sums

**(a) μ side, exact** (`E_beta.py`). Take L = Q^{0.79} and the 108 squarefree primary d with (d,6) = 1 and N(d) ≤ 400. Define

  S_d(R) = Σ_{x: d|x} β_x w_L(x)(λx/R)₂.

Only d with at least 2 columns in the window are used.

| Q | L = Q^{0.79} | #d | ⟨|S_d|²⟩_R / Σ_{d|x}|βw|²: mean (s.d.) [min, max] | max_R |S_d|²/Σ_{d|x}|βw|², for N(d) ≤ 43 (log #R) | untwisted |Σ_{d|x}βw|²/Σ_{d|x}|βw|²: median [max] | untwisted, all x |
|---|---|---|---|---|---|---|
| 6000 | 965 | 31 | 0.887 (0.046) [0.73, 0.96] | 1.2–6.1 (6.7) | 0.98 [2.0] | 0.11 |
| 12000 | 1669 | 52 | 0.891 (0.049) [0.74, 0.97] | 1.4–7.5 (7.4) | 0.90 [2.0] | 0.05 |
| 24000 | 2886 | 99 | 0.872 (0.071) [0.71, 0.99] | 1.5–6.6 (8.1) | 0.98 [2.2] | 0.35 |
| 48000 | 4990 | 106 | 0.888 (0.067) [0.74, 0.98] | 4.4–9.2 (8.8) | 0.71 [2.1] | 0.20 |

What this shows [R-computation / N]:
- **Random level for every modulus.** Restricting the cubic-theta columns to an arithmetic progression d | x does not create coherence with the quadratic characters (λx/R)₂ for any d. The R-mean square is slightly below the diagonal (0.8–0.95) because of the coprimality zeros.
- **No main term to subtract.** The untwisted progression sums have no theta main term. Heuristically [H], the angle character α(x) in β kills the θ₃ pole. There is nothing to subtract, and also no Type I main term that a decomposition could exploit.

**(b) m side** (`E_ls.py`, noisy). Rows are restricted to d | R, for the 12 smallest primes d (N(d) ≤ 43), with the same quadratic twists.

The summed ratio over all ν is 0.98–1.01 at every Q (for example 0.995–1.009 at Q = 48000 with 143–789 ideals per d). The median per-ν ratio is 0.87–0.96, the usual exponential-type skew. So π restricted to progressions is at the random level [N].

## 6. Type II: bilinear structure of π

**Test** (`E_typeII.py`). Take squarefree ideals R = a·b, with a prime and N(a) in a band, and b = R/a. Set M_c[a,b] = π(ζ₆^c a b).

Twisted multiplicativity

  π(ab) = f(a) g(b) ζ₆^{j(e(a,b)+e(b,a))},   e(a,b) = exponent of χ_a(b),

holds exactly when M̃ = M·ζ₆^{−j(e+e′)} has rank one.

Two statistics:
- **4-cycle statistic.** Γ_j/E_sym = Σ_{a≠a′, b≠b′} M̃[a,b]M̃̄[a,b′]M̃̄[a′,b]M̃[a′,b′] / Σ(P P + P P)/2. It is 1 for an exactly rank-one sequence, 0 for a random one, and ≈ φ² when a rank-one component carries a fraction φ of the energy. Each of the four entries uses a distinct scale, averaged over the 24 orders.
- **Spectral statistic.** The top eigenvalue of the unbiased cross-scale Gram M̃^H M̃, divided by its trace.

**Calibration [R].** ζ alone gives exactly 1 at j = 4 and ≤ 0.23 elsewhere. γ₁ gives 1 at j = 1. Random phases give |Γ/E| ≤ 0.007.

**Results for π** [N]. Γ/E at j = 4 (± random-sign null s.d.), for bands N(a) ∈ [7,60) / [60,400) / [400,3000) / [3000,20000):

| Q (data, window) | j = 4 (= ρ at j = 0) | largest other j |
|---|---|---|
| 1500 (K12, Y ≥ 2.7Q) | 0.16 ± 0.04 / 0.16 ± 0.07 / — / — | ≤ 0.10 (null s.d. 0.04–0.09) |
| 3000 (K12) | 0.086 ± 0.019 / 0.106 ± 0.029 / — / — | ≤ 0.08 |
| 6000 (K12, Y ≥ 1.3Q) | 0.083 ± 0.008 / 0.169 ± 0.022 / 0.045 ± 0.044 / — | ≤ 0.07 |
| 12000 (E, Y ≥ 4Q) | 0.052 ± 0.006 / 0.031 ± 0.012 / 0.070 ± 0.012 / — | 0.025 ± 0.006 (j = 0) |
| 24000 (E) | 0.044 ± 0.003 / 0.032 ± 0.005 / 0.061 ± 0.004 / — | 0.025 ± 0.004 (j = 1) |
| 48000 (E) | 0.060 ± 0.001 / 0.047 ± 0.001 / 0.047 ± 0.003 / 0.073 ± 0.004 | 0.017 ± 0.002 (j = 1) |
| 24000, conjugate pairs {a, ā} dropped | 0.070 ± 0.011 / 0.036 ± 0.012 / 0.051 ± 0.011 / — | ≤ 0.015 |
| 48000, conjugate pairs dropped | 0.060 ± 0.003 / 0.047 ± 0.004 / 0.052 ± 0.005 / 0.074 ± 0.012 | ≤ 0.010 |
| 24000, ρ's class means mod 36 removed | 0.025 ± 0.004 / 0.015 ± 0.004 / 0.048 ± 0.006 / — | 0.033 ± 0.007 (j = 0) |
| **48000, ρ's class means mod 36 removed** | **0.012 ± 0.003 / 0.004 ± 0.004 / 0.013 ± 0.007 / −0.011 ± 0.008** | 0.022 ± 0.002 (j = 0, band 1); 0.032 ± 0.010 (j = 0, band 4) |

**Spectral test, j = 4, band [7,60)** (top eigenvalue/trace):

| Q | π | null |
|---|---|---|
| 12000 | 0.085 | 0.058 ± 0.001 |
| 24000 | 0.089 | 0.060 ± 0.001 |
| 48000 | 0.081 | 0.059 ± 0.000 |

Other j sit within 0.01–0.02 of the null.

**Interpretation.**
- **Where the component comes from [R].** π's only factorisable component is the j = 4 one, and π(j) = ρ(j − 4) holds identically in every table, because ζ is exactly j = 4 twisted-multiplicative. So the component is ρ's *untwisted* (j = 0) correlation transported by ζ.
- **It is not a symmetry artifact [N].** The exact symmetry V_k(m̄) = conj V_k(m) makes the cycles (a, ā; b, b̄) contribute |π(ab)|²|π(ab̄)|² ≥ 0. Dropping one prime of each conjugate pair leaves the signal unchanged, so these cycles are not its source.
- **It is ρ's class-function part [N].** Removing the per-scale class means mod 36 of ρ before multiplying by ζ removes it at Q = 48000 (about 60 rows per class): 0.047–0.073 falls to ≤ 0.013. At Q = 24000 (about 30 rows per class) it is halved. The class means include ρ's genuine constant-term (Eisenstein-type) component and the estimator's deterministic small-Y component (§2.2); both are ray-class functions. Since a ray class function is a sum of multiplicative characters, its product with ζ is exactly of the j = 4 twisted-multiplicative type.
- **Earlier explanations superseded.** At Q = 6000 the K12 data use small-Y windows, and class projections there are too noisy. The decay 0.16 → 0.044 seen from Q = 1500 to 24000 does not continue to 48000 (0.047–0.073). The earlier description "concentrated on conjugate pairs" (`rowcorr` at Q = 6000) reflected only where the b-ranges overlap.
- **What is left.** Once ζ × (class function) is removed, π shows no twisted-multiplicative structure at any j: |Γ/E| ≤ 0.03, and mostly ≤ 0.013.

**For the other attacks [H].**
- A Vaughan or Heath-Brown-identity style decomposition of the m-sum has nothing to grip. Apart from ζ × (ray class function), π(m₁m₂) is not close to f(m₁)g(m₂) under any of the six reciprocity twists.
- The exactly factorisable objects are ζ itself (a cubic Gauss sum times an angle and ray character) and ρ's constant-term component.
- Splitting π = ζ·(ρ's class-function part) + π′ is natural.
  - **The first piece.** It is a cubic Gauss sum times a ray character. Its quadratic twists are at the random level up to Q = 3.8·10⁵ (§3.1), and they are amenable to cubic-theta (GL(2)-metaplectic) technology.
  - **π′.** It carries all the remaining difficulty and has no detectable structure.

## 7. Exact ζ twist scaling and the LP value at stake

**ζ alone at small conductor:** see §3.1. There is no visible degree-2 dual up to Q = 3.84·10⁵.

**LP** (`theta_split_model`/`lp_const_saving` machinery, with P²/Y dropped by (LoA6)).
- (V′_LS) at the *measured* rate, η = 1 − log L/log Q for all L ≤ Q, is the saving s(M′) = (1 − M′)/2. This is the model's "random" mode.
- It would give **θ = 0.851918** at (l_x, l_y, ℓ) = (0.1978, 0.7158, 0.0863), M = 0.914, or 0.854755 with marks ("randommarked").
- The "LSmarked" mode reproduces the old 0.86624.
- These are **conditional** values [H]. The numerics say the hypothesis behind them is consistent with the data up to Q ≈ 10⁵. They say nothing about whether it can be proved.

## 8. What the data can and cannot say

| question | answer | strength |
|---|---|---|
| Is the data pipeline exact? | yes: reproduces K12 to 4·10⁻¹⁷; Gauss sums to 2·10⁻¹⁵; ζ to 1.6·10⁻¹⁴ | [R] |
| R(L,Q) ~ L/Q (random), ~1 (saturation), or in between? | **random**: R_LS = (0.64–0.80)·#ν(L)/Q; b = 1 − γ ∈ [0.94, 1.01]; saturation excluded by ×30–50 at L = Q^{0.79} | [N], strong, Q ≤ 9.6·10⁴ |
| Any large-sieve-extremal ν (coherent twist)? | none; maxima at null level; no cross-Q reproducibility | [N], strong |
| (V′_β) with Patterson's β? | random level, ≈ #rows^{−1/2} of Cauchy–Schwarz | [N] (noisy: s.e. 0.1–0.35 per window) |
| ζ × cubic-theta columns coherent? | no (exact, Q ≤ 9.6·10⁴) | [N], noise-free |
| Type I (μ side, exact; m side) | random level for every d with N(d) ≤ 400 (μ) and ≤ 43 (m) | [N] |
| Type II / factorisation of π | only ζ × (ρ's class-function / constant-term part), at j = 4: Γ/E ≈ 0.05 (|corr| ≈ 0.22), and ≤ 0.013 once class means mod 36 are removed; nothing at other j | [N] |
| Conductor of π's quadratic twists (L²/Q vs L³/Q) | undeterminable: even exact cubic-Gauss-sum (degree-2) sequences show no small-conductor smallness up to Q = 3.8·10⁵, while degree 1 does. ρ's deficit at Q ≥ 24000 is most likely, at least partly, its degree-1 constant-term component | [N] + [H] |
| Is (V′_LS) *provable*? | the data cannot say | — |

**Gaps** [G].
1. **Ray factor and κ₀.** The fixed ray factor and κ₀ are handled only through per-class sums: unit classes, and R mod 4 for (V′_β). A κ_F with conductor beyond 36·4 is not tested. The same limitation applies in A2_numerics.
2. **Non-squarefree rows** are excluded (6.6%), because their ζ is not implemented.
   - **Columns** are restricted to squarefree ν with (ν,6) = 1. The full family ν = λ⁴μ ∈ O∖0 adds ν₀b² (the same character up to a principal factor, summable over b) and finitely many unit, λ and 2-power ray twists. That changes S(L) by a bounded factor [R-].
3. **Coupled weight.** The reflected weight V♯(q_μZ/q_{c(m)}²) couples μ and m. It is modelled by separable dyadic windows, which cost Z^ε via Mellin, but this is not tested.
4. **Additive character.** ϑ(λ⁴μ) is emulated by additive characters mod 9. The true modulus depends on the level L.
5. **Estimator systematics.** ±20% at L ≤ Q^{0.6} (§2.3), caused by ρ's deterministic small-Y component, which affects the row-level normalisation.
6. **Scale coverage.** At Q = 96000 only 3 scales have Y ≥ 4Q.
7. **Range of Q.** All statements are for Q ≤ 9.6·10⁴. A conductor effect of size N(ν)^k/Q with a large constant C (C ≳ 10⁴ is suggested by the cubic numerics) would only show at much larger Q.
8. **What random-level numerics prove.** Random-level numerics are what *any* generic sequence produces. In particular they are blind to the zero-detecting Möbius signal 𝒮 of `V_gaps.md` §1.6. That signal is absent if the relevant zeros are absent, and then any true-but-unprovable (V′) looks exactly like this.

## 9. Files

| purpose | file |
|---|---|
| primes and Gauss sums to 2^22 | `vls/E_primes.py` → `data/primes_X4194304.npz` (log `E_primes_4194304.log`) |
| sextic additive polynomial at large Q | `vls/E_build.py` → `data/A_Q{Q}_E21.npy`, `m_Q{Q}_E21.npy`, `At_Q{Q}_E21.npy` (t-components), `slist_X4194304.npz` |
| validation | `vls/E_validate.py` (against K12), `vls/E_validate_small.py` (direct) |
| library | `vls/e_lib.py` (RowsE with sparse incidence and Lemma Z ζ, quadratic twist matrix, columns, layout) |
| (V′_LS) | `vls/E_ls.py` → `data/ls_Q{Q}_E_y{ymin}.{log,json,npz}`; fits `vls/E_summary.py`; pairwise-scale diagnostic `vls/E_pairs.py` → `data/pairs_Q{48000,96000}*.{log,json}` |
| (V′_β) and Type I (μ side) | `vls/E_beta.py` → `data/beta_Q{Q}_E.{log,json}` |
| exact ζ × cubic theta | `vls/E_zeta_beta.py` → `data/zeta_beta.json` |
| Type II | `vls/E_typeII.py` (arguments: Q src ymin [F_proj] [noconj]) → `data/typeII_Q{Q}_E.log`, `typeII_Q{Q}_E_p36.log` (class means mod 36 removed), `typeII_Q{Q}_E_noconj.log`, `typeII_Q{1500,3000,6000}_K12*.json` |
| exact ζ twist scaling | `vls/E_zeta_twist.py` → `data/zeta_twist.{log,json}` |
| ρ's small-conductor deficit and class projection | `vls/E_rho_classproj.py` → `data/rho_classproj_Q{24000,48000}.log` |
| LP value at the random rate | `theta_split_model.py` modes "random" and "randommarked" through `wf/lp_const_saving.py`'s `th()` with P²/Y dropped (`vls/E_lp_random.py`) |
| optional faster Gauss sum | `vls/fastgauss2.py` (chunked, same output, about 2× faster under load) |

---

## Referee report

*Independent skeptical review, 9 Oct 2026.*

*Sources checked:*
- *this note, plus `e_lib.py`, `E_ls.py`, `E_build.py`, `E_beta.py`, `E_zeta_twist.py`, `E_summary.py` and `E_lp_random.py`, with the logs in `vls/data/`;*
- *`wf2/a2lib.py`, `wf2/A2_numerics.md` and its referee report;*
- *`TARGET_VLS.md`, `wf/V_gaps.md` §1.6–1.7 and `wf/AUDIT_chain.md` §3–4;*
- *`paper78.tex` ll. 3340–3410 (b_*, ξ, c(s)).*

*New scripts, all in `vls/`:*
- *`referee_E_twist.py`: independent Jacobi-symbol check of the twist matrices;*
- *`referee_E_zeta_mu.py` and `referee_E_gamma_cube.py`: the Möbius factor in ζ;*
- *`referee_E_rho_beta.py`: the ζ-free control in (V′_β);*
- *`referee_E_lp_floor.py`: LP at large constant saving.*

### Verdict: NOT REFUTED

Nothing is claimed to be proved, and nothing is. The numerical core reproduces exactly.

Three supporting statements are wrong:
- ζ is called a GL(2)-metaplectic sequence (§2 below);
- the ζ-free control in (V′_β) is said to be at the random level (§3);
- Gap 2 is labelled [R-] (§4).

None of them changes any statement about π. Two of them sharpen the "provability" caveat.

### 1. What I verified

- **Reproducibility [R as computation].** I re-ran `E_ls.py` at Q = 12000, Y ≥ 4Q, with output redirected to scratch.
  - Every π statistic is identical to `ls_Q12000_E_y4.log`: R_rand(Q^{0.79}) = 0.971 ± 0.024 and R_LS = 0.02779, per-bin Rat, extremes and Type I in m.
  - Only the random-sign null differs, because it uses a different number of draws.
- **Twist matrices [R].** For rows R = a + bω of squarefree rational norm n, O/R ≅ Z/n via ω ↦ −a·b⁻¹ mod n, so (ν/R)₂ is a Jacobi symbol. At Q = 12000:
  - `quad_twist`: 0 mismatches in 80,000 (row, ν) pairs;
  - the β-columns (λx/R)₂: 0 mismatches in 60,000 pairs.
- **ζ_m [R].**
  - |ζ| = 1 to 5·10⁻¹⁵.
  - ζ is identical on the six associates, so the unit classes absorb only m-side unit factors.
  - The product-formula check in `RowsE._zeta_fast` is the same as in a2lib.
- **LP [R as computation].** `E_lp_random.py` reproduces 0.851918 / 0.854755 / 0.866240.

### 2. Error: ζ is not of GL(2)-metaplectic (cubic-theta) type. It carries μ(m)

**The claim.** §0.7, §3.1 (the "decisive noise-free test") and §6 list ζ = ψ₁γ̄₂⁴ together with γ₂ and γ̄₂ᾱ as "sequences of GL(2)-metaplectic (cubic-theta) type with a Voronoi dual of length C N(ν)²/Q".

**Why it is wrong.**
- γ₂(p)³ = −α(p) on all 60 primes tested. So γ₂³ = μα on squarefree arguments, since the cubed pair phases are 1. Hence
  ζ = ψ₁ γ̄₂⁴ = ψ₁ · μ · ᾱ · γ̄₂,
  verified on all 4728 squarefree rows at Q = 6000 to 7·10⁻¹⁵. Dropping μ changes 52% of the rows.
- Patterson's coefficient g̃₃ (with ě) equals γ₂ times a cube root of unity at each prime (a ray character mod 9). So γ̄₂ᾱψ₁ is theta-type (up to a Grössen- and ray twist), while ζ is **μ(m) times** that sequence.
- μ is not a character: no ray character is −1 at every prime. So ζ has no Voronoi formula. It is also the object that `V_gaps.md` §1.6 identifies as "the Möbius information carried by the row".

**Consequences.**
- **§3.1.** The ζ column is not a degree-2 calibration. The conclusion "degree-2 duals are invisible at accessible Q" now rests on γ₂ and γ̄₂ᾱ only. Those columns do show the same behaviour, so the conclusion survives.
- **§6, "For the other attacks".** The claim that ζ·(ray class function) has quadratic twists that are "amenable to cubic-theta (GL(2)-metaplectic) technology" is **false**.
  - Its twists are μ-weighted theta sums Σ_m μ(m) τ(m) χ(m) (ν/m)₂.
  - Equivalently, via γ₁γ₂ = μαG, they are sums of a (twisted-)multiplicative function whose prime values are sextic Gauss sums. That is Heath-Brown–Patterson territory, not Voronoi.
  - So the split π = ζ·h + π′ does not isolate an easy piece.
- **§8, gap 8, sharpened.** π = μ · (ρψ₁ᾱγ̄₂) explicitly. (V′_LS) and (V′_β) are therefore Möbius-randomness statements for a sextic × cubic sequence.

### 3. Error: the ζ-free control in (V′_β) is not at the random level

**The claim.** §4 says "The same form with ρ in place of π … is also at the random level."

**E's own logs disagree.**
- At L = 3000, ρ:rand = 1.19 / 1.20 / 1.87 / 2.22 at Q = 6000 / 12000 / 24000 / 48000 (averaged over 8 additive characters).
- At L = 2886 (Q = 24000) it is 1.51.

**My check, with the trivial additive character and a leave-one-scale-out jackknife** (`referee_E_rho_beta.py`). At L = 3000, ρ gives 2.00 ± 0.35 / 2.95 ± 0.59 / 3.96 ± 0.95 at Q = 12000 / 24000 / 48000, growing with Q at fixed L. The source:

| Q | ρ, all x | ρ, b = 1 columns only | ρ at x = p⁴ | π at x = p⁴ |
|---|---|---|---|---|
| 12000 | 2.00 ± 0.35 | 0.94 ± 0.18 | 11.3 ± 1.9 | 1.08 ± 0.23 |
| 24000 | 2.95 ± 0.59 | 1.52 ± 0.20 | 20.1 ± 5.3 | 0.94 ± 0.24 |
| 48000 | 3.96 ± 0.95 | 0.56 ± 0.07 | 31.5 ± 8.6 | 1.03 ± 0.07 |

- **The excess is entirely the columns x = p⁴ and p̄⁴, with N(p) = 7.** Here x = nb³ with n = b = p, and |a(x)| = √7.
- **These are principal-class ("square") columns.** λx = λ(p²)², so (λx/R)₂ = (λ/R)₂ · 1_{p∤R}, which is a ray character mod 4. So the column sum is ρ's untwisted class sum, i.e. ρ's constant-term coherence.
- **The other cube-type columns are at the random level.** These are 7p² and b³ with N(b) = 13: 0.6–1.0.
- **π is at the random level at the same columns.** The μ in ζ destroys the coherence.

**What to change.** Correct the sentence in §4. The finding is also useful:
- **A positive control for (V′_β).** The noisy β-pipeline detects one coherent column carrying about 5% of the window's column energy.
- **π's principal-column sums are random-level for Q ≤ 4.8·10⁴.**

### 4. Gap 2 ("bounded factor [R-]") is not routine: the principal-class columns

**The columns.** The full family ν = λ⁴μ ∈ O∖0 contains ν = uλ^k b². For these, ψ_ν(m) = (uλ^k/m)₂ · 1_{(b,m)=1} is a ray character. There are about L^{1/2} such columns with N(ν) ≤ L. Each has F′(ν) ≈ Σ_m Ω π(m) χ(m), the **untwisted** (Möbius-type) sum of π.

**What (V′_LS) then demands.** (V′_LS) as stated in `TARGET_VLS.md` (unweighted, over all q_μ ≍ L) therefore needs

  |Σ_m Ω π χ|² ≪ Q^{1−η} L^{−1/2} ‖π‖².

At L = Q^{0.79} that is a saving of Q^{−0.2−η/2} in |Σ| over Cauchy–Schwarz for the untwisted π-sum. This is a non-trivial analytic statement about Σ π N(m)^{−z}, the A₂ residue series, which is [H] (μ₁₂ ⊄ Q(ω)). It is not a routine reduction.
- **Numerically it holds.** π's untwisted class sums are random (§2.2), and so is the principal-column test in §3.
- **The repair [R-].** Apply Cauchy–Schwarz with weights matched to |β_μ|², or restrict (V′_LS) to non-principal columns.
  - Inside supp β, the principal-class columns are x = b⁴ (and rarer n⁴c⁶), with |β_x|² = N(b)^{−3}.
  - Their total contribution to T is trivially ≪ Σ_{N(b)≍L^{1/4}} N(b)^{−3/2} Q^{1/2}‖π‖ ≪ Q^{11/12} L^{−1/8}.
- **Correct label for E's [R-]:** "[N] for the unweighted family; [R-] only after the principal columns are moved to the β-weighted bound".

### 5. LP bookkeeping (numbers correct; framing imprecise)

**Where the optimum is.** The random-mode optimum is at M = 0.9137 and b = 0.518, where b/6 = 0.0863 = 1 − M. The theta piece is cancelled exactly by s = (1 − M)/2 = 0.043. So:
- 0.8519 needs (V′_LS) at L = Q^{0.906} with η = 0.094. It does **not** need it at Q^{0.79} with s₀ = 0.087, and M = 0.914 lies outside the window [0.80, 0.86] of `TARGET_VLS.md`.
- It coincides with the constant-saving model: s₀ = 0.043 gives θ = 0.851987 (my run).
- The summary's sentence "s₀ = (1−M)/2 ≈ 0.087 at M = 0.826 … the LP then gives θ = 0.8519" mixes two geometries.

**For reference** (`referee_E_lp_floor.py`, P²/Y dropped), constant saving s₀:

| s₀ | θ |
|---|---|
| 0.06 | 0.844948 |
| ≥ 0.1 (theta piece gone; M → 1) | 0.836753 |

**Marks.** "0.8548 with marks" comes from a model with extra constraint pieces, so by construction it is at least the unmarked value.

### 6. Statistics and estimator

- **The fit for b is not independent evidence.** R_LS = R_rand · D(L)/(Q ΣΩ|π|²), and D ∝ #ν(L) ∝ L is deterministic. So b ≈ 1 is equivalent to "R_rand is flat". The ±0.007 comes from 45 *nested* points and is not an error bar. The γ fit over dyadic bins is the real statistic.
- **Jackknife errors with K = 3–7 scales are too small.**
  - At Q = 12000, R_rand = 0.962 ± 0.008 at both L = Q^{0.85} and Q^{0.9}. That is 4–5 jackknife σ below 1, and it does not recur at other Q.
  - At Q = 96000 (3 scales), R_rand = 1.22 ± 0.08 at Q^{0.5} against a null of 0.92. At Q^{0.79} it is 1.039 ± 0.011 in the 2.5Q window.
  - E's systematic envelopes (±5% for L ≥ Q^{0.7}, ±20% for L ≤ Q^{0.6}) are the right reading. Within them every π value is consistent with 1.
  - The Q = 96000 small-L excess is unexplained. Its attribution to systematics rests on non-replication only.
- **Correlated noise across scales.** The off-diagonal scale Gram (0.72–1.20 at Q = 48000; ~0.03–0.07 expected from independent row noise) shows that the premise "unbiased if the w_k are uncorrelated" fails at the 20–30% level for row-diagonal quantities.
  - R_rand self-normalises (numerator and denominator see the same pair weights).
  - Absolute normalisations (‖ρ‖², κ = 0.64–0.80, the "Q^{0.82}" growth) carry this systematic.
  - A genuine deviation of π from random would be diluted by up to about 25%. That could not hide saturation, which is a factor of 30–50.

### 7. Modelling gaps not in E's list

- **ρ_c versus ρ.** c(s) = 1_{s∈σ} χ_s(b_*)/(τξ(s)) twists the s side. For s prime to S, c(s) g̃(s,−m) = g̃(s, −m b_*⁵)/(τξ(s)). So ρ_c is a ray-twisted residue evaluated at m b_*⁵. Only c ≡ 1 is tested; the six t-components resolve only ray classes of s of small modulus.
- **ξ(m).** It has conductor b_*, the product of *all* primes of S, including every prime of norm ≤ P₀. This goes well beyond the 36·4 of Gap 1. Untested.
- **v ≠ 0.** Untested; [H] harmless.
- **Constants.**
  - The true column norms are N(ν) ≍ 81 q_{c_F}² Q²/Z. At Q ≤ 10⁵ this exceeds Q, so the actual problem never reaches its L < Q regime at accessible Q.
  - The tests use the stripped model L = Q^{0.79}. Exponents are unaffected, but "saturation excluded by 30–50×" is a statement about the stripped model.

### 8. Provability

I agree with E §8: the data cannot say. Two additions.

1. **What the data can and cannot show.**
   - By §2, π = μ·(…), and T contains the zero-detecting term 𝒮 (`V_gaps.md` §1.6; `AUDIT_chain.md` §3).
   - Any proof of (V′_LS) or (V′_β) with η > 0 is therefore itself a Möbius-randomness theorem, precisely the kind that fails under the zero hypothesis the probe is built to detect.
   - Random-level data at Q ≤ 10⁵ are expected whether or not (V′) is provable.
2. **No available tool saves.**
   - The Goldmakher–Louvel quadratic large sieve gives (L + Q)‖π‖², with no saving for L < Q.
   - Heath-Brown's cubic large sieve needs a cubic kernel. The kernel here is quadratic, and the row carries μ.
   - The A₂ WMDS route is [H] (μ₁₂ ⊄ Q(ω)) and does not save (`A2_theory.md` with its referee report).
   - So the "random" LP value 0.8519 (equivalently s₀ = 0.043) is a hypothetical target, not a reduction.

### 9. Strongest surviving claim

- **Proved.** Nothing toward (V′).
- **Numerics [N]** for Q ≤ 9.6·10⁴, in the stripped model (c ≡ 1, ξ ≡ 1, v = 0, squarefree rows, squarefree columns prime to 6), with an exact pipeline (twists, Gauss sums and ζ re-verified):
  - π's quadratic twists are at the random level for all L ≤ Q, with R_rand = 1 within ±5% for L ≥ Q^{0.7} and ±20% below that. Saturation is excluded by a factor of 30–50 at L = Q^{0.79}.
  - The (V′_β) bilinear form with Patterson's coefficients is at the random level for π.
  - π's principal-class (untwisted) sums are at the random level.
  - Degree-2 (cubic-theta) duals are numerically invisible at these Q, so the conductor exponent of π's twists cannot be determined.
- **Conditional [H].** If (V′_LS) held at the random rate at L = Q^{0.906}, i.e. s₀ = 0.043 at M = 0.914, the LP would give θ = 0.8519 (reproduced).

### Required edits

1. §0.7, §3.1 and §6: ζ = ψ₁μᾱγ̄₂ is μ times a theta-type sequence. Remove "GL(2)-metaplectic" for ζ, and remove "amenable to cubic-theta technology" for ζ·h.
2. §4: the ζ-free control is not at the random level at L ≈ 3000. Report the square-column effect (§3 above).
3. §8, Gap 2: relabel as in §4 above, and add the β-weighted repair.
4. Summary and §7: state the LP geometry (M = 0.914, L = Q^{0.906}, s₀ = 0.043), and say that it equals the constant-s₀ value.
5. §3: drop b = 1.008 ± 0.007 as independent evidence.
6. Gaps: add ρ_c versus ρ, ξ of conductor b_*, v, and the constant in L.
