# Numerical probe of π(m) = ρ(m)·ζ_m (the sequence behind (V′))

*9 Oct 2026. Scripts and logs are in `gram_split/wf2/`. The library is `a2lib.py`.*

*Data: `D_numerics/data/A_Q{1500,3000,6000}_K12.npy`, 12 scales Y_i = 256·2^i. We use the 8 scales Y = 4096 … 524288, all with Y > Q, which avoids the deterministic small-Y Poisson bias noted in `V_derivation.md` §8.3.*

Labels:
- **[R]** exact identity or exact implementation, verified to machine precision;
- **[N]** statistical numerical observation, with its noise level stated;
- **[H]** heuristic interpretation;
- **[G]** gap or untested assumption.

---

## 0. Bottom line

1. **ζ_m is implemented exactly** for all squarefree rows at Q = 1500, 3000 and 6000 [R].
   - Lemma Z's closed form ψ₁(R)·conj(γ₂(R))⁴ agrees with the paper's product formula (eq. `reflection-row-phase`, c_F = 1) to 1.2·10⁻¹⁴ on all 8280 squarefree rows.
   - It also matches `check_row_phase.py` at primes: on all 432 classes mod 72, w/ζ^Z is a constant sixth root of unity, while w itself is not constant.
2. **π is not "simpler" than ρ in any way the data can see. It is the opposite: π has lost every automorphic fingerprint that ρ shows.** [N]
   - **ρ fingerprints.** ρ has a reproducible coherent character spectrum mod 36 at k = 0. This is the θ₆ constant term: the cross-Q correlation of the 648-character profile is 0.74, and the top character grows 86 → 94 → 118. ρ is also coherent under the sextic twist j = 1 (the θ₆ multiplier).
   - **π has none of these.** For every Hecke angle α^k (|k| ≤ 6), with or without μ(m), and for every character mod 36, π's profiles are uncorrelated across Q (|corr| ≤ 0.12). No sextic twist j is coherent.
   - **The cubic-type degeneration is excluded.** The cubic analogue θ₃·conj(G) ≡ const would give coherence ratios of about 300 per unit class. π gives 1.3–2.4 (random = 1).
   - **Phase concentration is bounded.** The mean resultant length of the phase of π(m)α(m)^k inside any unit class is ≲ 0.05 on all squarefree rows and ≲ 0.12 on prime rows.
3. **Quadratic twists F′(ν) = Σ_m Ω π(m)(ν/m)₂ are at the random level** for every accessible N(ν) ∈ [7, 4Q] (N(ν)²/Q from 0.04 to 4·10⁴). They are neither coherent nor structurally small, and the same holds for ρ [N].
   - **The data cannot determine the conductor exponent k.** Exact calibrations show that at Q ≤ 6000 even degree-1 sequences (Hecke characters, the constant 1) give twist ratios of about 1 (0.3–0.7 at best for N(ν) < 128).
   - The asymptotic regime C·N(ν)^k < Q is out of reach for every k ≥ 1, so "dual length L²/Q" cannot be told apart from "L³/Q or L⁴/Q" here.
4. **The unreflected pairing ⟨ρ, B⟩ with the true B is at the random level** [N]: R_rand = 0.2–1.7 ± 0.3–0.7, against 1 for random. The ζ̄-inserted "(V)-model" pairing is indistinguishable from it at these sizes.
5. **What this means for (V′)** [H].
   - The numerics rule out the scenario that would kill (V′) outright: π, or the pairing, carrying a coherent component. They also rule out the cubic-type miracle (π a ray function).
   - They show π behaving like a "generic" sequence: no visible Eisenstein or theta constant term, no multiplier, no pole in the tested A₂-residue twists.
   - They give **no evidence about the conductor exponent**. In particular they neither support nor refute dual length L²/Q.

---

## 1. Step 1: exact ζ_m (`check_zeta.py`, `check_zeta_primes.py`)

**Implementation (`a2lib.Rows`)** [R].
- Each data row m (all six associates are separate rows) is factored into primary primes, split and inert, giving m = u·R with R primary.
- γ_j(R) is computed by twisted multiplicativity, γ_j(ab) = (χ_a(b)χ_b(a))^j γ_j(a)γ_j(b).
- Two forms of ζ are computed:
  - **Lemma Z form:** ζ^Z(R) = ψ₁(R)·conj(γ₂(R))⁴, with ψ₁(R) = χ_R(λ)⁻¹·conj(χ_R(4))·χ_R(−1) (c_F = 1).
  - **Product form:** ζ^P(R) = ᾱ(R)² ∏_{p∣R} χ_p(σ_p)⁻² τ⁻_{p,1} τ⁺_{p,3} χ_p(ε_p)⁻³, with σ_p = λ²R/p and ε_p = −((λ³R/p)σ_p)⁻¹. This is eq. `reflection-row-phase` (l. 1832) with the constant −(i/81)φ̂(h₀)κ̄_F dropped.
- Non-squarefree rows are excluded: 6.4–6.7% of rows. Their ζ is §1.3 of `V_gaps.md`, with an unknown frozen constant.

**Checks** (Q = 1500 / 3000 / 6000):

| check | result |
|---|---|
| γ_j(R) from twisted multiplicativity vs brute-force Gauss sum over O/R (68–80 composite R, incl. inert factors and π·π̄), j = 1, 2, 3 | max error 2.7·10⁻¹⁵ |
| eq. `signal`, γ₁γ₂ = μα·conj(χ(4))γ₃ and γ₂³ = μα, on all squarefree rows | max error 7·10⁻¹⁵ |
| ζ^P/ζ^Z − 1 over all 1164 / 2388 / 4728 squarefree rows | max 1.2·10⁻¹⁴ |
| `check_row_phase.py` convention w = ᾱ²τ⁻_{p,1}τ⁺_{p,3}, 3220 split primes N ≤ 3·10⁴: is w/ζ^Z constant on classes mod 72? | constant on all 432 classes; values are sixth roots of unity. Spread of w itself: 1.98 |

**Not implemented [G].**
- κ₀ = −(i/81)φ̂(h₀)κ̄_F ᾱ(c_F)² is a function of the class of R mod M², M = λ¹²L⁴. All tests below are therefore made *class by class*: per unit class (f = 3), per class mod 12, or with all 648 characters mod 36. This absorbs ψ₁, the c_F choices χ_R(c_F)⁻² (cubic, conductor | 18), G, and any ray factor of conductor | 36.
- A κ_F with genuinely larger conductor is not tested.

## 2. Estimator and noise levels

V_i(m) = Y_i^{1/3}A_m(Y_i)/Ŵ = ρ(m) + w_i(m), where the w_i come from disjoint s-ranges.

**Estimators.**
- For any linear statistic T_i = Σ_m Ω x_i h, the estimate of |F|² is the inverse-noise-weighted distinct-scale product
  [|Σ w_i T_i|² − Σ w_i²|T_i|²] / [(Σw)² − Σw²],
  which is unbiased if the w_i are mutually uncorrelated.
- The denominators Σ Ω²|x|²|h|² are estimated the same way, row by row.
- Standard errors are leave-one-scale-out jackknife, plus empirical nulls: the same statistic with h replaced by h·ε, ε random ±1, 20–60 draws.

| Q | rows (sf / prime) | ‖ρ‖²_Ω | noise per scale, Y = 4096 → 524288 | SNR, top scale | SNR, 8 scales combined |
|---|---|---|---|---|---|
| 1500 | 1248 (1164 / 570) | 2.10 | 22.3 → 4.7 | 0.45 | 2.0 |
| 3000 | 2526 (2388 / 1062) | 4.09 | 46.6 → 8.1 | 0.50 | 1.95 |
| 6000 | 5052 (4728 / 2028) | 7.93 | 89.7 → 17.3 | 0.46 | 1.95 |

**What this resolution allows** [N].
- Per row, the signal is only about 2× the noise power. So individual ρ(p) (or their phases) **cannot** be read off. Only aggregate statistics are meaningful.
- **Aggregate class statistics:**
  - Rat(f) = mean over the characters mod f of the coherence ratio;
  - jackknife s.e. is about 0.5–1 for f = 1 or 3, about 0.15–0.3 for f = 12 and about 0.05–0.1 for f = 36.
- **Single-character ratios mod 36:** noise s.d. about 2–3, with heavy tails. The maximum over 648 characters is about 9–18 under the null.
- **Second-order (4-scale) statistics:** null s.d. about 2–3 at f = 1, about 1 at f = 3 and about 0.3 at f = 12, with heavy tails.
- Hence only **cross-Q replication** is accepted as evidence.

## 3. Step 2: is π simpler than ρ?

### 3.1 Hecke angle × class functions (`test_char.py`)

The statistic is Rat(x; α(m)^k, mod f), with |k| ≤ 6 and f ∈ {1, 3, 12, 36}, on squarefree rows and on prime rows.

**ρ** [N].
- k = 0 is coherent at every f, and it grows with Q. At f = 36: 1.35 / 1.85 / 2.67 (s.e. about 0.1–0.3). At f = 1: 86 / 94 / 118, but with jackknife s.e. 30–140 (Y-dependent deterministic components).
- For k ≠ 0, ρ sits *below* 1 at f = 36 (0.6–0.9). This is consistent with energy conservation: about 20% of ‖ρ‖² lies in the (k = 0, mod-36 class-function) subspace.

**π** [N].
- Every (k, f) is consistent with 1. For example, at f = 36, Q = 6000 the values are 0.79–1.23 with s.e. about 0.08.
- The untwisted per-unit-class value (k = 0, f = 3) is 1.94 / 2.44 / 1.28 (s.e. about 0.5–1) for Q = 1500 / 3000 / 6000, compared with ρ's 15 / 16.7 / 20.2 (from the j = 0 column in §3.3).
- **There is no growing main term.** If Σπ(m)N(m)^{-s} had a pole at s₀ > 5/12, the ratio would grow like Q^{2s₀−5/6}.

**The cubic-analogue benchmark.**
- If π were constant on classes, as θ₃·conj(G) is in the cubic case, the f = 3 ratio would be about (#rows per class)·(E|ρ|)²/E|ρ|² ≈ 788 × 0.7 ≈ 500 at Q = 6000.
- For comparison, the exact sequence "1" gives 428.
- The observed value is ≤ 2.4. Converting the 95% upper limits into a mean resultant length r of the phase of π·α^k inside a class gives r ≲ 0.05 (all squarefree rows, f = 3) and r ≲ 0.12 (prime rows, f = 3), for all |k| ≤ 6 [N].
- **Phase concentration of ρ(p)γ₁(p)ᾱ(p)² within ray classes is therefore absent.**

### 3.2 Reproducible character spectrum mod 36 (`test_chars36.py`, `test_chars36_sweep.py`, `test_a2pole.py`)

All 648 characters of (O/36)^× are built explicitly. For each case we compute the 648-vector Rat(ψ) at Q = 3000 and at Q = 6000, and record two quantities:
- **corr:** the correlation of the two 648-vectors;
- **#both:** the number of characters with Rat > 8 at both Q.

A genuine coherent direction must reproduce. Results [N]:

| x · h | corr / #both at k = 0 | max \|corr\| over k ≠ 0, \|k\| ≤ 6 |
|---|---|---|
| ρ · α^k (sf rows) | **0.74 / 15** | 0.17 |
| ρ · α^k (prime rows) | **0.56 / 7** | 0.11 |
| ρ·γ₃^{±1} · α^k (control: γ₃ is a ray function) | **0.47 / 8** | 0.13 |
| π · α^k (sf rows) | 0.02 / 0 | 0.09 |
| π · α^k (prime rows) | 0.07 / 0 | 0.12 |
| π · μ(m) · α^k | 0.07 / 0 | 0.11 (#both ≤ 2, chance level) |
| ρ · γ₁^{±1} · α^k (A₂ residue Σρ·g̃₆(m)) | ≤ 0.11 / 0 | 0.13 |
| ρ · γ₂^{±1} · α^k | ≤ 0.17 / 2 | 0.11 |

**Interpretation** [H].
- ρ's constant-term coherence reproduces, so the test works.
- π, and the "A₂ residue" products ρ·γ_j^{±1}, show no reproducible coherent character for any k.
- If Res_w of the sextic A₂ series had a z-pole with non-negligible residue in these twist families, the double residue at z = 2/3 would give a ratio of about Q^{1/2} × const, and we would expect it to show. We see none.
- **Caveats.** Ray conductor > 36, |k| > 6, or a tiny residue are not excluded [G].
- A weak single cell, ρ·conj(γ₂)α⁴ ≡ π·μ·conj(ψ₁)α⁵ at f = 36 (1.02 / 1.27 / 1.15), is not reproducible character by character (corr 0.01–0.06) and does not grow with Q. It is noise.

### 3.3 Sextic twists (m/ν)₆^j, ν prime: the automorphy fingerprint and positive control (`test_jtable.py`)

The entries are per-unit-class ratios, Q = 6000 (Q = 1500 and 3000 agree):

| N(ν) | ρ, j = 0 | ρ, j = 1 | ρ, j = 2..5 | π, j = 0 | π, j = 1..5 |
|---|---|---|---|---|---|
| 7–32 | 18.9 | **4.6 ± 1.4** | 0.76–0.92 | 1.23 | 0.25–1.27 |
| 32–128 | 20.0 | **2.2 ± 0.3** | 0.60–1.21 | 1.26 | 0.68–1.11 |
| 128–512 | 20.1 | **1.65 ± 0.10** | 0.72–1.01 | 1.28 | 0.94–1.15 |
| 512–2048 | 20.2 | 1.38 ± 0.06 | 0.97–1.10 | 1.28 | 0.98–1.05 |

- **ρ** [N] shows the θ₆ multiplier (j = 1 coherent, decaying in N(ν)) and the constant term (j = 0). This reproduces `V_derivation.md` §3 and is the positive control.
- **π** [N] has no coherent j. Its j = 0 column is flat at the random-level value of its untwisted sum.
- **Interpretation** [H]: there is no sign that π is the coefficient sequence of a form with a constant term and a residue-symbol multiplier at the cusps a/ν.

### 3.4 Second order at primes (`test_order2.py`)

**Motivation.** Bröker–Hoffstein Conj. 5.7 gives τ₆(p)² = ζ₆^{k(p)}·(cubic Gauss sum)·F(p), with F in a finite set (powers of 3, and elements of norm 7 or 13). Algebraically, π(p)² = −ρ(p)²·ψ₁²·ᾱ³·γ₂(p). So if ρ² ∝ conj(γ₂)·F·ζ₆^k, then π² ∝ ᾱ³·F·ζ₆^k, a Hecke character times a finite-valued phase.

**Statistic.** Σ_cl |Σ_p Ω x(p)² h|² / Σ Ω²|x|⁴ is estimated with 4 distinct scales, for:
- x = ρ with h = α^k γ₂^e, |k| ≤ 9, e ∈ {−1, 0, 1};
- x = π with h = α^k;
- f ∈ {1, 3, 12}, at all three Q.

**Result** [N].
- 17 of the 228 (x, e, k, f) cells have z > 3 at some Q (684 cell–Q pairs in all). This is consistent with the heavy-tailed null, whose maximum reaches about 10 at f = 1.
- **No cell has z > 2 at all three Q.**
- ρ(p)², k = 0, e = 0, f = 1 gives −1.8 / 6.4 / 13.7 (z = −1.5 / 2.0 / 4.5). That is suggestive for ρ only, unreplicated, and irrelevant to π.

**Detectability.**
- Within unit classes, a non-uniform distribution of the finite phase ζ₆^{k(p)}F(p) with mean resultant length r ≳ 0.2 would have been seen. Smaller biases, or a uniform k(p), are invisible.
- This test **cannot** check Conj. 5.7 itself; that needs the 6th power, which is far too noisy.

## 4. Step 3: quadratic twists, the structural test (`test_twist.py`, `summarize_twist.py`)

**Statistic.**

F′_c(ν) = Σ_{m ∈ class c, m squarefree} Ω π(m) (ν/m)₂,   ν = λ⁴μ squarefree primary, 7 ≤ N(ν) ≤ 4Q,

with 1695 / 3390 / 6766 columns for Q = 1500 / 3000 / 6000. Note (λ⁴μ/m)₂ = (ν/R_m)₂ exactly. The ratio is Rat(ν) = Σ_c |F′_c|² / Σ Ω²|π|²|(ν/m)₂|², with the same statistic for ρ and for exact calibration sequences.

**Note on f.** f = 1 is contaminated: the six associates of one ideal have identical twists and identical R-only sequences, so exact R-only sequences give about 6. This is why per-class statistics are used.

f = 3 (unit classes separately), mean over ν, with s.e. over ν and (Q = 6000) the jackknife-over-scales s.e. of the bin mean:

| N(ν) | N²/Q · N³/Q · N⁴/Q (Q = 6000) | ρ, Q = 6000 | π, Q = 6000 | π, Q = 3000 | π, Q = 1500 | "1" (deg 1) | ᾱ² (deg 1) | γ₂ (θ₃-type) | rand |
|---|---|---|---|---|---|---|---|---|---|
| 7–32 | 0.04 · 0.7 · 11 | 1.07 (0.29)[0.4] | 1.37 (0.08)[0.3] | 0.76 | 1.03 | 0.37 | 0.72 | 0.75 | 0.97 |
| 32–128 | 1 · 76 · 6·10³ | 1.13 (0.14)[0.3] | 1.05 (0.11)[0.15] | 1.16 | 0.98 | 0.47 | 0.78 | 0.81 | 0.90 |
| 128–512 | 15 · 4·10³ · 10⁶ | 0.90 (0.06)[0.1] | 1.02 (0.07)[0.06] | 1.11 | 0.94 | 0.98 | 1.3 | 0.91 | 0.96 |
| 512–2048 | 240 · 3·10⁵ · 4·10⁸ | 1.12 (0.03)[0.05] | 1.01 (0.03)[0.05] | 1.00 | 0.96 | 1.15 | 1.1 | 1.06 | 0.99 |
| 2048–8192 | 4·10³ · 2·10⁷ · 9·10¹⁰ | 0.96 (0.02)[0.02] | 0.99 (0.02)[0.01] | 0.99 | 1.01 | 0.96 | 1.0 | 0.96 | 1.00 |
| 8192–24000 | 4·10⁴ · 6·10⁸ · 10¹³ | 0.99 (0.01)[0.01] | 0.98 (0.01)[0.01] | 0.99 | — | 1.11 | 1.0 | 1.02 | 1.01 |

**Pooled over all Q, N(ν) < 128** (105 columns):

| sequence | f = 3 | f = 12 |
|---|---|---|
| ρ | 0.94 ± 0.06 | 1.01 ± 0.02 |
| π | 1.06 ± 0.05 | 0.99 ± 0.02 |
| "1" | 0.71 ± 0.09 | 0.90 ± 0.03 |
| ᾱ² | 1.04 ± 0.08 | 0.96 ± 0.03 |
| γ₂ | 0.90 ± 0.09 | 0.97 ± 0.03 |
| rand | 0.97 ± 0.04 | 0.99 ± 0.01 |

(The s.e. over ν understates the error at small N(ν), because all columns share the same per-scale noise; use the jackknife values in brackets.)

**Findings.**
- **π: random at every N(ν)** [N]. There is no coherence (≫ 1) and no structural smallness (≪ 1), at f = 1, 3 and 12. The ratio is flat in N(ν), and hence flat in N(ν)^k/Q for every k. ρ behaves the same.
- **The test has no power to see the conductor** [N, decisive for interpretation].
  - Exact degree-1 sequences, whose quadratic twists have conductor ≍ N(ν), are already at about 1 for N(ν) ≥ 128. Below that they reach only 0.3–0.7.
  - This holds even on *all* rows, without the squarefree restriction (`test_jtable.py` second part): "1" gives 0.43 / 0.48 / 0.30 / 0.36 for N(ν) in 7–128, then 1.2 / 1.3 / 0.5 / 1.2 / 0.85 / 0.96 above.
  - The cause [H]:
    - the 2-adic and 3-adic conductor of the reciprocity factor (f ≥ 12, so N(f) ≥ 144);
    - the narrow radial weight Ω (annulus width ≈ 0.3√Q), which make C·N(ν)^k/Q ≥ 1 for every ν ≥ 7 even when k = 1;
    - the squarefree restriction, which turns a would-be ≈ 0 into ≳ (C N(ν)/Q)^{1/2} through the Möbius tail Σ_{d²|m}.
- **Consequence** [H]. Reaching the regime where k = 2 differs from k = 3 or 4 (C·N(ν)²/Q < 1 at N(ν) ≈ 7 with C ≳ 10⁴, as in the cubic numerics of `V_derivation.md` §8.2) needs Q ≳ 10⁶. Computing ρ there needs scales Y ≫ Q with Gauss sums of norm up to about 10⁸, which is out of reach. **The data say nothing about whether (V′) has dual length L²/Q.**

## 5. Step 4: unreflected pairing with the true B (`test_pairing.py`)

**Set-up.**
- B_m = Σ_{c sf, N c ∈ [Z, 2Z]} W(Nc/Z)·γ₂(c)ᾱ(c)·conj(G(c))·ν(c)·χ_c(qm)/√Nc, with ν = conj(χ_c(−1)), for 10 fixed q, i.e. 10 ray characters ν, as in `exp_pairing2.py`.
- This B carries ζ_m implicitly, since its reflection is ζ_m Σ_μ β_μ ψ_m(μ).
- T = Σ_m Ω ρ B_m, R_rand = E|T|²/Σ Ω²|ρ|²|B|² and R_CS = E|T|²/(ΣΩ²|ρ|²)(ΣΩ²|B|²).

| Q | Z (L = Q²/Z) | true ⟨ρ, B⟩: mean R_rand over 10 ν (s.e.) | R_CS | ζ̄-inserted ⟨ρ, ζ̄B⟩ (sf rows) | 1/#rows |
|---|---|---|---|---|---|
| 6000 | 3000 (12000) | 1.72 (0.68) | 1.1·10⁻³ | 0.45 (0.25) | 2.0·10⁻⁴ |
| 6000 | 10000 (3600) | 0.21 (0.25) | 1.4·10⁻⁴ | 1.01 (0.45) | 2.0·10⁻⁴ |
| 6000 | 20000 (1800) | 0.76 (0.63) | 5.2·10⁻⁴ | 1.44 (0.59) | 2.0·10⁻⁴ |
| 3000 | 3000 (3000) | 1.20 (0.33) | 1.6·10⁻³ | 1.07 (0.42) | 4.0·10⁻⁴ |
| 3000 | 20000 (450) | 1.30 (0.32) | 1.7·10⁻³ | 0.73 (0.24) | 4.0·10⁻⁴ |

- **Result** [N]: random level. This reproduces `exp_pairing.py`'s conclusion with the K12 data and the weighted estimator. The jackknife s.e. of a single pairing is about 1–2 in R_rand units.
- **What it can and cannot distinguish** [H].
  - It excludes coherence. At CS level R_CS would be about 1, but it is 10⁻⁴–10⁻³.
  - It cannot tell "random" from "a few times below random".
  - Numerically |T| ≈ #rows^{−1/2}·‖ρ‖‖B‖, i.e. full square-root cancellation at these sizes. This is what *any* generic sequence does, so it is no evidence that a power saving is *provable*.

## 6. Summary of what the data can and cannot distinguish

| question | answer | strength |
|---|---|---|
| Is ζ_m (Lemma Z) implemented correctly for the data rows? | yes, to 10⁻¹⁴, both forms | [R] |
| Is π a ray-class function times a constant (cubic analogue)? | no: r ≲ 0.05 within unit classes, against r = 1 for the cubic case | [N], strong |
| Does π keep ρ's automorphic fingerprints (constant term, j = 1 multiplier)? | no: none reproducible in any tested family | [N], strong for the tested families |
| Does π correlate with α^k × (characters mod 36), μα^k, or γ_j^{±1}α^k? | no reproducible cell (\|k\| ≤ 6, conductor \| 36) | [N]; larger conductors and \|k\| untested [G] |
| Does Res_w of the sextic A₂ series (ρ·γ_j^{±1}) show a z-pole in these twist families? | not visible | [N] / [H] |
| Is π(p)² concentrated (BH-type finite phase) relative to α^k within classes? | no replicated signal; would detect r ≳ 0.2 | [N], weak test |
| Are quadratic twists of π coherent? | no, ratio about 1 at all N(ν) ≤ 4Q | [N], strong |
| Are quadratic twists of π structurally small (dual sum empty)? | not at accessible sizes; but neither are degree-1 calibrations | [N] (no power) |
| Conductor exponent k (dual L²/Q vs L³/Q or L⁴/Q)? | **undeterminable** at Q ≤ 6000 (needs Q ≳ 10⁶) | — |
| Is the true pairing ⟨ρ, B⟩ coherent? | no, random level | [N] |

**For the (V′) programme** [H].
- The data remove two possibilities:
  - **a cheap cubic-type proof:** π is not degenerate;
  - **a hidden obstruction:** neither π nor the pairing has a coherent part.
- π looks like a generic, non-automorphic-looking sequence at these sizes. That fits Lemma Z's description, ρ·γ₁·ᾱ²·(ray), with the ᾱ² archimedean twist placing it outside the spherical A₂ series.
- Any proof of (V′) must come from structure (an A₂ / WMDS functional equation in the m-variable with an archimedean K-type twist, and uniformity in μ). That structure is invisible to finite-size numerics of this kind.

## 7. Files

| purpose | files |
|---|---|
| library | `a2lib.py`: rows, factorisation, γ_j, ζ^Z / ζ^P, cross-scale estimators |
| step 1 | `check_zeta.py` (writes `zeta_Q*.npz`), `check_zeta_primes.py` |
| step 2 | `test_char.py` (+ `test_char_Q*.log`), `test_chars36.py`, `test_chars36_sweep.py`, `test_a2pole.py`, `test_corr.py`, `test_mobius.py`, `test_jtable.py` (+ logs), `test_order2.py` (+ log) |
| step 3 | `test_twist.py` (+ `test_twist_Q*.log`, `twist_Q*.npz`), `summarize_twist.py` (+ log) |
| step 4 | `test_pairing.py` (+ `test_pairing_Q*.log`) |

`lp_dual_route.py` in this directory is not part of this note; it was written by another process.

---

## Referee report

*Independent skeptical review, 9 Oct 2026. Sources checked: this note; `a2lib.py`, `test_twist.py`, `test_char.py`, `test_pairing.py` and their logs; `V_derivation.md` §8.3 and §9; `wf/V_gaps.md` §1 and its referee report; `wf/AUDIT_chain.md`; `wf/FE_LoA6_proof.md`; `paper78.tex` ll. 1680–1850 and 2088–2200.*

*New scripts, all in `wf2/`:*
- *`referee_A2_window.py`: scale-window dependence and per-scale class sums;*
- *`referee_A2_profile.py`: per-scale π sums and conjugation symmetry;*
- *`referee_A2_f72.py`: classes mod 72;*
- *`referee_A2_calib_lib.py` and `referee_A2_calib.py`: independent Jacobi-symbol twist, and degree-1 calibration at large Q.*

### Verdict: NOT REFUTED

The main conclusion stands, and my checks strengthen it. At Q ≤ 6000 the numerics neither support nor refute (V′), and they cannot determine the conductor exponent of the quadratic twists of π.

Several supporting statements are wrong or overstated (items 2–4). One of them is a concrete estimator-bias error, but it affects only the ρ-controls, not π. None of them changes the verdict.

### 1. What I verified

- **Quadratic twist matrix [R].** For m with squarefree norm n, O/m ≅ Z/n via ω ↦ r, so (ν/m)₂ is the Jacobi symbol ((ν_a + ν_b r) mod n | n). This independent computation agrees with a2lib's ∏_{p∣m} χ_p(ν)³ on all 151,200 (row, ν) pairs tested (Q = 6000, 4320 rows, 35 columns with N(ν) < 128): 0 mismatches.
- **ζ_m [R].** I did not re-implement it. Lemma Z and the product formula had already been confirmed by `check_row_phase.py` and `wf/referee_composite_phase.py`. I re-derived π(p)² = −ρ(p)²ψ₁²ᾱ³γ₂(p) by hand from γ₂³ = μα.
- **A symmetry the note does not mention [R].**
  - V_i(m̄) = conj V_i(m) holds to 3·10⁻¹⁵ for every row at Q = 3000 and 6000, and ζ(m̄) = ±conj ζ(m).
  - So the unit classes come in conjugate pairs (labels 0↔5 and 1↔3; 2 and 4 are self-conjugate), and the columns ν and ν̄ give identical Rat(ν).
  - The f = 3 statistics and the ν-averages therefore have about half their nominal degrees of freedom. The empirical nulls absorb this, but the "s.e. over ν" in §4 understates the error by about √2.
- **Classes mod 72 [N], which the note did not test.** Mod 72 is the natural level (Γ(72), and the modulus used by `check_row_phase.py`). Setup: π·α^k with |k| ≤ 6, squarefree and prime rows, Q = 3000 and 6000, scales Y ≥ 8192.
  - π: Rat ∈ [0.86, 1.09], s.e. 0.01–0.05, against a null of 1.00 ± 0.03.
  - With 2.1 rows per class at Q = 6000, "constant phase on classes mod 72" would give Rat ≳ 1.3. It is excluded.
  - ρ, by contrast, does show mod-72 structure (k = 0: 1.44 ± 0.04 on squarefree rows, Q = 6000).
  - So the "no ray-class degeneration" conclusion extends from mod 36 to mod 72.

### 2. Estimator bias for periodic test functions [concrete error; affects only the ρ-controls]

**The claim.** The header says the 8 scales Y = 4096 … 524288 are "all with Y > Q, which avoids the deterministic small-Y Poisson bias". This is false on two counts.

**(a)** Y = 4096 < Q = 6000.

**(b) The bias persists well above Y = Q.** The per-scale untwisted sums Σ_{sf m} Ω V_i(m) are:

| Q | Y = 4096 | 8192 | 16384 | 32768 | 65536 | 131072 | 262144 | 524288 | per-scale s.d. |
|---|---|---|---|---|---|---|---|---|---|
| 6000 | +21.0 | −45.4 | −39.8 | +8.8 | +34.0 | +43.7 | +55.4 | +46.9 | ≈ 4–8 |
| 3000 | −23.8 | −16.4 | +1.7 | +18.1 | +27.1 | +30.8 | +30.2 | +18.6 | ≈ 3–6 |

Q = 1500 shows the same pattern. This is the deterministic k = 1 Poisson term of `V_derivation.md` §8.3 (Q·Y^{−1/2}Ŵ(√(Q/Y))), and it is alive up to Y ≈ 3Q. For periodic h, the noise w_i is not orthogonal to class indicators, so the cross-scale estimator is biased.

**Consequences.**
- **ρ's k = 0 numbers depend on the scale window:**

  | window | f = 1 (Q = 1500 / 3000 / 6000) | f = 3 (Q = 1500 / 3000 / 6000) |
  |---|---|---|
  | Y ≥ 4096 (as reported) | 86 / 94 / 118 | 15.2 / 16.8 / 20.2 |
  | Y ≥ 8192 | 101 / 125 / 121 | 17.8 / 21.7 / 20.3 |
  | Y ≥ 16384 | 121 / 169 / 195 | 21.1 / 28.9 / 32.7 |

  So the quoted values, and the growth "86 → 94 → 118", are not measurements of a fixed |Σρ|². The large-Y plateau of the class sums (about 14, 30 and 50 at Q = 1500, 3000, 6000) shows that ρ's constant-term coherence is real and in fact stronger than reported. The positive control therefore survives qualitatively, but its reported magnitudes and its "growth with Q" are unreliable. The ρ, j = 0 column of the j-table (≈ 20) is the same contaminated quantity.
- **‖ρ‖²_Ω at Q = 6000 is inflated by the Y = 4096 scale.** sig2 drops from 7.93 to 7.20 with Y ≥ 8192. On squarefree rows with Y ≥ 8192 the norms are 2.09 / 3.78 / 6.84, i.e. ∝ Q^{0.85}, which matches G3's Q^{5/6}. The table's 2.10 / 4.09 / 7.93 (∝ Q^{0.96}) is biased. All Q = 6000 ratios move by ≲ 10%.
- **π is not affected.** ζ_m is not periodic, so the Poisson term does not survive in π's sums.
  - The per-scale π sums fluctuate at the noise level with no trend in Y.
  - π's f = 3, k = 0 ratio is stable across the three windows: 1.94 / 2.44 / 1.28 → 1.76 / 2.47 / 1.39 → 1.79 / 2.39 / 1.31.
  - The same holds for the ν-twisted sums. **No π conclusion changes.**

### 3. "No power to see the conductor" is slightly overstated, and the Q threshold is too pessimistic [N]

**The test does separate degree 1 at Q = 6000.** At small N(ν):

| sequence | N(ν) ∈ [7, 32) | N(ν) ∈ [32, 128) |
|---|---|---|
| "1" | 0.37 | 0.47 |
| π | 1.37 | 1.05 |
| ρ | 1.07 | 1.13 |

So π's quadratic twists are not degree-1-like. What the test cannot do is separate ρ (heuristic dual length N(ν)²/Q) from π. Hence L²/Q versus L³/Q stays open. The correct wording is "no power beyond degree 1".

**Exact degree-1 calibration at larger Q.** Rows with squarefree norm, f = 3, same Ω:

| Q | N(ν) ∈ [7, 32) | [32, 128) | [128, 512) |
|---|---|---|---|
| 6·10³ | 0.37 | 0.57 | 1.04 |
| 6·10⁴ | 0.06 | 0.33 | 0.53 |
| 6·10⁵ | 0.027 | 0.11 | 0.15 |
| 2.4·10⁶ | 0.010 | 0.044 | 0.074 |

- A degree-1 object is already clearly visible at Q ≈ 6·10⁴, not 10⁶.
- Separating k = 2 from k = 3 needs C·N(ν)²/Q < 1 with an unknown C. The report's 10⁶ may be right for *that* question, but its basis ("C ≳ 10⁴ as in the cubic numerics") is an extrapolation [H].
- Q = 6·10⁴ costs about 10× the rows and 10× the s-range of the Q = 6000 run: expensive, but not obviously impossible.

**A structural limit [N].** The squarefree-restriction (Möbius-tail) floor decays only like about (N(ν)/Q)^{1/2}: it is still 0.15 at Q = 6·10⁵ for N(ν) ≈ 250. For squarefree-only sequences like π, the dynamic range of "structurally small" therefore stays small even at large Q.

### 4. Interpretive claims that overreach [H]

1. **"π has lost every automorphic fingerprint" / "non-automorphic-looking".** These tests detect only non-cuspidal structure: constant terms, multiplier coherence at the cusps a/ν, and ray-class degeneracy.
   - A cusp form's coefficients, or an ᾱ²-twisted A₂ residue without a pole, would give exactly the observed random profile.
   - So the data cannot tell "automorphic with small conductor" from "generic".
   - Correct wording: "no Eisenstein- or theta-type (non-cuspidal) component in the tested families".
2. **"The data rule out a hidden coherent obstruction."**
   - By `V_gaps.md` §1.6 and `AUDIT_chain.md` (R2 and the circularity check), the true pairing is T = 𝒮 + T′. Here 𝒮 is the zero-detecting Möbius signal over A ≍ Z: a structural coherent term whose size is set by the zeros of L(s,η), and which is ≍ Z^{C(β*)} under the contradiction hypothesis.
   - If those zeros do not exist, any (V′) that is true but unprovable also produces random-level pairings. Random-level data therefore say nothing about *provability*, which is the actual obstruction.
   - The note's own §5 hedge is right. The summary and key findings should not say that an obstruction is ruled out.
   - What *is* ruled out is a zero-independent main term visible at Q ≤ 6000, for example a pole of Σπ N(m)^{−s} at s₀ > 5/12, which would make the ratio grow like Q^{2s₀−5/6}.
3. **Cubic-analogue benchmark.** §0 says "about 300"; §3.1 and the summary say "about 500".
   - With the Ω weights, the exact sequence "1" gives 428, so the constant-phase benchmark is about 428 × π/4 ≈ 340.
   - The phase-concentration bounds then become r ≲ 0.06 (squarefree rows) and r ≲ 0.15 (prime rows), not 0.05 and 0.12.
   - They also assume |ρ| is roughly Rayleigh-distributed; heavier tails would weaken them.
   - The qualitative exclusion stands.
4. **The BBFH hypothesis.** The "A₂ residue" and "WMDS functional equation" language (§3.2, §6) presumes the sextic A₂ WMDS over Q(ω).
   - Brubaker–Bump–Friedberg(–Hoffstein) assume μ_{2n} ⊂ F, i.e. μ₁₂ ⊂ F for n = 6. This fails for Q(ω), exactly as `FE_LoA6_proof.md` found for Brubaker–Bump.
   - The series would have to be built in Kubota's or Friedberg–Hoffstein's way, with Hilbert-symbol corrections at 2. These are ray characters mod 4, so the class tests are unaffected.
   - Its functional equations, and the pole at z = 2/3 assumed in §3.2 ("ratio ≈ Q^{1/2}"), are therefore [H], not a citable theorem.
   - The note does label them [H], so nothing breaks. But the structural route it recommends has no ready-made theorem behind it.
5. **Column and weight modelling (harmless).**
   - (V′) as stated uses squarefree μ, while the paper's d(μ) is supported on uλ^k n b³. The extra unit, λ-power and cube factors have quadratic symbols that are ray characters mod 12, which f = 12 absorbs.
   - The reflected weight V♯(q_μX/q_c²) couples μ and R; here it is modelled by Ω alone.
   - The pairing test at (Q, Z) = (6000, 3000) has L = 12000 > Q, outside the regime L < Q. It is fine as a toy.

### 5. Bottom line

**NOT REFUTED.** "Inconclusive" is the right verdict, and the corrections above reinforce it:
- the ρ positive control is quantitatively biased, though it survives;
- the twist test can only exclude degree-1 behaviour;
- the data are blind, by construction, to the obstruction that matters: whether a bound can be *proved* for a pairing that contains the zero-detecting signal.

Nothing here supports or refutes a power saving s₀ > 0 for (V′), and nothing determines the conductor exponent. The arithmetic of the dual lengths in the task statement checks (L = Z^{0.652}, L²/Q = Z^{0.478}, L³/Q = Z^{1.13}, L⁴/Q = Z^{1.78}).

**Required edits.**
1. Delete "all with Y > Q". Use Y ≥ 8192 at Q = 6000, or state that the k = 0 numbers for ρ depend on the scale window.
2. Change "no power" to "no power beyond degree 1", and give the degree-1 threshold as Q ≈ 6·10⁴.
3. Soften the "automorphic fingerprint" and "obstruction ruled out" sentences.
4. Fix the 300 versus 500 inconsistency and the r-bounds.
5. Add the μ₁₂ caveat for the WMDS route.
