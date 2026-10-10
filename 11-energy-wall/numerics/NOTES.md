# numerics: the marked reflected energy at the weighted geometry (energy-wall program, step 1)

Agent `numerics`, 10 Oct 2026. Labels: **[P]** exact identity or exact implementation, checked to machine precision; **[N]** numerical statistic, with its error; **[H]** heuristic, i.e. an extrapolation of a fitted law.

## 0. Verdict

1. **[N] W(w) holds numerically in the weighted-geometry regime, with the random-matrix law and nothing else.** The relative off-diagonal ρ = E/Diag − 1 behaves like ρ_rms ∝ rows^{−1/2}. It has no dependence on F = N/K and none on the mark length L. Its size equals that of random-phase controls with the same moduli, to 1%. This holds:
   - for the single-mark energies E_P, for 13 primes with N(P) = 7 … 1009, i.e. up to and beyond N^{1/2};
   - for the coherent multi-mark core sum 𝒩 with P78's actual coherent mark coefficients c_P = ᾱ(P)²γ₂(P)N(P)^{−1/2} (one or two slots, with the theta cocycle) and the exact mask 1_{(P,k)=1};
   - throughout F ∈ [1, 181], with L/F up to ≈ 1000 (weighted geometry: L ≫ N/K).

   The configurations that satisfy both weighted-geometry conditions (L² ≤ N and L/F ≥ 10) show the same law.
2. **[H] Margin at the weighted geometry.** The geometry is (M, ℓ, w) = (0.851, 0.1915, 0.043), the lossless W(0.043) optimum of wall4 §8. The fitted laws give ρ ≈ Z^{−0.41 ± 0.02} for the coherent 𝒩 and Z^{−0.43 ± 0.01} for single marks. W needs ρ ≪ Z^ε, so the **margin is ≈ Z^{0.41}**. Over the BRIEF range (M, ℓ) ∈ [0.82, 0.86] × [0.18, 0.20] the margin is 0.39–0.43.
3. **[N] Rank-one pairing T (wall4 §7).** T was tested against the sextic Gauss-sum probe vector A(Y), whose pole part is ρ, and against π = ρζ. The ρ-component was isolated with a cross-scale estimator. Result: |T(ρ)|² = (0.99 ± 0.03) × its random-model size, and |T(π)|² = (0.93 ± 0.03) ×. The Cauchy–Schwarz-with-W need grows like rows^{0.975±0.006}. So |T|/Q^{11/12} ≈ K^{−0.49} → **Z^{−0.415}**, which is wall4's random model (Q^{5/12} against Q^{11/12}). A planted alignment is detected (§5.3).
4. **What the numerics cannot exclude** (§7):
   - log factors;
   - secondary terms below ≈ 0.2% (F ≤ 4) or ≈ 1% (F ≤ 64) relative size;
   - effects whose onset lies beyond rows ≈ 3·10⁴ or L ≈ 4·10³;
   - the asymptotic exponent itself.

   The extrapolation to Z^{0.85} is the hypothesis that the random-matrix law persists. Within reach there is no trace of a Patterson (NL)^{2/3} term, neither on generic rows nor on the principal row. If realised on generic rows, that term would exceed the diagonal by (L²/N)^{1/3} ≈ 2–9 in our L² > N configurations, where ρ = 0 ± 0.05.

## 1. Re-validation (step 1)

| check | result | file |
|---|---|---|
| `edge-theta/test_lib.py` | 38/20000 LUT mismatches, **all at the inert prime 2** (norm 4), where the sextic symbol is undefined. Rows and columns are prime to 6, so this is harmless [P] | `test_lib.out`, `check_lut.py` |
| identity (∗) of WALL_SPECTRAL §1.1, rerun | 540 tests, 0 failures; without μ(n), 351 fail [P] | `identity_check_rerun.out` |
| independent direct check (`validate.py`): β(n) = γ₂(n)ᾱ(n) from the library vs Gauss sums summed over Z[ω]/n (ω(n) ≤ 3, inert factors included); Möbius partner a = μḠ = conj(γ₁β); (∗) for random k, P | β: 9.8e-15; partner: 1.3e-14; (∗): 42 tests (15 with inert factors), 0 failures, max error 1.4e-13 once h = k³P⁴ is reduced mod n exactly. With floating-point h the error is ~3e-3, which explains wall-spectral's 1e-4 tolerance [P] | `validate.out` |
| fast data: `cubic_data.c` (cB-numerics, copied). Its G̃ uses e(Tr x/n), while P78 uses e(Tr(x/(n√−3))), so **γ₂(n) = (λ/n)₃·G̃(n)** with λ = 1+2ω | 778 primes plus 580 composite columns, max error 2.5e-13 [P] | `conv_check.out` |
| C kernel for ψ_k(n) = (n/k)₂ (Legendre on A+Br, or on N(n) for inert and conjugate-pair factors) and the row-sum engine | equal to the direct sextic-symbol evaluation to 2.8e-13 on 398 rows (98 with inert or conjugate-pair factors); masked sums and diagonals to 1e-15 | `test_engine.out`, `test_mask.py` |
| probe kernel: A_{m,σ}(Y) and γ₁(s) = μ(s)G(s)β̄(s) | A agrees with the direct evaluation to 2e-16; γ₁ agrees with brute-force Gauss sums to 1.2e-14 on 91 s | `probe_N*.json`, `test_gamma1.py` |

## 2. Setup

- **Reflected sums.** S_k = Σ_n c(n)ψ_k(n). Columns: n squarefree, primary, odd, N₀ < N(n) ≤ 2N₀, weight w = sin²(π(N/N₀−1)). Rows: k squarefree, primary, odd, N(k) ≤ N₀.
- **Exact diagonal.** D_k = Σ_n|c|²1_{(n,k)=1}.
- **Windows.** K = N₀/F with F = 2^{j/2} ≥ 1. A window contains the rows K/2 < N(k) ≤ K, weighted by V = sin²(2π(N(k)/K − ½)); the principal row k = 1 is excluded. Then ρ = ΣV(|S|²−D)/ΣVD.
- **Sectors.** The rows of a window are split into 12 angular sectors of arg k, giving ρ_s; ρ_rms = rms_s ρ_s, as in edge-theta.
- **Fits.** log ρ_rms = a + b log(rows/sector) + c log F [+ d log L]. Five values of N₀ break the collinearity rows ∝ N₀/F. Each fit uses 65 windows per family (rows/sector ≥ 8).
- **Sizes.** N₀ = 12500, 25000, 50000, 10⁵, 2·10⁵, giving 3534 … 56470 columns and up to 56434 rows; a window has up to 28210 rows, i.e. 8 to 2350 rows per sector. Each run took 4–342 s (`logs/`). The largest window is about 3× edge-theta's (rows ≤ 9·10³).
- **Coefficients** (`run_energy.py`):
  - unmarked βw;
  - E_P: βχ_P^{−2}w, using only rows with (P,k)=1;
  - 𝒩: βX_c w with X_c = Σ_P c_Pχ_P^{−2}, and the **exact mask** applied by recomputing every row that has a slot-prime factor.

  Each of these is paired with a **random-phase control** with the same modulus: e^{iθ}w, e^{iθ}w1_{(n,P)=1}, or e^{iθ}|βX_c w|.
- **Mark designs** (all × N(P)^{−1/2}):
  - `p78` = ᾱ(P)²γ₂(P). This is Cor. marked-reflection with one active mark, η = 1, c_F = 1. Note χ_P(λ³)² ≡ 1.
  - `gauss` = γ₂(P), `ones` = 1, `kub` = conj γ₁(P), `rnd` = random.
  - Single slots [L₀, 2L₀) with L₀ = 30, 100, 300, 1000, i.e. 6/20/44/136 primes and geometric-mean L = 37, 149, 425, 1459.
  - A 2-slot design 31×61 (96 tuples, L ≈ 3600). Its weights are ᾱ(P₁P₂)²γ₂(P₁)γ₂(P₂)χ_{P₁}(P₂)²χ_{P₂}(P₁)² (the theta cocycle, mark-design §0), plus a random 2-slot control.
- **Symmetry [P].** For β and for every conjugation-symmetric design (p78, gauss, ones, kub, p78x2), |S_{k̄}| = |S_k| exactly (to 1e-13). So the family has K/2 independent members. This leaves the exponents unchanged, doubles Var ρ_full (z-rms √2), and was corrected for in the bias tests.

## 3. Single-mark energies E_P (step 2) [N]

| N(P) | 7 | 13 | 19 | 31 | 43 | 61 | 97 | 151 | 211 | 307 | 409 | 601 | 1009 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| b (θ) | −.528 | −.508 | −.518 | −.512 | −.506 | −.534 | −.469 | −.478 | −.512 | −.483 | −.465 | −.454 | −.494 |
| c (θ) | −.014 | +.002 | −.033 | −.037 | −.019 | +.010 | +.051 | +.040 | −.022 | −.009 | +.016 | +.043 | −.016 |
| θ/control rms | 1.03 | 0.96 | 1.02 | 1.08 | 0.99 | 1.01 | 0.99 | 0.97 | 0.90 | 1.01 | 0.95 | 0.99 | 0.96 |

Each entry has s.e. ≈ 0.03; the controls give b = −0.46 … −0.58 and c = −0.09 … +0.04. Pooled fits:

| fit | b | c (F) | d (N(P)) | n |
|---|---|---|---|---|
| θ, all P | −0.497 ± 0.009 | +0.001 ± 0.009 | −0.019 ± 0.005 | 845 |
| control | −0.510 ± 0.009 | −0.011 ± 0.008 | −0.005 ± 0.005 | 845 |
| θ, unproved window 1 ≤ F < N(P) only | −0.500 ± 0.009 | −0.006 ± 0.010 | −0.015 ± 0.006 | 721 |
| θ, N(P)/F ≥ 10 | −0.490 ± 0.011 | −0.001 ± 0.016 | −0.009 ± 0.011 | 429 |
| θ, N(P)/F ≥ 10 and N(P)² ≤ N | −0.514 ± 0.018 | −0.021 ± 0.028 | −0.026 ± 0.021 | 215 |
| per-mark E_P, all 44 primes in [300,600) | −0.502 ± 0.005 | −0.009 ± 0.005 | — | 2860 |

- The θ/control rms ratio over all P is 0.990 ± 0.012.
- The bias test z = ρ_full/(rms/√12) has mean +0.016 and rms 0.99 for θ, against −0.057 and 1.03 for the control.
- **The level-P² twist changes nothing.** The weak d ≈ −0.02 makes ρ *smaller* for larger P.

## 4. The coherent multi-mark core sum 𝒩 (step 3) [N]

Per slot and design, b ranges from −0.41 to −0.57 and c from −0.06 to +0.08, with s.e. ≈ 0.04–0.05 (`fits.out`). The θ/control ratio is 0.86–1.03 per design and **1.003 ± 0.010** pooled. Pooled fits:

| fit | b | c (F) | d (L) | n |
|---|---|---|---|---|
| P78-coherent (p78 + p78x2, all slots) | −0.483 ± 0.020 | +0.006 ± 0.020 | +0.019 ± 0.011 | 325 |
| its control | −0.489 ± 0.013 | +0.010 ± 0.013 | +0.006 ± 0.007 | 325 |
| P78-coherent, L/F ≥ 10 | −0.484 ± 0.021 | +0.002 ± 0.024 | +0.019 ± 0.015 | 241 |
| P78-coherent, L/F ≥ 10 and L² ≤ N (weighted-like) | −0.540 ± 0.049 | −0.039 ± 0.072 | −0.039 ± 0.058 | 63 |
| `ones` / `gauss` / `kub` / `rnd` pooled | −0.50 / −0.45 / −0.49 / −0.51 (± 0.02) | +0.00 / +0.04 / +0.02 / −0.02 | ≈ 0 | 260 each |
| unmarked β | −0.458 ± 0.036 | −0.000 ± 0.036 | — | 65 |
| random (4 seeds) | −0.497 ± 0.015 | −0.014 ± 0.015 | — | 260 |

311 of the 325 P78-coherent points lie in the unproved window 1 ≤ F < L.

Raw values at N₀ = 2·10⁵ (`rep_table.out`):

| F | rows | 𝒩 s300 p78: ρ_full / ρ_rms | its control | 𝒩 2-slot p78 |
|---|---|---|---|---|
| 1 | 28210 | −0.0098 / 0.022 | +0.0015 / 0.025 | +0.0089 / 0.032 |
| 4 | 7058 | +0.0038 / 0.055 | −0.0004 / 0.064 | −0.0432 / 0.064 |
| 16 | 1760 | +0.0132 / 0.089 | −0.0033 / 0.134 | −0.0590 / 0.115 |
| 64 | 443 | −0.0389 / 0.169 | −0.0354 / 0.182 | +0.2006 / 0.321 |

The large-sieve worst case at these windows is ρ ≈ F.

### 4.1 Bounds on a deterministic off-diagonal component

These are inverse-variance means of ρ_full (± 1 s.e.):

| family | F ∈ [1, 4] | F ∈ [5.7, 64] | F ≥ 90 |
|---|---|---|---|
| single marks | +0.0004 ± 0.0008 | +0.0016 ± 0.0018 | −0.018 ± 0.009 |
| their controls | −0.0000 ± 0.0008 | −0.0025 ± 0.0018 | −0.023 ± 0.010 |
| P78 designs | −0.0011 ± 0.0016 | −0.0065 ± 0.0036 | −0.015 ± 0.017 |
| all designs | −0.0012 ± 0.0007 | −0.0033 ± 0.0016 | −0.013 ± 0.008 |
| design controls | −0.0000 ± 0.0006 | +0.0005 ± 0.0013 | −0.011 ± 0.007 |

- The bias test for the designs, corrected for the k̄-symmetry, has z-rms 0.97–1.01. Its mean is −0.07 (nominal t = −2.5, but adjacent windows share half their rows and the designs share columns). Any offset is slightly negative, which is the favourable direction.
- The unmarked β at F ≥ 90 has a genuine negative offset, −0.18 ± 0.04: E lies below its diagonal in the reflection-provable regime.

### 4.2 (H2) Gram on the 44-mark slot (`h2.out`)

- **Cross-correlations.** The normalised cross-mark correlations satisfy rms·√n_eff = 0.99–1.10 for F ≥ 8 and 1.26–1.35 at F = 1. The excess is the finite-aspect-ratio term 1 + O(1/F) that random columns share.
- **Rayleigh quotient** of the P78 coherent vector (= 1 + ρ of that slot's 𝒩): 0.81–1.12, inside the spread of random weight vectors (1.00 ± 0.01–0.09).
- **Top eigenvalue:** at the Marchenko–Pastur edge.

### 4.3 Principal row k = 1 (`principal.out`)

|S₁|²/d₁ ≤ 5.9 over 22 designs × 5 values of N₀, which is the size of an exponential maximum. There is no Patterson-bias main term for β = γ₂ᾱ with any of these marks at these sizes.

### 4.4 Sector-count slope

The rms over S ∈ {6, 12, 24, 48} sectors at a fixed window gives b = −0.546 (single marks), −0.575 (P78 designs), −0.54 (controls) and −0.65 (unmarked β). This estimator has a known Jensen/kurtosis bias of about −0.04, visible in the controls. The θ families agree with their controls.

## 5. The rank-one pairing T (step 4)

### 5.1 Why ρ is not computed per row

ρ(m) = Res_{w=2/3} Σ_s c(s)g̃(s,−m)q_s^{−w} is not computable per row. In 04/vls/E, even s up to 2²² left the per-row noise above the signal. So I used the object itself. The probe vector is

A_{m,σ}(Y) = Σ_{s sf, s ∈ σ} γ₁(s) χ̄_s(−m) W(N(s)/Y),

with σ the ray class of s mod 4, so 12 vectors. γ₁ was checked against brute force. By Corollary D its pole part is Y^{−1/3}Ŵρ_σ(m). On the window rows its rank-one content is f = |A|²/noise − 1 = 0.15–0.28 at Y = 2¹⁸ (Y/K = 5–170), growing with Y as the pole law predicts.

### 5.2 Statistics

- **The pairing.** T(u) = Σ_k V u(k) S_k is the reflected form ⟨Ψ*u, Π_V b⟩ of wall4 §7.2.
- **Row classes.** It is evaluated in the 12 row ray classes mod 4, which absorb unknown ξ and ψ₁ factors.
- **R(u)** = Σ_τ|T_τ|² divided by its random-model size; R = 1 if u is unrelated to S.
- **Need.** The Cauchy–Schwarz-with-W level is R_CS = n_eff/C (the adversary u = S̄ attains it).
- **Cross-scale estimator.** λ = Re Σ T(A(Y₁))·conj T(A(Y₂)) / √(f₁f₂E₁E₂), over pairs Y₁ ≠ Y₂ in {2¹⁶, 2¹⁷, 2¹⁸}. These s-ranges are disjoint, so their noise is independent. λ estimates |T(ρ)|²/(random size) for ρ, and for π = ρζ when ζ = μᾱγ̄₂ (Lemma Z up to a ray character) is inserted.
- **Runs.** N₀ = 5·10⁴, 10⁵, 2·10⁵, with the same rows N(k) ≤ 5·10⁴ and F = 1 … 128. Designs: the 5 coherent 𝒩 designs, E_97, E_409 and β (θ group); controls are XR (`probe_summary.out`).

| test vector u | θ designs (144 cases) | controls (36 cases) |
|---|---|---|
| A(Y = 2¹⁸) | R = 0.992 ± 0.009 | 1.004 ± 0.016 |
| A·ζ (π-proxy) | R = 1.011 ± 0.010 | 1.003 ± 0.015 |
| ζ, β, β̄, μ, 1, α^{±1}, ᾱ², α³ | R ≈ 1 (0.4–2.4 individually) | ≈ 1 |
| **λ(ρ)**, cross-scale | **0.986 ± 0.031** (max 2.24) | 0.943 ± 0.051 |
| **λ(π = ρζ)**, cross-scale | **0.931 ± 0.033** (max 2.14) | 0.926 ± 0.049 |
| need: R_CS (adversary) | 6.5 → 198 ∝ rows^{0.975 ± 0.006} | same |

### 5.3 Positive control (`probe_power.out`)

I planted S′ = S + ε‖S‖·Ā_{σ₀}(2¹⁵)/‖A‖. The plant's s-range is disjoint from the estimator's, so it correlates with them only through ρ. λ_{σ₀} rises by f₁₅ε²/(1+ε²)·n_eff/C, as calibrated: by +0.8 at ε = 0.3 and +5.9 at ε = 1, starting from 0.4 at ε = 0.

### 5.4 Conclusion

|T(ρ)|² and |T(π)|² are at their random-model size, while the need is ≈ K/12 times larger. Hence |T|/Q^{11/12} ≈ K^{−0.49}, which is **Z^{−0.415} at M = 0.851**. This matches wall4 §7.2's random model (Q^{5/12} against Q^{11/12}).

## 6. Extrapolation to the weighted geometry (step 5) [H]

ρ ≈ Z^{bM + cw + dℓ}; W needs ρ ≪ Z^ε, so the margin is −(bM + cw + dℓ), with propagated s.e. The table uses the W(0.043) optimum (M, ℓ, w) = (0.8511, 0.1915, 0.0426); the W(0.040) point agrees to 0.001.

| family | margin (Z-exponent) |
|---|---|
| single marks, pooled | 0.427 ± 0.008 |
| single, unproved window only | 0.428 ± 0.008 |
| single, N(P)/F ≥ 10 / also N(P)² ≤ N | 0.419 ± 0.010 / 0.443 ± 0.015 |
| **𝒩, P78-coherent, all slots** | **0.407 ± 0.018** |
| 𝒩, P78-coherent, L/F ≥ 10 / weighted-like | 0.408 ± 0.019 / 0.469 ± 0.040 |
| controls (design / random) | 0.415 ± 0.012 / 0.423 ± 0.013 |
| unmarked β | 0.390 ± 0.032 |
| rank-one T (|T|/need) | 0.415 |

- Over the BRIEF range M ∈ [0.82, 0.86], ℓ ∈ [0.18, 0.20], the P78-coherent 𝒩 has margin 0.39–0.41 and the single marks 0.41–0.43.
- Equivalently, the "need/truth" for ρ is K^{1/2}, i.e. Z^{0.41}, consistent with edge-theta's Z^{0.43}.
- The terms cw and dℓ contribute < 0.01. **At 2σ, c ≤ 0.05 even inside the unproved window**, so a slope-type excess F^f with f ≳ 0.05 is excluded at these sizes, unless its constant is below the bounds of §4.1.

## 7. What the numerics can and cannot exclude

**Excluded at the tested sizes** (rows ≤ 2.8·10⁴ per window, N ≤ 4·10⁵, F ≤ 181, L ≤ 3.6·10³):
- any growth of ρ with F inside 1 ≤ F < L: c = 0.00 ± 0.01–0.02;
- dependence on L or on the level P²: |d| ≤ 0.04 at 2σ;
- any θ-specific excess over random phases: ratio 1.00 ± 0.01;
- a deterministic off-diagonal component above ≈ 0.2% for F ≤ 4, ≈ 1% for F ≤ 64, and ≈ 5% for F ≥ 90;
- a Patterson (NL)^{2/3}-type term on generic rows: if realised, our L² > N configurations (s1000, 2-slot) would show a (L²/N)^{1/3} ≈ 2–9× excess, and none appears;
- a Patterson-bias main term on the principal row, for β = γ₂ᾱ and these marks;
- (H2)-type cross-mark coherence beyond Marchenko–Pastur;
- rank-one alignment of the reflected rows with ρ or π beyond ≈ 1.05× the random level (2σ), against a need of 6.5–198×.

**Not excluded:**
- **Log factors.** (log K)^A·K^{−1/2} is indistinguishable from K^{−1/2}. This is harmless for W, which tolerates Z^ε.
- **Small secondary terms**, e.g. ε₀F^f with ε₀ ≲ 10⁻³, or terms with a later onset. A main term with f > 0 at small amplitude would eventually violate W as a power; no such term is visible.
- **The asymptotic law itself.** Taking it from K ≤ 3·10⁴ to K = Z^{0.85} is a hypothesis.
- **The Patterson (NL)^{2/3} term on principal rows** at much larger L. This is wall-spectral's bias. It is invisible here; it would only cap W at 19/22 and is harmless at 41/47.
- **Simplifications of P78's object:**
  - η = 1 and c_F = 1, with no ray-class restrictions, punctures or ν;
  - the core dyad only: squarefree rows (O = 0), no dual cube part (ℓ_b = 0), sin² windows instead of P78's tests;
  - ρ only through probe vectors with 15–28% rank-one content;
  - ζ only up to its ray character.
- **Statistical caveat.** The OLS errors treat windows as independent. Adjacent windows share half their rows, so the true s.e. is about 1.4× larger. No conclusion depends on this.

## 8. Files (all in this directory)

| file | content |
|---|---|
| `validate.py`, `check_lut.py`, `conv_check.py`, `test_engine.py`, `test_mask.py`, `test_gamma1.py` (+ `.out`) | step-1 checks |
| `nlib.py`, `psi.c` → `libpsi.so`, `cubic_data.c` → `cubic_data` | exact arithmetic, the ψ_k kernel, fast γ₂ data |
| `run_energy.py`, `runs.sh`, `logs/run_N*.log`, `data/energy_N*.npz` | row sums S_k and exact diagonals D_k (steps 2–3) |
| `analyze.py` → `analysis_windows.json`; `fits.py` → `fits.out`, `fits.json` | window and sector statistics, fits, bias tests, extrapolation |
| `h2.py` → `h2.out/json`; `principal.py` → `principal.out/json`; `check_conj.py`; `rep_table.out` | Gram/(H2) check, principal row, k̄-symmetry, raw table |
| `probe.c` → `libprobe.so`, `probe.py`, `runs_probe.sh`, `probe_N*.json`, `probe_summary.py` → `.out`, `probe_power.py` → `.out` | rank-one pairing (step 4) and its positive control |
