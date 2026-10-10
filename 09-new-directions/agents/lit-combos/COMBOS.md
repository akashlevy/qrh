# lit-combos — (B) cheapest combinations of partial improvements (r4 model)

## Setup and run log

**Model.** `model.py` is an unchanged copy of `03-six-sevenths/routes/r4/model.py`.
- `model_x.py` is the same model plus one knob, `drop_avg`, which removes the energy branch (M′+4ℓ′−2)/3 − w.
- `combo.py` is a batch driver. It caches the row-count table R(δ,x;S) per count-knob family. This is exact: it reproduces 0.8749568251160069 bit-for-bit.
- `analyze.py` builds the iso-θ thresholds by piecewise-linear interpolation.

**Baseline for every run.**
- `caseA="refined"`: the proved refined powerful count (R4 §2).
- `delta0=1e-3`: the detector floor 1/2+η, argued provable in R4 §6. It matters only for the count levers.
- All other lemmas are as in P78. The baseline is θ₀ = 0.8749568.

**Run log.**
- Data: `runs/batch1.jsonl` (199 configs), `runs/batch2.jsonl` (37 configs), `runs/batch2x.jsonl` (5 configs).
- That is 241 evaluations. Because of the R-caching, they cost about 70 uncached single-run equivalents. The 4-core machine was shared with other agents, which inflated the per-run time.

**Interpolation checks** (the interpolated threshold, then the θ found when re-run there):

| re-run at | θ found | target | result |
|---|---|---|---|
| wE = 0.0638 | 0.8709679 | 27/31 = 0.8709677 | at target (within 3e-7) |
| wE = 0.1326 | 0.8666665 | 13/15 | ✓ |
| wE = 0.285 | 0.8571382 | 6/7 | ✓ |
| gram_b = 0.85, wE = 0.0361 | 0.8709667 | 27/31 | ✓ |
| gram_b = 0.7, wE = 0.025 | 0.8695589 | 20/23 | ✓ |
| gram_b = 0, wE = 0.0128 | 0.8571445 | 6/7 | 1.6e-6 above (needs 0.013) |
| κ_A = 0.8944 | 0.870950 | 27/31 | ✓ |
| κ_A = 0.853 | 0.869553 | 20/23 | ✓ |
| κ_A = 0.95, wE = 0.0339 | 0.8709659 | 27/31 | ✓ |
| κ_A = 0.9, wE = 0.0247 | 0.869568 | 20/23 | 3e-6 above |
| κ_A = 0.9, wE = 0.2176 | 0.8571415 | 6/7 | ✓ |
| w1 = 0.1, wE = 0.0395 | 0.8709658 | 27/31 | ✓ |

Thresholds below are good to about ±0.0005 in the knob.

## Global facts (from the runs)

1. **Wall slope.** θ(wE) ≈ θ₀ − wE/16, linear up to w = 0.3. At every optimum the wall binds exactly: M + ℓ − 1 = wE.
2. **Floor sensitivity.** The floor δ0 is irrelevant for the wall: wE = 0.0638 gives the same θ at δ0 = 1/50.
3. **Count levers saturate.** κ_A ≤ 0.7, w1 ≥ 0.45 and perfect DH counts all give θ = 0.866826 at δ0 = 1e-3. That tends to 13/15 as δ0 → 0 (0.866683 at δ0 = 1e-4) and equals 167/192 = 0.869792 at δ0 = 1/50. So counts alone never beat 13/15, and with the paper's floor 51/100 they never beat 20/23.
4. **No synergy.** Iso-θ curves are almost straight. Write each requirement as a fraction of its single-lever threshold; the fractions in a working combination add to 1.00–1.07. Example: the triple wE 0.02 + gram 0.9 + κ 0.95 has fractions summing to 1.10 and gives 0.871263, which fails 27/31.
5. **Gram does not combine with counts.** With gram_b ≤ 0.5 the count improvements are worth nothing, because the Gram geometry is floor-bound.
6. **The (NL)^{2/3} branch is never binding.** Dropping the (M′+4ℓ′−2)/3 − w branch with `model_x drop_avg` leaves θ identical (to 3e-6) at wE = 0, 0.06, 0.15 and 0.285.

## Table 1 — energy wall W(w) needed, by Gram exponent factor gram_b (paper counts)

| gram_b | s₀ (pairing saving, see §S2) | 27/31 | 20/23 | 13/15 | 6/7 |
|---|---|---|---|---|---|
| 1 (paper) | 0 | **0.0638** | **0.0862** | **0.1326** | **0.285** |
| 0.85 | 0.0025 | 0.0361 | 0.0581 | 0.1038 | 0.254 |
| 0.70 | 0.0057 | 0.0056 | 0.0249 | 0.0663 | 0.211 |
| 0.55 | ≈0.011 | 0 | 0 | 0.0263 | 0.162 |
| 0.40 | 0.015 | 0 | 0 | 0 | 0.117 |
| 0.25 | 0.018 | 0 | 0 | 0 | 0.075 |
| 0 | 0.024 | 0 | 0 | 0 | 0.013 |

gram_b alone (wE = 0) reaches 27/31 at **0.674**, 20/23 at **0.600** and 13/15 at **0.455**. It cannot reach 6/7: gram_b = 0 gives 0.858137.

## Table 2 — wall needed, by amplification strength κ_A (M_κ), gram_b = 1

| κ_A | 27/31 | 20/23 | 13/15 | 6/7 |
|---|---|---|---|---|
| 1 (paper) | 0.0638 | 0.0862 | 0.1326 | 0.285 |
| 0.95 | 0.0339 | 0.0560 | 0.1017 | ≈0.252 |
| 0.90 | 0.0030 | 0.0247 | 0.0697 | 0.2176 |
| 0.85 | 0 | 0 | 0.0428 | 0.1887 |
| 0.80 | 0 | 0 | 0.0226 | 0.1671 |
| ≤ 0.70 (= DH counts) | 0 | 0 | 0.0024 (→ 0 as δ0 → 0) | 0.1454 |
| GRH for rows | 0 | 0 | 0 | 0.1429 = 1/7 |

κ_A alone reaches 27/31 at **0.895** and 20/23 at **0.853**. It reaches 13/15 only in the limit δ0 → 0.

## Table 3 — wall needed, by marked-moment extension w1 (κ_eq = 1/(1+w1))

| w1 | 27/31 | 20/23 | 13/15 | 6/7 |
|---|---|---|---|---|
| 0.10 | 0.0395 | 0.0616 | 0.1073 | > 0.18 |
| 0.20 | 0.0121 | 0.0283 | 0.0624 | > 0.18 |
| 0.30 | 0 | 0 | 0.0151 | 0.159 |
| 0.45 | 0 | 0 | 0.0024 | 0.1454 |

w1 alone reaches 27/31 at 0.2215 and 20/23 at 0.251. **It is dominated by κ_A:** w1 = 0.2215 corresponds to κ_eq = 0.82, against κ_A = 0.895 for the same target.

## Table 4 — Gram factor needed, by κ_A (wE = 0)

| κ_A | 27/31 | 20/23 | 13/15 |
|---|---|---|---|
| 1 | 0.674 | 0.600 | 0.455 |
| 0.95 | 0.754 | 0.665 | 0.488 |
| 0.90 | 0.956 | 0.712 | 0.493 |
| ≤ 0.85 | any | any (κ_A = 0.85 alone gives 0.869462) | 0.493 |

6/7 needs gram_b = 0 together with κ_A ≤ 8/9 and δ0 → 0. The model gives 0.857157 at δ0 = 1e-4, which is R4's vertex B.

## Table 5 — "diagonal-quality" package (no beyond-the-edge claim)

DIAG = Lemma marked at its true diagonal (mk = 1, no 2r+8z cap) plus Lemma plain at its diagonal (cP = 1).

| configuration | θ |
|---|---|
| DIAG alone | 0.872550 |
| mk ≥ 0.75 alone | 0.873485 |
| cP = 1 alone | 0.873453 |
| relaxing 2r+8z ≤ 3 to 2r+4z ≤ 3, alone | 0.874957 (no effect) |
| DIAG + wE = 0.03 | 0.870134 |
| DIAG + wE = 0.06 | 0.867739 |

With DIAG, the wall needed falls to **0.0196** for 27/31 and **0.0371** for 20/23. DIAG does not help the Gram lever: DIAG + gram_b 0.6 gives 0.869560.

## Precise statements at the most attractive points

Geometry is (l_x, l_y, ℓ) at the optimum; b = l_y − l_x; M = l_x + l_y; h = 1 − l_x + ℓ.

### S1. W(w), the energy wall
- **P78 location:** Lemma probe-row-norm (l.8111, eq. old-eq:5.4), via Lemma reflected-energy (l.7747, eq. old-eq:4.6) and Lemma hybrid-energy (l.7510, eq. old-eq:4.6b).
- **Needed statement:** Σ_{q_m≪Z^{M′}} |B^J_{m,σ}(Z)|² ≪ Z^{M′+ε} for the d = 0 block whenever M′ + ℓ′ ≤ 1 + w. Equivalently, in (old-eq:4.6), replace max(H, v+ℓ_b) by H for retained dyads with v + ℓ_b ≤ H + w. Equivalently, in (old-eq:4.6b) the column term NB of (K+NB) may be dropped for NB ≤ K·Z^w, but only for the theta-reflected coefficients β(n,b).
- **Meaning:** Lindelöf on average over Z^{M} sextic twists of the prime-marked completed cubic-theta row of length Z^{1+ℓ}, when the reflected dual has length Z^{2M+ℓ−1}, longer than the family by Z^{w}.

| target | w | (l_x, l_y, ℓ) | M | dual length |
|---|---|---|---|---|
| 27/31 | 0.0638 | (0.366, 0.494, 0.204) | 0.860 | Z^{0.924} |
| 20/23 | 0.0862 | (0.370, 0.499, 0.217) | — | — |
| 13/15 | 0.133 | (0.378, 0.510, 0.244) | — | — |
| 6/7 | 0.285 | (0.406, 0.546, 0.333) | — | — |

The d = ℓ block (M−2ℓ, 0) is never near the wall.

### S2. Gram exponent, through the theta pairing (V′)
- **P78 location:** Prop probe-gram (l.8340, eq. gram-quantitative).
- **The bound itself cannot be improved.** P_a^{1/6} is a true main term: the energy of the sextic-theta vector.
- **What is needed instead:** use the 04-gram-split decomposition A = A^sig + A^rest, where ‖A^rest‖² ≪ (Q/Y)Z^ε is already derived there. Then bound the pairing |Σ_m A^sig_m B_m| (π = ρ·ζ_m) by (Cauchy–Schwarz bound)·Z^{−s₀}. The relation is s₀ = (1−gram_b)·b/12.

| target | gram_b | b | (l_x, l_y, ℓ) | s₀ |
|---|---|---|---|---|
| 27/31 | 0.674 | 0.241 | (0.301, 0.542, 0.157) | 0.0065 |
| 20/23 | 0.600 | 0.305 | (0.274, 0.578, 0.148) | 0.010 |
| 13/15 | 0.455 | 0.300 | (0.277, 0.577, 0.145) | 0.014 |
| 6/7 (with wE = 0.013) | 0 | 0.289 | — | 0.024 |

- **Comparison with known bounds:** in LoA language, any θ_π < 1.31 starts to help (04-gram-split). The large sieve gives 3/2.

### S3. M_κ (κ_A)
- **P78 location:** Lemma inverse-amplification (l.12362) and eq. raw-moment (l.12353).
- **Needed statement:** Σ_{0<q_u≤H} |D^{−1/2} Σ_n μ(n)ν(n)χ_n(u)^{±}W(q_n/D)|² ≪ H(HD)^ε for H ≥ D^{κ} (paper: H ≥ D^{1+c}), for D = U^r, 1 ≤ r ≤ 3/2. This gives e(r) = max(1, (1+5κr)/6).
- **Values:** 27/31 needs κ = 0.895 (geometry (0.409, 0.409, 0.183), so b = 0 and Gram is inactive). 20/23 needs κ = 0.853.
- It is a family of rows 10–15% (in log) smaller than the Möbius length. Prop O1 (03 R1) shows the *positive dual* form of this is false; the statement itself is open.

### S4. w1
- **P78 location:** Lemma marked (l.9250).
- **Needed statement:** replace r + 2z ≤ m − c₁ with r + 2z ≤ (1+w1)m − c₁. At z = 0 this is M_κ with κ = 1/(1+w1) on the marked family.

### S5. DIAG
- **Lemma marked:** r + z ≤ m − c₁, with no 2r + 8z condition.
- **Lemma plain (old-eq:3.9):** n₁+n₂+z ≤ M instead of n₁+n₂+6κz ≤ M.
- Both are at (not beyond) the large-sieve diagonal for Möbius/prime coefficients, so they are not covered by the move-calculus obstruction. They are still unproved, since the paper's recursion loses factors 2 and 4.5 in z.

## Cheapest combinations

| target | cheapest single lever | best two-lever combinations |
|---|---|---|
| 27/31 | W(0.064); or κ = 0.895; or (V′) with s₀ ≈ 0.0065 | W(0.036) + s₀ 0.0025; W(0.034) + κ 0.95; **W(0.020) + DIAG** |
| 20/23 | W(0.086); or κ = 0.853; or s₀ ≈ 0.010 | W(0.025) + κ 0.90; W(0.025) + s₀ 0.0057; W(0.037) + DIAG |
| 13/15 | W(0.133); or s₀ ≈ 0.014 | W(0.043) + κ 0.85; W(0.026) + s₀ 0.012 |
| 6/7 | W(0.285) | full (V′) (s₀ ≈ 0.024) + W(0.013); full (V′) + κ ≤ 8/9 (δ0 → 0); DH counts + W(0.145); GRH rows + W(1/7) |
