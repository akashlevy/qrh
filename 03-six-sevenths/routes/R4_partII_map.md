# R4 — Part II map: which lemma strengthenings give exactly 6/7

Route R4 working notes. Scripts and outputs: `wf67/r4/`.
Paper = 7/8 preprint (`.../September-30-2026/build/paper.tex`); companion = 11/12 paper (`paper2.tex`).

## 0. Bottom line

1. **No natural strengthening that is provable with the paper's tools gives 6/7 inside Part II.** I found and checked two legitimate improvements: a refined count of the powerful row parts, and lowering the detector floor to 1/2+η. Neither moves the optimum, which stays at 0.874957. Every configuration that reaches 6/7 needs at least one mean-value input "below the diagonal", which nobody can currently prove.
2. **Structural fact (exact).** Part II always gives θ ≥ 1 − h/6, where h = 1 − l_x + ℓ is the log-length of the Poisson row range. Equality holds when the low side is loss-free. So **θ = 6/7 ⟺ the argument closes with row range h = 6/7.** The paper closes at h = 13/16. With the paper's low-side lemmas, h is capped at 4/5 even under GRH for the rows, which gives 13/15.
3. **The exact-6/7 vertices.** Every place where 6/7 appears as a clean rational optimum has h = 6/7:
   - (A) energy wall moved by w = 4/7, everything else exactly as in the paper — exact certificate;
   - (B) Gram bound without the P_a^{1/6} term, plus raw inverse moment for H ≥ D^{8/9} (a sharp threshold), plus floor 1/2+η;
   - (C) wall moved by 1/7, plus the refined powerful count, plus raw inverse moment for H ≥ D^{2/3} (which makes the counts perfect), plus floor 1/2+η;
   - GRH-for-rows versions of (B) and (C) — exact rational LP certificates.
4. **Weakest hypothesis found.** The refined count is provable. Combined with it, W(w) for w ≥ w* ≈ 0.285 already gives θ ≤ 6/7, but w* is not a clean value.
5. **Likeliest source of a colleague's exact 6/7 is not Part II.** It is more plausibly the companion's prime extraction θ = max((6+5κ)/12, 1 − κ/6). That formula equals 6/7 at κ = 6/7, and κ = 6/7 is its minimum (prior analysis). Its 1 − κ/6 branch is the same "1 − (amplifier length)" law as Part II's 1 − h/6.

## 1. Model (r4/model.py) and validation
LP per supply value S = ℓ/h. For fixed S, every constraint is linear in (l_x, l_y, σ0); I then minimize over S.

**Low side.** Prop probe-low, Prop probe-gram, Lemma probe-row-norm:

    σ0 ≥ 5/6 + l_x/6 − ℓ/6 − d + g(d)/2 + e(d)/2,   d ∈ {0, ℓ}

- Every piece is affine in d, so the two endpoints suffice.
- g = max(0, b/6, 2b − l_y + d). For b ≤ 0, Lemma balanced-additive-norm gives g = 0.
- e = the Lemma reflected-energy exponent, maximized over block data (§2).

**High side.** Eq. common-high-exponent at d = h, with the slope argument for d < h:

    σ0 ≥ (1+δ)(1−l_y)/2 − ℓ/2 + xδℓ + h(R(δ,x) + δ/2 − 1/6)

- R(δ,x) = min_t max_r min{amplified inverse e(r) − δr, marked 1 − δr − 2qz_I, plain x2 with slots, plain x2}.
- Capacities are supply-limited (the model handles S < 7/37 directly).

**Floor bin.** δ0 = 1/50, R = 1, q ≤ δ0/2.

The intermediate-row check (d ≤ 1/2, t = 1) is dropped. It is a convenience, not a constraint: for d < h, selecting slots gives R_d ≤ R_h, and the slope R + δ/2 − 17/50 > 0. The older `sensitivity.py` imposed this check plus the closed-form supply condition, which is why its "Gram exponent 0" value 0.86398 differs from mine (0.85980).

**Validation.**
- Baseline 0.8749568, matching the prior barrier 0.874957.
- GRH-for-rows gives 13/15 exactly.
- R matches the paper's closed form R_* to 1e-9.
- Exact checker `r4/exact_cert.py` (Fractions, branch and bound) passes at the paper geometry for σ0 = 0.8755.

## 2. Energy LP: exact closed form (r4/energy_lp3.py, verified at 6000 random points)

    e(M',ℓ';w) = max(0, CA, M'+ℓ'−1−w, (M'+4ℓ'−2)/3 − w)

- Paper count of frozen powerful row parts, Z^{O/2}: CA = (1+3ℓ'−2M')/4.
- Refined count Z^{O/6 + 2A_nz/3 − A_4/3}: CA = (1+2ℓ'−2M')/3.

**The case-A extremal configuration is an artifact of the crude count.** With the paper count, the case-A maximum has A_0 = 0, meaning no active non-slot prime in the powerful row part. Then every powerful prime has row valuation ≡ 0 (mod 6), so k_pow is a sixth power. There are only Z^{O/6} of those, not Z^{O/2}.

**Refined count lemma** (elementary, proved here). Write the powerful part as k_pow = λ^6 κ with κ sixth-power-free and κ = κ1 κ2² κ3³ κ4⁴ κ5⁵. Since k_pow is powerful, rad(κ1) | λ. Put A = log-norm of rad(κ) (these primes have j ≠ 0, so they are active and A ≤ A_0) and A_4 = log-norm of κ4. Then

    #{k_pow : q ≍ Z^O, type (A, A_4)} ≪ Z^{A + (O − 2(A−A_4) − 4A_4)/6 + ε} = Z^{O/6 + 2A/3 − A_4/3 + ε}

with O((log Z)^5) dyadic types. Insert this in place of O/2 in the aggregate step of Lemma reflected-energy (last paragraph of its proof).

**Payoff.** Zero on its own: case A is not active at the paper optimum. It becomes decisive once the wall is relaxed.

## 3. Single replacements (float optimum; singles_out.json)

| replacement | θ | provable? |
|---|---|---|
| none (paper) | 0.874957 | — |
| refined powerful count (case A → (1+2ℓ'−2M')/3) | 0.874957 | YES (§2) |
| detector floor 1/2+η | 0.874957 | YES (region one of Lemma local-euler needs only x_r > 1/2; detector and saturation fine) |
| plain capacity c_P = 7/30, 1/3, ∞ | 0.874935, 0.874708, 0.873453 | no / no |
| Lemma marked width r+2z < 1+w1, w1 = 0.1, 0.3 | 0.873488, 0.869792 | no (below diagonal at z=0) |
| amplification with M_κ, κ_A = 6/7 | 0.869792 = 167/192 (floor-bound) | no |
| perfect counts R = 1−δ | 167/192; with floor 1/2+η → 13/15 | no (DH for rows) |
| GRH for rows | 13/15 | no |
| Gram without the 2b−l_y term | 0.874957 | — |
| Gram without the P_a^{1/6} term | 0.859796 | no (genuine Gram main term, §6) |
| wall W(w): w = 0.1, 1/7, 0.2, 0.3 | 0.868705, 0.866025, 0.863333, 0.861667 | no |
| **wall W(4/7)** | **6/7 exactly** | no |
| wall W(∞) | 0.849583 | no |

For w ≥ 0.2 the paper-count W(w) value is exactly 13/15 − w/60. This is the low-side-only barrier. Multipliers 17/30, 2/5, 1/30 on [low, g=b/6], [low, b/6 + caseA], [low, b/6 + wall] at d = 0; valid for every geometry.

## 4. Combination sweep (sweep1_out.json)
Rows: count variants. Entries: θ at wE = 0, 1/14, 1/10, 1/7, 1/5, 1/4, 2/7, 3/10, 1/3, 1/2, ∞.
Marks: `*` = 6/7 to 3e-6; `<` = below 6/7.

    paper counts      pap gb=1: .87496 .87049 .86870 .86603 .86333 .86250 .86190 .86167 .86111 .85833 .84958<
                      ref gb=1: .87496 .87049 .86870 .86603 .86245 .85933 .85709< .85620< .85412< .84373< .83750<
                      gb=0     : .85980 .85479< ...                                           (case A irrelevant)
    GRH for rows      pap gb=1: .86667 .86548 .86500 .86429 .86333 .86250 .86190 .86167 .86111 .85833 .84958<   (= 13/15 − w/60)
                      ref gb=1: .86667 .86190 .86000 .85714* .85333< ...                       (= 13/15 − w/15 for w ≤ 1/3)
                      gb=0     : .85714* .85204< ...
    perfect, δ0=1e-4: like GRH up to O(δ0);  kA=6/7, δ0=1e-4: gb=0 → .85716 (→ 6/7), ref w=1/7 → .86037

Targeted runs (targeted1/2_out.json), all with floor δ0 → 0 unless noted:
- **Gram without P_a^{1/6}**, by amplification strength κ_A:

  | κ_A | θ |
  |---|---|
  | 1 (paper) | 0.858137 |
  | 0.95 | 0.857453 |
  | 0.90 | 0.857153 |
  | 8/9 | 0.8571429 (= 6/7) |
  | ≤ 8/9 | 6/7 |

  So the threshold is **κ_A = 8/9 exactly**. Reason: at x = 0 the count slope near δ = 0 is c(0) = 2/(3κ_A), and vertex B needs c(0) ≥ 3/4.
- **W(1/7) + refined**: perfect counts give 6/7; κ_A = 0.7 gives 0.85724; κ_A = 0.6 gives 6/7. The threshold is κ_A ≤ 2/3: then t = 3/2 forces r ≥ 1 and R = 1 − δ exactly.
- **W(4/7), paper case A, paper counts, floor 1/50**: 0.857142857 at (4/7, 4/7, 3/7).

## 5. The exact 6/7 fingerprints (all have h = 6/7, θ = 1 − h/6)

**Structural identity.** Under the geometry relation h = 1 − l_x + ℓ, the d = 0 low constraint reads

    σ0 ≥ 1 − h/6 + (g + e)/2

so θ ≥ 1 − h/6 always. Under GRH for rows, the high side reads

    σ0 ≥ (1 − l_y)/2 − ℓ/2 + 5h/6

and the two balance at h = (1 + l_y + ℓ)/2. The low-side lemmas (the wall M+ℓ ≤ 1+w, Gram, case A) limit how large h can be:

| low-side lemma set | max h | θ = 1 − h/6 |
|---|---|---|
| paper | 4/5 | 13/15 |
| Gram without b/6 (2b ≤ l_y) | 6/7 | 6/7 |
| refined count, wall at 1+w | 4/5 + 2w/5 (= 6/7 at w = 1/7) | 6/7 at w = 1/7 |
| paper case A, wall at 1+w | 4/5 + w/10 (= 6/7 at w = 4/7) | 6/7 at w = 4/7 |

| tag | lemma changes | vertex (l_x, l_y, ℓ) | status |
|---|---|---|---|
| A | W(4/7) only (Lemma reflected-energy case B: no loss while T_d ≤ H + 4/7) | (4/7, 4/7, 3/7) | EXACT certificate: exact_cert.py, 109 boxes, all PASS at σ0 = 6/7. High side slack ≈ 0.09 with the paper's own counts. Lower bound 13/15 − w/60 holds for every geometry. Control w = 3/7 fails. |
| B | Gram without P_a^{1/6}; raw inverse moment for H ≥ D^{8/9}; floor 1/2+η | (2/7, 4/7, 1/7) | GRH version: exact LP certificate (2/3, 2/21, 1/21, 4/21). Non-GRH: θ → 6/7 as η → 0. Control with κ_A = 1 fails at δ ≈ 0.013, x = 0. |
| C | W(1/7); refined count; raw inverse moment for H ≥ D^{2/3}; floor 1/2+η | (3/7, 3/7, 2/7) | GRH version: exact LP certificate (4/15, 2/15, 2/5, 1/5). Non-GRH: analytic — R = 1 − δ gives E(h) = 6/7 − σ0 + δ(2x−1)/7 ≤ 6/7 − σ0, floor = 6/7 − σ0 + 6δ0/7. Fails with the paper's case A (low = 6/7 + 1/56). |

The GRH certificates are in `r4/grh_cert_verify.py`, pure Fractions; the multipliers on [low g0+e0], [low g0+wall], [low gq+wall], [high GRH] (B) and [low g0+e0], [low g0+wall(8/7)], [low gb6+e0], [high GRH] (C) give exactly s ≥ 6/7.

## 6. Provability, at the exact step of the paper that produces each constraint

- **Case-A term.** Source: the aggregate count Z^{O/2} in the proof of Lemma reflected-energy ("Every powerful ideal is uniquely x²y³...") combined with the branch s_hyb ≠ z_a in Lemma probe-row-norm (eq. old-eq:5.3). This is a TECHNICAL loss and is removed by §2. The remaining (1+2ℓ'−2M')/3 comes from O = 0 with a pure-cube dual (n = 1, 3ℓ_b = T_d). There the active mark sum Σ_P a(P) q_P^{−1/2} does not oscillate (χ_P(b³)^{−2} = 1), so the hybrid sieve's L term is attained. It is plausibly genuine.
- **Wall M+ℓ ≤ 1.** Source: in Lemma probe-row-norm, max(H, v+ℓ_b) = H uses T_d ≤ H + (M'+ℓ'−1). The excess is the column term NB of the quadratic large sieve inside Lemma hybrid-energy (Goldmakher–Louvel / Heath-Brown K + N), which is sharp for general coefficients.
  - The unreflected sextic large sieve (BGL; K + N + (KN)^{2/3}) gives loss max(1+ℓ'−M', (2+2ℓ'−M')/3). That is no better unless 4M'+ℓ' > 5; it is equal at vertex A.
  - Lemma canonical-moment yields only the column count Z^{F_0}, not the row diagonal.
  - So W(w) for any w > 0 is a below-diagonal ("Lindelöf on average over fewer rows than columns") mean value for reflected cubic-theta sums. This is a GENUINE barrier for the method.
- **Gram P_a^{1/6}.** Source: the exceptional frequencies (k) = a_1^6 e in Prop probe-gram. For k = a^6 one has F(n1,n2;k) = R(n1,n2)·1_{(n1n2,a)=1}, with no oscillation in the residual columns. So the Gram sum genuinely contains ≍ (Q/Y')P_a^{1/6}, and the Gram bound is sharp. Avoiding it would need a non-Cauchy–Schwarz treatment of the sixth-power-frequency component of A_m. Nothing in either paper does this; open. Fingerprint B needs M_{8/9} in addition anyway.
- **Amplification e(r) and raw-moment.** Source: Lemma marked at z = 0 needs r < 1 (H ≥ D^{1+c}). Going below is M_κ, the "theorem in mean-square form" (prior analysis; companion's dual target false for H < D). GENUINE barrier.
- **Count slope near δ = 0.** The paper has c(0) = 2/3: inverse 2nd moment against plain 4th moment, Montgomery-type. Vertex B needs c(0) ≥ 3/4. This requires M_{8/9}, or a plain sixth moment for polynomials of length ≤ 1/3 (not in the paper; the sextic large sieve's (UN)^{2/3} term blocks it).
- **Floor 51/100 → 1/2+η.** PROVABLE: Lemma local-euler region one uses x_r ≥ 51/100 only through 3/2 − 3x_r < 0 and 4 − 6x_r − 6z_r < −1. The detector Gamma integral and the saturation constants depend only on η.

## 6b. Rigor of the refined energy bound (r4/energy_dual_cert.py)
For each of the six LP branches (u ∈ {v, z_a, (v+z_a)/3} × max ∈ {H, dual}), dual vectors y ≥ 0 with Aᵀy ≥ c were checked EXACTLY in Fractions. By weak duality each gives a bound for every (M′, ℓ′, w). The certified branch bounds (energy loss E_ref − M′) are:

| branch (u, max) | certified bound |
|---|---|
| (v, H) | (1+2ℓ′−2M′)/3 |
| (z_a, H) | 0 |
| (z_a, dual) | M′+ℓ′−1−w |
| (avg, dual) | (M′+4ℓ′−2)/3 − w |
| (avg, H) | (1−2M′+ℓ′)/3 ≤ (1+2ℓ′−2M′)/3 |
| (v, dual) | −M′/2 + ℓ′ − w, which is ≤ (1+2ℓ′−2M′)/3 when M′+2ℓ′ ≤ 2+6w and ≤ M′+ℓ′−1−w when M′ ≥ 2/3; so dominated whenever ℓ′ ≤ 2/3 |

Hence, rigorously: E_ref − M′ ≤ max(0, (1+2ℓ′−2M′)/3, M′+ℓ′−1−w, (M′+4ℓ′−2)/3 − w) for ℓ′ ≤ 2/3. With w = 0 this is the paper's energy statement, with case A improved.

## 6c. Further checks
- **Lemma marked relaxed to the true canonical diagonal** (slots counted once: z ≤ 1 − r): θ = 0.873485. Also relaxing 2r+8z ≤ 3 to 2r+4z ≤ 3 gives 0.873159. Combined with the Gram fix and floor → 1/2: 0.858137, the same as without it. Marked slots cannot help the x = 0 slope (q = 0), which is what vertex B needs.
- **Thresholds**, all with paper counts and floor 1/50 (wstar_*.json):
  - refined count + W(w): θ = 6/7 at **w\* = 0.28493** (geometry ≈ (0.4054, 0.5464, 0.3331), h = 0.928);
  - Gram fix + W(w): w\* = 0.03785.
- **Exact certificates** (exact_cert.py, all checks: low, floor, small rows, high-side branch and bound, d ≤ h + ζ):

  | hypotheses | geometry (l_x, l_y, ℓ) | σ0 | result |
  |---|---|---|---|
  | (A) W(4/7), paper case A, paper counts | (4/7, 4/7, 3/7) | 6/7 | ALL PASS (109 boxes) |
  | control: W(3/7), otherwise A | (4/7, 4/7, 3/7) | 6/7 | low fails (0.929) |
  | (B) Gram fix, κ_A = 8/9, δ0 = 1/1000 | (2/7, 4/7, 1/7) | 6/7 + 1/1000 | ALL PASS (76379 boxes) |
  | control: B with κ_A = 1 | (2/7, 4/7, 1/7) | 6/7 + 1/1000 | high fails at δ ≈ 0.013, x = 0 |
  | (C) W(1/7), refined count, κ_A = 2/3, δ0 = 1/1000 | (3/7, 3/7, 2/7) | 6/7 + 1/1000 | ALL PASS (95435 boxes) |
  | control: C with paper case A | (3/7, 3/7, 2/7) | 6/7 + 1/1000 | low fails (0.875 = 6/7 + 1/56) |
  | **W(3/10) + refined count** (provable), paper counts, floor 1/50 | (51/125, 11/20, 171/500) | **6/7** | ALL PASS (17555 boxes); low = 5137/6000 = 0.856167 |
  | control: W(3/10) with paper case A | (51/125, 11/20, 171/500) | 6/7 | low fails (0.8699) |
  | **W(2/7) + refined count** (provable), paper counts, floor 1/50, ζ = 10⁻⁶ | (11489/28328, 14840/27149, 32673/97961) | **6/7** | ALL PASS (80235 boxes); low = 0.857094 |

  The checker's sanity run passes at the paper geometry for σ0 = 0.8755 (74087 boxes). Its t-choice uses a candidate-exact float search, and each R_t is then recomputed in Fractions.

## 7. Verdict for R4: OBSTRUCTION, plus conditional results
- Unconditional: Part II with provable strengthenings stays at 0.874957. Both strengthenings found (the refined powerful count and floor 1/2+η) are inactive at the optimum.
- Conditional, exact: W(4/7) ⇒ Re s > 6/7 (vertex A).
- Weakest single hypothesis, exactly certified: **W(2/7)** (case-B wall moved to M′+ℓ′ ≤ 9/7), plus the PROVEN refined count, plus every other lemma as in the paper ⇒ Re s > 6/7. The float threshold is w* = 0.28493; at w = 2/7 the optimum is θ = 0.857094.
- Single hypothesis giving EXACTLY 6/7 as a clean optimum: W(4/7) with the paper's crude count. Its lower bound 13/15 − w/60 holds for every geometry.
- The exact-6/7 numerology in Part II is the law θ = 1 − h/6 with h = 6/7. Every way of reaching h = 6/7 needs a below-diagonal mean value (the theta energy past the wall, or M_κ with κ ≤ 8/9 together with a non-Cauchy–Schwarz Gram treatment).
