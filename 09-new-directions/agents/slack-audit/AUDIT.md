# slack-audit: within-the-edge slack in P78 Part II

Agent `slack-audit`, 10 Oct 2026. These are working notes; the final report is the hand-back message.

**Files**
- `model_sa.py`: copy of r4 `model.py`.
  - With default knobs it reproduces 0.8749568251160069 bit-for-bit.
  - Extra knobs:
    - `k3`: plain sixth moment;
    - `interm`: paper intermediate rows; `True` = literal (t = 1), `"R0"` = t optimized;
    - `tmin`, `drop_dl`, `side`;
    - `theta_sc`: self-consistent cP = 1/(6κ), κ = max(κ_min, 2θ − 1), and dmax = 2θ − 1;
    - labelled LP duals (`_duals`).
- `batch.py`: parallel driver. Results are in `runs/b1–b4.jsonl`, with a digest in `runs/SUMMARY.txt`.
- `diag_closed.py`: closed form of the DIAG optimum.
- `frontier.py`, `frontier2.py`, `slackmap.py`: bin-by-bin slack at a fixed geometry. Outputs are `runs/frontier*.jsonl`.

**Notation**
- DIAG = Lemma marked at r + z ≤ 1 with no 2r + 8z cap (mk = 1, w2scale = 10), plus Lemma plain at 2m + z ≤ 1 (cP = 1).
- Unless stated otherwise, every run also includes the provable refined case-A count and δ0 = 10⁻³.
- In lemma language: knob mk means r + z/mk ≤ 1 (paper mk = 1/2), and knob cP means 2m + z/cP ≤ 1 (paper cP = 2/9, i.e. 2m + 4.5z ≤ 1).

## 0. Bottom line

**No.** A package of type-(A) improvements cannot reach 27/31 in P78 Part II.
- The best type-(A)-only value is **θ_A = 0.8725498**.
  - Fine grid (dgrid 641, xgrid 21): 0.8725499.
  - Closed form: 0.8725499.
- DIAG attains it, provided the intermediate-row fix of §2.1 is made.
- Every other within-the-edge item has payoff 0, both alone and on top of DIAG; see the table.

**Why DIAG stops at 0.87255.**
- The DIAG optimum is (l_x, l_y, ℓ) = (0.41177, 0.41177, 0.17647), with h = 0.76470 and S = 0.2308.
- Low side: b = 0 and M + ℓ = 1.0000. Hence g = e = 0, and the low side is the structural law θ = 1 − h/6, which is loss-free.
- High side: the single binding bin is δ = 0.329, **x = 0**.
- LP duals: low 0.730, high (δ 0.329, x 0) 0.175, wall 0.095. Nothing else binds.

**Slack in R at the DIAG geometry** (`slackmap.py`):

| x | 0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|---|
| slack | 0 | 0.0030 | 0.0058 | 0.0085 | 0.0109 | 0.0126 |

So DIAG over-delivers on rows with prime amplitude (§3). The binding rows are those whose slot sums cancel to square-root size; for them no slot lemma helps.

**The count at x = 0.** It is the slot-free count
R0(δ) = 1 − 2δα/(3α − δ), α = 5/6
(the long-inverse closed form). It has three inputs, each already at its natural diagonal range:
- the raw inverse moment for r < 1;
- sixth-power amplification for r > 1, with the minimal family H = D, which is the optimal choice within the diagonal;
- the plain fourth moment, 2m ≤ 1, at any lengths.

Within the witness calculus nothing diagonal beats it. At the crossing (r, m) = (2t/3, t/3), every sub-product M^a S^b of length ≤ 1 has gain ≤ 2δt/3:
- S³ requires 3m = t ≤ 1;
- M·S requires r + m = t ≤ 1;
- the adversary escapes both by taking t slightly above 1.
- The pure long witness (r = t, m = 0) admits only copies.

So beating R0 is "Ingham below 7/10" (B), as `long-inverse` found.

**Closed form** (`diag_closed.py`):
- θ_A = min over ℓ of max( 1 − (1+3ℓ)/12 , max over δ of [ (1+δ)(1+ℓ)/4 − ℓ/2 + ((1+3ℓ)/2)(R0(δ) + δ/2 − 1/6) ] ).
- This equals 0.8725499 at ℓ = 0.176467, matching the LP to 1e-7.
- Reaching 27/31 needs h = 24/31 (ℓ = 17/93). There the x = 0 bin gives 0.87857 > 27/31.

## 1. Audit table

Notes on the columns:
- Proved = what P78 proves. Natural = what the same method gives without lossy steps (the diagonal).
- Class: (A) within the diagonal or technical; (B) beyond the edge; (C) already optimal for the method.
- Payoff = model θ for the item alone, all else as in the paper (baseline 0.8749568), and the same with DIAG (0.8725498).

| # | item (P78 line) | proved | natural | class | alone | with DIAG |
|---|---|---|---|---|---|---|
| 1 | Lemma marked (9250) | r+2z ≤ 1, 2r+8z ≤ 3 | r+z ≤ 1 | A (DIAG) | mk=1: 0.873485 | part of DIAG |
| 2 | Lemma plain (12531, eq. 3.9) | 2m+6κz ≤ 1, κ = 3/4 | 2m+z ≤ 1 | A (DIAG) | cP=1: 0.873453 | part of DIAG |
| 3 | Intermediate rows d ≤ 1/2 (eq. adaptive-intermediate, 16161) | t = 1, no slots, R = 1 − 2δ/3, q = δ/2 | t optimized, R = R0(δ); this is Prop detector-counts with zero slots, no supply needed | A (trivial) | 0 | **required**: literal 0.873815; fixed 0.8725498 |
| 4 | Case A, Lemma probe-row-norm (8111, l. 8299) | Z^{O/2} | refined Z^{O/6+2A/3}, proved in R4 | A (proved) | 0 (≤ 1e-8) | 0 (inactive) |
| 5 | Floor 51/100 (eq. local-region-one 4055; floor-bound 15846) | δ0 = 1/50 | 1/2 + η, proved in R4 | A (proved) | 0 | 0 |
| 6 | Plain sixth moment (new; three copies of S_m) | none | 3m+z ≤ 1 | A (unproved) | 0 | 0 |
| 7 | Dynamic κ in Lemma plain | κ ∈ [3/4, 1] | κ = 2θ − 1 < 3/4 | A (model) | < 1e-7 | 0 (cP = 1) |
| 8 | Bin ceiling δ ≤ 2β* − 1 (15096) | model uses 5/6 | 2θ − 1 | model fix | < 1e-7 | < 1e-7 |
| 9 | Cubic large sieve (NL)^{2/3} term (7491) | HB | N + L | A (open conjecture) | 0 | 0 (branch never binds) |
| 10 | Gram P_a²/Y′ term (8340) | present | absent | A | 0 | 0 |
| 11 | Detector t ∈ [1, 3/2] (4510) | — | t < 1 and t > 3/2 useless | C | tmin = 0.5: 0 | 0 |
| 12 | Detector m ≤ 1/2 (eq. old-6.3) | functional equation | same | C | — | — |
| 13 | Zero bins, Euler regions (4281, 4010) | buffered loss O(e) | same | C | — | — |
|  | Region two (x_r ≥ 7/8) | — | extends trivially to x_r ≥ θ: closest exponent 4 − 6x − 6z < −1 for x > 0.67 | C | — | — |
| 14 | Bin-prime bound (15015) | q ≤ δ/2, from the zero-free strip | same | C | — | — |
| 15 | Reflected numerator U^{δ/2} at w = 1 − a | — | x + w ≥ 1 is forced by region one; w = 1 − a is optimal when U < Y² | C | — | — |
| 16 | Raw inverse moment, eq. raw-moment (12353) and P2 Prop thm:ms | H ≥ D^{1+c} | diagonal | C | — | — |
| 17 | Sixth-power amplification (12362) | e = (1+5r)/6 | best within diagonal (H = D minimal, P = (H/U)^{1/6} copies) | C (B below) | — | — |
| 18 | Local compensation, error slots (eq. old-5.13e, 8883) | error slots carry no amplitude | same | C | — | — |
| 19 | Low tuple sum f(d) = −d (eq. old-5.5, 8632) and the d = ℓ branch | — | — | C | never binds | never binds |
| 20 | Energy wall M+ℓ ≤ 1 (7747/8111) | wall | — | B | W(0.0638) → 27/31 | W(0.0196) → 27/31 |
| 21 | Gram P_a^{1/6} (8340) | genuine main term | — | B (V′) | gram_b 0.674 | 0.673 (DIAG adds nothing) |
| 22 | Amplification below diagonal (M_κ) | κ = 1 | — | B | κ_A 0.895 | **κ_A ≈ 0.964** |
| 23 | Marked moment past r < 1 at z = 0 (w1) | r < 1 | — | B | w1 0.2215 | ≈ 0.222 (DIAG adds nothing) |
| 24 | DH counts | — | — | B | 0.866826 (δ0 1e-3) | 0.866826 |
| 25 | Cauchy–Schwarz in m (8564) | — | — | B (V′) | — | — |

## 2. Model versus paper (discrepancies)

1. **Intermediate rows (looser than the paper as written; matters under DIAG).**
   - What r4 does: it drops eq. (adaptive-intermediate), arguing via R_d ≤ R_h and a positive slope.
   - Baseline: the paper's literal constraint (t = 1, no slots, q = δ/2, z-contour 17/50) is slack (θ unchanged).
   - DIAG geometry: l_y = 0.41, and the literal constraint gives 0.66845 + 0.299δ, which exceeds 0.8725 for δ ≳ 0.68. The DIAG value then becomes 0.873815.
   - Fix: at d ≤ 1/2, use Prop detector-counts with zero selected slots at t = α/(α − δ/3). This gives R0(δ) with no supply hypothesis.
   - With the fix, θ is restored to 0.8725498 exactly (`DIAG+intermR0`).
   - Alternatively, select slots at d ≤ 1/2 with a mesh refined by 1/d_min.
2. **dmax.** The model uses 5/6; the true ceiling is 2θ − 1. The model is stricter, but by less than 1e-7.
3. **cP.** The model uses 2/9. Lemma plain with κ = 2θ − 1 < 3/4 would need the lemma below its stated κ ≥ 3/4; the proof seems to use only 6κ − 1 ≥ 7/2 margins. Payoff is below 1e-7.
4. **Count R.**
   - The model takes the exact max over r of the min over all lemma bounds. This is at most the paper's crossing bound and equals the paper's closed form R_* at the baseline (R4 check).
   - Supply-limited capacities are more general than the paper's 7/37 hypothesis.
   - The pieces fA, fI and fP, the range x ∈ [0, 1/2], and m = t − r are exact. The adversary's best is uniform amplitude and minimal m.
5. **Low side.**
   - Two branches d ∈ {0, ℓ} suffice, since every piece is a max of affine functions of d.
   - The side conditions never bind (no LP duals).
   - For b ≤ 0 the balanced lemma (3588) gives g = 0, as encoded.
   - Small-row margin −l_y/2 + h(17/50 − 1/6) + 2d_min < 0: this holds at every optimum found (−0.053 at DIAG).
6. **Grid.** dgrid 161 and xgrid 11 understate θ by about 1e-7.
7. **Verdict.** No modelling error understates the paper's strength by more than 1e-7.

## 3. Cheapest partial DIAG reaching 0.8725498

These are min-slack values at the DIAG geometry (`frontier2.py`), confirmed by full LPs: P1 = (0.7, 0.75, cap kept) and P2 = (0.8, 0.6, no cap) both give 0.8725498.

| Lemma marked | Lemma plain needed |
|---|---|
| r + 1.33z ≤ 1 (mk 0.75), 2r+8z ≤ 3 cap kept | 2m + 1.43z ≤ 1 (cP 0.70) |
| r + 1.43z ≤ 1 (mk 0.70), cap kept | 2m + 1.33z ≤ 1 (cP 0.75) |
| r + 1.25z ≤ 1 (mk 0.80), no cap | 2m + 1.82z ≤ 1 (cP 0.55) |
| r + 1.29z ≤ 1 (mk 0.775), no cap | 2m + 1.67z ≤ 1 (cP 0.60) |
| r + 1.33z ≤ 1 (mk 0.75), no cap | 2m + 1.54z ≤ 1 (cP 0.65) |

With mk ≤ 0.675 or cP ≤ 0.5, the DIAG value is not reached.

## 4. Residual (B) requirement for 27/31 on top of DIAG (any one)

| lever | needed with DIAG | alone (for comparison) |
|---|---|---|
| M_κ: Σ_{q_u≤H}\|M_u(D)\|² ≪ H^{1+ε} for H ≥ D^κ | κ ≈ 0.964 (0.9627 gives 0.870928) | 0.895 |
| energy wall W(w) | w = 0.0196 (gives 0.8709686) | 0.0638 |
| Gram factor (V′) | 0.673 (0.674 gives 0.870982) | 0.674 |
| marked moment past r < 1 (w1) | ≈ 0.222 | 0.2215 |

Notes on M_κ:
- Binding data for DIAG + κ_A = 0.9627: geometry (0.4085, 0.4085, 0.1830); binding bin x = 0, δ = 0.27, t* ≈ 1.17.
- So the needed instance is D = U^{1.17}, H = D^{0.963} ≈ U^{1.13}.
- This is still beyond the edge, but it is the smallest such requirement found: the family is 3.6% (in log) shorter than the polynomial.
- With partial DIAG P1, κ_A = 0.9627 gives 0.870997, so 27/31 needs κ_A ≈ 0.962.
