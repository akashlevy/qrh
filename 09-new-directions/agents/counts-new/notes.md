# counts-new: new count mechanisms for sextic-twist rows at sigma in [0.63, 0.71]

Agent `counts-new`, 10 Oct 2026. Outcome: **negative**. No mechanism found that lowers the row count
R(delta, x) at delta in [0.27, 0.42] by proved means. Below: the precise obstruction, a sharper statement of what
*would* suffice (five beyond-edge targets, any one enough, with exact thresholds), and the payoff/status of every idea.

Files (all in this directory):
- `counts_lab.py`: count laboratory. It is an exact breakpoint version of model_sa's R_batch with optional new pieces
  (GM / Huxley / HM large values with copies, Moebius 2k-th moments, beyond-edge plain moments) and a brute-force
  `joint_R` with the detector's witness freedom. Installing it into model_sa reproduces 0.8749568251160069 bit-for-bit.
- `count_curves.py` -> `count_curves.out`, `.json`: R(delta, x) for each piece at x = 0 and x = 1/2.
- `joint_check.py`, `joint_x05.py`, `plaincube_check.py`, `ext_check.py`: brute-force checks of the closed forms.
- `diag_cf.py` -> `diag_cf.out`: DIAG closed form for an arbitrary x = 0 count. It matches the LP to 1e-7 for M_kappa.
- `batch.py`, `runs/c*.jsonl`: full LP runs. c1 GM runs were killed (OOM, before chunking); c3 supersedes them.

Notation: U = row range; witnesses (r, m) are from Prop. detector-witness. alpha = 5/6 is the copies slope, and
R0(delta; a) = 1 - 2 delta a/(3a - delta) is the long-inverse closed form.

## 0. Bottom line

1. **Large-values methods.** None helps where it matters (sec. 2).
   - HM, Huxley, Heath-Brown, Jutila and Bourgain help only for delta >= 1/2.
   - Guth-Maynard (GM) with copies helps only for delta >= 0.41 (x = 0) or delta >= 0.42 (x = 1/2).
   - GM beats the MVT on the r < 1 (balanced) witness only for delta > 2/5 (sigma > 0.7, as for zeta).
   - On the long witness, GM + copies beats MVT + copies only for delta > 17/30 (sigma > 0.78). Copies help MVT
     more than they help GM, because GM's T-term grows with the family size H.
   - LP payoffs (sec. 6): 0 to 1e-9, both with the paper lemmas and with DIAG.
2. **More or better copies.** Copies are rigid (sec. 3).
   - Kummer theory gives exactly u(F^x)^6, plus O(1) from conjugation, units and ray classes.
   - Subfamilies, cosets, sectors and base change leave the copies/family ratio invariant.
   - Lemma B does not shorten H.
   - Higher Moebius moments are strictly worse: count 1/6 + k r(5/6 - delta). Model hk=2 gives 0.8749568, payoff 0.
3. **Theta large values.** No-go (sec. 4).
   - Every duality/HM inequality for a sparse row set sees mu only through |mu|^2 (positivity).
   - The theta identity needs a smooth lattice sum over rows.
   - Even GRH-size off-diagonals need delta > 1/2.
4. **Detector.** The count for the pure long witness *is* Ingham-with-copies (sec. 5).
   - Turan loses U^{O(T0)}.
   - Two-length/joint arguments need joint moments of length > 1, which are worse.
   - Structural point: the detector only gives r + m >= t, not = t. This "witness freedom" is why every
     within-the-edge joint moment is useless.
5. **Precise obstruction ("Ingham-with-copies world", PROVED at exponent level, sec. 1.3).**
   - Take U^{R0(delta)} rows with:
     - M_r saturated at r = 2t*/3 and r = t*;
     - S_m saturated for m <= t*/3;
     - everything else unsaturated.
   - This is consistent with every proved input: the P2 mean square for all H >= D^{1+c}, Lemma plain,
     Lemma marked at x = 0, the generic sextic sieve, beta* <= 7/8, and the detector at every t in [1, 3/2].
   - So any improvement must be a statement at total length U^{t*} ~ U^{1.12-1.17}, beyond the edge, or a
     non-moment argument that sees a single row's joint saturation pattern.
6. **What would suffice with DIAG (closed form, any one):**
   - M_kappa with kappa = 0.9637;
   - a diagonal Moebius mean square sum_{u<=U}|M_u(U^r)|^2 << U^{1+eps} up to r = 1.043 (copies beyond);
   - a joint diagonal moment sum_u |M_r S_m|^2 << U up to r + m = 1.066;
   - a joint copies-type slope 0.616;
   - the plain sixth moment diagonal up to 3m = 1.236 (CF; the LP gives 1.22).

   All five are beyond the edge, and all are false in a "saturated-DH world". The truly needed statement is the
   count (large-values) version, which DH satisfies.

## 1. Structural facts

### 1.1 Witness freedom (PROVED; read off Prop. detector-witness, P78 l.4510-4715)
- The detector at t gives a dyadic pair with D <= D* = U^t, DN >= D* and DN <= U^{21}.
- Pointwise bounds then force m <= 1/2, hence r in [t - 1/2, t] and **s = r + m in [t, t + 1/2]**, not s = t.
- The long plain tail is unavoidable.
  - Truncating the k-sum at DN <= D* U^eta would require the plain tails sum_{k > K} psi(k) k^{-rho} to be
    negligible. Since L(rho) = 0, these equal minus the smoothed partial sums of L at rho, and those can be
    saturated for every K <= U^{1/2}.
  - The Gamma shift also needs Y* >> D*^{(1-a+c)/c} U^{O(1)}.
  - So pieces with s up to t + 1/2 always occur, and only the functional equation caps m, at 1/2.
- Consequence: any bound valid only for products of total length <= 1 never applies to the adversary's
  balanced witness (2s/3, s/3), which can be taken with s slightly > 1. This is why DIAG, the within-edge plain
  sixth moment (k3) and the within-edge joint moment all have zero effect on R(delta, 0).

### 1.2 Three-slope formula at x = 0 (PROVED as exponent algebra; checked by brute force)
Suppose the beyond-edge excess is linear with slope alpha (long Moebius M_r, r > 1), lamJ (joint M_r S_m,
s = r + m > 1) or lam6 (plain cubes S_m^3, 3m > 1). Then

  R(delta, 0) = 1 - 2 delta lam/(3 lam - delta),  lam = min(alpha, lamJ, lam6).

- If lam = alpha, the detector sits at t = alpha/(alpha - delta/3).
- If lam = lamJ or lam6, it sits at **t = 1**, where the pure long witness is diagonal (1 - delta).
- So the long inverse witness is *not* essential: a better joint or plain-cube count works with the detector at t = 1.
- Checks: `joint_check.py`, `plaincube_check.py`, error <= 5e-4 from the grid.
- Copies (P2-type mean squares) would give slope 5/6 for all three, hence R0 again. This is the long-inverse
  robustness statement, extended to plain cubes.
- At x = 1/2 a slot-free joint slope is much weaker than kA. For lamJ = 0.75 at delta = 0.39: 0.6623, against
  0.6381 for kA = 0.9 and 0.6654 for the paper (`joint_x05.py`). Without DIAG only the long-Moebius lever is efficient.

### 1.3 Ingham-with-copies world (consistency, exponent level)
- N = U^{R0} rows, each with:
  - M saturated at r1 = 2t*/3 and r2 = t*;
  - S saturated on m in [0, t*/3];
  - M unsaturated on (r1, t*).
- At each t the rows have a witness: (2t*/3, t*/3) for t <= t*, using the witness freedom, and (t*, t - t*) beyond.
- An even more economical adversary uses M saturated at {2t*/3} and on [t*, 3/2], and S saturated only at m = t*/3.
  Then no two plain scales are available to multiply, so short products of the form S_{m1} S_{m2}^2 are useless.
- Every proved bound holds with equality or slack:
  - 1 - delta r1 = R0;
  - e(t*) - delta t* = R0;
  - 1 - 2 delta m >= R0;
  - the generic sieve gives max(1, s, 2(1+s)/3) - delta s >> R0;
  - Lemma marked and Lemma plain give nothing extra at q = 0;
  - beta* <= 7/8 holds since delta <= 0.75.
- Whether such rows exist as L-functions is not decided. The left-of-rho part of 1/L, which is uncontrolled, must
  cancel the residue 1/L'(rho) on the range (r1, t*). That is the long-inverse "length flexibility" issue.
- Note: long-inverse's model G(D) (unsaturated for all D <= U) is *inconsistent* with the t = 1 detector, which
  forces some M_r with r <= 1 to be saturated. The corrected adversary is the configuration above. This does not
  change any count.

## 2. Item 1: large-values methods, with copies

- A coefficient-blind large-values (LV) bound F(D, V, H) on the family {w <= H} gives, through the copies
  u a^6, a primitive count of (U/H)^{1/6} F, minimized over H >= U.
- Closed forms (U-exponents, D = U^s, V^2 = D^{1+delta}):
  - MVT: (U/H)^{1/6}(H+D)D^{-delta}, giving max(1, (1+5s)/6) - delta s. This is the paper.
  - GM (|W| <= N^2V^-2 + N^{18/5}V^-4 + T N^{12/5}V^-4):
    max(1 + s(2/5 - 2delta), 1/6 + s(9/10 - 7delta/6), 1/6 + s(7/5 - 2delta)).
  - Huxley (|W| <= N^2V^-2 + T N^4 V^-6): max(1 + s(1 - 3delta), 1/6 + s(1 - 4delta/3)).
  - HM: copies *cancel*. The copy blocks of the Gram matrix are rank-one with norm (H/U)^{1/6} D, so HM over the
    copy set equals HM over primitive rows.
    - With off-diagonal g, HM gives N <= D^{1-delta} only if g <= D^delta.
    - GLH (g = D^{1/2}) needs delta > 1/2.
    - The a priori beta* (g = D^{7/8}) needs delta > 7/8 > d_max.
  - Plain witnesses: each estimate is applied to S_m^k, k <= 4 (coefficient-blind, copies allowed, optimistic).

| method | first delta with a gain, x = 0 / x = 1/2 | status |
|---|---|---|
| GM + copies | 0.41 / 0.42 | HEURISTIC: GM is not proved for this thin family. Gains appear only for sigma > 0.7 on r < 1 witnesses; on the long witness only for delta > 17/30 |
| Huxley + copies | 0.51 / 0.51 | HEURISTIC |
| HM (GLH off-diagonal) | 0.50 / 0.50 | CONDITIONAL on GLH; needs delta > 1/2 regardless |
| HM with beta* off-diagonal | none (needs delta > 7/8) | PROVED, useless |
| Heath-Brown / Jutila Q^2T-aspect, Bourgain | sigma >= 3/4 or 11/14 | not available for thin families; would be useless anyway |
| Turan power sums | sigma -> 1 | useless (sec. 5) |

At delta = 0.39 (x = 1/2), the long witness r = 1.12 gets:

| bound | value |
|---|---|
| paper (copies) | 0.663 |
| GM + copies | 0.861 (H = D^{6/5}) |
| Huxley + copies | 0.81 |
| MVT without copies | 0.683 |

LP payoffs (coarse grid 81x6; baselines on the same grid are paper 0.8749563 and DIAG 0.8725487): GM gives 0.8749563;
GM + Huxley + HM gives 0.8749563 (paper) and 0.8725487 (DIAG). **Payoff 0 to 1e-9.**

## 3. Item 2: more or better copies

(a) **Kummer rigidity (PROVED).**
- (u/n)_6 = (w/n)_6 for almost all n iff u w^5 is in F^{x6}. For sixth-power-free u, that means w = u a^6 (unit factors included).
- Other O(1) coincidences:
  - w = u-bar when nu∘sigma = nu-bar. Then |A_{u-bar}| = |A_u|, by the substitution n -> n-bar.
  - Unit multiples eps*u: the L-function of row (nu, u) reappears in row eps*u twisted by nu (eps/.)_6^{-1}, a fixed
    ray class, so it repeats only if the family also contains that twist. The family then grows by the same factor.
- So the copy multiplicity is (H/q_u)^{1/6} times at most 12. Exponent payoff 0.

(b) **Subfamilies (PROVED, exponent level).**
- Rows w ≡ w0 (mod c), sectors, annuli, or v a^6-structured families: Poisson's off-diagonal vanishes iff the
  family size H/N(c) >= D^{1+c'}.
- The copy proportion changes by at most 6^{omega(c)} N(c)/phi(c) = N(c)^{o(1)}.
- Per-target cosets (w0/u a sixth-power residue) partition the targets into 6^{omega(c)} classes.
- Payoff 0.

(c) **Base change K/F (PROVED, exponent level).** In the sextic family over K, the copies of the base-changed row are
u alpha^6 with alpha in O_K. Norms, lengths and conductors all scale by [K:F], so the count exponent is invariant.

(d) **Denser copies via quadratic twists (circular).**
- {nu chi_u lambda_d : d <= H} contains L(s, nu chi_u) at d = a^2, which is density 1/2.
- The Poisson dual (quadratic Gauss sums, no mu-absorption) has its own square class k = j^2. Its principal dual row
  is again a Moebius sum of a twist of the same row, contributing H^{1/2}|A_{u'}|^2, which exactly matches the copies.
- Without that circularity this would give a zero-free region 3/4.
- Sanity check of the framework: only mu-absorption (6 -> 3) breaks circularity, and it comes with density 1/6.

(e) **Bootstrap / Lemma B (PROVED, R2 sec. 2).**
- Lemma B bounds the copy defect |A_u - A_{u p^6}| << (D/Np)^{beta*+v}. That justifies using all a, including small
  ones, but it does not shorten H.
- M_kappa (kappa < 1) is false in a DH-world (R2 sec. 6), so no bootstrap from H(beta) can give it.
- The Ingham-with-copies world (sec. 1.3) satisfies H(7/8).

(f) **Higher moments (PROVED as algebra; model).**
- A diagonal 2k-th moment over the full family needs H >= D^k (the length of A^k).
- With copies the count is 1/6 + k r(5/6 - delta). That exceeds the k = 1 count by (k-1) r(5/6 - delta) > 0.
- The generic sextic sieve's (H D^k)^{2/3} term makes it worse still.
- Model `hk=2`: 0.8749568 (payoff 0); count curves show no gain at any delta.
- Only a *Lindelof-on-average* 2k-th moment (valid at H ~ D) would help, and that is GRH-strength.

(g) **Copies for products** (a P2-type mean square for mu*1 or plain-cube coefficients): slope 5/6, hence R0 (sec. 1.2).

## 4. Item 3: theta-based large values

(a) **Positivity.**
- Large-sieve duality and Halasz-Montgomery bound sum_{u in R}|<a, phi_u>|^2 by Cauchy-Schwarz with a positive
  weight omega(n).
- The Gram entries are therefore sum_n omega(n) chi_n(u u'^5) W, and depend on mu only through |mu|.
- The "mu-weighted Gram" sum_n mu(n) chi_n(u u'^5) W = A_{u u'^5}(D) (another row, conductor <= U^2) cannot arise
  in such an inequality.
- If it did: with GRH-size entries D^{1/2} HM still needs delta > 1/2, and with beta*-size entries D^{7/8} it needs
  delta > 7/8.

(b) **Theta needs smooth row sums.**
- P2's mechanism is Poisson in u over a lattice box, then mu gamma_{-1} ∝ gamma_2, then theta reflection, then the
  quadratic sieve.
- For a sparse R, Poisson gives a kernel R-hat(k/n) that does not separate. After mu-absorption the form
  sum_{n,k} gamma_2(n) chi_n(k) R-hat(k/n) needs >> |R| Mellin modes to separate. That is a hybrid full-family mean
  square: no gain over P2 + copies (HEURISTIC).

(c) **Exact Poisson identity.**
- S(H, D) = DH + (H/D) E_dual(D^2/H, D). The zero rows enter only through the total.
- At the optimum H = D every term (diagonal DH, dual column term (H/D) D^2, dual diagonal D^2) has size D^2.
- There is no subtraction available without a lower bound for the non-exceptional rows, which would need the fourth
  moment, beyond the edge.

(d) **Row-dependent amplifiers** (Sum_l chi_l(u)-bar chi_l(w)).
- This is Lemma marked with marks of length L. It needs H >= D L^{c_m}, with c_m = 2 in the paper and 1 under DIAG.
- It gives *individual* bounds |A_u|^2 <= D^{11/6} U^{1/6} L^{5c_m/6 - 1}, never counts.
- For the count it loses (5/6 - 2q)z (long-inverse), and it gives nothing at x = 0.

(e) **Target-averaged probe** (sum over targets u0 of the 7/8-probe twisted by chi_{u0}).
- The low side gains |R|^{-1/2} by a large sieve over targets, so |R| <= (low/signal)^2.
- This only counts zeros within o(1) of the individual zero-free line, and only for small conductors (copies). It is
  irrelevant at sigma ~ 0.67.

## 5. Item 4: the detector

(a) **The count for the pure long witness is exactly MVT + copies.**
- At D = U^{t*}, rows u <= U, the coefficient-blind bound is >= D^{1-delta} (random combinations of D^{1-delta}
  characters).
- P2 + copies gives U^{1/6} D^{5/6 - delta}. That already beats the coefficient-blind limit because it is
  mu-specific (P2).
- Improving it is the r1 / M_kappa target (sec. 6).

(b) **The pure long witness is harmless at t = 1** (count 1 - delta). The detector is pushed to t* > 1 only by the
balanced witness (sec. 1.2). The two binding objects can be traded against each other, but each needs a beyond-edge
input.

(c) **Two lengths / joint constraints.**
- Every zero row carries a witness at every t, but min_t max_w uses each t separately.
- Joint information about one row (M_{r1} and M_{r2}, or M and S at different t) needs joint moments of total length
  r1 + r2 > 1. With copies those give (1 + 5s)/6 - delta s at the larger s, which is worse.
- Hoelder between separate moments gives at best the geometric mean of the counts, which is >= the min.

(d) **Turan / power sums.**
- At s = rho the normalized Moebius block is Res_rho + sum over the other zeros rho' of c_{rho'} D^{rho' - rho}
  + (an integral left of rho).
- The residues 1/L'(rho') are unbounded, and >= (T0/pi) log U zeros compete.
- Turan's lemma (sup over I of |f| <= (C|I|/|J|)^{K-1} sup over J of |f|) loses U^{O(T0)}.
- Interpolating saturation across lengths therefore needs control of 1/L left of rho. The a priori beta* gives none:
  1/L is controlled only for Re s >= 7/8.

(e) **Residue structure.**
- Every dyadic block at s = rho contains Res = W-hat(0)/L'(rho), of size >= U^{-4e}. Saturation fails only where the
  left part cancels it.
- The adversary of sec. 1.3 needs that cancellation on the range (r1, t*), a "flat then drop" left part. This is
  consistent with an uncontrolled 1/L.
- No per-row implication was found (agrees with long-inverse).

## 6. Payoffs (theta; full LP unless marked CF = DIAG closed form, diag_cf.py)

| input | status | paper lemmas | with DIAG |
|---|---|---|---|
| none | - | 0.8749568 | 0.8725498 |
| GM + copies (coefficient-blind transplant) | HEURISTIC | 0.8749563 (= ref 0.8749563 on grid 81x6) | - |
| GM + Huxley + HM(GLH) | HEUR./COND. | 0.8749563 (payoff 0) | 0.8725487 (= ref, payoff 0) |
| Moebius 4th moment with copies (H >= D^2) | hypothetical | 0.8749568 | - |
| plain 6th moment, within edge (3m <= 1) | open (no gain) | 0.8749568 (slack-audit) | 0.8725498 |
| plain 6th moment, all m <= 1/2 (beyond edge) | open (zeta-6th-moment type) | 0.8727074 | 0.8694441 (CF 0.869351) |
| plain 8th moment, all m <= 1/2 | open | 0.8713661 | - |
| DH counts | open | - | 0.866827 (CF 0.866647) |

Thresholds for 27/31 *with DIAG* (CF; checked against brute force and the LP):

| beyond-edge input (any one) | needed |
|---|---|
| M_kappa: sum_{q_u <= H}|A_u(D)|^2 << H^{1+eps} D for H >= D^kappa | kappa = 0.9637 (LP: 0.9627 gives 0.870928) |
| Moebius diagonal over u <= U up to length U^{r1}, copies beyond | r1 = 1.043 |
| joint sum_u |M_r S_m|^2 << U for r + m <= 1 + w (no copies) | w = 0.066 |
| joint copies-type slope lamJ beyond s = 1 | lamJ = 0.616 |
| plain sixth moment sum_u |S_m|^6 << U for 3m <= 1 + w | w = 0.236 (CF); LP: m6 = 0.39/0.40/0.4118 gives 0.87251/0.87203/0.87051, so 3m6 = 1.22 |

- Note on the "slight improvement": M_kappa at kappa = 0.964 shifts the *edge* (it is diagonal up to r = 1/kappa)
  and also lowers the slope. In pure slope form the long count would need alpha = 0.616, not 0.80.
  The diagonal extension needed is small (4.3% past the edge) because the copies do the rest.
- Every target above is a moment statement that fails in a DH-world whose zero rows are saturated at all long
  lengths. Only the count version is consistent with DH.
- No method is known that proves a count (large-values) estimate beyond MVT + copies for sigma < 0.7. This is the
  classical Ingham barrier, transplanted.
