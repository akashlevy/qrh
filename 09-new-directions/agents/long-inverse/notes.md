# long-inverse: counting rows with long saturated inverse witnesses (r in [1, 1.35])

Agent `long-inverse`, 10 Oct 2026. Outcome: **no new mechanism**. The obstruction is precise, and it
is the family analogue of the classical sigma < 7/10 (Ingham) barrier. Files:
`model_li.py` (modified count, knobs `mix`, `mixfull`, `lam`, `r1`), `robust_check.py`, `run_*.json`.

## 0. Exact anatomy of the binding count (PROVED, elementary algebra; checked numerically)
Suppose a long object of length s > 1 is counted by e(s) - delta*s, with e(s) = 1 + alpha(s - 1).
Then at x = 0 the P78 count min_t max_{(r,m)} is exactly
  R(delta,0) = 1 - 2*delta*alpha/(3*alpha - delta)      (crossing t* = alpha/(alpha - delta/3)).
- alpha = 1 (the plain mean value theorem, U + N) gives 3(1-delta)/(3-delta). This is
  **Ingham's 3(1-sigma)/(2-sigma)**; the model reproduces it to 1e-5.
- alpha = 5/6 (sixth-power copies) is P78: 0.69194 at delta = 0.39.
- DH needs alpha <= delta.
General x: excess over DH = (alpha - delta) k / (2(alpha - delta + k)), with k = delta P_x/D_x.
For zeta, Ingham's exponent for sigma < 7/10 has not been improved; Guth-Maynard helps only for
sigma > 7/10 (delta > 0.4). The binding delta = 0.39 sits just below that edge.

## 1. What is needed (model, model_li.py)
| long count used for r >= 1 | theta |
|---|---|
| paper e = (1+5r)/6 | 0.874957 |
| slope only, e = 1 + lam(r-1): lam = 0.75 / 0.7 / 0.65 / 0.6 / 0.5 | 0.87451 / 0.87417 / 0.87377 / 0.87328 / 0.87223 |
| diagonal (e = 1) up to r1, then slope 5/6: r1 = 1.05 / 1.10 / 1.15 / 1.20 | 0.87342 / 0.87195 / **0.87058** / 0.86979 |
| diagonal for all r (lam = 0) | 0.86979 (= the perfect-count value) |

Slope reductions are useless. 27/31 needs a **diagonal Moebius mean square
sum_{u~U} |M_u(U^r)|^2 << U^{1+eps} up to r ~ 1.14**. That is M_kappa with kappa ~ 0.88: exactly
the blocked input (OBSTRUCTIONS section 3).

## 2. Ideas tried and the precise reason each fails
1. **Detector redesign forcing r <= 1.**
   - *Conservation law (PROVED).* Let a detector have coefficients c_n = sum_{l|n} mu(l) w(l)
     vanishing on 1 < n <= X. This forces w = 1 on [1, X]. For X < n <= 2X the only divisor above X
     is n itself, so c_n = -mu(n)(1 - w(n)) exactly: a pure Moebius block of length X with m = 0.
     It satisfies every pointwise bound, so it is an admissible saturated witness. Every
     mollifier-type detector at t > 1 therefore admits the case-(a) witness (t, 0).
   - Power damping w(l) = (X_1/l)^kappa kills blocks with r > t_1 + eta, but then the detector
     length becomes t_1. Smooth or Selberg weights (P(log(X/l)/log X)) change only constants. With a
     slow transition, the tail beyond X_2 again has Moebius blocks up to t_2, and the Lambda pieces
     give prime witnesses. Lemma marked counts prime witnesses only for z < 3/8 (2r + 8z <= 3), so
     their count is >= 1 - 3 delta/8.
   - Heath-Brown identity, products of detectors, L^k mollifiers: these split the length among more
     factors. The available moments see one inverse factor or two plain factors, so the count gets
     worse. A k-fold HB detector gives 1 - 2 delta t/(3k). A product (1 - LM_U)(1 - LM_{U^{t-1}})
     gives 1 - 2 delta/3, because the adversary puts the length t - 1 in a short pure inverse factor.
2. **Length flexibility (fails as a consequence of the bin data).** All a proof can use is: 1/L is
   analytic and U^{o(1)}-bounded on Re s >= a + 2e (heights <= T_1); there is a pole at rho; and
   |L'(rho)| <= U^{4e} (Phragmen-Lindelof, then Cauchy). The configuration
   G(D) = c(1 - Phi(log D/log U)) has the right properties: Phi = 1 on [0,1] and 0 beyond t, and
   its Mellin transform is entire and O(1/e) on the line. It is unsaturated at every D <= U and
   saturated beyond U^t. So no per-row implication from long to short saturation exists.
   - Laplace/Tauberian versions run into two problems. At shift kappa >> e the n = 1 term
     dominates. At kappa ~ e the residue (size U^{-O(e)}) and the losses U^{O(e)} are of the same
     order, while the saturation loss eps satisfies eps >> e. Any argument needs a power gain
     independent of e, which requires control of 1/L left of rho.
   - Doing this on average over rows is a Moebius mean square below the edge, i.e. M_kappa.
3. **Prime slots inside the amplification.** You need H >= D Z^{2z}, so the exponent changes by
   (5/3 - 2q) z > 0 (q <= delta/2). Even without Lemma marked's doubling it is (5/6 - 2q) z > 0.
   This always loses.
4. **A priori beta\*.** It gives |M_v(D)|^2 <= D^{2beta\*-1} for every row, hence a mean square
   U^{1 + 0.75 r}. That is worse than (1+5r)/6 for every r. For a < beta\* there is no gain in the
   pointwise bound. beta\* is already used in Lemma plain through kappa.
5. **Turan / Halasz-Montgomery (with near-orthogonality as in the cubic Gauss-sum paper, dual.tex).**
   - The BHM diagonal alone gives U^{r(1-delta)} = U^{0.70} (r = 1.15). That is worse than 0.6765.
   - The off-diagonal needs delta > 1/2.
   - Turan's method loses e^{O(K)} and is useful only as sigma -> 1.
   - Hybrid (t-aspect) families give slope 1 (Ingham). Combined with the copies, the exponent is
     (1 + 5r + lambda)/6 for T = U^lambda, which is worse.
   - **Only the sixth-power copies beat slope 1.**
6. **Robustness.** Suppose a full-family diagonal mixed moment existed. The copy identity also holds
   for plain polynomials: S_u(N) = sum_{j | a^inf} psi_u(j) q_j^{-1/2} S_{ua^6}(N/q_j). So products
   M_r S_m could be amplified, with count (1+5s)/6 - delta s at s = r + m. R(delta,0) is unchanged:
   the adversary moves to s = t*. `robust_check.py` confirms this to grid accuracy.

## 3. Coordinator's question: joint counting (M_r, S_m) with S_m as a slot
- (a) Lemma marked needs product-form prime marks. The marks are 1_{p|A} on the theta index, a
  bounded number of disjoint prime lists, and Cor. marked-reflection treats each marked prime as an
  active or inactive local branch. The paper says product form "not merely the pointwise divisor
  bound" is used (l. ~7320), and "arbitrary bounded weights are permitted only as the stated product
  of individual prime-slot coefficients" (l. ~9300).
  - M_r S_m has mark sum_{k|c, k~U^m} mu(k): all squarefree k, unbounded omega(k), tiny primes.
  - Splitting L = L_{<y} L_{>=y} makes the rough part slot-compatible: it is a bounded sum of
    products of at most m/eta prime slots.
  - But on the bin line L_{>=y} = 1 + O(U^{-c e eta}) inside the eps-slack. Saturation can
    therefore sit entirely in the y-smooth part, which has neither slot form nor a functional
    equation.
- (b) Model results:
  - **Full slot** (r + 2m <= 1, 2r + 8m <= 3; `mixfull`): never applicable, because r + m >= t >= 1.
    R and theta are unchanged at 0.874957.
  - **Fractional slot** (`mix`; it presupposes a saturated sub-factor of S_m, which a single plain
    polynomial does not provide): R(0.39; x = 0, 1/4, 1/2) goes from 0.6919, 0.6803, 0.6654 to
    0.6677, 0.6666, 0.6654. **theta = 0.874957, unchanged.** The binding x = 1/2 point is the pure
    long inverse witness (r = t = 1.125, m = 0), where no plain factor exists.
  - With kA = 0.9 or w1 = 0.3 the result equals that knob alone (0.87116 / 0.86979).
- (c) Powers S_m^k lengthen the slot to km, against the cap (1 - r)/2. Useless.

## 4. Status
- PROVED: the closed form in section 0, the conservation law, the robustness statement, and the
  negative model payoffs.
- HEURISTIC/MODEL: the length-flexibility counter-configuration. It is a Mellin-level model; no
  actual L-function is constructed.
- Conclusion: the long-inverse lever is equivalent in strength to "Ingham below 7/10 for this
  family". It needs a diagonal Moebius mean square to r ~ 1.14 (M_{0.88}); nothing weaker in the
  witness calculus suffices.
