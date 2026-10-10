# Part II high side at the slot-free geometry (ℓ = 0)

*Check of `paper78.tex` Part II (line 6807 onwards) and of the shared Part I analytic lemmas it invokes, at the geometry used by the slot-free certificate in `REPORT.md`. Written 9 Oct 2026.*

Labels:
- **[R]** rigorous, given the paper's stated lemmas;
- **[R-]** routine, but the written text has to be restated;
- **[H]** heuristic;
- **[G]** gap.

---

## 0. Verdict

**The high side is sound at the slot-free geometry, with routine restatements only. Nothing the model omits is binding.**

- **(a) [R]** The high exponent E(d) generalises *exactly* as the R4 model and `exact_cert_theta.py` encode it. It is re-derived in §2 directly from the proof of Lemma `bin-contour-accounting`, which is already written for general (l_x, l_y, ℓ).
- **(b) [R]** The zero-slot counts R_0(t) = max(L(t), 1 − 2δt/3), for every t ∈ [1, 3/2], follow from three stated results of the paper:
  - Prop `detector-witness`, which is stated for every t;
  - Eq `no-slot-inverse-count`, which the Remark after Prop `detector-counts` asserts for each fixed t;
  - Lemma `plain`, case z = 0.

  Minimising over t gives the closed form R_0(δ) = (15 − 16δ)/(15 − 6δ). It matches the checker's exact `R_t` at t* = 5/(5 − 2δ) with zero discrepancy.
- **(c) [R]** Slots supply nothing else. At K = 0:
  - the compensated probe *is* the base probe I_η(Z^{l_x}, Z^{l_y}, Z);
  - the full correction (eq `holomorphic-selected-tuple`) is the base correction 𝓗_{η,u};
  - the normaliser is A_T = 1, with residue error 𝓡 ≡ 0.

  Definition `stage-high-data` explicitly admits this case (line 5652): "For the base probe … ℓ = 0, 𝔥 = 𝓗_{η,u}". Compensation exists only to neutralise the high-side multipliers created by *marked* slot primes. With no marks there is nothing to cancel, and Lemma `local-euler` already gives the scalar Euler structure.
- **Hard-wired constant that must change [R-].** The paper hard-wires σ0 = 7/8 in four shared lemmas (hypothesis σ0 ∈ [7/8, 1)) and in the second Euler region (x_r ≥ 7/8). Both must be lowered to σ0 = 1733/2000. This is harmless: every local exponent still has a large margin (§5). It is a real textual gap in the sense that the paper as written does not cover σ0 < 7/8. It is not a mathematical obstruction.
- **No missed constraint.**
  - Exact certificate rerun with `theta=1 dropq=1 wE=0` (so no energy-wall relaxation): **ALL CHECKS PASS** (3055 boxes).
  - The LP optimum at S = 0 is exactly **123/142 ≈ 0.866197**, at l_x = 14/71, l_y = 134/213.
  - Three constraints bind: the low side, the detector floor bin, and the (V) cubic-large-sieve branch 7M/3 − 2.
  - The row counts are **not** binding: replacing R_0 by the ideal R = 1 − δ gives the same 123/142.

---

## 1. Parameters (exact)

| quantity | value |
|---|---|
| l_x | 989/5000 = 0.1978 |
| l_y | 3141/5000 = 0.6282 |
| ℓ (slots) | 0 (K = 0 slots) |
| h = 1 − l_x + ℓ | 4011/5000 = 0.8022 |
| M = l_x + l_y | 0.826 (so M + ℓ = 0.826, not 1) |
| b = l_y − l_x | 0.4304 |
| σ0 | 1733/2000 = 0.8665 |
| detector floor | a = 51/100 (δ0 = 1/50) |
| contour | z0 = 17/50 |
| row ranges | d_min = 1/100, ζ = 1/10000 |
| C(s) = s − 5/6 + l_x/3 + ℓ/6 | s − 11511/15000 |
| C(σ0) | 2973/30000 = 0.0991 |

The formula for C(s) is Def `stage-high-data`, eq `stage-mellin-exponent` **[R]**. Paper check: at (17/48, 23/48, 1/6) it gives s − 11/16.

**Low-side compatibility.** The split low bound |I_η| ≪ Z^{l_x/2+ε} = Z^{2967/30000+ε} sits below Z^{C(σ0)} by Z^{−1/5000}. So eq `low-probe-contract` holds with any ω > 0. This is exactly the model's low constraint σ0 ≥ 5/6 + l_x/6 at g = e = 0.

---

## 2. (a) The high exponent for general (l_x, l_y, ℓ)

Source: Lemma `bin-contour-accounting` (line 6024ff). Its proof never uses specific geometry; only its σ0-range hypothesis does (§5).

On the central contours Re s = a + 16e, Re w = 1 − a − 6e, Re z = z0, the outside factor 𝒲 · q_u^{−z0} contributes

  l_x(1/2 − z0) + a + z0 − 1 − a l_y − d z0 + (16 − 6 l_y)e.

The central row bound, eq `stage-central-bound`, is U^{R+δ/2} Z^{ℓ(z0−1/2)+g}. Adding it and subtracting C(σ0) = σ0 + l_x/2 − 1 + h/6 gives **[R]**

  E(d) = a − σ0 + h(z0 − 1/6) − a·l_y − ℓ/2 + g + d(R + δ/2 − z0),  with g = qℓ = xδℓ.

This is eq `stage-high-exponent` verbatim. At d = h with a = (1+δ)/2,

  E(h) = (1+δ)(1−l_y)/2 − σ0 − ℓ/2 + xδℓ + h(R + δ/2 − 1/6).

That is exactly the R4 / checker high constraint. The contour constant z0 cancels at d = h.

**Sanity check.** At the paper's geometry, the constant term is 25/96 − 7/8 − 1/12 − 13/96 + 13/16 = −1/48 = C_0, and the δ-coefficient is 25/96 + 13/32 = 2/3. Both match eq `common-high-exponent`.

**Where d ≠ h.** E is affine in d. The checker evaluates both endpoints d_min and h + ζ, so no slope-sign assumption is needed (unlike the paper, which relies on positive slopes).

**At ℓ = 0.**
- g = 0 and x is inert.
- The amplitude bins of §`amplitudes` are a single trivial bin (q = 0).
- Lemma `bin-prime-bound` is not invoked.
- The central bound is the Part I one: L^S(w, χ_•(u)) ≪ U^{δ/2+12e} (Lemma `buffered-bins`), and 𝓗_{η,u} ≪ U^ε on region one (Lemma `local-euler`).

**Floor bin [R].** R = 1 and q = 0 give E_floor(d) = 0.51(1 − l_y) − σ0 + h(z0 − 1/6) + d(1 + 1/100 − z0). The slope is positive, so the worst case is d = h + ζ:

| floor computation | value |
|---|---|
| exact | −2.93·10⁻⁴ |
| checker (uses the conservative slope 1 + δ0/2 − 1/6 on the ζ-extension) | −2.76·10⁻⁴ |

**Small rows [R].** Source: Lemma `outer-row-tails`, eq `small-row-sum`. The exponent is h(z0 − 1/6) − l_y/2 + (63/50)d_min:

| version | value |
|---|---|
| with 63/50·d_min | −0.1625 |
| with the paper's cruder 2·d_min | −0.1551 |

**Large rows [R].** The bound B_0 + (h+ζ)(1+ε) − ζ v, with B_0 = l_x/2 + 1 + l_y, tends to −∞ as v grows, for any geometry.

**Principal margins [R].**

| margin | value |
|---|---|
| m_w = l_y/20 | 0.0314 |
| m_z = h/600 | 1.34·10⁻³ |
| μ (residue error) | not needed, since 𝓡 ≡ 0 |

---

## 3. (b) Zero-slot row counts

**Inputs, all for a fixed bin with δ ∈ (1/50, 5/6] and any fixed t ∈ [1, 3/2]:**

1. **Prop `detector-witness` [R].** It is stated "for every t ∈ [1, 3/2]" and has no slot hypothesis. It gives dyadic r, m with
   - r ≤ t + o(1),
   - r + m ≥ t − o(1),
   - m ≤ 1/2 + O(ε),

   and rowwise lower bounds |M_r|² ≥ U^{δr−ε} and |S_m|² ≥ U^{δm−ε}.
2. **Eq `no-slot-inverse-count` [R].** It gives #𝓑 ≪ U^{e(r)−δr+ε}, with e(r) = max(1, (1+5r)/6). This comes from Lemma `inverse-amplification`, whose rows have valuations ≤ 5 and the common profile W_base. The Remark at line 15449 states it explicitly *for each fixed t ∈ [1, 3/2]*, with no supply hypothesis.
3. **Lemma `plain`, case z = 0 [R].** It has no length restriction and no κ / zero-free hypothesis. Applied to two copies of S_m on the rows 𝓡_0, with effective width 1 (ν fixed, no moving radical), it gives #𝓑 ≪ U^{1−2δm+ε}. The retained rows induce non-principal characters by the ramification paragraph (line 4259).

**Combination [R].** On each witness subdivision, #𝓑 ≤ U^{min(e(r)−δr, 1−2δm)+ε}, and m ≥ t − r − o(1).
- **Case r ≤ 1.** The minimum of 1 − δr and 1 − 2δ(t−r) is at most 1 − 2δt/3, with the crossing at r = 2t/3 ≤ 1.
- **Case 1 ≤ r ≤ t.** Then 1 − α + (α−δ)r ≤ L(t), using δ ≤ α = 5/6.

Hence R_0(t) = max(L(t), 1 − 2δt/3) + ε, exactly the paper's case analysis with selection gain x = 0. The choice of t is per bin and made before the target, which the high-bin lemma permits (R may depend on d, a, β*).

**Optimal t.** t* = 5/(5 − 2δ) ∈ [1, 3/2] for δ ≤ 5/6, giving

  **R_0(δ) = (15 − 16δ)/(15 − 6δ).**

Spot values: R_0(0.3) = 0.7727 and R_0(5/6) = 1/6 = 1 − δ.

| check | result |
|---|---|
| `scratchpad/p2/check_slotfree.py`: exact Fractions, checker's `R_t` at t* against the closed form, δ = 0.01 … 0.83 | max discrepancy 0 |
| brute-force minimisation over t at δ = 0.3 | agrees |
| R4 model at S = 0 (fI = 1 − δr, fP = fP0 = 1 − 2δm) | identical rule |

**Is the general-t count needed?**
- Yes, as an input: with t = 1 only (the paper's literal zero-slot clause, R = 1 − 2δ/3), the LP gives **0.877911**.
- But it is not binding: the non-floor high side has slack ≈ 0.009 at its worst point (δ ≈ 0.114, d = h + ζ). The ideal count R = 1 − δ gives the same optimum 123/142.

**Caveat [R-].** The paper's Prop `detector-counts` *statement* gives the zero-slot bound only at t = 1. The general-t version is a one-paragraph restatement of the proof (the "Cases requiring no selected primes" and "Small capacities" paragraphs, both slot-free). It should be written as a separate corollary.

---

## 4. (c) What the slots provide, and what disappears at K = 0

| slot ingredient | role in the paper | status at K = 0 |
|---|---|---|
| marks 1_{D \| cn³}, rescaled terms (eq `compensated-probe-definition`) | lengthen the completed row to Z^{1+ℓ'}; give M + ℓ = 1 on the low side | absent; the probe is I_η(X,Y,Z) **[R]** |
| compensation (marked minus rescaled) | turns each marked prime's high multiplier η̄(p)Q^{x+z−1}P_p*/P_p into −χ̄_p(u)-type sums Q_i plus errors D_i | nothing to compensate; the correction is ∏H_p = 𝓗_{η,u} (eq `holomorphic-selected-tuple` with an empty tuple) **[R]** |
| G_p, 𝓑_p, Prop `probe-errors`, Lemma `compensated-absolute-local-bounds` | local errors of the two-term operation | vacuous; small and large rows use Lemma `local-euler` directly, as in Part I's "first-stage outer ranges" **[R]** |
| normaliser A_T(Z) = (−1)^K Z^{−ℓ/6}∏S_i, κ_P, even K | principal residue with slots | A_η = 1, 𝓡 ≡ 0 ("The base correction" paragraph, line 6314) **[R]** |
| amplitudes q, Lemma `bin-prime-bound` | selection gain 2qz in counts; +qℓ in E | q = 0; the gain is irrelevant at ℓ = 0 |
| supply 8/39 > 7/37, mesh, b_round, ζ < 5ℓ − h | feasibility of selected capacities | not needed: the zero-slot counts have no supply hypothesis |
| Lemma `plain` with z > 0 (κ = 2β* − 1, Δ/4 capacity comparison, Θ restriction) | plain capacity | only z = 0 is used: no κ, no zero-free hypothesis |
| Lemma `marked` with slots, canonical recursion | inverse capacity | only the no-slot case → eq `raw-moment` → Lemma `inverse-amplification` is used; the statement covers the empty slot list |

The paper's remark that "the subtraction is designed to cancel the scalar prime contribution on the high side" (line 6913) refers only to the multipliers created by marked primes. Without marks, the high series of the base probe is the scalar Euler product of Lemma `local-euler` (eq `probe-high`). The paper uses it as such in Part I.

**The split is an analysis device on the low side only.** The physical expression bounded on the low side and expanded on the high side is the same I_η(Z^{l_x}, Z^{l_y}, Z), so Prop `continuation-criterion` applies to one function J = I_η/c_S **[R]**.

---

## 5. Every fixed-geometry use in Part II (and the shared lemmas it calls)

| place (line) | paper's fixed value | slot-free status |
|---|---|---|
| bootstrap eq `part-II-bootstrap` (6817) | Δ = β* − 7/8 ≤ 1/24, κ ≤ 5/6 | Δ0 = β* − σ0 ∈ (0, 11/12 − σ0]. κ = 2β* − 1 ≤ 5/6 still holds by Thm 11/12. Δ ≤ 1/24 itself was used only in the slotted capacity comparison and in m_hi, so it is not needed **[R]** |
| bin ceiling eq `part-II-bin-ceiling` (6825) | δ ≤ κ ≤ 5/6 | unchanged; the checker uses DMAX = 5/6 **[R]** |
| C_II(s) = s − 11/16 (6835, 15518, 8641) | 3/16 at 7/8 | C(s) = s − 5/6 + l_x/3; C(σ0) = 0.0991 **[R]** |
| geometry eq `old-eq:5.1` (6860) | (17/48, 23/48, 1/6), h = 13/16, M + ℓ = 1 | (989/5000, 3141/5000, 0), h = 4011/5000, M = 0.826 |
| slot definition, compensated probe (6866–6917) | K ≥ 1 | K = 0, empty **[R]** |
| Lemma `ray-prime-normalizer` (7007) | slot normaliser | unused |
| Lemma `probe-row-norm` (8112) | uses M' + ℓ' − 1 = −3d | General form, from the displayed T_d bound before substitution: E ≤ M' + max(0, (1 + 3ℓ' − 2M')/4), valid when M' + ℓ' ≤ 1. At (0.826, 0) this is ≤ M' **[R-]**. This is low side, recorded for completeness. |
| Prop `probe-low` conditions (8596) | l_x ≥ 3/16, l_y ≥ 5/16, l_y − 11b/6 ≥ 1/12 | The last condition fails: l_y − 11b/6 = −0.161. So P_a²/Y' ≫ P_a^{1/6}, and the paper's Gram bound would force g ≥ 2b − l_y = 0.2326 (σ0 ≳ 0.98). This is precisely why (LoA6) must remove the P²/Y term. With (D1) only (P²/Y kept), the slot-free LP gives 0.876322. Low side, not part of this check. |
| eq `stage-mellin-exponent`, Def `stage-high-data` (5590ff) | general (l_x, l_y, ℓ) | admits the base probe with ℓ = 0 explicitly **[R]** |
| **Lemma `local-euler` region two** (eq `local-region-two`) | x_r ≥ 7/8 | **must become x_r ≥ σ0 [R-]**, §5.1 |
| **Lemma `fixed-bin-contour`** (5811) | σ0 ∈ [7/8, 1) | not used in the proof (moves stay in D_1 because a ≤ β*); restate as σ0 ∈ (51/100, 1) **[R-]** |
| **Lemma `bin-contour-accounting`** (6024) | σ0 ∈ [7/8, 1) | only inherited; the algebra is general **[R-]** |
| **Lemma `principal-residue-interface`** (6178) | σ0 ∈ [7/8, 1); paths in D_2 | needs D_2(σ0) **[R-]**, §5.1 |
| **Lemma `outer-row-tails`** (6345) | σ0 ∈ [7/8, 1); small rows in D_1(3/8) | β* + e + 1/2 > 1.3665, so use D_1(1/3) **[R-]** |
| eq `shared-principal-nonvanishing` (H_η contraction) | sup over Re s > 7/8 | needed on Re s > σ0; P_0 exists, §5.1 **[R-]** |
| contour Re z = 17/50; regions D_1 | tied to the 51/100 floor, not to geometry | unchanged **[R]** |
| Lemma `high-bin` (15700), eq `common-high-exponent` | C_0 = −1/48 etc. | general E(d) of §2 **[R]** |
| floor eq `floor-bound` (15848) | −7/1200 | −2.9·10⁻⁴; this is the binding high constraint **[R]** |
| small rows eq `small-frequency-margin` (15877) | −79/800, m_small = 63/800 | −0.1625 **[R]** |
| large rows (15900) | B_0 = 53/32 | B_0 = l_x/2 + 1 + l_y = 1.7271; any ζ > 0 **[R]** |
| principal normaliser §`principal-signal`, m_w, m_z, κ_P | 23/960, 13/9600, 7/(96K) | 0.0314, 0.00134; κ_P unused **[R]** |
| Prop `detector-counts`, eq `plain-zero-loss-capacity`, eq `capacity-comparison`, 8/39 vs 7/37 | slotted counts | replaced by R_0(t), §3 **[R]** (statement restatement **[R-]**) |
| endpoint §`adaptive-endpoint`, Lemma `balanced-endpoint` (49/440640), eq `adaptive-intermediate` (R = 76/75 − 2δ/3 on d ≤ 1/2), "frequency ranges" | paper's hand certificate | replaced by the exact branch-and-bound certificate over d ∈ [d_min, h+ζ], δ ∈ [1/50, 5/6] **[R]** (exact arithmetic) |
| Prop `parameter-order` (16197) | margins ∝ Δ, mesh, K | margins are now absolute and target-independent (§6); losses chosen below them **[R-]** |
| Lemma `late-height-closure`, Prop `continuation-criterion`, Prop `hecke-dirichlet-transfer` | σ0 ∈ (1/2, 1) / [1/2, 1) | already general **[R]** |

### 5.1 Lowering the second Euler region to x_r ≥ σ0 [R-]

The proof of Lemma `local-euler` bounds H_p − 1 by the largest exponent in its table. At x_r ≥ σ0 = 0.8665, w_r ≥ 19/20, z_r ≥ 33/200:

| primes | exponents | largest | paper's value |
|---|---|---|---|
| good p ∤ u | R ~ Q^{4−6x−6z}, Q^{1−x−w−6z}, DW ~ Q^{−x−w}, DV ~ Q^{−x−6z}, VW ~ Q^{−6z−w} | −1.8065 | −363/200 |
| ramified p \| u | 1 − x − w, 3/2 − 3x, 2 − 3x − w, 2 − 4x, 5/2 − 4x − w, 3 − 6x | −0.8165 | −33/40 |

The good-prime bound needs only x_r > (5 − 6z_r)/6 = 0.668 (from the R-term). So D_2(σ0) is valid for every σ0 > 0.67.

**The principal function.** For H_η(s) = 𝓗_{η,1}(s, 1, 1/6), the defect is O(q_p^{−1−Re s}) = O(q_p^{−1.8665}) on Re s ≥ σ0. The P_0-tail is O(P_0^{−0.8065}), so c_H = 4/5 still holds. One fixed P_0, chosen before the target, gives sup_{Re s > σ0}|H_η − 1| ≤ 1/2.

**Uses of D_2(σ0).**
- The principal paths (Re s = β* + e ≥ σ0, 19/20 ≤ Re w ≤ 1 + e, 33/200 ≤ Re z ≤ 1/6 + e) lie in D_2(σ0).
- The large-row line (2, 2, v) lies in D_2(σ0).
- Nothing else uses D_2.

---

## 6. Numbers

**Exact certificate.** Command:

`exact_cert_theta.py lx=989/5000 ly=3141/5000 ell=0 sigma0=1733/2000 theta=1 dropq=1 wE=0`

| check | result |
|---|---|
| low side | 8663/10000 ≤ σ0 (margin 2·10⁻⁴) |
| small rows | −0.162452 |
| floor E(h+ζ) | −2.757·10⁻⁴ |
| high side | E < 0 on all 3055 boxes; worst box upper bound −5.03·10⁻⁴ (the exact maximum over non-floor bins is −9.0·10⁻³) |
| overall | **ALL CHECKS PASS** (30 s) |

`wE` is irrelevant here: every energy-wall piece is negative at M' = 0.826, ℓ' = 0.

**LP**, `theta_split_model.py` mode `LSmarked` with the P²/Y piece dropped (copy in `scratchpad/p2/tsm_dropq.py`). Every row except the "S > 0" row is at S = 0:

| variant | θ |
|---|---|
| S = 0 | **0.866197 = 123/142**, at (l_x, l_y) = (14/71, 134/213) |
| S > 0 | increases monotonically: 0.866620 (S = 0.005), 0.870307 (S = 0.05) |
| ideal counts R = 1 − δ | 0.866197 (counts not binding) |
| t = 1 only | 0.877911 |
| floor δ0 = 1/1000 | 0.863295 (matches the REPORT's 0.8633) |
| P²/Y kept (no LoA6) | 0.876322 |
| no theta split | 0.894737 |

The REPORT's "0.86624" is slightly above the exact slot-free optimum 123/142. This does not matter, since the certificate is at 0.8665.

**Closed form of the slot-free optimum.** Three constraints bind:
- **low:** σ0 ≥ 5/6 + l_x/6;
- **(V) cubic large sieve:** b/6 + 7M/3 − 2 ≤ 0, i.e. (13/6)l_x + (5/2)l_y ≤ 2;
- **floor bin at d = h:** σ0 ≥ (51/100)(1 − l_y) + (1 − l_x)(253/300).

Solving gives l_x = 14/71, l_y = 134/213, σ0 = 123/142. At the certified geometry the (V) slack is b/6 + 7M/3 − 2 = −7/7500.

**Loss budget.** The target-independent margins are absolute (they do not scale with Δ):

| margin | value |
|---|---|
| m_lo | 2·10⁻⁴ (+ Δ0) |
| m_hi | ≥ 2.76·10⁻⁴ (floor) |
| m_small | 0.162 |
| m_w | 0.031 |
| m_z | 1.3·10⁻³ |

The adjustable losses must total below ~10⁻⁴:
- (16 − 6l_y)e + 12ed + O(ε) ≈ 22e;
- the reciprocal U^{ε_d}, ε_c, the ζ-extension, and the normaliser power.

So e ≲ 4·10⁻⁶. This is admissible: e is any fixed number in (0, 10⁻³) and the bin count is O_e(1). It should be stated explicitly in a write-up **[R-]**.

---

## 7. Remaining gaps (high side only)

1. **[R-]** Restate Lemmas `fixed-bin-contour`, `bin-contour-accounting`, `principal-residue-interface` and `outer-row-tails` for σ0 ∈ (2/3, 1). Replace eq `local-region-two` by D_2(σ0) and D_1(3/8) by D_1(1/3). The exponent checks are in §5.1; no new idea is needed.
2. **[R-]** State the zero-slot count of Prop `detector-counts` for every t ∈ [1, 3/2]. The proof is already in the paper (§3).
3. **[R-]** Record the general-geometry form of Lemma `probe-row-norm` for M' + ℓ' ≤ 1 (low side; §5 table).
4. **[R-]** Explicit loss allocation, since the margins are about 2·10⁻⁴.
5. **Outside this check.** The low side rests on (V) and (LoA6), the hypotheses of `REPORT.md` §1, and on the theta-pairing bookkeeping. At this geometry (LoA6) is *essential*, not cosmetic: keeping P²/Y moves the slot-free optimum to 0.8763 (§5 table).

No fatal issue was found, and no constraint missed by the R4 model / exact checker was found on the high side at ℓ = 0.

---

## Referee report

*Independent skeptical check, 9 Oct 2026. Sources re-read: `paper78.tex` lines 3313–3860 (base probe, low separation, Poisson identity), 4005–4200 (Lemma `local-euler`), 4200–4760 (detector, Lemma `buffered-bins`, Prop `detector-witness`), 5505–6560 (shared principal data, Def `stage-high-data`, Lemmas `fixed-bin-contour`, `bin-contour-accounting`, `principal-residue-interface`, `outer-row-tails`), 6807–6930, 12362–12600 (Lemmas `inverse-amplification`, `plain`), 15060–15500 (Prop `detector-counts` and Remark), 15500–16460 (Part II conclusion, Lemma `high-bin`, endpoint, frequency ranges, Prop `parameter-order`). I also re-ran every script.*

### Verdict: not refuted

Conclusions (a), (b) and (c) stand. I found no constraint that the R4 model or the exact checker misses on the high side at ℓ = 0. The corrections below are minor. One of them (R3) works *in favour of* the claim.

### What I re-derived or re-ran

1. **E(d) [R].** I redid the computation in the proof of Lemma `bin-contour-accounting` (line ~6140).
   - Outside factor on (a+16e, 1−a−6e, z0): l_x(1/2−z0) + a + z0 − 1 − a·l_y − d·z0.
   - Adding d(R+δ/2) + ℓ(z0−1/2) + g and subtracting C(σ0) = σ0 + l_x/2 − 1 + h/6 gives eq `stage-high-exponent` for general (l_x, l_y, ℓ).
   - At d = h with a = (1+δ)/2, this is exactly the checker's high constraint.
   - At the paper's geometry: constant −1/48, δ-coefficient 25/96 + 13/32 = 2/3, q-coefficient ℓ = 1/6.
2. **Base-probe data [R].**
   - The Poisson identity (eq `probe-poisson-identity`, line 3831) is an exact identity for independent positive X, Y, Z.
   - The sentence after Def `stage-high-data` admits I_η(Z^{l_x}, Z^{l_y}, Z) with ℓ = 0 and 𝔥 = 𝓗_{η,u}.
   - The base residue identity 𝓗_{η,1}(s,1,1/6) = H_η(s), with A_η = 1 and 𝓡 ≡ 0, does not depend on geometry.
   - The only principal numerator row is u = 1 (line ~5505).
   - Principal raw exponent: l_x/2 + Re s − 1 + h·Re z + l_y(Re w − 1). This gives m_w = l_y/20 and m_z = h/600.
   - Small rows: z0(1−l_x) − h/6 − l_y/2 = h(z0 − 1/6) − l_y/2, because ℓ = 0.
   - Large rows: B_0 + h·v.
   - All of these agree with §2.
3. **Counts [R]; the statement needs restating [R-].**
   - Prop `detector-witness` is stated "for every t ∈ [1,3/2]". It needs only a > 51/100 and δ ≥ 1/50, and has no slot hypothesis.
   - The Remark at line 15449 gives eq `no-slot-inverse-count` for each fixed t, for 1/50 < δ ≤ 1.
   - Lemma `plain`, case 1 (z = 0), has "no restriction on the bounded nonnegative lengths … and no lower bound on the conductor".
   - The case analysis with x = 0 gives R_0(t) = max(L(t), 1 − 2δt/3), and the minimum over t is (15−16δ)/(15−6δ).
   - Supporting evidence that a zero-slot count at t ≠ 1 is already used by the paper itself: eq `large-delta-endpoint` takes t = 3/2 with no selected primes.
4. **Region two at x_r ≥ σ0 [R-].** The defect exponents −1.8065 (good primes) and −0.8165 (p | u) are reproduced. For H_η = 𝓗_{η,1}(s,1,1/6) one has V = Q^{−1} and R = Q^{3−6x}, so the defect is O(Q^{−1−Re s}).
5. **Numerics.** All re-run.
   - `exact_cert_theta.py lx=989/5000 ly=3141/5000 ell=0 sigma0=1733/2000 theta=1 dropq=1 wE=0`: ALL CHECKS PASS on 3055 boxes. Low side 0.8663, floor −2.757·10⁻⁴, worst box −5.03·10⁻⁴, small rows −0.162452.
   - `scratchpad/p2/check_slotfree.py`, `run_lp.py` and `run_lp2.py` reproduce every number in §6: 123/142; 0.866620 at S = 0.005; 0.877911 with t = 1 only; 0.876322 with P²/Y kept; 0.863295 at δ0 = 1/1000; 123/142 with ideal counts.
   - The closed-form solution of the three binding constraints, (14/71, 134/213, 123/142), checks out exactly.

### Corrections and remarks

- **R1 (gap 3 is moot) [R].** At K = 0 there is no compensated row, so Lemma `probe-row-norm` is not needed.
  - The base-probe row norm is Part I's Lemma `unmarked-completed-moment` (line 3139) with (M, N) = (l_x + l_y, 1).
  - It gives Σ|B|² ≪ Z^{O/2 + max(M−O, 2M−O−1) + ε} = Z^{M − O/2 + ε} ≤ Z^{M+ε}, because M = 0.826 < 1.
  - So no restatement is required. Gap 3 in §7 can be deleted.
- **R2 (absolute margins are not needed).** Every comparison is against C(β*) = C(σ0) + Δ0:
  - low side: ω < Δ0;
  - high side: eq `saving-criterion` has right side (β* − σ0) − ε_hi.

  Prop `parameter-order` chooses e, ζ, the amplitude width and all losses *after* Δ ("All these choices depend only on Δ"). Δ0 = β* − σ0 is a single global constant (β* is the global supremum), so it is not target-dependent. Two consequences:
  - The "e ≲ 4·10⁻⁶" requirement in §6 and gap 4 in §7 are not real constraints. Losses below Δ0/2 suffice.
  - Separately, applying Lemma `bin-contour-accounting` literally (ε_c ≥ 12e entering as (1+d)ε_c) gives an e-coefficient of about (16 − 6l_y) + 12(1+d) ≈ 34, not 22. This is immaterial, both for the reason above and because Part I notes that the Z^{ε_c} enlargement is optional.
- **R3 (sharpening, in the claim's favour) [R-, under the same (V) and (LoA6)].** By R2, the slot-free argument closes at σ0 = 123/142 ≈ 0.866197 exactly, at geometry (14/71, 134/213, 0).
  - Exact check at that point:

    | quantity | value |
    |---|---|
    | (low − σ0) | 0 |
    | (V) slack b/6 + 7M/3 − 2 | 0 |
    | floor E(h) | 0 |
    | max non-floor E(h) | ≈ −0.00876 (at δ ≈ 0.113) |
    | small rows | −0.1628 |
    | m_w | 0.0315 |
    | m_z | 0.00134 |

  - Take ζ below a fixed fraction of Δ0, as the paper does (ζ < m_high/16), and losses below Δ0/2. The continuation criterion then contradicts β* > 123/142.
  - The exact checker requires strict negativity, so it cannot certify this point as written.
  - The REPORT's "0.8665 (exactly certified)" is therefore conservative. Its "0.86624" LP figure should read 123/142 = 0.866197.
- **R4 (range nit).** D_2(σ0) needs the good-prime R-term exponent 4 − 6x_r − 6z_r < −1 at z_r = 33/200, i.e. σ0 > 401/600 ≈ 0.6683. Gap 1's "(2/3, 1)" should read "(401/600, 1)". This is irrelevant at 0.8665.
- **R5 (κ outside Lemma `plain`'s stated range).** The slot-free bootstrap allows β* ∈ (σ0, 7/8), so κ = 2β* − 1 ∈ (0.733, 0.75). That is below the stated hypothesis 3/4 ≤ κ ≤ 1 of Lemma `plain`.
  - This is harmless: only case z = 0 is used, and its bound does not involve κ. The lemma can be applied with κ = 1, which needs no zero-free hypothesis.
  - A write-up should say this explicitly. §5's "no κ" covers it only implicitly.
- **R6 (REPORT inconsistency).** REPORT §0 says "t = 1 alone … 0.8751". The current scripts (LSmarked, P²/Y dropped, S = 0, tmax = 1) give 0.877911, as stated here. Both are above Liu, so the conclusion is unchanged, but the REPORT figure should be updated.

### Attempts to break the high side that failed

- **Any use of h = 13/16, M + ℓ = 1, 17/48 or 23/48 in the base-probe high side.** None. The identity is exact for any (X, Y, Z), and E(d), the outer rows and the principal margins are all general.
- **Floor bin.** O(U) physical rows give R = 1 with no witness, unchanged at any geometry.
  - Bins just above the floor (δ = 1/50 + 2e) have actual witnesses and use R_0(δ) ≈ 0.9866. The Remark requires 1/50 < δ, which holds.
- **Buffered moderate range.** It needs d ≤ h + ζ < 1 (line ~16395). Here h + ζ = 0.8023.
- **Detector height ceiling.** It needs T_1 ≤ U^{1/100} on d ≥ d_min. This is unchanged.
- **D_1 contours.** Re w = 1 − a − 6e ≥ 1/12 − 6e > −1/100 because a ≤ β* ≤ 11/12; this is geometry-free.
- **Small rows.** The small-row line lies in D_1(1/3), since β* + e + 1/2 > 4/3.
- **Final steps.** The continuation criterion (σ0 ∈ (1/2,1)) and Prop `hecke-dirichlet-transfer` (σ0 ∈ [1/2,1)) are already general.

### Not checked here

The low-side hypotheses (V) and (LoA6), and the theta-pairing bookkeeping, are outside the scope of this check. At the optimum, the (V) cubic-large-sieve condition binds with zero slack (−7/7500 at the certified point). So the improvement past Liu rests entirely on those hypotheses, as `REPORT.md` §1 states. Nothing on the high side stands in the way.
