# Agent `combo`: one consolidated exponent model of the weighted high side (10 Oct 2026)

Everything here is under `agents/combo/`. Labels:
- **[P]** proved (relative to P78, the 08 lemma, Collet's PR #9, and the session-10 agent lemmas);
- **[M]** an exact statement about the exponent model;
- **[C]** conditional;
- **[H]** heuristic.

## 0. Verdict

1. **Model.** `combo_model.py` contains every knob found this session, in one model.
   - Its defaults are the proved inputs.
   - Its count LP is solved exactly. It matches the closed-form R of 08 and of Collet to within 10⁻⁹.
   - It reproduces the 08 certificate (−1.29·10⁻⁴) and Collet's boundary b₀ (to 10⁻⁹ at her geometry).
   - It reproduces every agent number, up to the agents' grid artefacts (§2).
2. **Best θ with proved inputs.** It is **θ = 0.8748866**, exactly Collet's b₀ = 0.874886647 [M].
   - κ = 2θ − 1 (variable κ) is worth only 3·10⁻⁷.
   - The other proved additions are worth exactly 0: slot powers (unlimited marked supply), the sextic large sieve for short plain witnesses, and the per-bin choice of E_old vs E_wt. Lemma D has no model effect on its own.
   - The floor bin and the dyads d < h never bind.
   - The gap to 41/47 is **Δθ = 0.0025462**.
3. **Frontier.** The levers are almost exactly additive (§5).
   - Pairs of power-saving levers lie on Σf = 1. Here f is the fraction of a lever's standalone 41/47 requirement.
   - Pairs involving α, η or g are sub-additive, with Σf up to 1.75.
   - No pair is super-additive by more than 3%. **There is no cheap synergy.**
   - The cheapest spread uses all seven active levers at **≈ 17 % of their standalone requirement each** (sM, eP, eL, sLS, w, g, α). The five additive levers alone need ≈ 20 % each.
   - Equivalently, a common power saving **ε ≈ 0.0099** is needed simultaneously in five lemmas: lem:marked cross terms, lem:plain edge, amplified inverse moment, restricted large sieve, and wall.
4. **Provability.** **None of the partial savings is provable by known techniques** (§6).
   - Even a hypothetical ε = 10⁻³ in the most efficient lemma (restricted large sieve on the bad rows) gives Δθ = 9·10⁻⁵, which is 3.6 % of the gap.

## 1. The model (`combo_model.py`)

### Geometry and low side

- Parameters: slot length ℓ; M = l_x + l_y; b = l_y − l_x; h = 1 + ℓ − l_x; U = Z^d with d = h; supply L = ℓ/h.
- Low side:

      B_low = 5/6 + M/12 − ℓ/6 − b(1−g)/12 + (M+ℓ−1−w)₊/6.

  We take the relaxed wall M = 1 + w − ℓ. Then B_low = 11/12 + w/12 − ℓ/4 − b(1−g)/12.
- Wall optimality [M] (`check_wall.py`): on every line B_low = T, moving M below the relaxed wall is worse. This holds for the proved inputs, for g = 0.74, for w = 0.036 and for sLS = 0.027.
- Side conditions, taken from the 08 certificate (`side_ok`):
  - ℓ ≤ 1/5; l_x, l_y > ℓ;
  - Gram domination l_y − ℓ − (2 − g/6)b > 0;
  - 5ℓ > h; 2l_y > h; l_y < 1/2;
  - small rows;
  - κ ≥ 7/10.

### High side

For a bin (δ, x), with a = (1+δ)/2 and q = δx:

    E + B = a + h(z_c − 1/6) − ℓ/2 + qℓ − l_y/2 + d(R + σ − z_c),   z_c = 17/50,

where σ is the smallest of the following bounds:

| term | formula | status |
|---|---|---|
| σ_pw | δn₀/2 | [P], ≡ E_old |
| σ_4 | (1−R)/4 − q z_N/2 | [P], the 08 weighted numerator |
| σ_6 | (1+η−R)/6 | sixth moment, lever η |
| σ_LS | (max{R, l_y/d, (1+R)/2 − q z_N − sLS} − R)/2 | restricted large sieve, lever sLS |

Here n₀ = max(0, 1 − l_y/d) and z_N = min(L, c_N(1 − 2n₀ + 6eP)).

The **floor bin** (δ = 1/50, R = 1, q ≤ 1/100) is included; it never binds. We also checked the dyads d = 0.5h … h (`check_d.py`): d = h binds in every configuration, with margins ≥ 0.013 at d = 0.95h (Collet: convex in d).

### Count

    R = max(1−δ, min_{t∈[1,3/2]} max_{r+m=t, m≤1/2} min(W_wit(r), W_pl(m))).

The witness lines are:

- **inverse-marked, r < 1:** 1 − δr − 2q z_M, with
  z_M = min(supply_M, mk(1 + s_M,eff − r), (3 + 4sM − 2r)/8 [cubic branch]),
  where s_M,eff = max(sM, w/2). The wall feeds lem:marked as r + 2z ≤ 1 + w/2. A cross-mark saving sM lowers the P78 excess e = max(0, r+2z−1, (2r+8z−3)/4) by sM (xmark §1).
- **long inverse, r ≥ 1:** max(1, 1 + α(r−1) − eL) − δr.
- **plain:** 1 − 2δm − 2q z_P, with z_P = min(L, max(c_P(1 − 2m + 6eP), 1/2 − 2m [sextic LS, P])).

The rates are c_P = c_N = 1/(6κ) with κ = min(3/4, 2θ − 1) (Collet's variable κ ∈ [7/10, 3/4], applied dynamically by a fixed-point iteration).

**Exact solution.** The inner maximum (the adversary's split) is computed exactly by bisection at the crossing of a non-increasing and a non-decreasing line. Monotonicity holds for mk ≤ 1 and c_P ≤ 1, and is asserted in the code. The outer minimum over t uses a grid plus golden refinement; any t is admissible, so this can only overestimate R.

The result agrees with the closed form R = 1 − δ + (α−δ)δP_x/(2J) (08/Collet) to within 2·10⁻⁹ on all bins (`test_count.py`). It also agrees for α = 0.7 where δ < α; where δ > α the closed form is invalid, and the LP gives the DH floor.

**Agent grids.** The agents' grid LPs, with 26–201 t-points, overestimate R by up to 7·10⁻⁴. All corrections in §2 come from this.

### Knobs and defaults

The default `PROVED` = `Lev()` has:

- mk = 1/2 and the cubic branch kept;
- `mark_supply='unlimited'` (slot-power lemma, amplitude Prop 3.3 [P]);
- c_P = c_N = 1/(6κ), joint; eP = 0;
- `plainLS=True` [P, plain-cap P2];
- α = 5/6, eL = 0;
- η off, sLS = 0, E_old on;
- w = 0, g = 1, Gram domination enforced;
- κ dynamic.

`PAPER08` is the 08 setting: κ = 3/4, physical supply, no LS.

### Search

- `theta_star(lev)`: for each b, take ℓ_T(b) from B_low = T. Then compute min_b B_high(ℓ_T(b), b) − T, using a grid plus golden search in b. B_high increases in ℓ, as checked. Bisect in T, with an outer κ fixed point.
- `closes41(lev)`: the same test at T = 41/47 with κ = 35/47. A negative value means θ < 41/47.
- Bins:
  - search: δ step 0.01 over [0.02, 0.75], x ∈ {0, .1, .2, .3, .4, .45, .5};
  - mid: δ step 0.0025;
  - fine: δ step 0.0025, x step 0.025.

## 2. Validation

| check | target | this model |
|---|---|---|
| 08 certificate, ℓ = 0.167, b = 0.086, κ = 3/4 (`validate08.py`) | margin ≈ −1.3·10⁻⁴ | **−1.29·10⁻⁴** at δ = 0.470, x = 1/2, R = 0.58970 (all bin sets) |
| Collet b₀ at her (ℓ₀, x₀, y₀, κ₀) (`validate_collet.py`) | 0.8748866470 at δ₀ = 0.469345 | **0.8748866479** at δ = 0.469355 |
| proved θ* (`proved2.py`) | ≈ 0.87489 | **0.87488661**: b = 0.087, ℓ = 0.167120, κ = 0.749773 |

Agent values (`validate_agents.py`, `validate_thresholds.py`; κ = 3/4, ℓ = 0.1773, B = 41/47, exact LP). Agent values are in parentheses:

| quantity | this model (agent value) |
|---|---|
| paper | **+0.01109** (+0.01119 amplitude, +0.01124 marked-cap, +0.0117 coarse); closed form on a fine grid also gives +0.01109 |
| mk = 1, physical supply, cap off | +0.00166 (+0.00167 / +0.00187) |
| mk .8 + cP .4 [cN fixed] | +0.00017 (+0.00032) |
| mk .8 + cP .45 | −0.00038 (−0.00018) |
| mk sat + cP .5 | −0.00102 (−0.00085) |
| mk .75 + α .75 | −0.00119 (−0.00088) |
| mk sat + α .78 / .76 | −0.00022 / −0.00093 (−0.00004 / −0.00070) |
| mk sat + cP = 1 | −0.00870 (−0.00828) |
| r+2z ≤ 1.13 / 1.045 | +0.00201 / +0.00786 (+0.00231 / +0.00794) |
| α = 0.7 | +0.00569 (+0.0057) |
| DH | −0.02837 (−0.028) |
| mk = 1, unlimited supply, cap kept | +0.00412 (+0.0043) |
| η = 0 | −0.0124 (−0.0116) |
| ideal sLS | −0.02866 (−0.029) |

Thresholds:

| threshold | this model (agent value) |
|---|---|
| mk* (unlimited supply, cap off) | **0.741** (0.744) |
| joint c* = c_P = c_N | **0.655** (plain-cap 0.73) |
| c* with mk .75 | **0.327** (0.37) |
| η* | **0.097** (numerator 0.0885) |
| w* alone | 0.0407 (0.041–0.042) |
| w* with marked feed | 0.0359 (0.037–0.038) |
| α* | 0.583 |
| sLS* | 0.0275 |
| cross-mark δ* | 0.155 (xmark 0.159) |

**Corrections.** Two thresholds were overstated because of the t-grid:
- **c\* = 0.655, not 0.73.** Checked by an independent 2001×2001 brute force in `check_c.py`.
- **η\* = 0.097, not 0.088.**

The qualitative conclusions of the agents are unchanged.

## 3. Best θ with PROVED inputs (`proved2.py`, `singles.py mid`)

| inputs | θ* |
|---|---|
| 08 inputs, κ = 3/4 fixed, exact optimum (no 10⁻⁴ margin) | 0.8748869 |
| + κ = 2θ−1 (variable κ) | 0.8748866 |
| Collet model (weighted branch only, κ dynamic) | 0.8748866 |
| ALL PROVED (+ slot powers, sextic LS short plain, E_old per bin, floor) | **0.8748866** (b = 0.087, ℓ = 0.16712) |

- The improvement of PR #9 over 10499/12000 is almost entirely the exact optimisation: ℓ goes from 0.167 to 0.16712 with zero margin. Variable κ adds 3·10⁻⁷.
- No proved item from this session changes θ. The binding bin is δ ≈ 0.469 with x = 1/2. There:
  - the counts are capacity-limited, not supply-limited;
  - the plain witness is long (2m ≈ 0.83 > 5/14);
  - E_wt < E_old.

## 4. Single small hypothetical improvements (θ*, mid bins; gap to 41/47 = 0.0025462)

| lever (proved value) | θ* at small values | rate dθ/dv near proved | standalone 41/47 threshold |
|---|---|---|---|
| sM, cross-mark saving in lem:marked (0) | .001: .8748702; .005: .8748047; .01: .8747229; .02: .8745590; .05: .8740667 | 0.0164 | **0.155** |
| mk, marked rate, cap off (0.5) | .55: .8744464; .60: .8739009 | 0.0089 | 0.740 |
| eP, power saving at lem:plain edge, joint (0) | .001: .8748440; .005: .8746735; .01: .8744605; .02: .8740352 | 0.0427 | **0.060** |
| c = c_P = c_N rate (1/(6κ) = 0.2238) | .25: .8747451; .30: .8744870 | 0.0052 | 0.655 |
| α, long-witness slope (5/6) | −.001: .8748792; −.01: .8748117; −.02: .8747344; −.05: .8744871 | 0.0075 | 0.583 |
| eL, saving beyond the diagonal, U^{max(1,(1+5r)/6−eL)} (0) | .001: .8748414; .005: .8746620; .01: .8744394; .02: .8739998 | 0.0449 | **0.059** |
| η, sixth-moment excess (0.34 provable) | ≥ .19: no effect; .18: .8747301; .15: .8739282; .10: .8724261 | 0 above 0.194 | **0.097** |
| sLS, restricted LS saving on bad rows (0) | .001: .8747951; .005: .8744282; .01: .8739670; .02: .8730361 | **0.0918** | **0.0274** |
| w, lossless wall, + r+2z ≤ 1+w/2 (0) | .001: .8748157; .005: .8745319; .01: .8741769; .02: .8734663 | 0.0710 | **0.0358** |
| g, Gram factor (1) | .99: .8748287; .95: .8744848; .90: .8739472 | 0.006–0.009 | 0.7385 |

Notes on the table:
- **Linearity.** The five power-saving levers (sM, eP, eL, sLS, w) are linear in θ up to their thresholds: rate × threshold = 0.00254–0.00265.
- **Small savings.** A saving ε = 10⁻³ in any single lemma gains at most 9·10⁻⁵ (sLS).
- **Unit conversions.**
  - sLS ≡ a uniform saving sLS/2 in the bad-row numerator exponent σ. 41/47 needs Σ_{u∈C}|S̃_u| to be U^{0.0137} below the Hölder/4th-moment bound.
  - η only acts below its activation point 0.194, where σ₆ < σ₄ at the binding bin. Above it η is inert, so improvements of the provable 0.34 down to 0.2 are worthless.
  - eP at the binding plain split (1 − 2m ≈ 0.16) multiplies the plain capacity by 1 + 6eP/0.16, and the numerator capacity (1 − 2n₀ ≈ 0.08) by 1 + 6eP/0.08.

## 5. Pareto frontier of minimal joint requirements (θ < 41/47, κ = 35/47)

Definitions:
- **f-units.** The value of lever i is p_i + f_i(t_i − p_i), where p_i is the proved value and t_i is the standalone threshold of §4.
- **Scales.** For sM the scale is 0.155. For η the anchor is the activation point 0.194.
- **Scripts.** `frontier.py pairs|subsets` and `epsilon.py`.

### 5.1 Pairs

For each pair, f_j is the minimum needed when f_i = 0.25 / 0.50 / 0.75. Σ is f_i + f_j; Σ = 1 means additive.

| pair (i, j) | f_j at f_i = .25 / .50 / .75 | Σ |
|---|---|---|
| sM, eP | .750 / .500 / .250 | **1.00** |
| sM, sLS | .750 / .500 / .250 | **1.00** |
| eP, sLS | .750 / .500 / .250 | **1.00** |
| sLS, w | .750 / .500 / .250 | **1.00** |
| eP, w | .745 / .494 / .245 | **0.995** |
| w, g | .721 / .503 / .277 | 0.97–1.03 |
| sM, w | .847 / .563 / .280 | 1.03–1.10 |
| sM, g | .766 / .602 / .374 | 1.02–1.12 |
| sM, η | .783 / .571 / .336 | 1.03–1.09 |
| sLS, g | .792 / .584 / .337 | 1.04–1.09 |
| α, sLS | .816 / .578 / .290 | 1.04–1.08 |
| α, w | .805 / .560 / .278 | 1.03–1.06 |
| eP, α | .815 / .602 / .350 | 1.07–1.10 |
| sM, α | .835 / .646 / .407 | 1.09–1.16 |
| eP, g | .830 / .644 / .428 | 1.08–1.18 |
| η, w | .837 / .568 / .284 | 1.03–1.09 |
| eP, η | .927 / .854 / .781 | 1.18–1.53 |
| α, η | .903 / .807 / .708 | 1.15–1.46 |
| α, g | .886 / .793 / .719 | 1.14–1.47 |
| η, sLS | 1 / 1 / 1 | **no joint gain** (alternative numerator bounds: σ = min(σ₆, σ_LS)) |
| η, g | 1 / 1 / .997 | **no joint gain** (η binds at b = 0, g needs b ≈ 0.15) |

Pairs with eL (`pairsA3.out`):

| pair (i, j) | f_j at f_i = .25 / .50 / .75 | Σ |
|---|---|---|
| sM, eL | .750 / .500 / .250 | 1.00 |
| eL, sLS | .750 / .500 / .250 | 1.00 |
| eL, eP | .750 / .500 / .250 | 1.00 |
| eL, w | .740 / .487 / .240 | 0.99 |
| eL, g | .779 / .651 / .512 | 1.03–1.26 |

The sM pairs above use the cubic branch relaxed by 4sM (`pairsA2.out`). With the cubic branch kept fixed, the sM pairs at f_sM = 0.75 lose about 0.1 (`pairs_A.json`).

**Structure** [M]:

1. **Additive group.** The count capacities (sM, eP), the numerator saving (sLS), the wall (w) and eL act on the same binding bin (δ ≈ 0.47, x = 1/2) through linear terms. They add exactly.
2. **Sub-additive levers.**
   - α changes which witness binds (the long/short balance), so its effect saturates against the other count levers.
   - g only pays at large b, where the other levers are weaker.
   - η and sLS are alternative bounds for the same σ.
3. **No super-additivity.** No pair beats Σf = 0.97. The earlier "mk sat + small cP" closures were artefacts of the physical-supply cap. With slot powers (proved), the marked lever is linear like the others.

### 5.2 Uniform fractions over subsets

f* is the smallest f such that every lever in the subset at fraction f closes 41/47. The complete lists are in `subsets_*.json`.

Cheapest subsets of each size. Values are in natural units; f* has precision ±0.002, from re-verification on fine bins.

| size | cheapest subsets (f*) | lever values at f* |
|---|---|---|
| 1 | sLS (1), w (1), … | sLS 0.0274 / w 0.0358 / eP 0.060 / eL 0.059 / sM 0.155 / η 0.097 / α 0.583 / g 0.739 |
| 2 | eP+w **0.497**; eP+sLS, sLS+w, sM+eP, sM+sLS 0.500 | e.g. sLS 0.0137 + w 0.0179 |
| 3 | eP+sLS+w **0.333**; sM+eP+sLS 0.334; sM+eP+w, sM+sLS+w 0.345 | eP 0.020, sLS 0.0091, w 0.0119 |
| 4 | eP+eL+sLS+w **0.249**; sM+eP+sLS+w 0.257; eP+α+sLS+w 0.267 | eP 0.0149, eL 0.0147, sLS 0.0068, w 0.0089 |
| 5 | sM+eP+eL+sLS+w **0.204**; sM+eP+α+sLS+w 0.219 | sM 0.032, eP 0.012, eL 0.012, sLS 0.0056, w 0.0073 |
| 6 | sM+eP+eL+sLS+w+g **0.185**; sM+eP+α+sLS+w+g 0.199 | … g 0.952 |
| 7 | sM+eP+eL+sLS+w+g+α **≈ 0.167** | sM 0.026, eP 0.010, eL 0.010, sLS 0.0045, w 0.0059, g 0.957, α 0.792 |
| (η) | any subset containing η: η is inert above 0.194, so it either adds nothing (η at f ≲ 0.3 in large subsets) or makes f* large (α+η 0.72, eP+η 0.77, η+g 0.92, η+sLS 0.96) | — |

Findings:
- f* ≈ 1/k for k levers from the additive group {sM, eP, eL, sLS, w}.
- Adding the sub-additive levers α and g still helps, but less than proportionally.
- **No combination of partial savings below ≈ 17 % of each standalone requirement closes 41/47.**

### 5.3 A common power saving ε in natural units

The levers sM, eP, eL, sLS and w are all power savings, in U- or Z-exponents. The smallest common ε that closes 41/47:

(`epsilon.py`; search bins; the final value was re-verified on fine bins.)

| levers sharing ε | ε* | linear prediction 0.0025462 / Σ rates |
|---|---|---|
| sLS | 0.0274 | 0.0277 |
| w | 0.0358 | 0.0359 |
| eL / eP | 0.0591 / 0.0600 | 0.0567 / 0.0596 |
| sM | 0.155 | 0.155 |
| sLS + w | 0.0155 | 0.0157 |
| eP + sLS / eL + sLS | 0.0188 / 0.0187 | 0.0189 / 0.0186 |
| sM + sLS | 0.0233 | 0.0235 |
| eP + w / eL + w | 0.0223 / 0.0220 | 0.0224 / 0.0220 |
| eP + eL | 0.0298 | 0.0291 |
| sM + w / sM + eP | 0.0321 / 0.0433 | 0.0291 / 0.0431 |
| eP + sLS + w / eL + sLS + w | 0.0123 / 0.0122 | 0.0124 / 0.0123 |
| sM + sLS + w | 0.0148 | 0.0142 |
| eP + eL + sLS + w | 0.0102 | 0.0102 |
| sM + eP + sLS + w | 0.0118 | 0.0114 |
| **all five** | **0.0099** (fine bins: 0.00984 → +2·10⁻⁷, 0.0100 → −1.8·10⁻⁴) | 0.0095 |

θ-rates at the proved point, per unit of saving: sLS 0.0918, w 0.0710, eL 0.0449, eP 0.0427, sM 0.0164.

**So a common saving of U^{0.01} (Z^{0.01} for the wall) is needed at once in all five lemmas.** These are:
- lem:marked cross terms;
- lem:plain at its edge;
- the amplified inverse moment beyond the diagonal;
- the restricted large sieve on the bad rows;
- the reflected-energy column term beyond the wall.

The frontier is the hyperplane Σ_i rate_i·ε_i ≥ 0.00255, to within 3 %. In raw exponent units the "cheapest" single input is the restricted large sieve: sLS = 0.0274, i.e. a saving of U^{0.0137} in Σ_C|S̃_u|.

## 6. Which partial savings are plausibly provable? (honest assessment)

The cheapest spreads need every lever at ≈ 17–25 % of its standalone requirement. In absolute terms that is:
- sM ≈ 0.026–0.04;
- eP ≈ 0.010–0.015;
- eL ≈ 0.010–0.015;
- sLS ≈ 0.0045–0.007;
- w ≈ 0.006–0.009;
- α ≈ 0.77–0.79;
- g ≈ 0.95.

For each one, the nearest known technique and the exact place it fails:

| lever | what a positive saving means | nearest technique | why it gives no ε > 0 | verdict |
|---|---|---|---|---|
| sM | Cross-mark (H2-mk) decorrelation: the off-tuple part E_off of lem:marked's canonical energy is ≪ Z^{F+e−ε}. | Cubic large sieve (Heath-Brown), de Faveri, Heath-Brown–Patterson / Dunn–Radziwiłł (Gauss sums at primes). | P78's bound already equals the sum of single-mark bounds, because the cubic sieve returns exactly the P-diagonal (marked-cap §2.5). A coefficient-blind saving is false: the Gauss-sum counterexample lives in E_off. The pair forms land inside the quadratic sieve beyond its edge (xmark §3–4). Numerically the cross terms are random-size, with no main term to exploit. | not provable |
| eP | lem:plain beyond Theorem G: a twisted fourth moment of self-dual sextic sums with a power-length prime twist at the edge (WR(λ < κ)). | Twisted fourth moments of quadratic families (Li; Shen–Stucky). | Theorem G: the move calculus {reflect, transform, remove, terminal} is exactly optimal. Any ε needs a new move. Even for quadratic families only log savings are known at the edge. | not provable |
| eL / α | Möbius mean square over the sextic family beyond the diagonal, with a power saving (M_κ type, "Lindelöf on average for 1/L beyond the edge"). | Sixth-power amplification; large-values methods (Huxley, Guth–Maynard, Halász–Montgomery, Heath-Brown). | Copies are exactly u·(F^×)⁶ (Kummer), so slope 5/6 is forced (09 §28). All large-values methods were killed for this family in session 09. | not provable |
| sLS | Σ_{u∈C}\|S̃_u\|² ≪ U^{(1+R)/2−qz_N−ε} on the bad rows C. | Duality / Halász–Montgomery; Burgess on pair characters; sparse large sieves. | Every black-box input admits a flat configuration that saturates the 4th-moment bound (numerator §2). Pair sums have length U^{0.42–0.46}, which is below Burgess' range U^{1/2} for conductor U². Even granting Lindelöf for every pair, Halász gives σ ≤ n₀/4 > σ_wt. A saving must use the zero structure of C, which is a zero-density-type input. | not provable |
| w | A lossless wall offset W(w): the reflected-energy column term stays at Z^F up to M + ℓ ≤ 1 + w. | Power-saving asymptotics in affine moment Dirichlet series; GL2 quadratic-twist second moments. | A *relative* saving beyond the wall pays nothing unless η ≳ 0.85 (09 §14). A lossless offset is a Conrey–Rodgers 2-swap statement for a thin family; only log savings are known (Li 2024, Shen–Stucky). | not provable |
| η | Sixth-moment Lindelöf on average at length 3n₀ ≈ 1.3–1.4. | Sextic large sieve; pointwise β*. | Provable 0.34; a perfect large sieve gives 3n₀ − 1 ≈ 0.26–0.39. Both lie above the activation point 0.194, so small savings from what is provable are inert. A proof needs the d₃ structure. | not provable (and inert) |
| g | Shrinking the Gram P_a^{1/6} term. | Theta split / (V′). | The Gram term carries the signal's own Mellin coefficient Ŵ₁(1)Ŵ₀(5/6). Reducing it is the (V′) pairing problem, which is signal-circular (04/07). It is probably not even true for this probe. | not provable |
| κ | Already dynamic (κ = 2θ − 1). | — | A smaller κ in lem:plain needs a zero-free region beyond β* on the transformed families, which is circular. | — |

**Conclusion.**
- No lever admits a provable ε > 0 with known techniques.
- Even granting a hypothetical ε = 10⁻³ in one lemma, θ moves by ≤ 9·10⁻⁵.
- Because the frontier is additive, tiny savings cannot be pooled to reach 41/47. The total "θ budget" Σ_i rate_i·ε_i ≥ 0.00255 must be paid in full. With every lever in play, this needs ε ≈ 0.01 in five independent beyond-the-edge statements at once (§5.3).

**Distance to the conjectured truth.** This is heuristic, but it says where the needed saving is a small fraction of what is expected to be true:
- sLS needs about 14 % of the ideal Bessel bound (sLS_ideal ≈ 0.20 at the binding bin);
- w needs w = 0.036, while the 09 numerics show W(0.064) holding with a margin of about Z^{0.43};
- eP needs about 43 % of WR(0);
- sM needs about 50 % of square-root cross-mark cancellation;
- η needs 50 % of its active range.

So the restricted large sieve and the wall are the inputs that are "least far from the truth", but they are exactly the beyond-the-edge statements catalogued in docs/OBSTRUCTIONS.

## 7. Files

| file | content |
|---|---|
| `combo_model.py` | the consolidated model: `Lev`, `count_R`, `B_low`, `B_high`, `margin_T`, `theta_star`, `closes41`, `side_ok`, `dcheck` |
| `levers.py` | lever builders on top of `PROVED` |
| `test_count.py`, `validate08.py`, `validate_collet.py`, `validate_agents.py`, `validate_thresholds.py`(.out), `check_c.py` | validation |
| `proved2.py`, `proved_theta.py` | proved θ* |
| `singles.py` → `singles_mid.out/json` | §4 |
| `thresholds41.py` → `thresholds41.json` | standalone thresholds |
| `frontier.py` → `pairs_*.json`, `subsets_*.json`, `pairs*.out`, `subsets*.out`; `subsets_eL.py` → `subsets_eL.json` | §5.1–5.2 |
| `epsilon.py` → `epsilon*.out/json` | §5.3 |
| `check_wall.py`, `check_d.py`, `verify_final.py` (+ `vf2.json`, `vf3.json` → `verify_final.out`) | wall optimality, d-range and floor checks, final re-verification on fine bins |
