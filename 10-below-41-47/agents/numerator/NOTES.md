# Agent `numerator` — can the bad-row numerator beat the fourth-moment Hölder bound?

Session of 10 Oct 2026 (afternoon). Scope: the weighted-numerator high side of `BRIEF.md`.
Scripts in this directory: `gmodel.py`/`wmodel.py` (copies), `jmodel.py` (joint witness/numerator
level-set model with hypothetical inputs), `scan_cN.py`, `scan_eta6.py`, `eta_thresh*.py`, `scan_np3.py`,
`denmodel.py`, `barrier_check.py`.

**Verdict.** No rigorous improvement found. I prove (at the exponent level; §2) that the current
numerator exponent σ_wt is *optimal* for every argument built from the available black-box inputs
(lem:plain in any admissible form, joint moments with the plain witness, lem:marked/amplification for
the witnesses or the denominator, the bin's pointwise bounds, the sextic large sieve, Hölder/level sets),
because an explicit "flat" configuration satisfies all of them with σ = σ_wt. Every route listed in the
mission reduces to an input of degree ≥ 6 in plain polynomials at total length 3n₀ ≈ 1.34–1.39, i.e.
beyond the large-sieve edge. The sixth-moment threshold is **η ≤ 0.088** (not 0.10 as in BRIEF).

---------------------------------------------------------------------------------------------------

## 0. Exponent dictionary (U-units, U = Z^d, d = h)

Binding point (ℓ = 0.1773, B = 41/47): b ∈ [0, 0.075], δ ≈ 0.45–0.55, x = 1/2, q = δx ≈ 0.22–0.28,
R ≈ 0.51–0.61, n₀ = 1 − l_y/d ≈ 0.448 (b = 0.05) / 0.463 (b = 0), L = ℓ/h ≈ 0.224,
z_N = (1 − 2n₀)/(6κ) ≈ 0.023, κ = 3/4.

For a witness class C (count U^{R_s}), write λ for the exponent of |S_u| := |S_{χ̄_u}(n₀)| on C.
Inputs, each a linear constraint on (count, λ):
- (P1) lem:plain, zero slots: Σ_u |S(n₀)S(n')|² ≪ U for every n' (reflection makes it length-free).
- (P2) lem:plain, slots: Σ_u |S(n₁)S(n₂)Q|² ≪ U if n₁+n₂+6κz ≤ 1. Used with (n₀,n₀): z_N; with (n₀,0): z₁=(1−n₀)/(6κ).
- (PW) joint numerator × plain witness S_m (|S_m|² ≥ U^{δm} on C), capacity (1−n₀−m)/(6κ).
- (B) bin pointwise: |S_ψ(n)|² ≤ U^{δ min(n,1−n)} (eq. old-6.3), so λ ≤ δn₀/2.
- (M) lem:marked / lem:inverse-amplification for Möbius polynomials (witnesses, or the denominator).
- (LS) sextic large sieve: Σ_u|Σ_{k≤K} c_k ψ_u(k)|² ≪ (U + K + (UK)^{2/3})‖c‖².

Level-set (dyadic in λ) + Hölder: Σ_C|S_u| ≤ max_λ U^{λ + min(R_s, 1−4λ−2qz_N, 1−2λ−2qz₁, …)}.
At the binding classes the 4th-moment line binds: σ_wt = (1−R)/4 − q z_N/2.

## 1. Model numbers (jmodel.py; per witness split, adversary maximises R_s + σ, we choose t)

| input (status) | effect at ℓ = 0.1773 (worst min(E_old,E_wt), b ∈ [0,0.15]) |
|---|---|
| paper (num4 + pointwise) — PROVED | +0.0122 (b=0.05; finer δ-grid: +0.0132 at b=0, δ=0.55) |
| N6: Σ_u|S(n₀)|⁶ ≪ U^{1+η} — OPEN | needs **η ≤ 0.0885** (b=0), 0.072 (b=0.05), 0.063 (b=0.075); η=0: −0.0116 |
| NP3: Σ_u|S(n₀)|⁴|S_m|² ≪ U^{1+η'} (same presentation) — OPEN | needs η' ≤ 0.086 (b=0) |
| NM3: Σ_u|S(n₀)|⁴|M_r|² ≪ U^{1+η''} — OPEN | needs η'' ≤ 0.25 (b=0); total length 2n₀+r up to 2.4, trivial excess ≥ r |
| NM: Σ_u|S(n₀)M_r Q|² ≪ U for n₀+r+cz ≤ 1 (c=2 or 1) — OPEN | **zero gain** (adversary's split has r ≥ 0.62 > 1−n₀ ≈ 0.55 at every t) |
| NM extended past the diagonal by sixth-power amplification | zero gain |
| numerator prime-capacity rate c_N (paper 2/9) | c_N = 1: +0.0032; needs c_N ≈ 1.4 (gmodel) |
| denominator bounded by lem:inverse-amplification moment instead of pointwise Z^a (PROVABLE type) | zero gain at binding bins (e(1/d) − δ/d ≥ R there: 0.651 vs 0.611 at b=0.05, δ=0.45) |

Correction to BRIEF: the `N6` knob in `model/gmodel.py` sets z_N = 0 and σ = (1−R)/6, i.e. it assumes
**η = 0**; with η = 0.10 the worst exponent is +0.0012 (not enough). Threshold η* ≈ 0.0885 at b = 0
(bisection, `eta_thresh2.py`), where n₀ = 0.463, so the sixth moment is needed at length 3n₀ ≈ 1.39.

## 2. PROVED (exponent level): saturation of the black-box calculus

**Proposition (flat configuration).** Fix a bin (δ,x), detector parameter t and witness split (r,m)
with witness count R_s (from the paper's LP). Let λ* be the level-set LP optimum over the
numerator-only lines {R_s, 1−4λ−2qz_N, 1−2λ−2qz₁, λ ≤ δn₀/2} and c* the corresponding count. Consider:
#C = U^{c*}; on C, |S_{χ̄_u}(n)| = U^{λ*} for n near n₀ and its reflection 1−n₀, all other plain
polynomials of χ̄_u-presentations of size U^{o(1)}, witnesses exactly at their lower bounds
U^{δr/2}, U^{δm/2}, slots at amplitude q; off C everything of size U^{o(1)}. Then every inequality
(P1), (P2) (any lengths/capacities), (PW), (B), (M), (LS) holds at the exponent level, and
Σ_C|S_u| = U^{c*+λ*}. Hence **no Hölder/level-set combination of these inputs gives a numerator
exponent below σ_wt** (the LP value; = σ_wt at the binding classes).

*Proof.* (P1): c*+4λ* ≤ 1 and c*+2λ* ≤ 1 by the LP. (P2) with (n₀,n₀) and (n₀,0): the LP lines.
(PW): with c* ≤ R_s ≤ 1 − 2δm − 2qz_p (plain witness count, z_p = (1−2m)/(6κ)) and
λ* ≤ (1−c*)/4 − qz_N/2,
c* + 2λ* + δm + 2q z_j ≤ (1+c*)/2 − qz_N + δm + 2qz_j ≤ 1 + q(2z_j − z_N − z_p),
and 2z_j − z_N − z_p = (2(1−n₀−m) − (1−2n₀) − (1−2m))/(6κ) = 0. So (PW) holds with equality in the
plain branch: the joint numerator–plain-witness moment is **exactly log-linear**. (B): λ* ≤ δn₀/2.
(M): involves only Möbius witnesses/slots, whose count is the witness LP. (LS): slack ≥ 0.20
(numerator×inverse witness) and ≥ 0.27 (S³) on the whole grid. `barrier_check.py` verifies all
slacks ≥ 0 for δ ∈ [0.3,0.6], x ∈ [0,0.5], t ∈ [1,1.5], b ∈ {0,0.05,0.075,0.1}. ∎

This is a statement about what the listed lemmas can imply, not about L-functions (the configuration is
presumably false: a sixth moment would exclude it).

## 3. Routes checked (all fail; reasons)

### (i) Sixth moment / third factor
- **Pointwise × 4th**: Σ_C|S|⁶ ≤ U^{δn₀}Σ|S|⁴ gives σ = ((1−R)+δn₀)/6, a convex combination of σ_pw and
  σ_wt (log-linear). Over all rows, β* ≤ 7/8 gives Σ|S|⁶ ≤ U^{1+0.75 n₀} (η ≈ 0.34).
- **Large sieve on S³** (length K = U^{3n₀}): BGL gives (UK)^{2/3} ⇒ η = (2+6n₀)/3 − 1 ≈ 0.55–0.59.
  Even an optimal large sieve U + K gives η = 3n₀ − 1 ≈ 0.34–0.39, and for general-coefficient
  families this is the truth (full character group mod q: Σ_χ|Σ_{n≤N}χ(n)|⁶ ≍ qN³ + N⁶ from
  n₁n₂n₃ ≡ m₁m₂m₃ mod q). So any proof must use the d₃-structure of S³: it is a sixth-moment
  Lindelöf-on-average statement. We know of no unconditional bound of this strength even for the
  quadratic family (Heath-Brown's quadratic large sieve gives only X + N³ for S(N)³).
- **lem:plain on S·T, T a piece of S²**: lem:plain is a two-plain-factor statement. Its length-freeness
  comes from reflecting each factor longer than M/2; here every factor (n₀ < 1/2) is already short, so
  three factors at total 3n₀ > 1 cannot be shortened. Rewriting S² by the hyperbola identity,
  S(N)² = 2Σ_{a≤N}χ(a)S(N²/a) − D(N²), produces the divisor sum D of length 2n₀ whose functional equation
  (L², conductor U²) reflects it to 2 − 2n₀ ≈ 1.1 (longer); the first term is again three plain factors.
- **Copies** (S_{ua⁶} ≈ S_u up to Euler factors at a): the 4th moment at row range UP⁶ gives capacity
  +p/κ at cost +5p (P = U^p): net 5p − 2qp/κ > 0, a loss. For the 6th moment amplification is worse.

### (ii) Duality / Parseval with the probe variable
The Poisson/Parseval link Σ_u|S_u|² ↔ ‖A‖² is the (degree-2) large sieve: it gives σ = (1−R)/2.
**Halász–Montgomery over C** (duality): (Σ_C|S|^k)² ≤ |C|N^k + N^{k/2}Σ_{u≠v∈C}|S^{(k)}_{χ̄_uχ_v}|, pair
characters of conductor ≍ U². Even assuming Lindelöf for *every* pair character, this proves only
σ ≤ max(n₀/2 − R/(2k), n₀/4) ≥ n₀/4 ≈ 0.112 > σ_wt ≈ 0.095. Halász cannot help because N = U^{0.45}
is shorter than |C| = U^{0.6} (PROVED, exponent level). Unconditionally the pair sums have no
cancellation at all (length U^{0.45} ≪ conductor^{1/2}).

### (iii) Zero-free region of the numerator's own L-function (same family 𝒳_u)
- It yields exactly (B), σ_pw = δn₀/2. The reflected-side convexity interpolation between the 4th moment
  on Re = 1/2 and (B) at Re = a has slope n₀ − (1−R)/(2δ) ≈ +0.004 (δ=.47,R=.59) to +0.016
  (δ=.45,R=.611,b=.05) at the binding bins (affine: an
  endpoint is optimal; consistent with NOTES §5 of the session).
- HEURISTIC: |S_u(n₀)| ≈ U^{0.1} needs arg L(1/2+it, χ̄_u) to rotate at rate ≈ −n₀ log U on a unit
  t-window, i.e. |ΔS(t)| ≈ c·log U (a zero-count deficit or ~c log U displaced zeros). Backlund/Jensen
  with zero-free-right-of-a and convexity bounds S(t) by ≈ (δ/2)·log U/ log 2, which permits this. Zero
  density counts zeros with multiplicity, so "c log U zeros per row" saves only a logarithm, never a power.

### (iv) Non-log-linear combinations with witnesses, marks, copies
- Numerator × plain witness: exactly log-linear (§2).
- Numerator × inverse witness (a *new* mixed plain–Möbius moment, not in P78): even granted diagonal up to
  n₀ + r ≤ 1, zero gain — at every t the adversary's split has r ∈ [0.62, 0.74] (t ≤ 1.15) or a long
  witness, both beyond 1 − n₀. With sixth-power amplification past the diagonal, e(n₀+r) = (1+5(n₀+r))/6
  gives σ ≥ 0.1 (worse).
- Marks: generalized Hölder between the 2nd moment of S·Q (capacity z₁ ≈ 0.12) and the 4th moment of S²Q
  (z_N ≈ 0.02) is linear in the interpolation parameter; the 4th-moment endpoint is optimal. Putting all
  marks inside any moment needs 2n₀ + 6κL ≤ 1 (impossible) or, with DIAG, a second moment of S·Q
  (σ = (1−R)/2 − qL ≈ 0.15, worse).
- Unit associates εu: 𝒳_{εu} = 𝒳_u, numerators differ by Θ-twists; orthogonality over the six units
  gives only a constant.
- Slot resonance at x = 1/2 (slots share the numerator character): saturation |Q_i| ≈ P^{δ/2} only means
  a bias P^{a−1} ≈ P^{−0.27} of χ̄_u on slot primes; it imposes no structure on S(n₀).

### (v) Y-weights and contours
- The reflected length C_u/Y is profile-independent; multi-scale weights are sums of single-scale probes;
  shorter reflected lengths only for non-squarefree u (rad(u) < q_u), never the worst case.
- Sieve/mollifier s-weights (1*λ)(s): the numerator becomes L(w,χ_u)M_u(w), which reflects only L, so the
  reflected sums get longer (n₀ + log_U d); it also changes the low side.
- Per-amplitude-class min(E_old, E_wt) is legitimate (the w-integral of each row is contour-independent,
  amplitude classes are formed at fixed (s,z)); it is already what the model uses.

### (vi) Denominator moment (new check)
The s-integral is a Möbius sum of length Z = U^{1/d} (log-normal profile from e^{(s+z−1)²}). Replacing
the pointwise bin bound Z^a by lem:inverse-amplification (e(r) = (1+5r)/6) in a joint level-set LP with
the numerator gives zero gain at the binding bins because e(1/d) − δ/d ≥ R there (`denmodel.py`).
(It would also require integrating s fully through the coupled correction H(s,w,z).)

## 4. Precise obstruction

Within the inputs of P78 + the 08 lemma, the bad-row numerator obeys exactly
  Σ_{u∈C}|S_u| ≤ U^{R + min((1−R)/4 − q z_N/2, δn₀/2) + o(1)},
and this is sharp for the calculus (§2). To reach 41/47 one needs, on top of it, a mean value of degree
≥ 6 in plain polynomials of one sextic row character at total length ≈ 1.39 with excess ≤ 0.088
(e.g. Σ_u|S(n₀)|⁶ ≪ U^{1.088}, or Σ_u|S(n₀)²S_m|² ≪ U^{1.086}); the best provable excess is ≈ 0.34
(pointwise β*) and the general-coefficient truth is 3n₀ − 1 ≈ 0.39, so the d₃ (sixth-moment Lindelöf)
structure must be exploited — beyond the large-sieve edge, like the other open inputs in BRIEF.
