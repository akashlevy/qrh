# 09 — New directions below 7/8 (session of 10 Oct 2026)

Goal: a *meaningful* improvement of the QRH exponent (target ≤ 27/31 ≈ 0.870968, ideally 6/7),
not a re-tuning of the 7/8 lemma set. Working log; newest entries at the bottom.

## 0. State of the art at session start (scan 10 Oct 2026, ~04:50 UTC)

- Tracker (GTimG/QuasiRiemannTracker `main`): best framework-verified entry is Nielstron's
  874957019420098946128604623/10^27 (the κ-feedback cubic root θ* ≈ 0.874957019420).
- Open tracker PRs: #2 (Nielstron), #3 (Liu), #4 (Argonaut), #5 (Nielstron, catalogue
  housekeeping), #6 (ours, 10499/12000 ≈ 0.8749167), #7 (Gao, "numerical slack retuning").
- No public result is meaningfully below 7/8. Everything clusters in [0.87491, 0.875].
- arxiv.org is blocked by this environment's network policy (also doi.org, Springer, AMS,
  Semantic Scholar); GitHub and web search work.

## 1. Model of the architecture (my own re-derivation)

Probe I = Q^{-1/2} Σ_m A_m B_m with Q = XY:
- A_m = Y^{-3/2} Σ_{s~Y} c(s) g(χ_s, −m) (sextic Gauss sums, additive polynomial);
- B_m = completed cubic-theta row of length Z twisted by χ_A(m) (marked by slot primes).

Low side: Cauchy–Schwarz, ‖A‖² ≪ (Q/Y)(1 + P_a^{1/6} + P_a²/Y), ‖B‖² ≪ Q (energy, diagonal
below the wall M+ℓ ≤ 1). Bound X^{1/2} P_a^{1/12}.

High side: Poisson in m gives rows u (sixth-power-free frequency parts):
I = ∫∫∫ X^{1/2−z} Z^{x+z−1} Y^{w−1} Σ_u q_u^{−z} ζ(6z) L(w,χ_u)/L(x,ηχ̄_u) · (marks) · H.

**Derivation check.** Row u with rightmost zero at a = (1+δ)/2, contour x = a, w = 1−a
(reflected numerator, cost U^{δ/2}), marks |Σ_p χ̄_p(u) p^{z−1}| = P^{z−1/2+xδ}, count U^R,
U up to ZP/X. The z-dependence cancels and one gets exactly the paper's
eq. `common-high-exponent`:
(1+δ)(1−l_y)/2 − β − ℓ/2 + xδℓ + h(R + δ/2 − 1/6) < 0, h = 1+ℓ−l_x.
The signal is (number of copies a⁶ ≤ Z/X)^{} × (copy size) × (slot amplification P^{1/6}).

## 2. Lever sensitivities (r4 model, all other lemmas as in the paper)

| change | θ |
|---|---|
| none | 0.874957 |
| energy wall relaxed by w (W(w)): w = 0.05 / 0.1 / 0.2 | 0.87183 / **0.86870** / 0.86333 |
| Gram exceptional exponent ×0.5 / ×0.75 / ×0.9 | **0.86949** / 0.87232 / 0.87394 |
| Gram exceptional term removed | 0.85980 |
| amplification slope kA = 0.8 (needs M_κ, κ = 0.8) | 0.86979 |
| perfect row counts (DH for rows) | 0.86979 |
| plain capacity c_P = 1/3 | 0.87471 |

So 27/31 needs a *low-side* input (W(w) with w ≈ 0.07, or a ~35% cut of the Gram exponent),
or a below-diagonal Möbius mean square. Each of these is a "beyond the edge" mean value
(see docs/OBSTRUCTIONS.md).

## 3. Dead end: rescaled marks (κ-marks)

Idea: mark with Z → Z p^κ (paper: κ = 1) and compensate with p^{−3/2} I(X/p, Y/p, Z p^{κ−1}).
For κ = 2 the reflected dual length Q²P²/(ZP²) = Q²/Z no longer grows with the marks,
so the energy wall M+ℓ ≤ 1 seems to disappear, and the signal amplification becomes
P^{1+(κ−1)β−5κ/6} > P^{1/6}.

Result (`kappa_marks_model.py`): θ = 0.8749568 for κ = 0.6, 0.8, 0.9, 1.0 (identical to 1e−9).
Under GRH for rows every κ gives exactly 13/15 (closed form: β = (8/9 + ℓ(8κ/9 − 1))/(1+(κ−1)ℓ)
at the binding ℓ = 1/(6−κ)).

Reason: all slot primes satisfy q_p ≍ P_i, so Z p^κ = (Z P^{κ−1})·p·O(1). The κ-probe is
the paper's probe at base scale Z̃ = Z P^{κ−1}: a reparametrization of the geometry. The
energy wall reappears as M ≤ 1 in the old scale. **Lesson:** any change that only rescales
lengths by powers of the slot length is exponent-equivalent to the paper.

## 4. Anatomy of the 13/15 GRH floor (closed form)

With normalised sizes (low side X^{1/2}; signal S = X^{1/3}Z^{β−5/6}P^{1/6} = copies × copy size;
GRH rows R = F^{1/2}Y^{−1/2} where F = ZP/X is the number of Poisson frequencies, i.e. rows
bounded WITHOUT cancellation lose exactly F^{1/2} against the "random" size Y^{−1/2}):
- high: β > 4/3 − 5l_x/6 − l_y/2 + ℓ/3,
- low: β > 5/6 + l_x/6 − ℓ/6 + b⁺/12,
- wall: l_x + l_y + ℓ ≤ 1.
In the regime b ≥ 0 the low side becomes 11/12 − ℓ/4 (independent of l_x: the Gram excess
exactly compensates), the high side at b = 0 becomes 2/3 + ℓ, hence ℓ = 1/5, θ = 13/15.
The "true" (random-model) size of the probe is Y^{−1/2}: C–S loses Q^{1/2}, the row-by-row
high side loses F^{1/2}. Marks multiply rows by P^{1/2} (more frequencies, P^{−1/2} each) and
the signal by P^{1/6}; they pay off only because the low side does not grow with P below the wall.

## 5. Dead end: twisting the m-sum by χ_p(m) instead of marking A

Inserting χ_p(m) with weights a(p) = conj(γ₁(p)) gives a coherent principal factor (Σ_p 1)
and the row factor Σ_p χ̄_p(u); but after the Poisson normalisation p^{−1/2} the signal and
the low side both gain P^{1/2}: no net gain, rows suppressed only relatively, and the
frequency range still grows by P. Strictly worse than the paper's marks.

## 6. Poisson duality of the marked energy (no free lunch)

Unreflected marked energy Σ_{m≤Q}|B_m|² (theta length ZL, marks L): Poisson in m plus
μ-absorption turns the off-diagonal into (Q/(ZL)²)·[Σ_{H≤(ZL)²/Q} |Π_H|²|M_H(Z)|² − R Z L],
i.e. the MARKED MÖBIUS MEAN SQUARE (Lemma `marked`'s object) minus its diagonal. Lemma
`marked` applies there under M ≤ 1 (not M+ℓ ≤ 1), but it only gives an upper bound with a
Z^ε loss, which returns the column term ZL. A diagonal energy beyond the wall is EQUIVALENT
to a power-saving ASYMPTOTIC for the marked Möbius mean square with relative error Q/(ZL).
Same move-calculus obstruction.

## 7. Row-count anatomy at the binding point (`count_anatomy.py`)

δ = 0.39: R ≈ 0.665 (x = 1/2) vs DH 0.61. The worst witness is either the pure long inverse
witness r = t ≈ 1.12 > 1 (counted by sixth-power amplification e(r) = (1+5r)/6 − δr), or a
balanced pair (r ≈ 0.7, m ≈ 0.4). The detector forces m ≤ 1/2 (plain polynomials longer
than U^{1/2} are controlled by the functional equation), so t = 3/2 would give DH counts
if the adversary could not use r > 1. Model: marked inverse moment valid to r < 1.3 (w1=0.3)
→ 0.86979; amplification slope kA = 0.9 → 0.87116, kA = 0.85 → 0.86979.

## 8. Infrastructure

- Lean 4.34.1 installed manually from the GitHub release into `~/.elan/toolchains/leanprover--lean4---v4.34.1`
  (elan's default mirror release.lean-lang.org is blocked). The Mathlib cache hosts
  (lakecache.blob.core.windows.net, mathlib4.lean-lang.org) are blocked: a final build must
  compile Mathlib `d13f23b` + PNT+ + the OAI closure from source (4 cores, 15 GB RAM, ~28 GB free).
- Five research agents launched 05:45 UTC on: energy wall, long inverse witnesses, theta
  pairing, new architectures, literature + combination table (`agents/<name>/`).

## 9. More quick dead ends (coordinator, 06:00–06:20 UTC)

- **Marks on the s-variable** (p | s with Y → Yp): signal gains only a log (Σ p^{-1}),
  the frequency range is unchanged, rows are suppressed like a longer Y: equivalent to the
  paper's probe with Y' = YP and special s-weights; the signal/Gram-excess ratio is invariant.
- **Marks on the cube part** (p | n in A = cn³): the principal local factor restricted to
  t = v_p(A) ≥ 3 is O(Q^{4−6x−6z}) at u = 1: no amplification.
- **Averaging over Mellin height t** (targets η N^{-it}, |t| ≤ T₀): the hybrid large sieve makes
  the marked energy diagonal per unit t once T₀ ≥ ZL/Q, but the signal lives on an O(1)
  t-interval, so one loses T₀^{1/2} = (ZL/Q)^{1/2}: exactly the column-term loss again.
- **Single-mark energy**: the wall M+ℓ ≤ 1 already holds for ONE marked prime P (the theta row
  restricted to P | A has reflection conductor (QP)² against length ZP). Summing over P
  cannot remove it; a beyond-wall energy bound is a statement about a single level-P² theta.
- **Lemma `marked` slots must be primes in product form** (P78 l.9309–9314: "no additional
  residual coefficient b(n)"); slots cost 2z of row length, so a plain witness used as a slot
  at t = 1 leaves no room (r + 2m ≤ 1 with r + m = 1 forces m = 0).

## 10. Literature check: unconditional beyond-the-edge moments for quadratic twists

X. Li (Invent. Math. 2024, arXiv 2208.07343) proved the second moment of quadratic twists of a
modular L-function unconditionally (Soundararajan–Young needed GRH); Shen–Stucky (arXiv
2402.01497) the fourth moment of quadratic Dirichlet L-functions (affine D4 type) with four
main terms. Key input: "optimal" large-sieve-type bounds for quadratic characters. BUT their
error terms are only logarithmic (Li: O(X(log X)^{1/2+ε})). Our levers need POWER savings
beyond the edge (W(w), M_κ, (V′)). This supports the MDS heuristic in R1: affine-type
moments admit asymptotics with log savings, not power savings, with current technology.
Also the dual of the Möbius family after μ-absorption is the cubic theta (twisted
multiplicative, no Euler product), so the SY/Li contour-shifting of dual L-functions does
not apply directly.

## 11. Binding bin and two quick checks (coordinator, 06:30–06:50 UTC)

- **Binding constraint at the optimum** (`coord/binding.py`): a single narrow bin
  δ ∈ [0.37, 0.41] at mark amplitude x = 1/2 (R ≈ 0.664 vs DH 0.61). Everything else has slack.
- **Mark amplitude is not a lever** (`coord/xcap.py`): capping x at 0.45/0.4/0.3/0 gives
  0.87485/0.87470/0.87426/0.87255. Even x ≡ 0 misses 27/31.
- **L²-in-Z continuation (dead end).** The continuation criterion (P78 l.400) only needs the
  Mellin integral of f to converge, so mean-square bounds over Z in dyadic (or [Z0, Z0^{1+ε}])
  ranges suffice. On the high side this turns Σ_rows|row| into an L² sum over the zeros' heights
  γ (Montgomery–Vaughan / Gallagher). The gain is only (T0·L)^{1/2}, with T0 the Mellin height range
  of the kernel. With T0 = Z^τ the hybrid (u, t) zero count grows like T0^{R}, R ≈ 0.66 > 1/2
  (even DH has 1−δ = 0.61 > 1/2), and a rough kernel also moves the energy wall by 2τ. Net loss.
  Square-root cancellation across rows would need counts with R < 1/2 at the binding bin.

## 12. Agent report: long-inverse (finished 06:55 UTC) — NEGATIVE, obstruction sharpened

Files: `agents/long-inverse/{notes.md, model_li.py, robust_check.py, run_*.json}`.
- **Closed form of the binding count.** Suppose objects of length s > 1 are counted with slope α, e(s) = 1 + α(s−1).
  Then P78's count at x = 0 is exactly R(δ,0) = 1 − 2δα/(3α−δ), at t* = α/(α−δ/3).
  - α = 1 is **Ingham's 3(1−σ)/(2−σ)**. For ζ this is still the record below σ = 7/10. The binding bin sits at σ ≈ 0.695.
  - The sixth-power copies give α = 5/6, which is 0.6919 at δ = 0.39.
  - DH needs α ≤ δ.
- **Detector redesign is impossible** (proved): any mollifier detector vanishing on (1, X] has a pure Möbius block on (X, 2X],
  so the pure long witness (t, 0) is always available to the adversary when t > 1. Heath-Brown identities, products of
  detectors and L^k mollifiers only give the adversary more freedom.
- **What would suffice.** The marked inverse moment must be diagonal up to r1 ≈ 1.14 (r1 = 1.15 → 0.87058). That is M_κ with κ ≈ 0.88.
  Lowering the slope alone (slope 0.5 → 0.87223) never reaches 27/31.
- Other failures: prime slots inside the amplification lose (5/3 − 2q)z; the a priori β* gives a mean square of U^{1+0.75r}, worse
  than the paper; Turán / Halász–Montgomery need δ > 1/2 even with near-orthogonality; the joint M_r·S_m-as-slot idea
  leaves θ at 0.874957 because the binding point is the pure long inverse witness, with m = 0.
- **Verdict:** improving the counts is "Ingham below 7/10" for the sextic family. That is a classical open problem
  (Guth–Maynard improves only σ ≥ 7/10).

## 13. Agent report: new-architectures (finished ~06:27 UTC) — NEGATIVE

Files: `agents/new-architectures/{NOTES.md, template_lp.py, scenarios.py, scenarios_output.txt, device.py, aside_marks.py}`.
- **13/15 is pinned twice:** by wall ∩ case A on the low side and by wall ∩ rows.
  - θ_GRH = 1 − c(7+W)/10 for W ≥ 1, and 1 − 2c(1+W)/5 for W ≤ 1, with c = 1/6 the copy exponent.
  - Better row suppression and removing case A both still give 13/15.
- **Device criterion:** a new structure that touches neither the wall nor case A lowers the GRH floor iff its
  row cost per unit length is < 4× its low gain (ratio 3 → 41/48, ratio 2 → 5/6). P78 marks have ratio 2 but consume the wall.
- **Higher-rank thetas:** μ-absorption is a Gauss-sum fact (Jacobi J ≡ −1). It collapses only for 6→3, 4→4, 3→6, 10↔5.
  Kazhdan–Patterson thetas (GL_3^{(3)}, GL_5^{(6)}, …) realize the same duals at degree r ≥ 3, so the wall becomes
  (r−1)(M+ℓ) ≤ 1. Optimistic GRH values are 11/12, 9/10, 17/18 and ≥ 1. P2-shape: θ = 1/2 + (1−1/N)(1−1/r).
- **Trilinear probes:** A·A′·B equals P78 with Y ↦ YY′ (twisted multiplicativity of Gauss sums); B·B̄′ gives 11/12;
  A·B·B′ has device ratio ≈ 11–14.
- **Denser copies (c = 1/N):** these would be very valuable (c = 1/3 → 11/15 under GRH), but they need a degree-2 object with order-4/6
  Gauss-sum coefficients, which Kazhdan–Patterson non-uniqueness rules out.
- **Mark placement:** the wall QP ≤ Z appears wherever the amplifier sits on the theta side. On the Gauss side the
  Gram excess (YP²/X)^{1/6} cancels P^{1/6} exactly (8/9). χ_p(m) with γ₆ weights → BGL wall → 8/9. A sup-norm on
  B → Z ≥ Q² → 14/15.

## 14. The energy wall needs a LOSSLESS range, not a power saving (coordinator, `coord/releta.py`)

Replace the excess (M'+l'−1) beyond the wall by (1−η)(M'+l'−1), i.e. a relative power saving η on the
column term:

| η | 0.1–0.8 | 0.9 | 1.0 (no wall) |
|---|---|---|---|
| θ | 0.874957 | 0.870689 | 0.863636 |

The optimum sits exactly on the wall (M+ℓ = 1.000). Every step beyond it costs (1−η)/2 on the low side and
gains less than that, unless η ≳ 0.85. So a cubic-first-moment-style "tiny ineffective saving beyond the sharp
sieve" is useless here. W(w) needs the energy to be **diagonal up to Z^ε** on a range of width w ≈ 0.07
beyond the wall (Lindelöf-on-average strength for a thin family of twisted theta rows).

## 15. Agent report: theta-pairing (finished ~06:50 UTC) — NEGATIVE, obstruction quantified

Files: `agents/theta-pairing/{THETA_PAIRING.md, model_ext.py (knobs gram_s, gram_sQ, gram_ls), model_summary.txt, *.log, kp_forced_zero_density.py}`.
- **Every LP optimum has M + ℓ = 1.0000.** This holds for the paper, Gram halved, Gram removed, and uniform savings up to 0.0127. The marked reflected dual length is
  L = Q·Z^{M+ℓ−1}, so at the binding geometry **L = Q**: the pairing T is a square Q×Q quadratic-symbol bilinear form exactly
  at the large-sieve edge.
- **Uniform saving s on T at the wall (Z-units):** s = 0.002/0.004/0.005 → 0.87333/0.87171/0.87090. The **threshold for 27/31 is s ≈ 0.0049**.
  A mean-square saving that vanishes on the wall (a fraction η of the random rate) needs η > 0.15 to help at all and η ≈ 0.277 for 27/31.
- **Obstructions:**
  - The Gram excess is a fixed finite-rank positive form in c, while the signal is linear in c. No s-weighting gains a power.
  - Inequalities that see only |π|, |β| and the kernel are capped at Q^{11/12}.
  - A restricted operator norm stays ≥ L = Q for every row set, so sparsity arguments are dead. Kazhdan–Patterson forced zeros cover ≤ 5.95% of rows.
  - Petersson projection onto θ₆ is Bessel = Cauchy–Schwarz. The z = 5/6 residue *is* the signal.
- **Precise obstruction:** a saving must come from the joint phases of π = ρ·conj(γ₂)⁴ and Patterson's column coefficients.
  That is a first moment of the A₂/4-node series. Evaluating it exposes the Möbius signal (circular), and its functional-equation images do no better than Cauchy–Schwarz.

## 16. Functional-equation group of the (unmarked) energy MDS is AFFINE (coordinator)

Z(s₁,s₂,w) = Σ_m m^{−w} Σ_{A₁,A₂} τ(A₁)τ̄(A₂)χ_{A₁}(m)χ̄_{A₂}(m)A₁^{−s₁}A₂^{−s₂} has three functional equations.
Write x_i = s_i − 1/2 and y = w − 1/2.
- σ_i: theta twisted by χ_·(m), conductor N(m)². It maps x_i ↦ −x_i, y ↦ y + 2x_i.
- σ₃: L(w, χ_{A₁Ā₂}), conductor A₁A₂. It maps x_i ↦ x_i + y, y ↦ −y.

Roots α₁ = (1,0,−1), α₂ = (0,1,−1), α₃ = (−½,−½,1); coroots 2x₁, 2x₂, 2y.
Cartan matrix [[2,0,−1],[0,2,−1],[−2,−2,2]], **det = 0: affine (C̃₂-like)**.

Lindelöf-on-average beyond the edge (W(w)) is therefore equivalent to a power-saving asymptotic in an affine MDS
(Poisson in the family plus subtraction of the dual diagonal: §6 and the dual-diagonal argument). This is the same class
as power-saving second moments of quadratic twists of GL2 forms at the edge. Li (2024) and Shen–Stucky reach that class
only with log savings, and power savings there are open. BFKMM-type power savings at the edge exist for the family of
ALL characters mod q (exact orthogonality and shifted convolutions), not for quadratic families.

## 17. Agent report: lit-combos (finished ~06:55 UTC) — SOTA unchanged, combination map

Files: `agents/lit-combos/{LIT.md, COMBOS.md, combo.py, runs/}`.
- **SOTA:** tracker PRs #1–#7 only; best 0.874957019420; openai/math HEAD fd4aeeb (Oct 7). Gao now claims 0.8749571516.
  Nothing below 0.8749 anywhere.
- **Literature:** no beyond-the-edge input for any lever (full table in LIT.md). DDHL (1.11) is the only published
  "fewer rows than conductor" mean value, at large-sieve strength 3/2.
- **Levers trade off linearly.** As fractions of their single-knob thresholds, the requirements in a working combination sum to about 1, so combining gives no bonus.
  - Single knobs for 27/31: W(0.0638), κ_A = 0.895, s₀ ≈ 0.0065 on the Gram pairing, w1 = 0.2215.
  - For 6/7: full (V′) + W(0.013), or W(0.285).
- **DIAG package:** Lemma marked at r + z ≤ m with no 2r + 8z cap (mk = 1), and Lemma plain at n₁ + n₂ + z ≤ M (cP = 1).
  - Alone: **0.87255**. With DIAG, 27/31 needs only W(0.0196), and 20/23 needs W(0.037).
  - Partial versions: mk ≥ 0.75 alone gives 0.87349; cP = 1 alone gives 0.87345.
  - Caveat (coordinator): the factor 2 in z comes from the same mechanism as the wall: a mark p costs conductor p² for length p in the
    reflected recursion (P78 invariant 3M + 6z₀ ≤ 4F₀). In the primal picture DIAG is "within the diagonal"
    (product length ≤ family), but there the sextic large sieve's theta term U^{1/3}N^{5/6} exceeds U for N > U^{4/5}.
    So DIAG needs a coefficient-sensitive sextic sieve for Möbius × prime-slot coefficients, which is unproved.

## 18. Agent report: energy-wall (finished ~07:00 UTC) — NEGATIVE, sharp obstruction plus numerics

Files: `agents/energy-wall/{ENERGY_WALL.md, dual_energy.py (+json), exchange_rate.py (+out), model_wall.py, mw_*.json}`.
- θ(wE) ≈ 0.874957 − wE/16 on [0, 0.1]. 27/31 needs wE ≥ 0.0638, at (l_x, l_y, ℓ) = (0.366, 0.494, 0.204).
- **Only an offset helps.** If the excess grows at rate f past the wall, there is no gain unless f ≲ 0.15, and 27/31 needs f ≲ 0.103.
  This agrees with §14. Every transfer mechanism found gives f ≥ 1/2.
- **Obstructions:**
  - O1 (proved): phase-blind arguments cannot beat the column term. This kills Hölder with slot moments and dual.tex's moments-of-different-lengths applied to X.
  - O2: exchange-rate bookkeeping. Every Cauchy–Schwarz transfer of a variable into the rows leaves excess ≥ w when ℓ ≥ w.
  - O3: the "rough" part (n prime, or a large prime times something small) needs Type-I twisted theta sums at conductor ≫ length.
  - O4: β* does not help.
  - O5: every mark redesign fails.
- **Numerics:** per-row energy divided by its diagonal stays between 0.7 and 1.17 up to N/K ≈ 5700, with Gaussian-like maxima. W is numerically true.
- **The missing estimate:** Σ_{k≤K}|Σ_{n~N} c(n)ψ_k(n)|² ≪ (K + N^{1−η})Z^ε with η ≥ 0.068, for c = β·X.
  It needs a spectral input (amplification plus Kuznetsov on the cubic cover, uniform, upgraded to an average).
  The round-2 agent `wall-spectral` is pursuing this.

## 19. Agent report: diag (round 2, finished ~07:35 UTC) — NEGATIVE

Files: `agents/diag/{DIAG_NOTES.md, model.py (+wz), verify.jsonl, sweep_*.jsonl, diag_anatomy.py, xcap0.py, theta_slot_test.py, lp_bhm.py}`.
- **Re-verified** (refined case A, δ0 = 10⁻³): DIAG 0.8725498; mk = 0.75 0.8734852; cP = 1 0.8734526; DIAG + wE 0.0196 0.8709686.
  Partial versions: mk 0.6/0.7 give 0.87399/0.87349; cP 0.3/0.4/0.5 give 0.87479/0.87452/0.87418.
- **Ceiling for all slot levers: 0.8725498.** At the DIAG optimum every binding bin has mark amplitude x = 0
  (DIAG, mk = cP = 5 and x ≡ 0 give bit-identical θ).
- **Where the losses are.**
  - Lemma marked: the second z is paid once, at the final marked theta reflection. A marked prime enters the conductor squared (l.7374–7383, old-eq:4.4),
    and the quadratic sieve needs rows·marks ≤ columns (old-eq:4.10, l.9910), i.e. r + 2z ≤ m. The 2r+8z cap is the cubic (NL)^{2/3} branch.
    This is the energy-wall mechanism again.
  - Lemma plain: no conductor effect. The slots are bounded pointwise through β*, |Q|² ≤ Z^{κz} with κ = 2β*−1 (l.12876), and the slot removal
    exchange rate is 6:1 (l.14236–14276). So cP = 1/(6(2β*−1)); cP = 1 would need β* = 7/12 (circular).
- **Rigorous negative result.** DIAG stated uniformly over all bounded slot coefficients is FALSE for z > m/2. Take a(p) = conj(sextic Gauss sum);
  then Σ_u|Q_u|² ≫ U^{5/6}P^{1/3}, attaining the sieve's theta term. The paper's own hypothesis r + 2z ≤ m keeps z < m/2, so the paper is unaffected.
  A DIAG proof would have to use the Hecke structure of the actual slots.
- de Faveri's sieve gives Lemma marked for r + z ≤ m/2 (zero payoff).
- OpenAI's moments-of-different-lengths technique cannot work for n = 6. Theta profiles satisfy every power and BHM constraint;
  near-orthogonality excludes them only when 1/n > 1/4, i.e. only for n = 3.
- The β* bootstrap costs κ = 3/4 per unit of slot length against a gain ≤ 0.39.

## 20. Why every lever is thin-family (coordinator synthesis, 07:40 UTC)

- **Marks are intrinsically thin.** By twisted multiplicativity, τ(pA) = τ(p)τ(A)(p/A)₃^{±}, so a first-power mark is a twist of the
  cubic theta by a cubic character of conductor p. That multiplies the reflection conductor by p² for length gain p. The P^{1/6}
  amplification and the conductor cost come from the same cubic symbol.
  - Cube marks p³ (Kazhdan–Patterson periodicity, no symbol) cost nothing but also amplify nothing (§9).
  - Coherent amplification needs the slot sum inside the square, hence a family of fixed size Q with conductor (QP)².
  - Summing the slots outside the square loses the coherence factor P (checked).
- **One character per modulus.** The probe family {(·/m)_6 : m ~ Q} and its dual {(·/k)₂ : k ≤ K} have one character per modulus,
  as in quadratic-twist families.
  - For full families of all characters mod q, exact orthogonality turns beyond-the-edge energies into shifted convolutions
    Σ τ(A)τ̄(A+hq). These are Fourier coefficients of |θ|², which is an ORDINARY (non-metaplectic) automorphic function on a Bianchi
    congruence group, since the Kubota multipliers cancel. Ordinary spectral theory would then give power savings (BFKMM mechanism).
  - The μ-absorption, however, needs the varying-modulus Gauss sums g(χ_A, ·), so a fixed-modulus family cannot carry the signal.
- **The Poisson dual of the energy** (via g(χ)g(χ²) = J(χ,χ²)g(χ³), with J multiplicative on squarefree A) is the second moment of
  sextic twists of 1/L(s, ξ) for a Grössencharacter ξ, in the diagonal regime. W(w) needs its asymptotic with relative error
  Z^{−w}: an "asymptotic large sieve" (cf. Conrey–Iwaniec–Soundararajan) for a thin (one-character-per-modulus) family.

## 21. Agent report: wall-spectral (round 2, finished ~07:45 UTC) — NEGATIVE, conditional theorem

Files: `agents/wall-spectral/{WALL_SPECTRAL.md, identity_check.py, fe_group.py, fe_group_compare.py, mds_reach.py, opnorm.py, model_ws.py}`.
- **Identity (proved, 540/540 checks):** β(n)ψ_k(n)χ_P(n)^{−2} = μ(n)·G(n)·𝔯(n,k)·conj γ₁(n; k³P⁴). The Grössencharacter cancels but μ remains.
  So E is a second moment of Möbius-twisted sextic Kubota sums on the thin index set h = k³P⁴. That is P78's own dual (§6).
- **Rankin–Selberg / MDS:** the group is infinite with invariant form (x₁+x₂+y)², affine C̃₂. It is the same as Soundararajan–Young / Li for GL2
  quadratic twists. Contrast: the David–de Faveri–Dunn–Stucky cubic-family second moment has a finite group (order 24/48).
  The group fixes I = x₁+x₂+y, so Bochner hulls never go below the initial inf I = 1/2. The cross-mark terms need I ≤ 0.263 (proved no-go).
- **Shifted convolutions (proved no-go).** A full character group containing every ψ_k (k ≤ K) has modulus ≥ lcm(k ≤ K) ≫ N, so
  congruences become equalities and no shifted sums appear. Twist-averaging moves the wall the wrong way.
- **Kuznetsov / amplification:** even Lindelöf for every individual twist gives E ≪ KLN, which is worse than the large sieve (KL/N = Z^{1−M} > 1).
- **Conditional theorem:**
  - (H1) Σ_{k≤K}|D_P(½+it,k)|² ≪ K^{1+ε}L^ε uniformly in the level L²;
  - (H2) decorrelation across distinct marks, Σ_k D_{P₁}·conj D_{P₂} ≪ K^{1+ε}/L.
  - Together they give E ≪ KN·Z^ε for all N ≤ KL, hence θ = 19/22 = 0.863636.
  - 27/31 needs (H2) with saving ≥ L^{−0.93}.
- **CORRECTION to earlier model values.** The Patterson bias of unweighted slots produces the (NL)^{2/3} branch of Lemma
  hybrid-energy, of excess (M+4ℓ−2)/3. It is harmless up to wE ≈ 0.064 (0.8709554 either way). If it is kept, any W(w)
  **saturates at 19/22 = 0.863636**, so energy-wall's large-wE values (0.86333, 0.85) were optimistic.
  Numerics (opnorm): W holds uniformly in the slot coefficients (top singular value/K′D tracks the random control).

## 22. The structural chain behind every obstruction (coordinator synthesis, 07:50 UTC)

1. Zero detection needs a Möbius-twisted sum (1/L, not L).
2. μ-absorption works only for order-6 twists, via μ·γ₋₁ ∝ γ₂ (6→3). So the theta row must be twisted by the SEXTIC symbol χ_A(m).
3. Theta reflection absorbs the cubic part of the twist into the theta index. The QUADRATIC part χ_k³ survives on the dual.
   In Kubota-index form (wall-spectral identity) the rows enter only as cubes h = k³P⁴.
4. A family indexed by cubes k³ is thin: one quadratic character per row. Its second moment has an affine functional-equation group
   (C̃₂, invariant (x₁+x₂+y)², the same as Soundararajan–Young / Li for GL2 quadratic twists).
   - Contrast: full-index families Σ_h |Kubota(h)|² are finite type (DdFDS: group of order 24/48), where power savings exist.
   - With cubic twists (j = 2, 4) the dual is a full-index family, but the probe would then detect L(s,η), not 1/L(s,η).
5. Coherent prime amplification (marks) adds conductor P² at fixed family size, which makes the family thinner still (wall M+ℓ ≤ 1, Lemma marked's 2z).
6. Hence every lever (W(w), DIAG/marked, (V′), M_κ/counts) is a power-saving moment at or beyond the edge of an affine family.
   This is open even for GL2 quadratic twists.

Any meaningful improvement must break this chain somewhere. The candidates:
- (a) a zero-detecting identity that does not go through μ-absorption of order 6;
- (b) a way to make the quadratic residual twist "full" (e.g. average it into a full family without losing the signal);
- (c) solve the affine edge problem (power-saving second moments of quadratic twists) for theta-type coefficients;
- (d) a high-side (count) breakthrough, i.e. "Ingham below 7/10" for the sextic family.

## 23. SOTA scan 07:57 UTC (coordinator)

- Tracker refs: main 40844a7 (unchanged), PRs #1–#7 only.
- openai/math HEAD: fd4aeeb (unchanged).
- Web: only Liu 2610.12234 (0.874957069799) and spin-offs, e.g. DrMurphyIsIn/qrh-debruijn-newman (Λ ≤ 9/32, irrelevant here).
- Nothing below 0.8749. Next scan no earlier than ~10:00 UTC.

## 24. Agent report: slack-audit (round 2, finished ~08:05 UTC) — no within-edge path to 27/31

Files: `agents/slack-audit/{AUDIT.md, model_sa.py, batch.py, diag_closed.py, frontier*.py, slackmap.py, runs/SUMMARY.txt}`.
- **Best type-(A) (within-edge) value is 0.8725498**, attained by DIAG. Closed form agrees with the LP to 1e-7.
  - At the DIAG optimum the low side is loss-free (b = 0, M+ℓ = 1, θ = 1 − h/6), at l_x = l_y = 0.4118, ℓ = 0.1765.
  - The only binding high bin is δ ≈ 0.33, x = 0, where the count is the slot-free Ingham-type R0 = 1 − 2δα/(3α−δ), α = 5/6.
  - All its inputs are at the diagonal, so beating it is "Ingham below 7/10" (type B).
- **Model discrepancy found.** r4 drops the paper's intermediate-row constraint (eq. adaptive-intermediate, l.16161).
  - It binds at the DIAG geometry: DIAG with the literal paper treatment gives 0.873815.
  - Fix: apply detector-counts with zero slots at t = α/(α−δ/3) for d ≤ 1/2. This restores 0.8725498.
  - No modelling error that understates the paper was found.
- **Zero-payoff (A) items:** refined case A, floor 1/2+η, plain sixth moment, optimal cubic sieve, dropping P_a²/Y, t < 1, dynamic κ, bin ceiling.
  **Already optimal:** detector t ∈ [1,3/2], m ≤ 1/2, zero bins, x ≤ 1/2, w = 1−a contour, H ≥ D^{1+c}, amplification, d = ℓ branch.
- **Cheapest partial DIAG reaching 0.87255:** Lemma marked at r + 1.33z ≤ 1 together with Lemma plain at 2m + 1.43z ≤ 1, among other pairs.
- **Residual for 27/31 on top of DIAG** (any one):
  - M_κ with κ ≈ 0.964 (D = U^{1.17}, H ≈ U^{1.13}: family 3.6% shorter than the polynomial);
  - W(0.0196);
  - Gram factor 0.674;
  - w1 = 0.222.
  - DIAG + DH counts gives 0.866826.
  - Coordinator check: M_{0.966} at σ ≈ 0.635 is marginally FALSE in a DH-world (U^{2.31} > D^{1+κ} = U^{2.30}). It is stronger than needed: only counts are needed there.

**Round 2 summary:** all three agents negative. Cheapest known requirement for 27/31: DIAG (coefficient-sensitive sextic sieve for Möbius × Hecke-structured prime slots; a 1.33z partial version suffices) plus a count improvement at δ ≈ 0.27–0.33 below Ingham-with-copies. Both are open.

## 25. Agent report: edge-theta (round 3, finished ~08:40 UTC) — truth has big margin, no proof route

Files: `agents/edge-theta/{EDGE_THETA.md, edge_lib.py, primal_energy.py, dual_mobius.py, pairs_mobius.py, bias_check.py, analyze.py, model_et.py}`.
- **NUMERICAL truth** (exact Z[ω] computations, 12 angular sectors for repeated samples):
  - Primal relative off-diagonal ~ rows^{−0.50±0.08} (theta) and rows^{−0.57±0.07} (marked), with no F = N/K dependence over F = 1.2–312.
  - (H2) cross-mark correlation 1.43·rows^{−0.500±0.012}, the same as random.
  - Dual Möbius-sextic family ~ (#h)^{−0.53±0.02}, beyond the wall too.
  - Extrapolated: W(0.064) holds with margin about Z^{0.43}; (H2) truth ~Z^{−0.45} against a need of Z^{−0.27}.
  - **The needed statements are true with lots of room.** If true as measured: θ = 19/22.
- **Proposition U (CONDITIONAL on P78's energy lemma):** the UNMARKED Möbius-sextic family has a power-saving asymptotic for
  Z ≤ H ≤ Z^{3/2}, with relative error Z^{1+ε}/H. Generic coefficients give none there, so theta beats "affine ⇒ logs" without marks.
- **Precise obstruction.** W sits at F_m ≈ L², deep inside the diagonal regime, not at the SY/Li edge fixed point.
  - For GL2 Hecke-eigenvalue coefficients this regime is unconditional (Perron + edge moments of L-values).
  - For theta coefficients the only tool is the reflection R, which stops at the wall.
  - β* and GRH reproduce exactly the large sieve. Theta structure (twisted multiplicativity, Kazhdan–Patterson periodicity, Kubota residue) generates nothing beyond R.
  - Mark choices (coherent, random, Kubota-aligned) all give random-size cross-correlations.
- **Lemma M(w)** (open, Euler-product form, as hard as (H2)): for P₁ ≠ P₂ ~ L and Z^{1−w}L³ ≤ H ≤ ZL³,
  Σ_{h~H} ψ_{P₁}(h)conj ψ_{P₂}(h)|Σ_{n~Z} a(n)(h/n)₆|² = C(P₁,P₂) + O(Z^{2+ε}). It implies W(w).
- Function fields: affine moments have power-saving truth (DGH secondary terms), but proofs need family-level inputs (DPP functional equation, homological stability).

## 26. Agent report: counts-new (round 3, finished ~08:50 UTC) — NEGATIVE, thresholds sharpened

Files: `agents/counts-new/{notes.md, counts_lab.py, diag_cf.py, count_curves.out, runs/}`.
- **Large-values methods with copies** (thresholds proved by exponent algebra; validity for this family heuristic):
  - Guth–Maynard + copies wins only for δ ≥ 0.41 (x = 0) and δ ≥ 0.42 (x = 1/2).
  - Huxley needs δ ≥ 0.51.
  - Halász–Montgomery: the copies form a rank-one Gram block, so it needs δ > 1/2.
  - Heath-Brown, Jutila, Bourgain: only σ ≥ 3/4.
  - Model payoff 0 (to 1e-9).
- **Copies are exactly u·(F^×)^6** (Kummer theory). Conjugation, units and ray classes add at most a factor 12; subfamilies and base change do not change the ratio.
  - Quadratic-twist copies would have density 1/2, but their Poisson dual is circular (θ = 3/4 sanity check).
  - Higher Möbius moments are worse: count 1/6 + kr(5/6−δ).
- **Theta-based large values: no-go.** Over sparse row sets, C–S and duality see μ only through |μ|², and the theta kernel needs smooth lattice sums.
- **Detector:** the adversary can place r+m anywhere in [t, t+1/2]. At x = 0, R = 1 − 2δλ/(3λ−δ), with
  λ = min(α_long-Möbius, λ_joint(M·S), λ_plain-cubes). So a better joint or plain-cube count would also work at t = 1.
  - Obstruction configuration: U^{R0} rows with Möbius saturated at U^{2t*/3} and U^{t*}, and plain blocks saturated at U^{t*/3}.
    It is consistent with all proved inputs.
- **What 27/31 needs with DIAG** (any one):
  - M_κ with κ = 0.9637;
  - a diagonal Möbius mean square up to length U^{1.043};
  - a joint Σ|M_r S_m|² ≪ U up to r+m = 1.066;
  - a plain sixth moment up to 3m ≈ 1.22.
  In slope form the long-witness slope must drop from 5/6 to 0.616.
- **Conditional payoffs** (open sixth/eighth-moment Lindelöf-on-average for the sextic family, beyond the edge):
  - plain sixth moment for all m ≤ 1/2: **0.8727074** (paper), **0.8694441** (with DIAG);
  - plain eighth moment: 0.8713661 (paper).

## 27. Agent report: alt-detect (round 3, finished ~09:05 UTC) — FIRST NEW STRUCTURE: the cubic θ₃⊗θ₃ architecture

Files: `agents/alt-detect/{NOTES.md, p2_orbit.py, p78_cubic_lp.py, model_cubic.py, conductor_test*.py, check_gauss.py, runs/}`.
- **The quadratic residual is forced (PROVED).**
  - Stickelberger gives μ·Hecke = γ_aγ_bγ_c for a+b+c ≡ 0 (mod N).
  - For single-Gauss-sum rows, Poisson supplies γ_a and the row γ_b, so γ_c must be a ray quantity. Then c ≡ N/2 (c ≡ 0 detects L, not 1/L), and the residual twist has order exactly 2, for every N, field and Kazhdan–Patterson theta.
  - Finite type (conductor exponent 1) is incompatible with μ-absorption for single-Gauss-sum rows. The Cartan matrix is [[2,0,−1],[0,2,−1],[−c_B,−c_B,2]], with det = 8 − 4c_B.
- **New: double-Gauss-sum rows.** The only new case is N = 3, triple (1,1,1).
  - The row coefficients are γ₁², the coefficients of θ₃⊗θ₃ (equivalently the cubic A₂ WMDS residue, or by Davenport–Hasse a sextic Gauss sum).
  - Numerically checked: G̃(π)³ = −π/|π|, so μ = λ̄γ³ and μγ̄ = λ̄γ².
  - The twist is absorbed into both theta factors, so the low-side dual is FULL-INDEX (Kubota index h, no residual twist), with copies u = a³ (copy exponent 1/3).
- **The decisive open number is c_B, the conductor exponent in h of the dual**
  B_h(s) = Σ_{n sf,(n,h)=1} μ(n)ν(n) g̃(χ_n) χ̄_n(h) N(n)^{−s}
  (the coprime part of a twisted sextic Kubota series at index h⁴).
  - Kubota's rank-one theory proves c_B ≤ 4. Heuristics give 3; c_B = 2 is "optimistic".
  - Numerics up to X ≤ 6·10⁴ were inconclusive.

  | c_B | 1 | 2 | 3 | 4 |
  |---|---|---|---|---|
  | P2 shape | 5/6 | 37/42 | 23/24 | 1 |
  | P78 shape, GRH-rows floor | — | 13/17 = 0.765 | 0.8596 | 10/11 |
  | P78 shape, r4 counts transplanted (HEUR) | — | **0.797** | 10/11 | 10/11 |

  - At c_B = 2, sensitivity: without the prime-gain lemmas 0.798; with h_min = 1.25 0.819; with h_min = 1.5 0.872.
- **Normaliser trap.** Small-conductor completions equal B·L(3s−1,(αν)³); used naively, this gives a bootstrap fixed point of exactly 7/8 in P78 shape.
  It is bypassed because B is holomorphic in Re s > 1/2 (sextic representation), but that only matters if c_B ≤ 2.
- **Dead (with proofs or structural arguments):** quartic/quintic/sextic double rows; Λ / log L / L′/L; |θ|², θ³ and Shintani zeta
  (full families give congruences, never γ(A)).
- **Single sufficient new input:** a conductor-N(h)² functional equation (with polynomial growth) for h ↦ B_h(s).
  Then cubic analogues of the count lemmas are needed.

**Coordinator decision (09:10):** settle c_B with two round-4 agents: theory (exact FE, p-part polynomials, orientation)
and numerics (fast Z[ω] Gauss sums, AFE-consistency conductor test).

## 28. Coordinator prediction on c_B (09:45 UTC, before the round-4 agents report)

- Using γ³ = μλ (from G̃(π)³ = −λ(π)), B_h(s) = Σ_n λ(n)ν(n)·conj(g̃(h,n))²·N(n)^{−s}. This is the square of the h-shifted Gauss sum:
  "θ_h ⊗ θ_h", both factors shifted, as alt-detect says.
- **Degree heuristic.** θ₃ is degree 2 over F (P78's reflection kernel R(t) = Π_± Γ(1+t±1/6)/Γ(1−t±1/6) has two Γ_C factors, and the dual length
  involves q_c² for a twist of conductor c). θ⊗θ is degree 4, and a cubic twist of conductor h on a degree-4 object gives conductor h⁴.
  In the sextic picture the Kubota series of order k has gamma factor of degree k−1. That would suggest degree up to 5,
  but Kubota's rank-one FE gives N(h)⁴.
- So the expected answer is **c_B = 4** (dead), unless the coprime-part local polynomials absorb a factor (unlikely: on the critical line
  each P_p ≈ 1 ± N(p)^{−1/2}, so dividing by them does not change the AFE length).
- **Dichotomy (if confirmed).** Single-Gauss-sum rows give degree 2 but a thin (quadratic-residual) family. Double-Gauss-sum rows give a full family
  but degree 4. Both lose against the family/conductor balance. The metaplectic μ-absorption method would then appear capped near 7/8
  with current mean-value technology.
