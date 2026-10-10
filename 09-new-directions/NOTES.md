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
