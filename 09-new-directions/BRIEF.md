# Shared brief for research agents (read fully before starting)

## Goal
Prove a zero-free half-plane Re s > θ for every finite-order Hecke L-function over F = Q(√−3)
(hence ζ and all Dirichlet L-functions) with θ **meaningfully below 7/8**: at least
27/31 ≈ 0.870968 (bare minimum), ideally 6/7 ≈ 0.857. Every public result is ≥ 0.87491
(re-tunings of OpenAI's 7/8 proof). A re-tuning of constants is NOT acceptable; we need a
genuinely new input (new lemma / new mechanism). Eventually the proof must be formalized in
Lean, so it must be rigorous.

## Sources (local)
- OpenAI 7/8 paper (P78): `/home/user/openai/math/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex`
  (16.7k lines). Key labels/lines: continuation criterion `lem:continuation-criterion` (l.400);
  probe `eq:general-probe` (l.3397); low separation `eq:low-separated` (l.3470); Poisson/high
  identity `eq:probe-poisson-identity` (l.3831) and Euler factorization `lem:local-euler` (l.4010);
  detector `prop:detector-witness` (l.4510); compensated probe (l.6857–6916); reflected energy
  `lem:reflected-energy` (l.7747), hybrid energy `lem:hybrid-energy` (l.7510); compensated
  low estimate (l.8061–8650: `lem:probe-row-norm`, `prop:probe-gram`, `prop:probe-low`);
  local compensation (l.8651–9200); marked inverse moment `lem:marked` (l.9250);
  sixth-power amplification `lem:inverse-amplification` (l.12362); plain fourth moment
  `lem:plain` (l.12531); row counts `prop:detector-counts` (l.15185); high exponent
  `eq:common-high-exponent` (l.15726).
- OpenAI alternate 11/12 paper (P2): `/home/user/openai/math/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex`
  (cleanest form of the mechanism; read the outline §2, l.206–660).
- OpenAI "An unconditional first moment for cubic Gauss sums" (Patterson's conjecture, proved
  unconditionally): `/home/user/openai/math/preprints/An-unconditional-first-moment-for-cubic-Gauss-sums-September-25-2026/build/`
  — `sections/dual.tex` has a power saving BEYOND the sharp cubic large sieve for
  prime-convolution coefficients ("moments of different lengths" + near-orthogonality).
- Prior project work (READ the synthesis first): `/home/user/qrh/docs/{OBSTRUCTIONS,STRATEGIES,RESULTS,OPEN_LEADS}.md`,
  `/home/user/qrh/01-seven-eighths-lemma-set/QRH-analysis.md`, `/home/user/qrh/03-six-sevenths/SIX_SEVENTHS.md`
  and `routes/R1..R5`, `/home/user/qrh/04-gram-split/REPORT.md`, `/home/user/qrh/07-method-barriers/`.
- This session's notes: `/home/user/qrh/09-new-directions/NOTES.md`.
- Exponent model of P78 Part II: `/home/user/qrh/03-six-sevenths/routes/r4/model.py`
  (`python3 model.py '{"wE":0.1}'` etc.; ~1 min per run; knobs: wE wall relaxation, gram_b
  Gram exceptional exponent, kA amplification slope, cP plain capacity, w1 extends marked
  inverse moment past r<1, perfect=True gives DH counts, grh=True GRH for rows).
- Network: arxiv.org and most academic sites are BLOCKED. GitHub (raw.githubusercontent.com,
  git clone) and WebSearch work. Python with numpy/scipy/sympy/mpmath is installed.

## The architecture in one page (my re-derivation; verified against the paper's formulas)
Probe (Part II) I = Q^{-1/2} Σ_{m∈O, q_m≍Q} A_m B_m, Q = XY, with
- A_m = Y^{-3/2} Σ_{s≍Y} c(s) g(χ_s,−m) (sextic Gauss sums; additive polynomial);
- B_m = (prime-marked) completed cubic-theta row Σ_{A=cn³≍Z} c_θ(A) η(A) χ_A(m) (twisted by
  the sextic symbol χ_A(m) = (m/A)_6), marks 1_{p|A} on K disjoint prime slots of total length
  P = Z^ℓ, each marked term rescaled Z→Zp and compensated by p^{-3/2} I(X/p, Y/p, Z).
Geometry: X = Z^{l_x}, Y = Z^{l_y}, b = l_y − l_x, M = l_x + l_y, row range h = 1 + ℓ − l_x.

LOW side (Cauchy–Schwarz in m): ‖A‖² ≪ (Q/Y)(1 + P_a^{1/6} + P_a²/Y) (P_a = Y/X; the
P_a^{1/6} "Gram exceptional term" is the energy of a rank-one sextic-theta vector, a genuine
main term = the shadow of the signal); ‖B‖² ≪ Q (energy; diagonal only below the WALL
M + ℓ ≤ 1, beyond which the quadratic large sieve after theta reflection has a column term).
Result: |I| ≪ X^{1/2} P_a^{1/12} Z^ε.

HIGH side (Poisson in m): I = ∫∫∫ X^{1/2−z} Z^{x+z−1} Y^{w−1} Σ_u q_u^{−z} ζ(6z) L(w,χ_u)/L(x,ηχ̄_u)·(marks)·H.
Principal row u=1 gives the signal X^{1/3} Z^{β−5/6} P^{1/6}/L (copies a⁶ of the principal
frequency, plus slot amplification). A row u with rightmost zero at a = (1+δ)/2 contributes,
after zero detection and counting (count Z^{hR}), the paper's condition
(1+δ)(1−l_y)/2 − β − ℓ/2 + xδℓ + h(R + δ/2 − 1/6) < 0 (x = prime amplitude ∈ [0,1/2]).
Counts R(δ,x) come from saturated witnesses (inverse Möbius polynomial M_r of length U^r,
plain polynomial S_m of length U^m, r+m ≥ t ∈ [1,3/2], m ≤ 1/2) bounded by: marked inverse
moment (r<1, diagonal regime), plain fourth moment, and sixth-power amplification for r ≥ 1
(e(r) = (1+5r)/6, uses the P2 mean square over ALL rows ≤ H for H ≥ D^{1+c}).

## Key numbers (model)
- Paper lemmas re-optimized: 0.874957 (Liu). 13/15 = 0.8667 if GRH for all rows (floor of the
  architecture with unchanged low side). Perfect (DH) counts: 0.869792.
- Energy wall relaxed by w (W(w)): w=0.05 → 0.87183; **w=0.1 → 0.86870 (< 27/31)**; w=0.2 → 0.86333.
- Gram exceptional exponent ×0.5 → **0.86949**; removed → 0.85980.
- Amplification slope kA=0.9 → 0.87116; kA ≤ 0.85 → 0.86979 (needs below-diagonal Möbius mean square M_κ).
- Marked inverse moment extended to r < 1.3 (w1=0.3) → 0.86979.

## Known obstructions (do not rediscover; see docs/OBSTRUCTIONS.md)
Every lever above is a "beyond the edge" mean value (family smaller than length / column term
of a large sieve). Poisson (P) and theta reflection (R) are involutions preserving
|ε−δ| = 1−κ (move calculus) — no sequence of P, R, large sieves proves a below-diagonal bound.
Coefficient-blind large-values bounds cannot beat the sixth-power copy multiplicity. The
Gram excess is the signal's shadow (same Mellin coefficient). Rescaling marks Z→Zp^κ is a
reparametrization (θ invariant in κ). The quartic transplant over Q(i) gives θ = 1.
New input must come from somewhere other than {Poisson, reflection, large sieve, positivity}:
e.g. arithmetic structure of the coefficients (primes, Möbius, theta Hecke relations),
analytic information (a priori QRH β* ≤ 0.8749, zero-free rectangles), spectral theory,
or a different probe.

## Working rules
- Be rigorous and skeptical. Separate PROVED / CONDITIONAL / HEURISTIC. Check every claimed
  exponent gain with the model (or a modified copy of it in your own directory).
- Write all notes/scripts under `/home/user/qrh/09-new-directions/agents/<your-name>/`.
  Do not modify files outside that directory. Do not commit or push.
- Final message: ≤ 600 words: what you tried, what you proved (with precise statements),
  the exponent payoff, and the precise obstruction if it failed. Point to your files.
