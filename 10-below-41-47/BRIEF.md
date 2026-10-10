# 10 — Brief: θ < 41/47 ≈ 0.8723404 (session of 10 Oct 2026, afternoon)

## Goal
Prove a zero-free half-plane Re s > θ with **θ < 41/47 = 0.87234043** for every finite-order
Hecke L-function over F = Q(√−3) (hence ζ and all Dirichlet L-functions), rigorously enough to be
formalized in Lean 4 on top of the existing OpenAI 7/8 development and our 10499/12000 package
(`08-weighted-numerator-10499-12000/`). A re-tuning of constants cannot reach this; a genuinely new,
PROVABLE input is required. Separate PROVED / CONDITIONAL / HEURISTIC rigorously.

State of the art: ours 10499/12000 = 0.8749167 (pending, tracker PR #6); best pending public
0.874887 (Hailey Collet, tracker PR #9 = our weighted numerator + variable κ ∈ [7/10,3/4] in
lem:plain + completed zero-slot moment). Nothing public below 0.8748.

## Sources (local)
- P78 (OpenAI 7/8): `/home/user/openai/math/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex`
  (key labels: lem:continuation-criterion l.400; probe; lem:marked l.9218ff; lem:inverse-amplification l.12343;
  lem:plain l.12477ff; prime amplitudes l.14995; prop:detector-counts l.15180ff; lem:high-bin & eq:common-high-exponent l.15700ff).
- P2 (OpenAI 11/12): `.../The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex` (outline l.206–660).
- Our weighted numerator: `/home/user/qrh/08-weighted-numerator-10499-12000/paper/weighted_numerator_{proof,lemma}.tex`,
  `audits/independent_count.txt`.
- Hailey Collet PR #9 paper: `/home/user/ext/hailey/public_proofs_paper.md` (public; may be read and cited).
- Prior synthesis (READ FIRST): `/home/user/qrh/docs/{OBSTRUCTIONS,STRATEGIES,RESULTS,OPEN_LEADS}.md`,
  `/home/user/qrh/09-new-directions/{SUMMARY,NOTES}.md` (38 sections, 18 agents: every lever of the old
  pointwise-numerator model was shown to need a beyond-the-edge input).
- PROJECT RULE: do not fetch or copy Li Xiang's (lixiang90) work.

## Exponent model of the WEIGHTED high side (this session; `model/gmodel.py`, `model/wmodel.py`)
Geometry: l_x, l_y, slot length ℓ, M = l_x + l_y, b = l_y − l_x, h = 1 + ℓ − l_x (row range), L = ℓ/h.
Low side (proved): B = 5/6 + M/12 − ℓ/6 (+ excess beyond wall M+ℓ ≤ 1); on the wall B = 11/12 − ℓ/4.
**41/47 needs ℓ ≥ 0.1773 on the wall** (currently the high side only closes at ℓ ≈ 0.167).
High side at row scale U = Z^d, d = h, zero bin a = (1+δ)/2, prime amplitude q = δx, x ∈ [0,1/2]:
- E_old = a − B + h(z_c − 1/6) − a·l_y − ℓ/2 + qℓ + d(R + δ/2 − z_c)       (P78 pointwise numerator)
- E_wt  = a − B + h(z_c − 1/6) − l_y/2 − ℓ/2 + qℓ + d((3R+1)/4 − z_c) − q(2l_y − d)/(12κ)  (weighted)
- need min(E_old, E_wt) < 0 for all (δ, x); z_c = 17/50; κ = 3/4.
- R(δ,x) = 1 − δ + (α−δ)δP_x/(2J), α = 5/6, c = 2/9, D_x = 3 − (1+4c)x, P_x = (2−4cx)(1−x), J = (α−δ)D_x + δP_x
  (max over adversary witness split r+m = t of: inverse-marked 1 − δ(r + 2q·(1−r)/2 /δ), plain 1 − δ(2m) − 2q(1−2m)c,
  long inverse 1 − δ + (α−δ)(r−1); min over t ∈ [1,3/2]).
- Numerator: Σ_{u∈C}|S_u| ≤ U^{R+σ}, σ_wt = (1−R)/4 − q z_N/2, z_N = (1 − 2n₀)/(6κ), n₀ = 1 − l_y/d.

**Deficit at ℓ = 0.1773, b ≈ 0.075:** worst min(E_old,E_wt) ≈ +0.011, at δ ∈ [0.3, 0.6], x ≥ 0.2 (max at
δ ≈ 0.45–0.5, x = 1/2, R ≈ 0.59). All bins with x = 0 are fine. At δ = 0.47, x = 1/2 the allowed numerator
exponent is σ ≤ 0.0857 vs σ_wt = 0.0995 (need −0.014 in U-units), or equivalently R ≤ 0.571 vs 0.590.

**Single hypothetical inputs that WOULD reach 41/47 in the weighted model (all currently unproved):**
| input | effect at ℓ = 0.1773 |
|---|---|
| numerator 6th moment Σ_u|S_u(n₀)|⁶ ≪ U^{1+0.10} (n₀ ≈ 0.44) | −0.012 (enough) |
| DIAG: lem:marked at r+z ≤ 1 and lem:plain at n₁+n₂+z ≤ 1 | −0.008 (enough) |
| lem:marked capacity 0.8(1−r) + lem:plain capacity 0.4(1−2m) | −0.0001 (just enough) |
| lem:marked capacity 0.75(1−r) + amplification slope α ≤ 0.75 | −0.0001 |
| lossless energy wall W(0.045) (M+ℓ = 1.045 at no low-side cost) | ≈ −0.001 |
| DH row counts | −0.028 |
| lem:marked capacity (1−r) alone | +0.003 (not enough) |
| lem:plain capacity 1 alone | +0.0044 (not enough) |

## Known obstructions (do not rediscover; see docs/OBSTRUCTIONS.md and 09 NOTES)
Poisson/reflection/large-sieve "move calculus" cannot cross the diagonal; coefficient-blind large values
cannot beat the a⁶ copy multiplicity; Gram P_a^{1/6} term is the signal's shadow; marks cost conductor p²
(wall M+ℓ ≤ 1 and the factor 2 in lem:marked); DIAG false for general slot coefficients (a(p) = sextic Gauss sum);
lem:plain's prime capacity 1/(6κ) comes from bounding slots pointwise by β*; β* ≈ 0.87 is too weak for Halász–
Montgomery (needs off-diagonal P^{1/2}); L²-in-Z averaging loses (count grows faster than T₀^{1/2}); κ-rescaled marks
are reparametrizations; multi-Poisson and higher-rank theta splits are worse; mollified/structured s-weights must
keep the w=1 pole of the u=1 numerator.

## Working rules for agents
- Be rigorous and skeptical; check every claimed gain with the model (copy `model/gmodel.py` into your directory).
- Write notes/scripts only under `/home/user/qrh/10-below-41-47/agents/<your-name>/`. Do not commit or push.
- Final message ≤ 600 words: what you tried; what you PROVED (precise statements, with proofs written in your notes);
  exponent payoff in the weighted model; the precise obstruction if it failed.
