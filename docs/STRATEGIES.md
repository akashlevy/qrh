# Strategies tried

Every route attempted, in rough chronological order, with its outcome and the reason. The outcome labels are:
- ✅ worked
- ⚠️ conditional only
- ❌ failed or obstructed
- ↩️ retracted
- ⏸ not pursued

## Pushing the 7/8 architecture (Oct 8)

| Route | Outcome | Why | Where |
|---|---|---|---|
| Re-optimize P78's fixed geometry (b, ℓ) | ✅ 34999/40000 | Only the geometry changes; the endpoint gets a Bernstein certificate. It is essentially the lemma-set optimum (Liu: 0.874957069799; floor 0.8749570). | [01](../01-seven-eighths-lemma-set) |
| Improve counting lemmas (c_P, c_M, α, perfect counts) | ❌ ≤ 0.008 | Even perfect counts give 167/192 or 13/15. The bottleneck is the low side. | [01](../01-seven-eighths-lemma-set), [02](../02-strengthening) |
| Relax the marked-energy wall W(w) | ↩️ then ⚠️ | The first draft's claims (6/7 from W(1/5)) were retracted after audit. Corrected: W(1/5) ⇒ 108/125; 6/7 needs w = 4/7, or w ≈ 0.285 with the refined count. | [02](../02-strengthening), [03](../03-six-sevenths) |
| Lower the detector floor 51/100 → 1/2 + η | ✅ but no gain | Legitimate, but it buys nothing alone | [02](../02-strengthening) |
| Remove the Gram P_a^{1/6} term (sensitivity) | ⚠️ 0.864 if possible | It is a genuine main term; see the Gram split below | [01](../01-seven-eighths-lemma-set) |
| Below-diagonal mean square M_κ | ⚠️ | κ = 6/7 ⇒ 6/7 exactly via either extraction. It is the cleanest single lemma, and it is open. | [02](../02-strengthening) |
| Quartic setup over Q(i), first heuristic | ↩️ | ≈ 0.81 via theta (needs Patterson's quartic conjecture). Later refuted. | [02](../02-strengthening) |

## The 6/7 workflow (Oct 8)

| Route | Outcome | Why | Where |
|---|---|---|---|
| R1: prove M_κ directly | ❌ | Tried: positivity, Φ(0) = 0 dyadic rows, double reflection, asymptotic evaluation, rewriting via γ₁(n; u⁵), Heath-Brown/Vaughan + Hölder, zero-free region + large sieve, other families, GRH for non-principal rows, asymptotic large sieve, D-averaging, mollifiers, prime-amplified weights. The column term D² reappears at the positivity step (Prop O1); moves cannot escape (move calculus). By-product: Theorem C1. | [03](../03-six-sevenths) R1 |
| R2: bootstrap H(β) | ❌ (by-products ✅) | Consistency barrier W_β. Rigorous by-products: Lemma B, Theorem C. | [03](../03-six-sevenths) R2 |
| R3: zero density | ❌ | Coefficient-blind impossibility; the multiplicity-beating LV hypothesis ⇔ the conclusion; 13/15 floor | [03](../03-six-sevenths) R3 |
| R4: map every lemma replacement in Part II | ⚠️ | Weakest single certified hypothesis: W(2/7) (fragile) or W(3/10) + the proved refined count | [03](../03-six-sevenths) R4 |
| R5: literature | ❌ | Nothing usable. de Faveri's optimal n-th order large sieve changes θ in neither paper. | [03](../03-six-sevenths) R5 |
| Prime amplifier in the mean square | ❌ | 11/12 under the proved walls; 6/7 needs the lossless marked moment (C4) | [03](../03-six-sevenths) R1 |
| Quartic transplant (Q(i), Kubota series, de Faveri n = 4) | ❌ refuted | The true conductor is N(h)³ and the dual twist is quartic, giving θ = 1. The earlier 0.854 came from the wrong conductor. | [03](../03-six-sevenths) `quartic-transplant/` |

## Gram / theta split (Oct 9 early)

| Route | Outcome | Why | Where |
|---|---|---|---|
| Theta split + Voronoi (V) + Heath-Brown cubic large sieve | ↩️ 0.8665 / 0.871 / 0.864 | F1: the pairing is against π = ρζ | [04](../04-gram-split) |
| Cancel the excess with signed probe combinations | ❌ | The signal has the same Mellin coefficient | [04](../04-gram-split) |
| Clean Gram via (FE)/(LoA6) alone | ⚠️ but no gain | Proved modulo 3 lemmas; gives exactly 0.874957 | [04](../04-gram-split) `wf/` |
| Unreflected route T = 𝒮 + T′ | ❌ | The T′ bound is worse than Cauchy–Schwarz; 𝒮 is the signal | [04](../04-gram-split) `wf/` |
| A₂ weighted multiple Dirichlet series | ❌ | Conductor N(μ)³, coupled dual: min(0.874957, 0.934) | [04](../04-gram-split) `wf2/` |
| Attacks A–E on (V′): functional-equation geometry, Gauss-sum structure, MDS/spectral, probe redesign, numerics | ❌ s₀ = 0 | See [OBSTRUCTIONS](OBSTRUCTIONS.md) §4 | [04](../04-gram-split) `vls/` |
| Kuznetsov on the dual form (Idea 2) | ❌ | No metaplectic spectral large sieve | [04](../04-gram-split) |
| Prove the energy wall W(w) directly (Idea 3) | ❌ | Beyond the edge on both sides | [04](../04-gram-split) |
| Other μ-absorbing families or fields (Idea 4) | ❌ | The sextic case is unique | [04](../04-gram-split) |

## Formalization (Oct 8–9)

| Route | Outcome | Why | Where |
|---|---|---|---|
| Lean 11/12 per paper2 | ⏸ dropped | A 23/24 proof already existed; 11/12 is a months-scale side goal | [05](../05-lean-23-24) |
| Declaration-level tree-shaking, computed in Lean | ✅ −30% | 173k → 121.5k code lines; Comparator-accepted | [05](../05-lean-23-24) |
| Paper-numbered spine + exposition paper | ✅ | 5 spine files, 46-page paper, numbering checked by script | [05](../05-lean-23-24) |
| In-place pruning of openai/math `DirichletL` | ✅ −25.7% | All 5 challenges pass; proof-term identity gotcha | [05](../05-lean-23-24) |
| Convert paper2 into a 23/24 paper by scripted replacements | ✅ | Reproducible byte-for-byte | [05](../05-lean-23-24) |
| One-tactic automation golf | ❌ 0.4% | The proofs depend on project lemmas | [05](../05-lean-23-24) |
| Merge α-equivalent duplicates | ❌ about 157 lines | Negligible | [05](../05-lean-23-24) |
| Paper-form descent rewrite | ⏸ declined | At most about 15% for weeks of work | [05](../05-lean-23-24) |
| Agent-driven module-by-module proof rewriting | ⏸ not started | Estimated 20–40% at tens of millions of tokens; pilot first | [05](../05-lean-23-24) |

## The GPT-6 Astra session and its completion (Oct 9–10)

| Route | Outcome | Why | Where |
|---|---|---|---|
| κ-feedback bootstrap (lem:plain for κ ≥ 13/18) | ✅ but superseded, not novel | Gives 43747851/50000000 and converges to the public cubic root. The κ-feedback mechanism and its limit were already public (Li Xiang; also Guo, Shunia/GPT-6 Astra); they are cited, not copied. | [06](../06-kappa-bootstrap) |
| Unweighted central L⁴ of the numerator | ❌ | Did not beat the public feedback limit | [08](../08-weighted-numerator-10499-12000) `audits/independent_count.txt` |
| Replace global κ by the parent zero-bin δ | ❌ invalid | Child rows need not inherit the bin | same |
| Hölder of the inverse and plain moments | ❌ | Gives only convex combinations | same |
| **Weighted original-numerator fourth moment** (full w-line, reflection, selected primes, Hölder per amplitude class, quotient-free Euler expansion) | ✅ **10499/12000** | Escapes the pointwise-numerator barrier on the high side | [08](../08-weighted-numerator-10499-12000) |
| Lean formalization: partial (GPT-6 Astra), then completion (Opus 5.5) | ✅ | The prior session stopped one compile error short (`NumeratorActualSource.lean:59`). The completion also fixed a stale-object import. | [08](../08-weighted-numerator-10499-12000) `history/` |
| Projection onto ρ; angular modulation; Gauss-phase cancellation | ⏸ inconclusive | A reformulation only / needs a new high identity / τ₆(m³) unknown | [07](../07-method-barriers) `independent_low.txt` |
| Signed dual of the P2 mean square | ⚠️ partial | The unit j-block gap | [07](../07-method-barriers) |

## Toward θ < 41/47 (Oct 10 afternoon)

Session 09 (`09-new-directions/`, Oct 10 morning) mapped the levers of the pointwise-numerator architecture; see its `SUMMARY.md`. Session 10 redid the analysis for the weighted architecture of 10499/12000. Everything here is in [10](../10-below-41-47); the conclusions are collected in [MISSING_ESTIMATE_41_47.md](MISSING_ESTIMATE_41_47.md).

| Route | Outcome | Why | Where |
|---|---|---|---|
| Weighted exponent model of the high side | ✅ tool | Reproduces the 10499/12000 certificate; locates the 41/47 deficit (+0.011) at δ ≈ 0.47 with maximal prime amplitude; the numerator, not the counts, is the bottleneck | `model/`, `NOTES.md` §2 |
| Numerator × plain witness; copies; sextic large sieve; β*-numerator; Halász–Montgomery; separate numerator/denominator bins | ❌ | Each reproduces the current bound or is worse | `NOTES.md` §3 |
| Numerator contour Re w = c with moment convexity; low-side Cauchy–Schwarz on bad rows; multi-scale numerators | ❌ | Affine in c (an endpoint is optimal); needs an M_κ-type large-values bound; additive | `NOTES.md` §5 |
| Probe redesign (s-weights, product numerators, annihilation, coupled probes) | ❌ | All worse or neutral; new open sufficient condition: a restricted large sieve on the bad rows | `agents/redesign/` |
| Beat the fourth-moment numerator bound | ❌ with a barrier | A flat configuration meets every available input at σ_wt; a sixth moment would need excess ≤ 0.097 (provable: 0.34) | `agents/numerator/` |
| Raise `lem:marked`'s capacity | ❌ (Lemma D ✅) | Reduces to cross-tuple terms, which need ≈ 60% of square-root cancellation across marks | `agents/marked-cap/` |
| Exploit the prime amplitudes (slot powers, amplitude windows) | ❌ (slot-power lemma ✅) | Unlimited prime supply, but no exponent effect at the paper's capacities | `agents/amplitude/` |
| Raise `lem:plain`'s capacity | ❌ (Theorem G ✅) | The region is optimal for the proof's moves; 41/47 needs a twisted fourth moment at the edge | `agents/plain-cap/` |
| Cross-mark cancellation | ❌ (Prop X1 ✅) | Needs a Möbius fourth moment beyond the diagonal; cross terms are random-size numerically | `agents/xmark/` |
| Combine small savings | ❌ | The levers are additive; best θ from proved inputs is 0.8748866 | `agents/combo/` |
| Relax the energy wall | ❌ (rank-one reduction ✅, conditional) | W(w) reduces to one core sum; Lemma-D split, large-values dichotomy, slot groups and the cube part all fail | `agents/wall4/` |
| Literature | ❌ | No published tool supplies any of the needed inputs; tracker PR #6 merged | `agents/lit2/` |
| Gram exceptional factor g < 1 | ⚠️ g ≤ 0.739 would give 41/47 | Equivalent to the signal-circular pairing problem of 04/07 | `NOTES.md` §6 |

## The energy-wall program (Oct 10 evening)

The lossless wall W(w) was judged the least hopeless input, so it got a three-step program. Everything is in [11](../11-energy-wall).

| Route | Outcome | Why | Where |
|---|---|---|---|
| Step 1: exact Z[ω] numerics at the weighted geometry | ✅ evidence only | W holds numerically: random-matrix off-diagonal, margin ≈ Z^0.41 | `numerics/` |
| Step 2: the single-mark case (extend Proposition U to one twist) | ❌ (Proposition U_P ✅) | The only lossy step is the quadratic large sieve on the reflected dual; Perron, Hecke, spectral, first-moment and Poisson-in-P routes each fail at a stated step | `single-mark/` |
| Step 3: function-field analogue over F_q[T] | ✅ evidence; ❌ for θ | W holds there with random-matrix margin and is provable at large q under hypotheses (Proposition A), but the mechanisms (Deligne, homological stability) have no number-field substitute; the C̃₂ functional-equation invariant explains why finite functional equations stop at the wall | `function-field/` |
