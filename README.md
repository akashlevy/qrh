# QRH: zero-free half-planes below 7/8

This is a research log, proofs, certificates and Lean formalizations from work on **quasi-Riemann hypothesis (QRH) bounds**, done 8–10 October 2026 by Akash Levy with AI assistants. A QRH bound is a zero-free half-plane Re s > θ, with θ < 1, for ζ, every Dirichlet L-function and every finite-order Hecke L-function over Q(√−3).

The starting point was OpenAI's 7/8 paper and Lean development ([openai/math](https://github.com/openai/math)). This repo collects everything we found, including:
- results;
- partial and conditional results;
- obstructions and dead ends;
- strategies tried, and why each succeeded or failed;
- prior work and online resources.

## Headline results

| θ | Result | Evidence | Where |
|---|---|---|---|
| **10499/12000 ≈ 0.8749167** | **β\* ≤ 10499/12000.** ζ, all Dirichlet and all Q(ζ₃) Hecke L-functions are nonzero for Re s > 10499/12000. This is about 4·10⁻⁵ below every other framework-verified bound (the cubic-root endpoint ≈ 0.874957019). Hailey Collet's pending tracker PR #9 (0.8748866) re-optimizes this construction. | **Lean**, standard axioms only. Isolated rebuild, Comparator, and a con-ron + NanoDa + Lean three-kernel PASS (author-side, unsandboxed). Merged into the tracker on 10 Oct as its best framework-verified record ([GTimG/QuasiRiemannTracker#6](https://github.com/GTimG/QuasiRiemannTracker/pull/6)). Rebuilt from source on Linux on 10 Oct (endpoints: standard axioms only). No human peer review yet. | [`08-weighted-numerator-10499-12000/`](08-weighted-numerator-10499-12000) |
| 23/24 | **The shortest verified QRH proof:** ζ(s) ≠ 0 for Re s > 23/24. It is 121.5k code lines, 72% smaller than the 7/8 proof, with a paper-numbered spine and a matching 46-page paper. A Dirichlet version and an in-place pruning of openai/math (−25.7%) are in [akashlevy/math#1](https://github.com/akashlevy/math/pull/1). | **Lean**, Comparator (all 5 QRH challenges after pruning) | [`05-lean-23-24/`](05-lean-23-24) |
| 34999/40000 | The 7/8 lemma set with re-optimized geometry. Liu found the same value independently the same day. | Exact rational certificate (paper-level) | [`01-seven-eighths-lemma-set/`](01-seven-eighths-lemma-set) |
| 43747851/50000000 | κ-feedback bootstrap. **Superseded and not novel.** | Exact certificates; partial Lean | [`06-kappa-bootstrap/`](06-kappa-bootstrap) |
| 6/7 (conditional) | Theorem C1: the below-diagonal mean square M_{6/7} ⇒ 6/7. Exact conditional certificates also give 6/7 under refined count + W(3/10), and under GRH-for-rows + Gram fix. | Proof + exact certificates | [`03-six-sevenths/`](03-six-sevenths) |
| 41/47 ≈ 0.8723 (conditional, exponent model) | Any one of about ten estimates past a large-sieve edge would give θ < 41/47. The cleanest is a **lossless energy wall** W(0.041): a second-moment bound for quadratic twists of twisted cubic theta series. It holds numerically (margin ≈ Z^0.41) over Q(ω) and in an F_q[T] model, but no known method reaches it. With proved inputs the method stops at 0.8748866. | Exact LP exponent model; numerics; written up for experts in [docs/MISSING_ESTIMATE_41_47.md](docs/MISSING_ESTIMATE_41_47.md) | [`10-below-41-47/`](10-below-41-47), [`11-energy-wall/`](11-energy-wall) |

### Key obstructions

The full list is in [docs/OBSTRUCTIONS.md](docs/OBSTRUCTIONS.md).

**Ceilings of the 7/8 method**
- **The 7/8 lemma set stops at 0.8749570.** The exact LP floor is 33216105/37963129.
- **κ-feedback stops at the cubic root θ\* ≈ 0.874957019420.** This is why every other public bound clusters there. The weighted numerator escapes it by changing the high-side numerator estimate.

**The 13/15 wall**
- **The low side caps the whole architecture at 13/15**, even assuming GRH for every auxiliary row. Anything below 13/15 needs a new low-side estimate.

**Routes ruled out**
- **6/7 needs a mean value below the large-sieve diagonal.** Coefficient-blind tools provably cannot supply one; the needed hypothesis is strictly stronger than its conclusion.
- **The Gram/theta split** only gives 0.874957: its claimed 0.8665 was retracted because the pairing is against π = ρ·ζ_m.
- **The quartic (Q(i)) transplant gives nothing (θ = 1).** The reflection conductor is N(h)³.
- **θ = 1/2 is out of reach** for this method without GRH-strength input.

**Below 41/47 (sessions 10–11)**
- **Small savings don't add up.** The five power-saving levers of the weighted argument are additive to within 3%, so 41/47 needs about a 0.01 saving in five lemmas at once (or one large saving).
- **The energy wall stops every known method.** Its loss is a single quadratic-large-sieve column term, sharp for general coefficients. Every finite functional equation preserves an invariant that the wall estimate needs to change. Function fields get past such walls with algebraic geometry, which has no number-field counterpart.

## Repository map

| Directory | Date | Topic |
|---|---|---|
| [`01-seven-eighths-lemma-set/`](01-seven-eighths-lemma-set) | Oct 8 | Re-optimize P78's geometry (34999/40000); exact LP floor; GRH-for-rows 13/15; sensitivity |
| [`02-strengthening/`](02-strengthening) | Oct 8 | What moves θ: W(w) energy wall (audited; several claims retracted), M_κ ⇒ 6/7, numerics |
| [`03-six-sevenths/`](03-six-sevenths) | Oct 8 | Five routes toward 6/7 (R1–R5): conditional theorems, exact certificates, obstructions, quartic refutation |
| [`04-gram-split/`](04-gram-split) | Oct 9 | Removing the P_a^{1/6} Gram loss: theta split, Lemma Z, retraction, five attacks on (V′) |
| [`05-lean-23-24/`](05-lean-23-24) | Oct 8–9 | Shortest verified QRH (23/24): tree-shaken Lean, spine, exposition paper, upstream PR, abandoned-routes study |
| [`06-kappa-bootstrap/`](06-kappa-bootstrap) | Oct 9 | κ ≥ 13/18 for `lem:plain`, then the feedback bootstrap (historical) |
| [`07-method-barriers/`](07-method-barriers) | Oct 9 | κ-feedback cubic barrier, 13/15 low-majorant barrier, low-side and signed-dual no-gos |
| [`08-weighted-numerator-10499-12000/`](08-weighted-numerator-10499-12000) | Oct 9–10 | **Headline result**: manuscript, Lean package, certificates, kernel evidence, audits, submission |
| [`09-new-directions/`](09-new-directions) | Oct 10 | 18 agents on every lever of the 7/8 architecture: Proposition U, edge-theta numerics, (H1)/(H2) ⇒ 19/22 (conditional), why each lever needs a beyond-the-edge input |
| [`10-below-41-47/`](10-below-41-47) | Oct 10 | Toward θ < 41/47: weighted exponent model, nine agents (all negative), proved-input frontier 0.8748866, single-input thresholds, Lemma D, slot powers, Theorem G |
| [`11-energy-wall/`](11-energy-wall) | Oct 10 | The energy-wall program: exact Z[ω] numerics (margin Z^0.41), single-mark reduction (Proposition U_P), function-field analogue |
| [`docs/`](docs) | — | Cross-cutting synthesis (below) |
| [`tools/`](tools) | — | `render_pdf.swift` (PDF page renders for review) |

**Synthesis documents**
- [docs/RESULTS.md](docs/RESULTS.md): every result by strength of evidence (Lean / certified / conditional / partial / numerical / **retracted**), plus errata.
- [docs/OBSTRUCTIONS.md](docs/OBSTRUCTIONS.md): every barrier, with precise statements and status.
- [docs/STRATEGIES.md](docs/STRATEGIES.md): every route tried and why it worked or failed.
- [docs/OPEN_LEADS.md](docs/OPEN_LEADS.md): ranked next steps for mathematics, formalization and verification.
- [docs/MISSING_ESTIMATE_41_47.md](docs/MISSING_ESTIMATE_41_47.md): a short note for experts: what would give θ < 41/47, the energy-wall estimate, the evidence, and why known methods stop.
- [docs/PRIOR_WORK.md](docs/PRIOR_WORK.md): OpenAI's papers, the competing-claims table, and literature by topic.
- [docs/RESOURCES.md](docs/RESOURCES.md): our repos and PRs, upstream sources, Lean pins, verification tools.
- [docs/TIMELINE.md](docs/TIMELINE.md): what happened when.

Each numbered directory has its own README with findings, layout and reproduction commands.

## Quick checks

```bash
python3 01-seven-eighths-lemma-set/certificate.py                  # 34999/40000: ALL CHECKS PASS
python3 01-seven-eighths-lemma-set/barriers.py                     # LP floor 0.8749570 and 13/15
python3 07-method-barriers/feedback_barrier_certificate.py         # κ-feedback cubic barrier
python3 07-method-barriers/method_barriers.py                      # 13/15 low-majorant, 6/7 extraction floor
python3 08-weighted-numerator-10499-12000/certificates/weighted_numerator_certificate.py --output /tmp/cert.json   # 10499/12000 exponent certificate
```

All five scripts above are stdlib-only and passed when this repo was assembled. Lean builds need Lean 4.34.1, Mathlib `d13f23b` and about 25 GB of disk; see `05-…/README.md` and `08-…/README.md`. For the other scripts, `pip install -r requirements.txt`.

## What is not in this repo

- **Build artifacts and clones,** about 50 GB: `.lake`, oleans, the openai/math clone, toolchains. The READMEs give the pins and commands to rebuild them.
- **Regenerable numerics caches,** about 190 MB of `.npy`/`.npz`/`.bin`. Each README gives the regeneration commands.
- **Third-party papers and text dumps:** linked in `docs/PRIOR_WORK.md` instead.
- **Raw search API dumps:** they contained third-party e-mails.
- **Li Xiang's work:** nothing of his was copied, per the project's standing rule.

Some scripts still hard-code the original workspace paths (`/Users/akashlevy/Documents/riemann/...`). The affected READMEs say how to relocate them.

## Attribution and licence

Author: Akash Levy. For the 10499/12000 result:

> GPT-6 Astra (Ultra) was used to discover the proof and perform partial formalization; Opus 5.5 Ultracode was used to complete the formalization and submit the result.

The other workstreams were done with Claude (Opus 5.5).

Lean and TeX derived from openai/math are Apache-2.0 ([NOTICE](NOTICE), [LICENSE-APACHE-2.0.txt](LICENSE-APACHE-2.0.txt)).
