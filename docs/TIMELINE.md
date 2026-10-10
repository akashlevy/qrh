# Timeline

Local times are PDT and come from file modification times; UTC is marked where used. "We" means Akash Levy working with AI assistants:
- Claude (Opus 5.5) sessions;
- one GPT-6 Astra (Ultra) research session on Oct 9.

## Background

| Date | Event |
|---|---|
| 2026-09-30 | OpenAI posts *The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane Re s > 7/8* (P78) and a Lean formalization in openai/math |
| 2026-10-05 | OpenAI posts the alternate 11/12 paper (P2) |
| 2026-10-06 | Yongxi Lin extracts the 7/8 Lean sources (CoolRmal/quasi-riemann-seven-eighths @ 6ff3bd5) |
| 2026-10-07 | openai/math `fd4aeeb2e`. Kintali posts a 47/48 short proof. |
| 2026-10-08 | Baiying Liu posts 34999/40000 and (1507−2√921)/1653 (arXiv:2610.12234). Argonaut posts 3499999/4000000. |
| 2026-10-08/09 | Several groups publish the cubic-root endpoint θ\* ≈ 0.874957019420 or rationals just above it: GPT-6 Astra & Shunia, Guo, Li Xiang (repo created Oct 9 02:43 UTC), ProofCouncil (Oct 9 11:37 UTC), Nielstron (Oct 9 17:21 UTC). The QRH Bounds tracker launches on Oct 9. |

## Our work

| When | Workstream | What happened |
|---|---|---|
| Oct 8 00:22–00:56 | [01](../01-seven-eighths-lemma-set) | Transcribed P78 Part II bookkeeping. Found 34999/40000 and wrote the exact certificate. Found the LP floor 33216105/37963129 and the GRH-for-rows 13/15 barrier. Ran sensitivity. |
| Oct 8 01:00–03:45 | [02](../02-strengthening) | Explored count-lemma payoffs, the marked-energy wall W(w), ℤ[ω] wall experiments and a first quartic model |
| Oct 8 12:00–14:07 | [02](../02-strengthening) | Multi-agent workflow plus adversarial audit. Three modeling errors in W(w) were found and the optimistic claims **retracted**. Corrected W(1/5) ⇒ 108/125. Identified M_κ (κ = 6/7) as the lemma that gives 6/7. Wrote `QRH-analysis.md` and `STRENGTHENING.md`. |
| Oct 8 17:36–20:54 | [03](../03-six-sevenths) | "6/7" workflow: routes R1–R5, referee passes on R4. Theorem C1, the refined powerful-row count, exact conditional certificates, obstructions. |
| Oct 8 21:10–23:17 | [03](../03-six-sevenths) | Focused check of the quartic transplant: conductor N(h)³, θ = 1, **refuted** |
| Oct 8 23:23 | [03](../03-six-sevenths) | `SIX_SEVENTHS.md`: no unconditional 6/7, so nothing to formalize as 6/7 |
| Oct 8 23:00 – Oct 9 02:04 | [05](../05-lean-23-24) | Lean environment set up. Found that openai/math already proves 23/24. Pivoted to the *shortest* verified QRH: tree-shaking (10 passes), spine, 46-page exposition, Comparator OK. |
| Oct 9 00:37–12:14 | [04](../04-gram-split) | Gram/theta split. Claimed 0.8665/0.871/0.864 (01:30), then the `wf/` audit found F1 and the claims were **retracted** (~03:00). A₂ lead (`wf2/`), five attacks on (V′) (`vls/`), final note (12:14). |
| Oct 9 02:30–13:30 | [05](../05-lean-23-24) | Fork akashlevy/math: paper2 converted to 23/24, `DirichletL` pruned in place (−25.7%), Dirichlet generalization. All 5 Comparator challenges pass. PR akashlevy/math#1. |
| Oct 9 13:00–13:25 | [05](../05-lean-23-24) | Abandoned-routes study of the 6,419 dead declarations |
| Oct 9 ~13:00–16:38 | [06](../06-kappa-bootstrap), [07](../07-method-barriers), [08](../08-weighted-numerator-10499-12000) | GPT-6 Astra session. κ-bootstrap note (13:23); method and feedback barriers (13:20–13:24); independent audits; public prior-work searches; the **weighted-numerator** lemma and proof (14:30–14:40); partial Lean formalization. Stopped by the user one compile error short (handoffs at 16:36–16:38). |
| Oct 9 16:44 – Oct 10 03:22 UTC | [08](../08-weighted-numerator-10499-12000) | Claude (Opus 5.5) session completed the formalization: fixed `NumeratorActualSource.lean:59` and a stale-object import, added `Final.lean`/`FinalAxioms.lean`, isolated rebuild, Comparator, public-source Lake build, final audits. Manuscript updated (18 pp). |
| Oct 10 01:53:26 UTC | [08](../08-weighted-numerator-10499-12000) | Three-kernel judge (con-ron, NanoDa, Lean) **PASS** |
| Oct 10 03:22 UTC | [08](../08-weighted-numerator-10499-12000) | Submitted [GTimG/QuasiRiemannTracker#6](https://github.com/GTimG/QuasiRiemannTracker/pull/6) |
| Oct 9 (evening, PDT) | — | This repository assembled from the working directories and session scratchpads |
