# Report: a verified Lean proof of a quasi-Riemann hypothesis, with a matching human proof

*9 Oct 2026.* The goal changed partway through. The brief asked for an 11/12 proof. On 8 Oct the user decided: "shortest possible QRH proof, any constant is OK." On 9 Oct the target became the best tradeoff between Lean size and human readability, with the formal proof and a written proof matching each other.

## Deliverables

| | What | Where |
|---|---|---|
| Formal proof | ζ(s) ≠ 0 for Re s > 23/24, Comparator-accepted, standard axioms only | `QRH.lean`, `QRH/` |
| Spine | each step of the written proof as a Lean theorem, in the paper's notation | `QRH/Spine/*.lean` |
| Human proof | 46-page exposition, ζ only, ϑ = 1/10, integer rows | `paper/qrh.pdf` (source `paper/src/*.tex`, built by `paper/assemble.py`) |
| Correspondence | paper label ↔ Lean name | paper Appendix C; docstrings in `QRH/Spine/Milestones.lean` |

## Result

```lean
-- QRH.lean
theorem QRH.riemannZeta_ne_zero_of_twentyThree_twentyFourths_lt_re
    {s : ℂ} (hs : (23 / 24 : ℝ) < s.re) : riemannZeta s ≠ 0
```

- **Axioms:** only `propext`, `Classical.choice` and `Quot.sound` (`#print axioms`).
- **Comparator:** accepts the theorem against the independently stated `ComparatorChallenges/QRH.lean`. The log is `verification/comparator-qrh.log`; it ends "Lean default kernel accepts the solution". Pinned tools: Comparator `d03acab1` and lean4export `076e8e57`, the same commits as lean-base. Caveat: on macOS the Landrun sandbox (Linux-only) is replaced by `.tools/landrun-shim`, so the build was not sandboxed. That is fine for code we produced; for an adversarial-grade run, use lean-base's Linux workflow.
- **Import audit:** `scripts/audit_imports.py` passes. The `QRH` closure contains no lean-base module and no module of the 7/8 chain.
- **lean-base:** unmodified (`git -C lean-base status` is clean).

## How the proof and the paper line up

The main theorem in `QRH.lean` is a three-line composition, identical to the paper's proof of Theorem 1.1:
- `meanSquare` is Proposition 3.1;
- `saving_of_meanSquare` is Proposition 3.2;
- `riemannZeta_ne_zero_of_saving` is Proposition 3.3.

`meanSquare` is itself `poissonReduction ∘ dualBound` (Propositions 4.5 and 5.1). `QRH/Spine/Milestones.lean` names the remaining results of §§4–7 and Appendix A in the paper's numbering, as aliases of library theorems on the proof path:
- Prop 5.2 (the completed mean square);
- Lemma 5.3 and Prop 5.4 (cube inversion, short cubes, the descent step, the transfer);
- Prop 6.2 (the theta reflection);
- Lemma 6.5 (the quadratic large sieve);
- Lemma 4.1 (Gauss–Jacobi);
- Lemma 4.2 (Poisson with exclusions);
- Lemma 7.2 (the second Poisson identity).

The paper's §§1–3 are written from the Lean proofs, and its §§4–7 follow paper2. Where the formal bookkeeping in §§4–7 differs, the paper says so (Remarks 4.6 and 5.5, Appendix C):
- test functions are expanded in height twists;
- excluded primes are carried as a mask;
- rows are weighted by a compactly supported majorant;
- the descent uses a fixed cube cutoff `D^{1/1000}`, at most 3000 levels, and an explicit loss budget.

**Decision: the descent was not rewritten.** The plan had been to rewrite the descent in the paper's form, with a stop condition. On inspection, the inherited descent is already paper2's descent with different bookkeeping:
- `outsideLogEnergy` is the paper's dual mean square;
- `canonicalLog_actual_step` is cube reduction plus transfer;
- `outsideLogEnergy_finite_rank` is the induction on scale.

The middle layer is about 34k code lines, of which perhaps 8–10k is bookkeeping (ranks, margins, caps, tails, bins). The transfer and Poisson machinery would have to be reimplemented in any form, and a paper-form version was estimated at 10–20k new lines. That is weeks of work for at most about 15% fewer lines overall. The stop condition was therefore applied before starting, and statement-level matching was chosen instead. A full rewrite remains possible as a separate, piloted project.

## Where the proof comes from

lean-base, the extraction of OpenAI's 7/8 formalization, already contains an unconditional theorem `OAI.ShortDraft.riemannZeta_ne_zero_of_re_gt` (23/24 < Re s, s ≠ 1).
- Its closure has 286 modules. It is paper2's argument at the fixed value ϑ = 1/10, which gives 11/12 + 5ϑ/12 = 23/24.
- It averages over integer rows b⁶ instead of prime rows, so it needs no prime counting.
- The 7/8 proof imports these modules as infrastructure but does not use this theorem.
- Neither paper states 23/24.

`MAP.md` gives the paper2 → Lean correspondence and the gaps found in paper2.

`QRH/*` is a **declaration-level tree-shaking** of that closure:
- `scripts/reach.py` computes, inside Lean, every constant reachable from the theorem through types and proof terms.
- `scripts/treeshake.py` copies the 286 modules under the `QRH.` prefix and deletes every declaration that is not reachable.
- Kept anyway: closed rfl lemmas and `simp` lemmas (`simp`/`dsimp` use them without leaving a trace in proof terms), notation, and `… in` prefixes of kept declarations.
- Also handled: `open NS (…)` and `open NS hiding …` lists, duplicate `open`s, and `example`s.
- The process converged after 10 build passes. Each failure exposed one construct class, which was then fixed generically (see `logs/build-qrh-pass*.log`).

## Sizes

`scripts/measure.py` reports the local, non-Mathlib import closure. "Code lines" are non-blank, non-comment lines; non-whitespace characters are included so that formatting cannot fake savings.

| | modules | lines | code lines | non-ws chars |
|---|---:|---:|---:|---:|
| 7/8 proof (`Nonvanishing`) | 2924 | 486,490 | 441,096 | 18.9M |
| 23/24 proof as shipped in lean-base | 286 | 190,290 | 172,968 | 7.04M |
| **QRH (tree-shaken 23/24 + spine)** | **292** | **136,217** | **121,523** | **5.06M** |

That is 30% fewer code lines than the shipped 23/24 proof and 72% fewer than the 7/8 proof.
- **Where it stands:** 107.3k lines are declarations reachable from the theorem. Only about 2.3k lines are conservatively kept but unreachable (notation helpers, rfl/simp lemmas). The rest is headers, `variable` lines and scopes.
- **Where the lines are (code lines):**
  - Eisenstein/ (cubic theta function, built from scratch as a residue of the metaplectic Eisenstein series with Fredholm theory): 36.6k (30%)
  - CubicSieve/ (descent and transfer): 20.4k
  - QuadraticSieve/: 16.2k
  - GaussSum/: 13.3k
  - MeanSquare/: 7.6k
  - RowCompletion/: 6.3k
  - everything else: 21k

## What did not help
- **Merging exact duplicate theorems:** only 22 classes of α-equivalent reachable theorems, about 157 lines (`audit/dup-theorems-2324.tsv`).
- **Blind automation golf:** `scripts/golf_probe.py` tried replacing each whole proof by one of `simp`, `fun_prop`, `positivity`, `omega`, `norm_num`, `linarith` or `grind`. On 254 sampled theorems (3.9k lines, 8 modules) only 4 succeeded, saving **0.4%**. The proofs rely on project lemmas that the automation does not know (`audit/golf/`).

## Options for going shorter (not started)
1. **Cosmetic** (merge `open` lines, drop empty sections and unused `variable`s): about 5k lines, almost no characters. Not a real saving.
2. **Agent-driven proof rewriting**, module by module, with each module compile-checked:
   - Plausibly 20–40% on the verbose machine-written proofs, but expensive: tens of millions of tokens over days.
   - Should be piloted on 2–3 modules first to measure yield and cost.
3. **Structural:**
   - specialize to ζ (χ = 1): small, because most χ code just passes χ through;
   - re-formalize paper2 economically from scratch: the only path to a dramatically shorter proof, but a months-scale project. The theta construction (30%) is the crux.

## Reproduce
```bash
bash scripts/setup-macos.sh                       # elan, Lean 4.34.1, deps, Mathlib cache
python3 scripts/reach.py --import OAI.NumberTheory.DirichletL.RowCompletion.ZeroFreeRegion \
  --root OAI.ShortDraft.riemannZeta_ne_zero_of_re_gt --out audit/reach2324.tsv
python3 scripts/treeshake.py --reach audit/reach2324.tsv \
  --root-module OAI.NumberTheory.DirichletL.RowCompletion.ZeroFreeRegion \
  --src-prefix OAI.NumberTheory.DirichletL --dst-prefix QRH --out-dir QRH
lake build QRH
python3 scripts/audit_imports.py
bash scripts/run-comparator-macos.sh ComparatorChallenges/QRH.json
python3 paper/assemble.py && (cd paper && tectonic qrh.tex)   # the human proof
```
The commands above regenerate `QRH/` from lean-base. Since 9 Oct `QRH/` is hand-maintained (see `QRH/HAND_EDITED`): `treeshake.py` refuses to overwrite it, and later pruning uses `--in-place`.

## Paper2
No mathematical errors were found. Minor gaps and items worth checking are listed in `MAP.md` §4. The most consequential is the Kubota-character convention at paper2 line 3153: it decides whether the dual twist is quadratic. Lean's version is internally consistent.
