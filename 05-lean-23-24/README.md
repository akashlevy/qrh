# 05 — Shortest verified QRH: a tree-shaken 23/24 Lean proof (Oct 8–9, 2026)

**Result.** `QRH.riemannZeta_ne_zero_of_twentyThree_twentyFourths_lt_re : 23/24 < s.re → riemannZeta s ≠ 0`. It uses only the standard axioms, and Comparator accepts it.
- 292 modules, **121,523 code lines** (5.06M non-whitespace characters).
- That is −30% against the 23/24 proof that shipped inside openai/math, and −72% against the 7/8 proof (441,096 code lines).

There is also a Dirichlet-general version, `OAI.DirichletCharacter.LFunction_ne_zero_of_twenty_three_twenty_fourths_lt_re`. It lives in the fork PR [akashlevy/math#1](https://github.com/akashlevy/math/pull/1) (branch `qrh-23-24-pruned` @ `c388e88e6`). That PR also prunes openai/math's `DirichletL` in place (−25.7%), and Comparator still accepts all 5 QRH challenges afterwards.

## How we got here

1. **Original goal (`notes/HANDOFF.md`).** A Lean proof of ζ ≠ 0 on Re s > 11/12 following paper2, with no 7/8-chain imports.
2. **M1 finding (`notes/MAP.md`).** openai/math already contains an unconditional 23/24 proof, `OAI.ShortDraft.riemannZeta_ne_zero_of_re_gt`. It is paper2's argument at ϑ = 1/10, using integer rows b⁶ with weight φ(𝔫)/N𝔫, so no prime ideal theorem is needed. Neither paper states 23/24.
3. **Pivot 1: shortest proof, any constant.**
   - `lean/scripts/reach.py` computes, inside Lean, every declaration reachable from the root theorem.
   - `lean/scripts/treeshake.py` deletes everything else at the declaration level. It took 10 build passes.
4. **Pivot 2: readable proof.**
   - `QRH/Spine/*.lean` states each step in the paper's numbering: Prop 3.1 `meanSquare`, Prop 3.2 `saving_of_meanSquare` (Y = D^{11/60}), Prop 3.3 `riemannZeta_ne_zero_of_saving`.
   - `paper/qrh.pdf` (46 pp) is the matching human proof.
   - `lean/scripts/check_numbering.py` keeps the numbers cited in the Lean docstrings in sync with the paper.
5. **Pivot 3: upstream PR.**
   - The work moved into a fork of openai/math, because openai/math does not accept outside PRs.
   - paper2 was converted *in place* into a 23/24 paper by `paper2-conversion/convert.py`. The result is byte-exact reproducible; the output is `paper2-conversion/paper2-23-24.pdf`.
   - See `upstream-pr/`.
6. **Abandoned-routes study (`abandoned-routes/`).** This is a forensic reading of the 6,419 dead declarations (about 99k lines) that the pruning deleted, explaining why each was abandoned.

## Layout

| Path | What |
|---|---|
| `lean/QRH.lean`, `lean/QRH/` | The standalone proof. **Hand-maintained** after tree-shaking (`QRH/HAND_EDITED`); keep it verbatim. |
| `lean/ComparatorChallenges/QRH.{lean,json}` | Comparator challenge statement |
| `lean/lakefile.toml`, `lake-manifest.json`, `lean-toolchain` | Lean v4.34.1, Mathlib `d13f23b`. The `OAI` stanza is removed, since lean-base is not vendored. |
| `lean/patches/` | Lean 4.34.1 compatibility patches for PrimeNumberTheoremAnd and rellich-kondrachov, from openai/math `lean/patches/` |
| `lean/scripts/` | `setup-macos.sh`, `prepare-dependencies.sh` (patch path now `./patches`), `measure.py`, `reach.py`, `treeshake.py`, `audit_imports.py`, `golf_probe.py`, `check_numbering.py` (now reads `../paper/qrh.aux`), `run-comparator-*.sh`, and `landrun-shim`, which replaces the Linux sandbox on macOS and runs **unsandboxed** |
| `lean/audit/` | Import-boundary policy, closure report, duplicate-theorem table, golf-probe results, treeshake stats |
| `lean/verification/` | Comparator logs, all ending "Your solution is okay!" |
| `paper/` | `src/*.tex` and `assemble.py` produce `qrh.tex`; compile with `tectonic -X compile qrh.tex --keep-intermediates`. Appendix C maps paper labels to Lean names. |
| `paper2-conversion/` | `convert.py` plus `pieces/`: scripted, exactly-once replacements turning paper2.tex (11/12) into the 23/24 Dirichlet paper. Output PDF included. |
| `upstream-pr/` | PR body, before/after closure and size reports, and the `workspace/` lakefile used to build the fork |
| `abandoned-routes/` | 7-page note `abandoned-routes.{tex,pdf}` and the five agent reports it condenses |
| `notes/` | `HANDOFF.md` (original 11/12 plan), `MAP.md` (paper2 → Lean map and the paper2 gaps), `REPORT.md` (decisions, sizes, reproduction) |

## Sizes (code lines = non-blank, non-comment; `scripts/measure.py`)

| Closure | Modules | Code lines |
|---|---:|---:|
| 7/8 ζ, openai/math as shipped | 2924 | 441,096 |
| 7/8 ζ after the PR's pruning | 2671 | 327,627 (−25.7%) |
| 23/24 as shipped (`RowCompletion.ZeroFreeRegion`) | 286 | 172,968 |
| 23/24 fork `DirichletL.QRH`, after the PR (Dirichlet-general) | 292 | 122,571 |
| **Standalone `QRH` (this folder, ζ only)** | 292 | **121,523** |

Code lines by directory: Eisenstein 36.6k (30%: the cubic theta function is built from scratch), CubicSieve 20.4k, QuadraticSieve 16.2k, GaussSum 13.3k, MeanSquare 7.6k, RowCompletion 6.3k, other about 21k.

## Key lessons and gotchas

- **Comparator compares proof terms of restated constants exactly.** Deleting an unrelated declaration earlier in a file renames shared auxiliary lemmas (`foo._proof_1`). The fix was to exempt from pruning the 5 modules that `HeckeSevenEighths` restates.
- **Pruning removed upstream's Dirichlet 23/24.** Commit 4 restored it for +23 lines.
- **Do not repeat blind automation golf.** One-tactic replacement (`simp`, `grind`, …) succeeded on 4 of 254 theorems, a 0.4% saving. Merging α-equivalent duplicates would save about 157 lines.
- **A paper-form rewrite of the descent was declined.** It would save at most about 15% for weeks of work.
- **Stale build objects.** Pruning left 253 stale `.olean` objects in the fork workspace build. Downstream work must reject imports absent from the fork commit; see `../08-weighted-numerator-10499-12000`.
- **The descent hard-codes 1/20, 1/40, 1/1000 and 3000** in 66 declarations. This blocks 11/12 without a parameterized descent.

## Open leads

See [../docs/OPEN_LEADS.md](../docs/OPEN_LEADS.md). The cheapest is a **15/16 short proof**: the formal mean square already holds for every σ > 1/20 (`outside_mean_square_of_canonical_density_at_radius`, `hσ : 1/20 < σ`). Only the extraction and Mellin steps hard-code ϑ = 1/10 (11/60, 11/10, 21/10, 137/120, 23/24). The estimate is 100–150 lines, not yet compiled.

## Reproduce

Requirements: Lean 4.34.1, about 25 GB of disk, and tectonic for the paper.

```bash
cd lean
bash scripts/setup-macos.sh             # elan + Lean 4.34.1 + patched deps + Mathlib cache
lake build QRH                           # ~4.4k jobs
bash scripts/run-comparator-macos.sh ComparatorChallenges/QRH.json   # needs comparator d03acab + lean4export 076e8e57
cd ../paper && python3 assemble.py && tectonic -X compile qrh.tex --keep-intermediates
cd ../lean && python3 scripts/check_numbering.py   # "numbering OK"
```

Some tools assume the original workspace, where lean-base (openai/math@adc7f124) sits at `lean-base/` and the fork at `math/`. These are `audit_imports.py`, `measure.py`/`importgraph.py` (which default to `lean-base/`; set `QRH_SRC_DIRS`), `reach.py`/`treeshake.py` for regenerating from lean-base, and `run-comparator-upstream.sh`. They are kept for provenance.

To rebuild the fork:

```bash
git clone --filter=blob:none -b qrh-23-24-pruned https://github.com/akashlevy/math
```

Then use `upstream-pr/workspace/lakefile.toml` (srcDir `../math/lean`).

## Provenance and licence

`QRH/` and the spine are **modified derivatives** of openai/math (Apache-2.0), commits `adc7f124` and `fd4aeeb2e`, via CoolRmal/quasi-riemann-seven-eighths@6ff3bd5. `paper/` and `paper2-conversion/` are derived from OpenAI's paper2. See [../NOTICE](../NOTICE).
