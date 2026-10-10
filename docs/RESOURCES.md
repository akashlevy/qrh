# Online resources, repositories, tools and pins

## Our repositories and pull requests

| Resource | What |
|---|---|
| [akashlevy/math](https://github.com/akashlevy/math), branch `qrh-23-24-pruned` @ `c388e88e6519d139ad2ca5571ef0c7fbf0e245b5` | Fork of openai/math. It holds the 23/24 Dirichlet proof, the paper2 → 23/24 conversion and the in-place pruning. Tag `qrh-weighted-10499-12000-base` points at `c388e88e6`, the base of the weighted proof. **This is not an OpenAI commit.** |
| [akashlevy/math#1](https://github.com/akashlevy/math/pull/1) | "QRH: a 23/24 Lean proof matched step by step to paper2 (converted in place), plus in-place pruning". Open, for discussion. It was opened on the fork because openai/math does not accept outside PRs. Upstream compare: <https://github.com/openai/math/compare/main...akashlevy:math:qrh-23-24-pruned> |
| [akashlevy/QuasiRiemannTracker](https://github.com/akashlevy/QuasiRiemannTracker), branch `akashlevy-weighted-numerator-20261009` @ `2fd60c0` | Tracker fork carrying the 10499/12000 package (`proofs/akashlevy-20261009-weighted-numerator/`, commit `ec626ec`) |
| [GTimG/QuasiRiemannTracker#6](https://github.com/GTimG/QuasiRiemannTracker/pull/6) | "Add 10499/12000 weighted-numerator QRH result (Akash Levy)", open |

## Upstream sources

- **openai/math:** <https://github.com/openai/math>, Apache-2.0. Commits `adc7f1241b42e322a6451854ab7e4b4c146bf78a` (initial) and `fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb` (Oct 7, merge of "Update manuscripts and Lean formalizations").
- **lean-base:** <https://github.com/CoolRmal/quasi-riemann-seven-eighths> @ `6ff3bd5`.
- **QRH Bounds tracker:** site <https://gtimg.github.io/QuasiRiemannTracker/>, repo <https://github.com/GTimG/QuasiRiemannTracker>.
  - The verification policy is `comparator-and-nanoda-v1`.
  - The Palomar kernel-dossier protocol (`PalomarSubmission` `d4e41c1`) was re-implemented for our three-kernel run.

## Lean toolchain and dependency pins

| Component | Pin |
|---|---|
| Lean | `leanprover/lean4:v4.34.1` |
| Mathlib | <https://github.com/leanprover-community/mathlib4> @ `d13f23b723b8a846827a245b89c10fc7d3f11612` |
| PrimeNumberTheoremAnd | <https://github.com/AlexKontorovich/PrimeNumberTheoremAnd> @ `c39a751132c88b6e8080b74c74023fd95b3d8be0`, plus the Lean 4.34.1 patch |
| rellich-kondrachov | <https://github.com/abenenson/rellich-kondrachov> @ `70f85d4c1bf99c6e7d61e8be4daa6f3664d08d23`, plus the Lean 4.34.1 patch |
| LeanArchitect | <https://github.com/hanwenzhu/LeanArchitect> @ `468e8f58` |
| checkdecls | <https://github.com/PatrickMassot/checkdecls> @ `3d425859` |
| leancert | <https://github.com/alerad/leancert> @ `7f91b6eb` |
| PrimeCert | <https://github.com/b-mehta/PrimeCert> @ `0803c2f6` |

Patches for the two patched dependencies are in `05-lean-23-24/lean/patches/`; the same files are in openai/math `lean/patches/`.

## Verification tools

| Tool | Pin / notes |
|---|---|
| Comparator | <https://github.com/leanprover/comparator> @ `d03acab154d269c06e60e4de7e4cc85deebff94b`, built locally with its `lean-toolchain` bumped to v4.34.1 |
| lean4export | <https://github.com/leanprover/lean4export> @ `076e8e57707e813375e8f9da8bf989799ace9680` |
| Landrun | Linux sandbox, <https://github.com/Zouuup/landrun> @ `811cfff5`. On macOS it was replaced by `05-lean-23-24/lean/scripts/landrun-shim`, which runs **unsandboxed**. |
| Three-kernel judge | `lake comparator` from the official Lean **v4.35.0-rc2** darwin_aarch64 release, with `nanoda_bin` (NanoDa) and `con-ron` in `bin/`. Driver: `08-…/lean/kernels/judge_three_kernels.py`. |
| tectonic | TeX engine for all papers |
| Python | 3.13. Most certificates are stdlib-only. Numerics need numpy ≥ 2, scipy, sympy and mpmath; see `requirements.txt`. |

## Data snapshots kept in this repo

- `docs/data/tracker_catalogue_main_20261009.json`: tracker `catalogue/results.json` on main (5 records).
- `docs/data/public_catalogue_snapshot_20261009.json`, `public_search_final_20261009.json`, `public_search_refresh_20261009.json`, `public_benchmarks_20261009.json`, `public_search_formalization_20261009.{json,txt}`: curated prior-work searches from the GPT-6 Astra session.

The raw GitHub, arXiv and Zenodo API dumps behind these snapshots were **not** kept. They contained third-party e-mail addresses and Li Xiang commit metadata.
