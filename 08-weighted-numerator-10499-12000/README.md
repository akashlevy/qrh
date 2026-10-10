# 08 — Weighted numerator moments: θ = 10499/12000 (Oct 9–10, 2026) — headline result

*GPT-6 Astra (Ultra) was used to discover the proof and perform partial formalization; Opus 5.5 Ultracode was used to complete the formalization and submit the result.*

## Theorem

Let β\* be the supremum of 1/2 and the real parts of the nontrivial zeros of all primitive finite-order Hecke L-functions over Q(√−3) = Q(ζ₃). Then

> **β\* ≤ 10499/12000 = 7/8 − 1/12000 = 0.874916666…**

Hence on the open half-plane Re s > 10499/12000 none of the following vanishes (the trivial character's pole at s = 1 is excluded):
- ζ;
- every Dirichlet L-function of every modulus;
- every finite-order Hecke L-function of Q(ζ₃) in OpenAI's `HeckeFamily` definition.

This is not GRH, and nothing is claimed for other fields.

The bound is about **4.04·10⁻⁵ below** the public cubic-root endpoint 0.874957019420… and below the tracker's best listed value. The Lean lemma `WeightedQRH.theta_lt_public` in `Parameters.lean` proves the comparison.

## Lean endpoints (only `propext`, `Classical.choice`, `Quot.sound`)

`lean/WeightedQRH/Final.lean`:
- `WeightedQRH.beta_le_theta : HeckeZeroSupremum.beta ≤ 10499/12000`
- `WeightedQRH.hecke_nonzero`
- `WeightedQRH.dirichlet_nonzero`
- `WeightedQRH.zeta_nonzero`

`lean/WeightedQRH/Statements.lean`, the Comparator form:
- `OAI.riemannZeta_ne_zero_of_10499_12000_lt_re`
- `OAI.DirichletCharacter.LFunction_ne_zero_of_10499_12000_lt_re`
- `OAI.SevenEighths.HeckeFamily.LFunction_ne_zero_of_10499_12000_lt_re`

The last step, `fullW_analytic_input := fullW_input_of_raw_source Numerator.actual_raw_fullW_input`, discharges the only analytic input of `HighRows/ConditionalFinal.lean`. The full chain is in `theorem-maps/final_assembly_theorem_map.txt`.

## The new idea

OpenAI's 7/8 high side bounds the numerator L(w, χ_u) **pointwise**, and uses fourth moments only to *count* zero-witness rows. The weighted argument changes this in three steps.

1. **Integrate the whole w-line.** Mellin inversion reconstructs the actual smoothed plain sum S_χ(Y/Nk), truncated exactly at Nk ≤ C·Y before reflection.
2. **Reflect it.** The Hecke functional equation produces a polynomial of length ≤ U/Y (U = Z^d). Since l_y < d < 2l_y, the reflected length is n₀ = 1 − l_y/d < 1/2.
   - The leftover capacity z_P = (2l_y/d − 1)/(6κ) in the original `lem:plain` carries selected prime factors.
   - Hölder with exponents (4, 4/3) on each amplitude class gives Σ_C|S| ≪ U^{(1+3R)/4 − q·z_P/2}.
3. **Handle row-dependent Euler-correction coefficients without putting them inside the moment.**
   - a weighted ℓ¹ coefficient norm at Re w = 1/2 − δ/4;
   - a quotient-free expansion E_p = G_p + χ̄_p(u)H_p;
   - a single conductor-deficit allocation for strict ramified labels;
   - dilation tails, paid because ρ − q/(6κ) ≥ 5δ/36.

New high exponent:

E_W(d) = a − θ + h(z₀ − 1/6) − l_y/2 − ℓ/2 + qℓ + d(3R/4 + 1/4 − z₀) − q(2l_y − d)/(12κ).

The old contour is used for δ ≤ 1/3, the weighted one for δ ≥ 1/3.

Parameters:
- **κ = 3/4**, unchanged; there is **no κ-feedback**;
- **geometry:** ℓ = 167/1000, b = 43/500, l_x = 747/2000, l_y = 919/2000;
- **low side:** θ_low = 11/12 − ℓ/4, which needs ℓ < 1/5;
- **Euler regions:** extended below 7/8.

Exact certificate: E ≤ −10⁻⁴ on 1 + 9 Bernstein rectangles plus 32 rational checks (`certificates/weighted_numerator_certificate.py`, stdlib).

## Verification

All of these are **author-side** runs. **No external human peer review** has taken place.

- **Source scan.** 231 files: no `sorry`, `admit`, `axiom`, `native_decide`, `extern`, … (`evidence/source-scan.txt`).
- **Isolated rebuild.** 231 modules, all exit 0 (`evidence/isolated-rebuild/`).
- **Lake build from public sources.** 11,847 jobs (`evidence/lake/`). The build directory was seeded, so the fetch/patch/cache steps of `setup.sh` were not exercised.
- **Comparator** (`d03acab1` + lean4export `076e8e57`, Lean 4.34.1) accepted all 3 challenges. A negative control (ζ at 1/2) was rejected. Landrun was replaced by an **unsandboxed** shim (`evidence/comparator/`).
- **Three-kernel judge** (`lean/kernels/judge_three_kernels.py`, re-implementing the Palomar protocol `d4e41c1`):
  - **con-ron** "accepted 123717 declarations (--verified)"; **NanoDa** and the **Lean** kernel also accepted.
  - **PASS** at 2026-10-10T01:53:26Z.
  - It used the tracker's three-target template with only the constant swapped, compiled against pristine openai/math `fd4aeeb2e`.
  - It ran **without a sandbox** on macOS (`evidence/kernels/`).

AI-assisted referee rounds are in `audits/`:
- `weighted_holder_referee.txt` found a real gap: a rowwise coefficient-mass bound cannot be swapped with the moment. It was fixed with the ℓ¹ supremum method.
- `weighted_numerator_referee.txt` required truncating before reflection.
- The assembly, complete and global referees found no missing hypothesis.
- `weighted_numerator_reflection_audit.txt` withdrew an earlier, wrong "zero-slot gap" claim.
- The final audits fixed 13 stale passages and the mislabeling of fork commit `c388e88e6` as an OpenAI commit.

## Submission

[GTimG/QuasiRiemannTracker#6](https://github.com/GTimG/QuasiRiemannTracker/pull/6) (open).
- Package id `akashlevy-20261009-weighted-numerator`, with proposed status `framework-verified` for maintainer review.
- Head: `akashlevy/QuasiRiemannTracker@2fd60c0`.
- Package commit `ec626ec`.
- The fork base is pinned by tag `qrh-weighted-10499-12000-base` → `akashlevy/math@c388e88e6` (openai/math `fd4aeeb2e` pruned in place; see `../05-lean-23-24`). **This is not an OpenAI commit.**

## Layout

| Path | What |
|---|---|
| `paper/weighted_numerator_proof.{tex,pdf}`, `weighted_numerator_lemma.tex` | The 18-page manuscript as submitted (sha256 of the PDF: 228cc5b0…) |
| `lean/` | The submitted formalization: 228-module `WeightedQRH/` closure, `ComparatorChallenges/`, `setup.sh`, `run_comparator.sh`, `kernels/` |
| `lean-extras/WeightedQRH/` | 24 audit and benchmark modules outside the shipped closure (`*Axioms`, `*Audit`, `CertificateBenchmark`, `FullWAnalyticAudit`, `HighRows/MixedCollected`) |
| `certificates/` | The exact certificate (`.py` + `.json`); `central_numerator_search.py`, an exploratory parameter search (numpy) |
| `evidence/` | Kernels, Comparator, Lake and isolated-rebuild logs, plus historical rebuild JSONs |
| `audits/` | Referee reports, independent derivations (`independent_count.txt`, `independent_weighted_numerator_audit.txt`), assumption/quantifier audits, final citation/statement audits |
| `theorem-maps/` | Manuscript ↔ Lean maps for each stage, port manifests, dev scripts (`build.py`, `generate_certificates.py`, `port_moments.py`) |
| `submission/` | Package README, NOTICE, LICENSE, PUBLICATION, transformation records, sha256 manifest, PR body, the public-download metadata, `tracker_requirements.md` (notes on the tracker's kernel pipeline), `check_upstream.py` + `lhash.py` (Lake trace-hash check against `c388e88e6`) |
| `history/` | `COMPLETION_20261009.txt` (completion record), `HANDOFF_ROOT_20261009.txt` and three agent handoffs from the stopped GPT-6 Astra session, the pre-submission verification manifest |

**Stale records to be aware of.**
- `history/weighted_numerator_verification_manifest.json` lists an earlier 16-page PDF hash. The authoritative hashes are in `submission/akashlevy-20261009-weighted-numerator.sha256.json`.
- `COMPLETION` §7 says NanoDa/con-ron "not run"; its submission addendum supersedes this.

## Reproduce

```bash
python3 certificates/weighted_numerator_certificate.py --output /tmp/cert.json   # stdlib; PASS
cd lean && bash setup.sh       # clones akashlevy/math@c388e88e6, fetches deps, applies patches, Mathlib cache, builds
export COMPARATOR=... COMPARATOR_LEAN4EXPORT=... COMPARATOR_LANDRUN=...
for c in Zeta Dirichlet Hecke; do bash run_comparator.sh $c; done
python3 kernels/judge_three_kernels.py --solution-workspace . --challenge-workspace <pristine fd4aeeb2e workspace> \
  --template <tracker>/public/proofs/palomar-20261009/qrh/src/Challenge.lean \
  --source-lean <lean-4.34.1> --judge-lean <lean-4.35.0-rc2> --lean4export <lean4export> --out <fresh>
```

Requirements:
- Lean 4.34.1 and Mathlib `d13f23b7`;
- the official Lean v4.35.0-rc2 release (with `nanoda_bin` and `con-ron`) for the three-kernel check;
- about 9 GB RAM and 2.2 GB of disk for the exports.

Use a Linux sandbox for a trustworthy rerun.

## Open directions

- **The weighted route alone appears saturated near θ ≈ 0.87489.** This comes from a coarse-grid pure-Python re-run of `central_numerator_search.py --weighted` (κ = 3/4, M = 1 − ℓ family), done while organizing this repo. It is *not* part of the submitted work. It found the margin −1.3·10⁻⁴ at ℓ = 0.167 and +9.6·10⁻⁴ (failing) at ℓ = 0.168.
- **The low side still caps everything at 13/15** (`../07-method-barriers`).

See [../docs/OPEN_LEADS.md](../docs/OPEN_LEADS.md).
