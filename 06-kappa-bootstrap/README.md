# 06 — κ-feedback bootstrap to 43747851/50000000 (Oct 9, 2026) — superseded, not novel

> **Status: historical.** This route is superseded by `../08-weighted-numerator-10499-12000` (θ = 10499/12000). It is **not new**:
> - stage 0 is Liu's rational 34999/40000;
> - stage 2 is exactly Jizhou Guo's published weaker rational 43747851/50000000;
> - iterating the feedback converges to the public cubic-root endpoint θ\* ≈ 0.874957019420.
>
> Shunia / GPT-6 Astra, Guo, Nielstron and Li Xiang had already published that endpoint.
>
> The κ-feedback mechanism and its limiting cubic were already described publicly by Li Xiang, and the note cites him. Per this project's standing rule, **none of Li Xiang's work was copied**. Only a citation of his public repository appears in `paper/kappa_bootstrap_proof.tex` (`\bibitem{lixiang}`).

## What it shows

1. **`lem:plain` holds for 13/18 ≤ κ ≤ 1**, not only 3/4 ≤ κ (`audits/plain_moment_extension.txt`). The proof uses only κz ≤ M/6 and 6κ − 1 ≥ 10/3.
2. **Feedback.** A proved bound B gives κ = 2B − 1, then capacity c = 1/(6κ) > 2/9, then a better row count; re-optimize the geometry and re-certify:

| Stage | Bound | κ used | Certificate |
|---|---|---|---|
| 0 | 34999/40000 | 3/4 | `../01-seven-eighths-lemma-set/certificate.py` |
| 1 | 17499141/20000000 = 0.87495705 | 14999/20000 | 24 Bernstein leaves, margin 10⁻⁸ |
| 2 | **43747851/50000000 = 0.87495702** | 7499141/10000000 | 26 leaves, margin 10⁻⁹ |

3. **Lean (`lean/`).** Only the energy-induction extension is formalized.
   - It is a 63-module overlay (6,861 lines) of OpenAI's `Energy/*` modules on fork prune commit `817ddb5e6`, with κ floor 13/18.
   - Endpoints `QRHResearch.KappaThirteenEighteenths.certified_bands` and `terminal_certificate` use only the standard axioms (`lean/logs/endpoint-axioms.log`).
   - **The improved zero-free theorem itself was never formalized.**
   - The overlay files are modified copies of openai/math modules (Apache-2.0) without per-file "adapted" headers; see [../NOTICE](../NOTICE).

The barrier for this route, that κ-feedback can never beat the cubic root θ\*, is in `../07-method-barriers`.

## Layout

| Path | What |
|---|---|
| `paper/kappa_bootstrap_proof.{tex,pdf}` | 11-page note. Its author line reads "prepared with OpenAI Codex", and it cites unpinned `openai/math/blob/main` URLs. |
| `certificates/bootstrap_certificate.py` | Exact stdlib certificate; `--stage all [--negative-control]` |
| `certificates/bootstrap_*.json` | Stage outputs |
| `audits/bootstrap_referee.txt` | AI referee report |
| `audits/plain_moment_extension.txt` | The κ ≥ 13/18 lemma |
| `audits/verification_manifest.json` | Its novelty field records the overlap with prior public work |
| `lean/` | `KappaFormal/` overlay sources, `KappaThirteenEighteenths.lean`, build/verify scripts, provenance, logs. Oleans are not kept. The scripts reference absolute paths of the original workspace. |

## Reproduce

```bash
python3 certificates/bootstrap_certificate.py --stage all
```

The Lean overlay needs the fork build (`../05-lean-23-24`); run `lean/build_overlay.py`, then `lean/verify_endpoint.py`.
