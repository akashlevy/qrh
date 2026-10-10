# 01 — Re-optimizing OpenAI's 7/8 lemma set, and where it stops (Oct 8, 2026)

**Question.** Can OpenAI's 7/8 zero-free half-plane be pushed lower using only that paper's own lemmas? If so, where does the lemma set provably stop?

## Findings

| Result | Status | File |
|---|---|---|
| **θ = 7/8 − 1/40000 = 34999/40000** for all Q(√−3) Hecke L-functions, all Dirichlet L-functions and ζ. Only the fixed geometry changes: b = 1/8, ℓ = 1/6 + 10⁻⁴, M + ℓ = 1 (l_x = 21247/60000, l_y = 28747/60000, h = 16253/20000). | Exact stdlib certificate, 24/24 checks, including a new Bernstein endpoint certificate (35 boxes, margin 9·10⁻⁵). Paper-level only; not Lean-checked, since 19 Lean files hard-code 17/48, 23/48, 13/16. | `certificate.py` |
| The checker is discriminating: it re-certifies the paper's margin 49/440640, fails for margin 3·10⁻⁴, and fails for ℓ = 1/6 + 1/5000. | Re-run Oct 9 | `tests/test_certificate_discriminates.py` |
| **Lemma-set floor:** no geometry gives θ < 33216105/37963129 ≈ 0.8749570. | Exact LP-dual certificate | `barriers.py` (A) |
| **GRH-for-rows floor 13/15:** even with GRH for every auxiliary row L-function, Part II gives θ ≥ 13/15. The multipliers are 4/5 (low) and 1/5 (high), and equality holds at l_x = l_y = 2/5, ℓ = 1/5. | Exact LP-dual certificate | `barriers.py` (B) |
| Single-ingredient sensitivity: improving counting or detection is worth ≤ 0.008. The bottleneck is the low side (Gram/energy). | Float model | `sensitivity.py` |

**Context.**
- Baiying Liu independently posted the same rational 34999/40000 the same day, together with the algebraic optimum (1507 − 2√921)/1653 ≈ 0.874957069799 ([arXiv:2610.12234](https://arxiv.org/abs/2610.12234)).
- Liu's optimum lies 8·10⁻⁸ above our LP floor. So the true optimum of this lemma set is in [0.8749569879, 0.874957069799], and 7/8 was within about 4.3·10⁻⁵ of it.

`QRH-analysis.md` is the full write-up:
- how the 7/8 proof works;
- the exponent bookkeeping table;
- the certificate;
- the LP duality argument;
- the sensitivity table;
- §6, the barriers to θ = 1/2: energy wall, sextic large sieve conjecture, mollifier-length consistency argument, and the 11/12 companion's κ wall.

Its §7 mentions a third party's unpublished lower θ. That passage is speculation about which lemma they might have changed; it contains none of their work.

## Reproduce

```bash
python3 certificate.py                          # ALL CHECKS PASS (0.04 s, stdlib)
python3 barriers.py                             # 33216105/37963129 and 13/15, "both certificates verified exactly"
python3 tests/test_certificate_discriminates.py # True / [False] / [False]
python3 sensitivity.py                          # ~45 s, needs numpy + scipy
```
