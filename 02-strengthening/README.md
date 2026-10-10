# 02 — Strengthening the 7/8 argument: what actually moves θ (Oct 8, 2026)

**Question.** Which single lemma would have to be replaced or strengthened to move θ substantially? In particular, does relaxing the marked-energy wall W(w) reach 6/7?

`STRENGTHENING.md` is the audited write-up. A multi-agent workflow first produced optimistic claims; an adversarial audit then found three modeling errors in how W(w) entered the lemmas, and those claims were **retracted** (§1 there).

## Findings

| Result | Status |
|---|---|
| **W(1/5) ⇒ Re s > 108/125 = 0.864**, at (l_x, l_y, ℓ) = (91/200, 93/200, 7/25) | Exact conditional certificate, 2181 boxes; re-run Oct 9. `certify.py 1/5 91/200 93/200 7/25 108/125` |
| **W(w) low-side barrier:** low ≥ 13/15 − w/60 for every geometry, so 6/7 needs w ≥ 4/7 with the paper's case-A term | Exact LP (`energy_lp.py`) |
| **6/7 comes from M_κ at κ = 6/7:** the twisted-Möbius mean square below the large-sieve diagonal, Σ_{N(u)≤D^κ}\|A_u(D)\|² ≪ D^{1+κ+ε}. Prime extraction gives max((6+5κ)/12, 1 − κ/6), whose minimum is 6/7. Integer extraction gives (6+5κ)/12. | Proved implication. Developed into Theorem C1 in `../03-six-sevenths` |
| M_κ numerics: per-row energy / diagonal ≈ 0.91–0.98 at κ = 6/7 for D = 3000 and 10⁴ | Numerical evidence; consistent with GRH, so not evidence of provability. `mobius_meansquare.py` |
| W(w) numerics: case-B energies stay diagonal to 8× (unmarked), 128× (one mark) and 100× (slot sums) past the wall | Numerical evidence. `wall_experiment*.py` |
| Quartic μ-absorption identity μγ₋₁ = Φ·ᾱ·γ₁ verified on 538 primes | Numerical. `quartic_absorption_check.py` |
| Detector floor 51/100 → 1/2 + η is legitimate, but buys nothing alone | Audited |

**Retracted claims** (kept for the record; do not cite):
- W(1/10) ⇒ 0.865;
- W(1/5) ⇒ 0.858 or 6/7;
- W(17/100) with δ₀ = 1/1000 ⇒ 6/7;
- "one wall controls three lemmas".

The errors:
- (a) the case-A energy (1+2ℓ′−2M′)/3 assumed no powerful row part; the true value is (1+3ℓ′−2M′)/4;
- (b) W(w) does not relax `marked` past the canonical diagonal r + z < 1;
- (c) W(w) does not relax the amplification slope α = 5/6.

The certificate *arithmetic* was sound; the *model* was wrong.

## Files

| File | Status |
|---|---|
| `conditional_certificate.py` + `certify.py` | Corrected exact W(w) certificate. **Caveat:** the built-in defaults of `conditional_certificate.py` still encode the retracted W(1/10) ⇒ 0.865 claim and print `SOME CHECK FAILED` when run with no arguments. Always run through `certify.py` with explicit arguments. Its energy term also lacks the (M′+4ℓ′−2)/3 − w piece, which is inactive at the W(1/5) geometry; `../03-six-sevenths/routes/r4/exact_cert.py` includes it. |
| `audit/optimize_corrected.py`, `audit/run_w02_0864.json` | The auditor's corrected float model and the W(1/5) ⇒ 108/125 log, rescued from the session scratchpad |
| `energy_lp.py`, `canonical_lp.py`, `row_counts.py` | LPs for the reflected-energy and canonical-moment lemmas; generalized row counts |
| `mobius_meansquare.py`, `eisenstein.py`, `wall_experiment*.py` | ℤ[ω] numerics (numpy) |
| `quartic_absorption_check.py` | stdlib |
| **Superseded:** `conditional_certificate_SUPERSEDED.py`, `optimize_w.py`, `wall_relaxation.py`, `payoff.py`, `count_lemma_payoffs.py`, `order_n_engine.py`, `quartic_routes.py` | Pre-audit models encoding retracted claims. The ≈ 0.81 quartic theta route was later refuted in `../03-six-sevenths` §5. |

`STRENGTHENING.md` line 20 speculates about a third party's unpublished result; it contains none of their work.

## Reproduce

```bash
python3 certify.py 1/5 91/200 93/200 7/25 108/125   # "ALL_PASS": true, ~50 s, stdlib
python3 energy_lp.py 0.2 91/200 93/200 7/25          # scipy
python3 mobius_meansquare.py 3000                    # numpy
```
