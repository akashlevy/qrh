# 03 — Toward θ = 6/7: conditional theorems, certificates and obstructions (Oct 8, 2026)

**Bottom line.** **No unconditional 6/7 proof** was found, so nothing was formalized. Every route to 6/7 needs a mean-value estimate *below the large-sieve diagonal*, and §3 of `SIX_SEVENTHS.md` explains why existing machinery cannot supply one.

How the work was done:
- Five independent research routes (R1–R5) ran in a multi-agent workflow.
- Positive claims were adversarially refereed (`routes/ref_R4*`).
- A focused workflow then checked and refuted the one surviving lead, the quartic transplant (`quartic-transplant/`).

## Results

**Conditional theorems**
- **Theorem C1** (complete conditional proof, `SIX_SEVENTHS.md` §1):
  - Hypothesis M_{6/7}: ∫_D^{2D} Σ_{N(u)≤D^{6/7}} |A_u(y)|² dy/y ≪ D^{13/7+ε}.
  - Conclusion: every Hecke L over Q(√−3), every Dirichlet L and ζ are zero-free on Re s > 6/7.
  - The hypothesis is *strictly stronger* than the conclusion: it goes beyond the density hypothesis for the sextic family.
- **Theorem C:** M_κ ⇒ β\* ≤ (6+5κ)/12 for every κ ∈ (0, 1]. This gives 6/7 at κ = 6/7, 3/4 at κ = 3/5, and 1/2 as κ → 0. It uses the new **Lemma B**, a bootstrapped prime extraction (`routes/R2_bootstrap.md`).

**Exact conditional certificates inside Part II** (`routes/r4/exact_cert.py`)

| Hypothesis | Result | Notes |
|---|---|---|
| W(4/7) | 6/7 | 109 boxes |
| **refined count + W(3/10)** | 6/7 | 17321 boxes, worst −8.3·10⁻⁹ — the robust headline |
| refined count + W(2/7) | 6/7, low side 0.857094 | Passes only with ZETA = 10⁻⁶ (80235 boxes, worst −4.7·10⁻⁹). It **fails** at the default ZETA = 10⁻⁴ (`routes/r4/certW27_ref.txt`), so treat it as fragile. |
| GRH rows + Gram bound without P_a^{1/6} | 6/7 | |
| GRH rows + W(1/7) + refined count | 6/7 | |
| Vertex B, vertex C | 6/7 + 1/1000 | |
| Lossless one-slot marked moment (C4) | 6/7 | |

**New proved lemmas**
- **Refined powerful-row count.** The case-A loss in Lemma `reflected-energy` improves from (1+3ℓ′−2M′)/4 to (1+2ℓ′−2M′)/3. Elementary; alone it does not move θ.
- **Detector floor.** 51/100 → 1/2 + η is legitimate; alone it does not move θ.
- **Lemma B** (bootstrapped prime extraction), as above.
- **Best provable density exponents for sextic rows:** A(6/7) = 7/3, A(7/8) = 16/7, A(11/12) = 2.

## Obstructions (proved; see `SIX_SEVENTHS.md` §3 and the route notes)

- **Paper2's dual target is false below the diagonal (Prop O1).** This is GRH-conditional in general and unconditional for κ < 1/4.
- **Move calculus.** |ε − δ| = 1 − κ is invariant under the available duality moves.
- **Coefficient-blind impossibility (R3, Theorem B).** Any large-values bound valid for every |c_n| ≤ 1 must allow ≥ H^{1/6} large rows. This rules out the large sieve, Halász–Montgomery, Huxley and Guth–Maynard. Theorem A: the multiplicity-beating large-values (LV) hypothesis is *equivalent* to the conclusion.
- **13/15 floor.** It holds for any row-count input, GRH for rows included (Theorem D).
- **Bootstrap consistency barrier (R2).** A counterfactual with one zero at β ≤ 13/15 is consistent with every lemma conclusion of both papers.
- **Prime-sixth-power extraction floor.** max((6+5κ)/12, 1 − κ/6) ≥ 6/7.
- **Quartic transplant refuted** (§5, `quartic-transplant/`):
  - the row twist enters as h³, so the reflection conductor is N(h)³;
  - the dual twist stays quartic, and there is no polar term;
  - the transplanted Part II therefore gives **θ = 1** (no half-plane);
  - R5's heuristic 0.854 came entirely from a wrong N(h)² conductor.

## Retracted or annotated

- `routes/R4_partII_map.md`: the claim "θ = 6/7 ⟺ the argument closes with h = 6/7" was **refuted** by the referee (`routes/ref_R4b/fixed_h_check.py`). Only θ ≥ 1 − h/6 holds.
- `routes/R5_literature.md`: its quartic lead was later **refuted**.
- `R4_partII_map.md` §0 item 5 speculates about a third party's exact 6/7; it contains none of their work.

## Layout

| Path | What |
|---|---|
| `SIX_SEVENTHS.md` | Summary: C1, R4 table, obstructions, quartic refutation, Lean status |
| `routes/R1_direct_meansquare.md` | R1: prove M_κ directly. Prop O1, move calculus, MDS signatures, C1–C4 |
| `routes/R2_bootstrap.md` | R2: bootstrap. Lemma B, Theorem C, BMS, the W_β counterfactual |
| `routes/R3_zero_density.md` | R3: zero density. Theorems A, B, D, D′; density exponents |
| `routes/R4_partII_map.md` | R4: every lemma replacement in Part II; exact certificates |
| `routes/R5_literature.md` | R5: literature survey (de Faveri, dFDH, BGL, Baier–Young, Gao–Zhao, …) |
| `routes/*.py` | `exact_certificates.py`, `r3_verify.py`, `wbeta_consistency.py` and `r1_bookkeeping.py` are exact stdlib scripts; `r1_fake_diagonal.py` and others are numerics |
| `routes/r4/` | Canonical Part II model `model.py`, generalized exact certificate `exact_cert.py`, energy LPs and dual certificates, GRH duals, certificate logs `cert*.txt` |
| `routes/ref_R4/`, `routes/ref_R4b/` | Referee checks. `chk_energy.py` and `fixed_h_check.py` point at the original scratchpad path `wf67/r4`; repoint them to `../r4`. |
| `routes/r5/` | Quartic heuristic from R5. **Superseded:** the lead was refuted. |
| `quartic-transplant/` | Rescued from the session scratchpad: Gaussian-integer numerics (`shift_check.py` h³ spread 10⁻¹³ against 2.00 for h; `polar*.py`; `meansquare_test.py`; C helpers `q*.c`) and the derived-geometry model `model/` (A: dF → 1.0; P2 route → 89/88, 65/64). There are no markdown notes; the narrative is `SIX_SEVENTHS.md` §5. |
| `lp_partII_relax.py` | GRH-rows LP with toggles |

Third-party paper text dumps are omitted; `lit/` held arXiv 0804.2233, 1112.1650, 2212.09241, 2607.07911, 2610.04045 and Liu's preprint. See [../docs/PRIOR_WORK.md](../docs/PRIOR_WORK.md).

## Reproduce

Requirements: Python 3; `numpy scipy sympy mpmath` for non-stdlib scripts. Some numerics scripts hard-code `/Users/akashlevy/Documents/riemann/strengthening` on `sys.path`; point it at `../02-strengthening`.

```bash
cd routes && python3 exact_certificates.py && python3 r3_verify.py && python3 wbeta_consistency.py && python3 r1_bookkeeping.py
cd r4 && python exact_cert.py                                                   # W(4/7) vertex A, 109 boxes
python exact_cert.py lx=51/125 ly=11/20 ell=171/500 sigma0=6/7 wE=3/10 caseA=refined          # ~1-3 min
python exact_cert.py lx=11489/28328 ly=14840/27149 ell=32673/97961 sigma0=6/7 wE=2/7 caseA=refined ZETA=1/1000000   # ~13 min
python3 grh_cert_verify.py && python energy_dual_cert.py
```

Quartic data must be rebuilt; the 36 MB and 8 MB data files were not kept:

```bash
cc -O2 -o qprimes qprimes.c -lm
cc -O2 -o qsqf2 qsqf2.c -lm
./qprimes 2000000 ...
./qsqf2 2000000 primes_2e6.txt sqf2_2e6.bin
```
