# 04 — The Gram / theta split: removing the P_a^{1/6} loss (Oct 9, 2026, 00:37–12:14)

**Goal.** Beat Liu's lemma-set optimum θ = (1507 − 2√921)/1653 ≈ 0.874957069799 by removing the **P_a^{1/6} Gram loss** in the low side of the 7/8 paper (Prop. `probe-gram`).

**Net status: no improvement on 0.874957 along this route.**
- The structural split survives.
- The claimed improvements to **0.8665 / 0.871 / 0.864 were retracted**.
- An improvement needs a power saving s₀ > 0 for a pairing (V′). (V′) is itself a Möbius-randomness, i.e. zero-free-region, statement in disguise.

## What survives

- **The theta split.** The P_a^{1/6} excess is exactly the energy of one rank-one "theta vector" A^sig = Y^{−1/3}Ŵ₁(1/6+iv)ρ_c, where ρ_c are sextic theta coefficients.
  - **Corollary D:** ‖A − A^sig‖² ≪ (Q/Y)Z^ε for all Y. This is proved modulo three unwritten lemmas of known type (`wf/FE_LoA6_proof.md`).
- **Lemma Z:** the reflected row carries the scalar ζ_m = κ₀ψ₁(m)·conj(γ₂(m))⁴. It is checked numerically to 10⁻¹⁴–10⁻¹⁵.
- **Further proved items:**
  - generalized Lemma Z (a trade-off theorem: zero detection and ζ-free rows are mutually exclusive);
  - signal non-cancellation;
  - the norm-only obstruction;
  - the 5/6 barrier for the residual-Cauchy–Schwarz template;
  - Lemmas C0, R and P;
  - zero-slot counts R₀(δ) = (15−16δ)/(15−6δ);
  - Weyl-groupoid computations;
  - the Coxeter (2,6,6) obstruction.
- **Bookkeeping-only exact certificates** (they certify "input ⇒ θ", not the input):
  - 1733/2000, 3055 boxes;
  - 108/125 with floor 1/200, 79,327 boxes;
  - 871/1000 slotted;
  - an independent pure-rational re-check (`audits/check_conditional_slotfree.py`).

## Why the headline claim was retracted (F1)

The Voronoi / large-sieve analysis bounded Σρ(m)(μ/m)₂. But the actual pairing is against **π = ρ·ζ_m**, with the reflection scalar of Lemma Z.
- The `wf/` audit found this (`wf/AUDIT_chain.md`).
- Three referees and `audits/gram_split_audit.txt` confirmed it.

## Conditional table

(V′) with saving s₀ in every block gives:

| s₀ | 0 | 0.002 | 0.005 | 0.01 | 0.02 | 0.04 |
|---|---|---|---|---|---|---|
| slot-free | 0.8947 | 0.8929 | 0.8902 | 0.8858 | 0.8775 | 0.8635 |
| slotted | 0.874957 | 0.8733 | 0.8709 | 0.8668 | 0.8615 | 0.8532 |

Other conditional values:
- The random-rate (V′_LS) would give 0.851918.
- A large-sieve-shaped saving gives 123/142.

## Five refereed attacks on (V′), all s₀ = 0 (`vls/`)

| Attack | Approach |
|---|---|
| A | Functional-equation geometry and Weyl groupoids |
| B | Gauss-sum and bilinear structure |
| C | Multiple-Dirichlet-series moments and spectral tools (Coxeter (2,6,6) is hyperbolic) |
| D | Probe redesign: 50 variants. Zero detection excludes shortening. |
| E | Numerics to Q = 9.6·10⁴: (V′_LS) sits at the random rate, but numerics cannot address provability |

## Layout

| Path | What |
|---|---|
| `REPORT.md` | Master status with correction banners; read the banners first |
| `note/theta_split_note.{tex,pdf}` | The final write-up |
| `TARGET_VLS.md` | Precise (V′), (V′_LS), (V′_β) targets |
| `D_derivation.md` | Lemma D1 and K₀ |
| `LIT_ideas234.md` | Literature scan for Ideas 2–4 |
| `IDEA1_NOTES.md`, `V_derivation.md`, `LoA6.md` | **Superseded** working notes; they contain retracted claims |
| `wf/` | Audit workflow: `AUDIT_chain.md` (finds F1), `FE_LoA6_proof.md`, `PARTII_slotfree.md`, `V_gaps.md`, plus checks and LPs |
| `wf2/` | The A₂ WMDS lead: conductor N(μ)³, no saving |
| `vls/` | Attacks A–E with referee reports and scripts. `vls/data/` keeps only json/log results. |
| `D_numerics/` | Lemma D1 numerics |
| `audits/gram_split_audit.txt`, `audits/check_conditional_slotfree.py` | Independent audit, from the later GPT-6 Astra session |
| `referee-scratch/` | Referee scripts that lived only in the session scratchpad: `p2/` (slot-free LP checks), `refB/`, `fc/`, `n2check.py`, `n6check.py` |
| `*.py` at top level | Libraries (`zomega.py`, `cubic_theta.py`), LP models (`theta_split_model.py`, `model_r4.py` = `../03-six-sevenths/routes/r4/model.py`), the certificate `exact_cert_theta.py`, and experiments `exp_*.py` |

**Known errata inside working files.** These come from referee reports; the full list is in [`../docs/RESULTS.md`](../docs/RESULTS.md) §F.
- `REPORT.md`'s "t = 1 → 0.8751" should be 0.877911, and "0.86624" should be 123/142.
- `exact_cert_theta.py` is not stdlib-only: it imports numpy/scipy to pick the witness t.
- `B_bilinear.md` overclaims (V′_β) ⇔ (V′_LS); only ⇐ holds.

## Not included (regenerable, about 190 MB)

`*.npy` / `*.npz` caches, which can be rebuilt:
- `vls/data/primes_X4194304.npz`: `cd vls && python E_primes.py` (≈ 48 min, 11 cores)
- `slist_*`, `A*_E21`: `python E_build.py 21 6000 12000 24000 48000 96000`
- `ls_*.npz`: `python E_ls.py Q E 4 1`
- `D_numerics/data`: `python expA.py 12 1500 3000 6000`
- `wf2/*.npz`: `python check_zeta.py Q`, `python test_twist.py Q`
- top-level `A_Q*` / `m_Q*`: `python exp_rank1.py 3000 256 9`
- `cubic_gauss_cache_400000.npz`: `python exp_cubic_scaling.py 400000`

## Reproduce

Scripts hard-code `/Users/akashlevy/Documents/riemann/gram_split` in `sys.path`, `os.chdir` and `np.load`. Relocate them first:

```bash
grep -rl '/Users/akashlevy/Documents/riemann/gram_split' --include='*.py' . | xargs sed -i.bak "s#/Users/akashlevy/Documents/riemann/gram_split#$PWD#g"
python theta_split_model.py
python wf/lp_const_saving.py      # the s0 table
python exact_cert_theta.py theta=1 dropq=1 lx=989/5000 ly=3141/5000 ell=0 sigma0=1733/2000   # ~30 s
python audits/check_conditional_slotfree.py
```

Requirements: Python 3.13, numpy ≥ 2, scipy, sympy, and pdflatex for the note.
