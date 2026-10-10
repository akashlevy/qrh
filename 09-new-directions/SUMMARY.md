# 09 — Session summary: pushing QRH below 7/8 (10 Oct 2026)

**Status: no improvement yet.** Best framework-verified public bound: ≈ 0.874957. Ours is 10499/12000 ≈ 0.874917 (tracker PR #6).
The new pending record is 0.874887 (tracker PR #9, Hailey Collet), a re-tuning that reuses our weighted-numerator construction. This session mapped *why* every lever of OpenAI's 7/8 architecture is blocked, and found exactly one
genuinely new architecture, which turned out to be dead (c_B = 4). Details are in `NOTES.md` (§1–§31), with agent work under `agents/`.

## What a meaningful improvement would need (r4 exponent model, all numbers re-verified)

| lever | single-lever requirement for 27/31 | type |
|---|---|---|
| energy wall W(w) | lossless offset w ≥ 0.0638 (a relative saving needs η ≳ 0.85; slope f ≲ 0.10) | beyond the edge |
| Gram / theta pairing (V′) | absolute saving s ≥ 0.0049 at the edge (circular: T carries the whole signal) | at the edge |
| row counts | long-witness slope 5/6 → 0.746 (M_κ, κ = 0.895), i.e. "Ingham below σ = 7/10" for the sextic family | beyond DH (moment form) |
| DIAG (Lemma marked r+z ≤ m, Lemma plain cP = 1) | alone only 0.87255; with W(0.0196) or M_{0.964} it reaches 27/31 | needs a coefficient-sensitive sextic sieve; false for general slot coefficients |

**Conditional results:**
- (H1)+(H2), Lindelöf on average plus mark decorrelation for quadratic twists of the level-L² cubic theta ⇒ **19/22 = 0.8636**;
- plain sixth-moment Lindelöf for the sextic family ⇒ 0.8727, or 0.8694 with DIAG;
- perfect (DH) row counts ⇒ 0.8698.

Numerically, W(w) and (H2) hold with a margin of about Z^{0.43} (agent `edge-theta`): the obstruction is methodological, not truth.

## The structural chain (NOTES §22, §27–§29)

1. Zero detection needs 1/L, i.e. Möbius.
2. μ-absorption works only through Jacobi/Stickelberger triples μ·Hecke = γ_aγ_bγ_c (a+b+c ≡ 0 mod N).
3. A single-Gauss-sum theta row forces a quadratic ray factor. So after reflection the rows carry a QUADRATIC residual twist
   (proved for every order N, every field and every Kazhdan–Patterson theta). Quadratic-twist families are thin and their moment MDS is affine
   (C̃₂: the same class as GL2 quadratic-twist second moments, where only log savings are known).
4. Coherent prime amplification (marks) costs conductor P² at fixed family size. This is the energy wall, and the factor 2 in Lemma marked.
5. Double-Gauss-sum rows (cubic θ₃⊗θ₃, triple (1,1,1)) make the dual full-index (copies a³, exponent 1/3). But they are degree 4,
   with conductor exponent c_B = 4 (proved via Kubota's rank-one theory: the coprime part at index r has conductor Π q_p^{v_p(r) mod n}).
   That gives ≥ 10/11. The c_B = 2 orientation contains no μ.

Dichotomy: **degree-2 rows ⇒ thin family; full family ⇒ degree-4 rows.** Both lose against the family/conductor balance.

## Routes killed this session (with reasons)
- **κ-rescaled marks:** a reparametrisation.
- **m-twists and s-marks:** no net amplification.
- **Mellin-height averaging; L²-in-Z continuation:** gain T0^{1/2} against count growth T0^{R}, with R > 1/2.
- **Relative savings on the wall:** useless unless η ≳ 0.85.
- **OpenAI's prime-convolution beyond-sieve technique:** the phase-blind barrier (O1), and near-orthogonality excludes theta profiles only for n = 3.
- **Spectral routes:**
  - Rankin–Selberg / MDS: the affine group fixes x₁+x₂+y, so Bochner hulls cannot reach the needed region.
  - Shifted convolutions: there is no full character group.
  - Kuznetsov: loses the coherent sum.
- **Large-values methods with copies:** Guth–Maynard, Huxley, Halász–Montgomery, Heath-Brown, Jutila, Bourgain.
- **Theta-based large values; higher moments.**
- **Higher-rank Kazhdan–Patterson thetas:** the wall becomes (r−1)(M+ℓ) ≤ 1.
- **Trilinear probes; denser copies; Λ / log L detection; |θ|², θ³, Shintani zeta; other fields and orders.**
- **Within-edge slack:** the best is DIAG at 0.87255, itself unprovable. Every other item has zero payoff.

## Infrastructure
- `/home/user/leanbuild` builds the full OAI Lean library from source (Mathlib cache hosts are blocked here): 6817 jobs, OK.
- Exponent models:
  - `03-six-sevenths/routes/r4/model.py`, with variants in the `agents/*/model_*.py` files;
  - `agents/new-architectures/template_lp.py` (GRH template);
  - `agents/alt-detect/p78_cubic_lp.py` (cubic template).
- Z[ω] numerics: `agents/edge-theta/edge_lib.py`, and `agents/cB-numerics/cubic_data.c` (exact cubic Gauss sums up to norm 10⁷).

## The space of μ-absorbing probes is now fully mapped (NOTES §32)
| split of the Jacobi triple | dual | GRH floor |
|---|---|---|
| P78: sextic, 1 Poisson + θ₃ row + quadratic ray | thin, degree 2 | 13/15 (7/8 with the proved counts) |
| cubic, 1 Poisson + θ₃⊗θ₃ row | full, degree 4 (c_B = 4) | 10/11 |
| cubic, 2 Poisson + θ₃ row | full, degree 2 | 13/12 (each Poisson costs Z^{-1/6}) |
| cubic, 3 Poisson + Hecke row | full, degree 1 | 4/3 |
| sextic multi-Poisson splits | — | > 1 |
| higher-rank Kazhdan–Patterson thetas | degree r | wall (r−1)(M+ℓ) ≤ 1 |

So P78's split is the unique viable member. Literature (NOTES §33): beyond-edge Lindelöf on average is known only for full
character groups on powerful moduli or for full spectral families; μ-absorption forces squarefree, one-character-per-modulus
families. (H1) ⇔ a subconvex bound C^{0.202} on the cubic cover, where no subconvexity is known. W(w) is a Conrey–Rodgers "2-swap"
statement, uncontrolled for any thin family.

## Rounds 5–6: the 19/22 conditional theorem dissected (NOTES §34–§37)
- (H1) can be replaced by its mark-average, which follows from a theta-free hybrid large sieve (HLS) for {ψ_kχ_P^{−2}}.
- (HLS) is numerically true (random-matrix size), but unproved at its edge N ≈ KL. The Poisson involutions in n and in k
  share the fixed point N = KL, and a sextic cocycle blocks Heath-Brown-style recursions. Alone it gains nothing.
- (H2) cross-mark decorrelation is pure random-matrix cancellation (no main term, no structure across 54 phases). Alone it gains nothing.
- Both reduce to one "freeness / uncertainty principle": ‖Π_VΠ_X‖² ≲ K/N for the quadratic-row and mark subspaces.
  - A partial sieve KL + N·L^a with perfect (H2) gives 27/31 at a ≈ 0.687 and 19/22 at a ≤ 1/3.
- Mark design cannot dodge this: amplification and the conductor P² come from the same cubic symbol (lemma L2).

## Open leads, ranked
1. A lossless 2-swap or energy bound beyond the wall (W(0.064)) exploiting the theta coefficients' self-duality: the closest lever to
   existing technology, though still open. Mark-design freedom is under study (agent `mark-design`).
2. (H1)+(H2) for quadratic twists of the marked cubic theta (⇒ 19/22): true numerically with a big margin; needs a new spectral or
   amplification idea on the cubic cover.
3. A count mechanism beating Ingham-with-copies for sextic twists at σ ≈ 0.63–0.70. Ingham is unbeaten below 7/10 even for ζ, so this is deprioritised.
