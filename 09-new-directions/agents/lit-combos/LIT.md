# lit-combos — (A) literature scan for "beyond the edge" inputs (10 Oct 2026)

Network: arxiv.org, msp.org, github.io, personal sites, Crossref/OpenAlex are all denied by the egress
policy (checked; not retried). Sources available: WebSearch snippets (abstract-level), GitHub
(git/raw), and the project's earlier full-text reads (03-six-sevenths/routes/R5_literature.md,
04-gram-split/LIT_ideas234.md, marked [V] there). Labels below:
[V] statement read in full text earlier in this project (R5 / LIT_ideas234);
[A] abstract / search-snippet level only (this session);
[I] my inference; [M] checked with the exponent model (files in this directory).

Levers (see BRIEF): W(w) energy wall; (V′)/Gram exponent gram_b; M_κ (kA, w1); row counts R(δ,x).

## 1. New QRH results (status 10 Oct 2026, ~06:00 UTC)

* Tracker GTimG/QuasiRiemannTracker: PR refs are #1–#7 only (`git ls-remote ... 'refs/pull/*/head'`),
  so there is no PR beyond #7. `main` = 40844a7, catalogue as_of 2026-10-09, best entry
  nielstron 0.874957019420099 (framework-verified). Snapshot: scratch/tracker_main_20261010.json.
* openai/math: HEAD = main = fd4aeeb (the Oct 7 merge of PR #1, `codex/update-10-7` = 3014888, an
  ancestor). There are no preprints newer than Oct 7.
* Competitor repos re-checked (HEADs): Gao moved to 8c99552 and now claims 2187392879/2500000000 =
  0.8749571516 (conditional on P78 inputs). aster2024 (Guo) is at 867b0dc, cubic root
  0.8749570194. tomoto0's "0.87451 bootstrap fixed point" and "13/15 density barrier" are model
  computations it labels optimistic, not proofs. Argonaut, Liu, CrossNox, Shunia and Creech are
  unchanged.
* Web search ("quasi-Riemann hypothesis" Oct 2026): only Liu arXiv:2610.12234 (0.874957069799),
  Kintali 47/48, and press items.
* **No claimed θ below 0.8749 anywhere.**

## 2. Items, exact statements as far as verifiable, lever, expected gain

| # | item | statement (verification) | lever it could feed | expected gain |
|---|---|---|---|---|
| 1 | de Faveri, *Optimal large sieve for fixed order characters*, arXiv:2610.04045 (2 Oct 2026) | [V via R5] Thm 1.1: for n ≥ 3 and K ⊇ μ_n, Θ_n(A,B) := sup_{‖λ‖=1} Σ♭_{Na≤A} \|Σ♭_{Nb≤B} λ_b χ_a(b)\|² ≪ (AB)^ε (A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}), where ♭ means n-th-power-free. Cor 1.5 gives zero density for the n-th-power-free family. [A] The abstract says the bound is "expected to be optimal in all ranges", with applications to the cubic second moment by conductor and to moments and zero density for all n ≥ 3. | (i) Part I sextic row count: θ unchanged. (ii) Energy wall, unreflected route: with K = Z^{M′} rows and N = Z^{1+ℓ′} columns, the loss is max(1+ℓ′−M′, (5+5ℓ′−4M′)/6), worse than the reflected wall loss M′+ℓ′−1 for every M′ > 2/3 [I, computed]. (iii) Gram: the n = 6 term A^{5/6}B^{1/3} reproduces P_a^{1/6} exactly. (iv) Row counts: Cor 1.5 is weaker than P78's own counts. | **0** (all four) |
| 2 | de Faveri–Dunn–Hoffstein, arXiv:2607.07911 (Jul 2026) | [V via R5] Unconditional lower bounds Ξ_3 ≫ (AB)^{−ε}(A+B+(AB)^{2/3}) and Ξ_4 ≫ (AB)^{−ε}(A+B+A^{3/4}B^{1/2}+A^{1/2}B^{3/4}); conjecture (1.8) for general n. The extremiser is the Gauss-sum (theta) vector. The method is a Zagier-regularised Rankin–Selberg second moment of Kubota series. | It is an **obstruction**: the Gram P_a^{1/6} and the cubic factor (NL)^{2/3} in Lemma hybrid-energy are sharp for general coefficients. The Rankin–Selberg technique is the natural tool for (V′). | 0 directly. It is a technique template for gram_b. |
| 3 | David–de Faveri–Dunn–Stucky, arXiv:2410.03048 | [V via LIT_ideas234] Lemma 6.6: Σ_{N(m)≤M} \|ψ^{(3)}(hm,σ+it)\|² ≪ N(h)^{1/2}M^{1+ε}(1+\|t\|)² for 1<σ≤5/4 (Lindelöf on average, cubic Kubota series). Lemma 6.7 is the polar term ∝ g̃(a)N(a)^{−1/6}B^{5/6}. The mollifier is shorter than the family, so the regime is above the diagonal. | (V′) would need the sextic, quadratically twisted version (conductor ≈ L³, family L) at σ near 1/2, not σ > 1. | 0 |
| 4 | Castillo–de Faveri–Dunn, arXiv:2604.01316 (Apr 2026, 65 pp) | [A] Unconditional positive proportion (≫ 1/log) of non-vanishing for quartic Hecke L over Q(i) via mollified first and second moments; the same Rankin–Selberg Kubota input in the quartic case. DIPP 2607.27131 improves its error term. | Same regime as #3 (rows ≥ length) | 0 |
| 5 | David–Dunn–Hamieh–Lin, arXiv:2306.11875 (ANT, v5 Jan 2026) | [V via LIT_ideas234] (1.11): Σ_{N(a)≤A} μ²(a)\|ψ̃^{(4)}(½+ε+it, a²)\|² ≪ A^{3/2+ε}, with conductor ≈ A². This is the **only** published mean value with fewer rows than the conductor, and it holds at large-sieve strength (3/2). The authors say A^{1+ε} "seems likely". | (V′) needs the sextic analogue with exponent θ_π < 1.31 (04-gram-split; LS = 3/2, LoA = 1). A sextic analogue at 3/2 gives nothing. | 0 at 3/2. At θ_π = 1 it would give the full (V′), θ ≈ 0.852 per 04-gram-split [H]. The table in §3 shows the partial gram_b. |
| 6 | Dunn, *Metaplectic cusp forms and the large sieve*, ANT 19 (2025) 1823–1880, arXiv:2403.13151 v4 | [A] A new bilinear large sieve whose **kernel is the cusp-form coefficient sequence** (Browning–Vishe circle method, Voronoi, Gauss–Ramanujan sums), and a linear-form estimate past level 2/3 (Voronoi plus the HB cubic sieve). Consequence: power saving for Σ_p ρ_f(p). | It is a "special-coefficients beyond the generic sieve" template, but for cusp forms (not θ) on the cubic cover, and for first moments. None of our levers has this shape. | 0 now |
| 7 | OpenAI, *Unconditional first moment for cubic Gauss sums* (Sep 25), `sections/dual.tex` | [V via R5] Σ_{p,q}\|A(p,q)\|² ≪ Y^{7/3−η}, beyond the sharp cubic sieve, for prime-convolution coefficients; η is ineffective and small. | It could only hit the (NL)^{2/3} factor of Lemma hybrid-energy (marks are prime), i.e. the branch (M′+4ℓ′−2)/3 − w of the energy LP. [M] That branch is **never binding**: dropping it entirely leaves θ unchanged at w = 0, 0.06, 0.15 (runs/batch2). | 0 |
| 8 | Diaconu–Ion–Paşol–Popa, arXiv:2607.27131 (29 Jul 2026) | [A] First and second twisted moments of r-th order Hecke L-functions with power saving, over global fields **containing μ_{2r}**, r ≥ 3, by multiple Dirichlet series; positive-proportion non-vanishing for the squarefree and r-th-power-free families, with secondary terms. | It needs μ_12 for r = 6, which Q(ω) lacks. For r = 3 over Q(ω) it is a second moment at the edge (length ≈ family). No Möbius below the diagonal. | 0 |
| 9 | Sawin–Whitehead, arXiv:2507.08662 | [A] Functional equations of axiomatic MDS (function fields) via Weyl groupoids and quantum algebra; meromorphic continuation for many series. | The MDS signature of M_κ is affine (det 0), "beyond the edge" (03 R1). A function-field-only framework gives no number-field bound. | 0 |
| 10 | Čech (ratios conjecture, quadratic, 2110.04409); Gao–Zhao cubic/quartic ratios 2311.08626 / 2402.17198; Gao–Zhao first moments 2308.10602; Gao–Zhao Möbius quadratic twists 2212.09241 | [A]/[V] Negative moments (1/L) are **GRH-conditional** by necessity. 2212.09241 assumes GRH and is above the diagonal only. | M_κ needs the unconditional opposite range. | 0 |
| 11 | Friedberg–Offen, arXiv:2507.06963 (Jul 2025) | [A] Fundamental lemma for Hecke correspondences, cubic Shimura lift PGL₃ ↔ cubic cover of SL₃ (relative trace formula). It is **not** a Kuznetsov formula or large sieve on the 6-fold cover of SL₂. | spectral (V′) / W(w) | 0. No metaplectic spectral large sieve over Q(ω) exists. Louvel's Kuznetsov (Forum Math 2014) is on the 3-fold cover only. |
| 12 | Pascadi–Thorner, arXiv:2508.14888 (*GL_n large sieves … via positive semi-definiteness*) | [A] Large sieves for GL_n families covering the coefficients of L, L^{−1} and log L, with density estimates and no Ramanujan input. | Möbius coefficients appear, but the bound is coefficient-blind (PSD/duality), so it keeps the column term. The coefficient-blind impossibility applies (OBSTRUCTIONS §3). | 0 |
| 13 | Hamdar, arXiv:2509.06152 (Math. Ann. 2026) | [A] First moment of cubic Hecke L over Q(√−3) at every s in the strip; transition at s = 1/3; secondary term Q^{5/6} at s = 1/2 from Patterson residues. | It confirms that the theta-polar "shadow" is a true main term (Gram), so it is not an input. | 0 |
| 14 | Zero density for thin families: BGL Cor 1.6; de Faveri Cor 1.5 (n = 6: exponents 0.727, 0.466 at σ = 3/4, 7/8 vs DH 0.5, 0.25) | [V via R5] | row counts | 0. P78's counts are already stronger. [M] Even **perfect DH counts** give only 0.86667 (= 13/15 + O(δ0)) with the provable floor, and 0.86979 with floor 1/50. |
| 15 | Sextic theta coefficients: Bröker–Hoffstein 1312.0568 (Math. Comp. 2016), Wellhausen 1996, CFH conjecture (squares only; proved over F_q(t)), Friedberg–Ginzburg descent | [V via LIT_ideas234] τ_6(π)² is numerically a cubic Gauss sum times sporadic norm-7/13 factors, with "no pattern". Nothing new in 2023–2026 found [A]. | Gram/(V′) would need the arithmetic of τ_6 at squarefree indices. | 0 |

**Bottom line (A).** No 2023–2026 result supplies a beyond-the-edge input for any lever. The two
items closest to usable are:

- **(i)** a sextic, quadratically twisted analogue of the DDHL fewer-rows-than-conductor mean value
  (1.11), *improved below the large-sieve exponent 3/2*. It would feed (V′), i.e. gram_b.
- **(ii)** the DFDH / DdFDS Rankin–Selberg (Zagier-regularised) Kubota machinery extended to the
  6-fold cover. That is the tool such a bound would need.

Both are research programmes, not citations.
