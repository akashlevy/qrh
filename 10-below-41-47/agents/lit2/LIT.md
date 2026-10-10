# lit2: literature and state-of-the-art scan for the 41/47 inputs (10 Oct 2026, about 14:30–16:00 UTC)

Scope: the six input classes in the brief, plus a re-check of the QRH tracker and the web for new bounds.
I did not repeat `09-new-directions/agents/lit-edge/LIT_EDGE.md` or `lit-combos/LIT.md`. I cite items from them only where
this session adds a new mapping to the **weighted** model.

**Access.** arxiv.org, ar5iv, alphaxiv, pith.science, ias.edu, episciences and every other non-GitHub host fail DNS (WebFetch) or are blocked (curl).
So my sources are WebSearch snippets, GitHub (git, raw files and github.com pages), and the local `openai/math` preprints.
No Li Xiang (lixiang90) material was fetched.

**Labels.**
- [F]: full text read this session.
- [C]: statement read on a third-party GitHub card (verbatim formula).
- [A]: abstract or snippet only.
- [I]: my inference.
- [M]: checked with the weighted model. Scripts are in this directory: `zd_knob.py`, `n6_knob.py`, and `gmodel.py` copied from `model/`.

---------------------------------------------------------------------------------------------------------------------

## 0. Bottom line

1. **No published or preprint tool, up to 10 Oct 2026, supplies any of the six inputs.** Every item found is one of the following:
   - at or inside the large-sieve edge;
   - GRH-conditional;
   - for a full family (a character group mod q, a coset, or a spectral family);
   - weaker than what P78's own lemmas already give.

   The estimated payoff in the weighted model is **0** for every item.
2. **New calibrations [M]** (ℓ = 0.1773, B = 41/47, b scanned):
   - **Family zero density.** A density bound N_fam(σ) ≪ U^{A(1−σ)} reaches 41/47 only if **A ≤ 2.125 on all of σ ∈ [0.59, 0.80]**. Ingham suffices below 0.588.
     - This would beat the best known exponent **for ζ itself** on [0.59, 0.76]: the sup there is 2.308, at σ = 0.7 (Ingham = Guth–Maynard).
     - Transplanting ζ's entire ANTEDB table to the sextic family still leaves the worst exponent at +0.0118, i.e. no change.
     - P78's own counts already beat Ingham and Guth–Maynard at the binding bin. At δ = 0.475, x = 1/2, R = 0.586 corresponds to A = 2.232, against GM 2.243 and Ingham 2.376.
   - **Numerator sixth moment.** Σ_u|S_u(n₀)|⁶ ≪ U^{1+λ}, with n₀ ≈ 0.445, closes 41/47 only for **λ ≲ 0.09**:

     | λ | worst exponent |
     |---|---|
     | 0 | −0.0118 |
     | 0.08 | −0.0016 |
     | 0.10 | +0.0009 |
     | 0.15 | +0.0073 |

     For comparison, three available routes give much larger λ:
     - a perfect sextic large sieve: λ = 3n₀ − 1 = 0.34;
     - de Faveri: λ = 1/3 + 5n₀/2 − 1 = 0.45;
     - fourth moment × Burgess pointwise (μ ≈ 0.185): λ ≈ 0.37.

     The required statement is Lindelöf on average 0.25 beyond the large-sieve diagonal. Its quadratic analogue is open, and is known only under GRH (§3.1).
3. **Removing the (KD)^{2/3}-type term "for prime-supported coefficients" cannot hold in the needed generality [A+I].**
   Dunn–Radziwiłł (Ann. Math. 2024) show, under GRH, that the cubic large sieve is sharp up to X^{o(1)}. This is a by-product of Patterson's bias for Gauss sums at primes, so the extremal vector can be taken prime-supported.
   Any usable improvement must therefore use the Möbius × slot structure; that is exactly DIAG, which the NOTES record as false when a(p) is a sextic Gauss sum.
4. **SOTA.** No QRH bound below 0.8748 has been claimed anywhere, before or after 14:00 UTC.
   **News: tracker PR #6 (ours, 10499/12000) was merged at 14:08 UTC (16:08 +0200) after the maintainer's independent verification (commit a2033c6).**
   The catalogue on `main` (0f78cd9) lists `akashlevy-20261009-weighted-numerator` as *framework-verified*; it is now the best verified record.
   Collet's PR #9 (0.87488664704) is still open and pending.

---------------------------------------------------------------------------------------------------------------------

## 1. Tracker and web re-check (after 14:00 UTC)

**Tracker** (`git ls-remote` plus a clone of GTimG/QuasiRiemannTracker; github.com pages).
- `main` = 0f78cd9, "Merge pull request #6 from akashlevy/...", 2026-10-10 16:08:20 +0200.
- Catalogue `as_of` 2026-10-10. Records:

  | id | θ |
  |---|---|
  | openai-baseline | 7/8 |
  | argonaut | 3499999/4000000 |
  | liu | 0.874957069799… |
  | proofcouncil | 0.874957019421 |
  | nielstron | 0.874957019420098946… |
  | **akashlevy-20261009-weighted-numerator** | **10499/12000** |

  All six are framework-verified.
- Pull requests:

  | PR | author | content | last commit | state |
  |---|---|---|---|---|
  | #7 | Gao | slack retuning | 11:08 UTC | open |
  | #8 | Collet | 29/33, non-record | 11:16 UTC | open |
  | #9 | Collet | quartic boundary 0.87488664704 | 12:04 UTC | open |
  | #10 | Niels | website | — | open |
  | #14 | Niels | "Track pending proof PRs and queue verification" (infrastructure, 14:26 UTC) | 14:26 UTC | open |
  | #5 | — | historical baseline | — | closed |

  Numbers #11–#13 are issues by nielstron, not PRs: "Repo blowup", "Proof viewer", "Automatically turning PRs into pending-verification entries".
- Forks: the HaileyStorm, Gaozhongpai and nielstron forks have no branches beyond those PRs. There is no `submission-status` branch yet.
- **Hailey Collet's site** (HaileyStorm.github.io): last commit 12:08 UTC, "Add the quartic QRH paper and proof links". It holds only the PR #8 and #9 papers. Her cited manuscript "[C] shared detector witnesses and signed Gram estimates" is not public.

**Other repositories.**
- `openai/math` HEAD is still fd4aeeb (Oct 7).
- Unchanged HEADs: Argonaut 5971383, Liu 7d10420, Gao 8c99552, Shunia c3abed0, Rubinstein e2c0a74, tomoto0 b83a73c, Creech b1ebd8d, CrossNox f677cc0, LW56e2 b83dbbd.
- aster2024 (Guo) moved to a371d87. The new commits are non-QRH: almost-prime pairs, hot spots.
- New QRH-adjacent repositories, none with a bound:
  - `DrMurphyIsIn/qrh-debruijn-newman`: Λ ≤ 9/32 for the de Bruijn–Newman constant, from QRH, in Lean;
  - `CiangCing14/The-Quasi-Riemann-Hypothesis-Open-AI-Deepseek-Chinese`: a translation;
  - `the-omega-institute/trureturing` PR #14952 (merged Oct 10): a bot "literature card" for de Faveri Thm 1.1 used in a Robin-inequality project. It claims no QRH bound and only gives a "conditional 19/36 envelope" for its own packet.
- GitHub repository search ("quasi riemann", "zero-free half-plane" created after Oct 5) and web search ("quasi-Riemann" with 0.8748, 0.8747, 41/47, 13/15; arXiv 2610) find only Liu 2610.12234, Kintali 47/48, press pieces and the above.
- **No claimed θ < 0.8748, and no new mechanism.**

---------------------------------------------------------------------------------------------------------------------

## 2. Model calibrations [M]

`zd_knob.py` replaces the count by R = min(R_model, max(1−δ, A(a)(1−a))), with a = (1+δ)/2, in U-units (DH ⇔ A = 2 ⇔ R = 1−δ). The copy multiplicity (U/q)^{1/6} is harmless here because A(1−σ) ≫ 1/6.

| A(σ) supplied | worst min(E_old, E_wt) |
|---|---|
| none / Ingham 3/(2−σ) / GM 15/(3+5σ) / ζ ANTEDB table / Chen–Gupta–Li 7/3 / de Faveri n=6 (≈2.9) | +0.0118 (unchanged; binding at δ = 0.475, x = 1/2) |
| constant 2.2 | +0.0092 |
| constant 2.15 | +0.0043 |
| 2.14 or 2.13, each combined with Ingham | +0.0018 / +0.0011 |
| 2.12 combined with Ingham | **−0.0009** |
| constant 2.1 | −0.0041 (the binding bin moves to δ = 0.25) |
| DH | −0.0284 |

`n6_knob.py` gives the sixth-moment loss threshold λ ≈ 0.09 (table in §0). A Hölder bound with the sixth moment replaces the fourth-moment σ by (1−R+λ)/6. Prime selection is unavailable because 3n₀ > 1.

---------------------------------------------------------------------------------------------------------------------

## 3. Items by input class

### 3.1 Sixth or higher moments; Lindelöf on average beyond the edge (cubic or sextic families, quadratic twists of cubic theta)

| item | statement | regime | payoff |
|---|---|---|---|
| Munsch–Toma, *Bounds for moments of quadratic character sums and theta functions*, arXiv 2503.05269 [A] | Σ_d \|Σ_{n≤Y}χ_d(n)\|^{2k}: an asymptotic for all integer k in a *restricted* range (multivariate Tauberian); unconditional lower bounds in a wide range; **sharp upper bounds for arbitrary length only under GRH** | the quadratic analogue of our numerator 6th moment; a short range means Y^k below the edge | 0 (GRH) |
| Gao, arXiv 2101.08483 [A, already in 09] | sharp unconditional upper bounds for quadratic L-moments only for 2k ≤ 4 | — | 0 |
| *Simultaneous non-vanishing of Dirichlet L-functions*, arXiv 2604.11941 [A] | sharp upper bound for a mollified 6th moment **under GRH**; "asymptotic … well beyond reach … unless additional averaging" | full family | 0 |
| Huxley 6th moment; *The sixth moment of Dirichlet L-functions at the central point*, arXiv 2409.01457 (CLMR line) [A] | asymptotic for the 6th moment over all primitive χ mod q, with averaging over q removed | full family, at the edge | 0 |
| Güloğlu, PAMS 153 (2025), arXiv 2306.09474; Güloğlu–Yeşilyurt 2301.10979 [A] | cubic over the Eisenstein field: first-moment asymptotic and second-moment upper bound; higher mollified moments **under GRH** | ≤ 2nd moment | 0 |
| Lin–Wong, arXiv 2511.08783 [A] | Keating–Snaith lower bounds for a thin cubic Hecke family over Z[ω] via twisted one-level density, **under GRH** | — | 0 |
| David–Devin–Fazzari–Waxman, arXiv 2608.06286 [A] | average analytic rank of y² = x³ − dx **under GRH**. It has a section "Lindelöf on average" for Dirichlet series of quartic Gauss sums; per the snippet, the analytic properties come from earlier work (DDHL), so there is no new LoA bound | quartic analogue of the Kubota-series mean value (lit-combos #5) | 0 |
| dFDH 2607.07911 [A] | uses an LoA *upper* bound for the second moment of Kubota series (DdFDS Lemma 6.6 type, σ > 1) as an input | above the diagonal | 0 |

No sixth moment, and no LoA beyond the edge, exists for any cubic, sextic or quadratic thin family unconditionally. The needed statement (λ ≲ 0.09, n₀ = 0.445) is 0.25 beyond the large-sieve diagonal 3n₀.

### 3.2 Large values and zero density beating Ingham at σ ≈ 0.65–0.80

| item | statement | payoff |
|---|---|---|
| ANTEDB (teorth/expdb @ c8eda5e, 26 Sep; no newer commit) [F, `scratch/expdb`] | best A(σ) for ζ: 3/(2−σ) on [1/2, 7/10] (Ingham 1940, unimproved); 15/(3+5σ) on [0.7, 0.76] (Guth–Maynard); Ivić/Bourgain pieces above. Contains ζ entries only, no character families | 0 (calibration) |
| Bin Chen, with Gupta–Li (merged), arXiv 2507.08296 [A] | Σ_{χ mod q} N(σ,T,χ) ≪ (qT)^{7(1−σ)/3+ε} uniformly, i.e. Guth–Maynard in the q-aspect; new input: GCD-twisted affine sums. A reviewer flagged an N^{σ/6} vs N^σ typo | 0: A = 7/3 > 2.125, and it needs the full group mod q |
| Cossaboom, arXiv 2609.27715 [A] | GL₂ Guth–Maynard in the t-aspect, 8/3 → 5/2; ℓ^q coefficient control | 0 |
| Darbar–Kerr–Munsch–Shparlinski, arXiv 2604.02960 (3 Apr 2026) [A] | zero density on average over **thin subgroups** of characters mod q (sets with small multiplicative doubling), via Heath-Brown 1979 | 0: needs subgroup or coset structure, which our family lacks (cf. lit-edge §1.1) |
| Maynard–Pratt, arXiv 2206.11729 [A] | 24/11 **assuming zeros lie on finitely many vertical lines** | 0 |
| Tao–Trudgian–Yang, arXiv 2501.16779 (Math. Comp. 2025) [A] | new exponent pairs and ζ zero density near σ = 1 | 0 |
| Heath-Brown 1995 (quadratic); Corrigan(–Zhao) 2211.00260, 2310.16518; de Faveri Cor. 1.5 [from 09] | Ingham form or worse in the family aspect (n = 6: A ≈ 2.9 at 3/4) | 0 |

The needed A ≤ 2.125 on [0.59, 0.80] would be a first even for ζ (§2).

### 3.3 n-th order large sieves (de Faveri and follow-ups)

- **de Faveri, *Optimal large sieve for fixed order characters*, arXiv 2610.04045 (2 Oct 2026)** [C: trureturing card, verbatim; A: abstract].
  - Hypotheses: n ≥ 3; K ⊇ μ_n; a finite set S of places containing those above n with O_S of class number one; fixed ray data; reciprocity χ_a(b)χ_b(a)^{−1} = α([a],[b]); **both a and b n-th-power-free**.
  - Statement: Θ_n(A,B) ≪_ε (AB)^ε (A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}).
  - Abstract: "expected to be optimal in all ranges", with applications to the cubic second moment by conductor and to moments and zero density for every n ≥ 3.
  - For n = 6 the bound is diagonal only for B ≤ A^{1/2} or B ≥ A².
  - **No follow-up found.** Searches on arXiv 2610, "sextic large sieve" and "fixed order" return nothing newer.
- **Why removing the extra terms is out of reach [A+I].**
  - dFDH 2607.07911 proves the cubic and quartic extra terms are sharp unconditionally (theta vector).
  - Dunn–Radziwiłł, Ann. Math. 200 (2024), prove under GRH that the cubic large sieve is "sharp up to X^{o(1)}", a by-product of the Patterson bias Σ_p g̃(p) ≍ X^{5/6}/log X. So prime support alone does not remove the term; the sextic analogue is expected.
  - Only a coefficient-specific statement (Möbius × slot polynomial with non-Gauss slot coefficients) could. That is the DIAG row of the BRIEF table (−0.008), which the NOTES record as false for a(p) = sextic Gauss sum.
- Other 2026 large-sieve papers: Baier 2503.18009 (square moduli); "weighted large sieve through Parseval" 2609.25885 (classical, constants only). Both irrelevant.
- **Payoff 0.**

### 3.4 Metaplectic spectral large sieve or Kuznetsov (3- or 6-fold cover, level aspect)

- Louvel (Forum Math 2014): fixed level, cubic cover.
- Dunn, ANT 19 (2025): a cubic Duke–Iwaniec bilinear large sieve with cusp-form coefficients and cubic Kloosterman sums, |K(m,n,c)| ≤ 2^{ω(c)}N((m,n,c))N(c)^{1/2} via Livné–Patterson.
- Kaneko, arXiv 2502.18366: a subconvex *metaplectic prime geodesic theorem*, Kubota character, soft input.
- Friedberg–Offen 2507.06963: a fundamental lemma on SL₃.
- Herman 1208.1705: spectral sums of a cubic Dirichlet series.
- **None has a level aspect, and nothing exists on the 6-fold cover** (all [A]). Searches on "Kuznetsov n-fold metaplectic cover", "spectral large sieve metaplectic" and "level aspect cubic cover" found nothing new.
- **Payoff 0.**

### 3.5 Power-saving second moments of quadratic twists at the edge

- With only log savings at the edge:
  - Li, Invent. 2024 (2208.07343);
  - Shen–Stucky 2402.01497 (four main terms);
  - Jiang–Shen–Tang 2603.21739 (L′);
  - Kumar–Mallesham–Sharma–Singh.
- Chen–Lin, arXiv 2608.25812: a power-saving second moment of L(f⊗χ)L(χ), but over the **full** group mod q and **conditional on Selberg's eigenvalue conjecture** (attribution inconsistent across indexes).
- Gao–Zhao: arXiv 2506.22256 (Ramanujan J. 70 (2026) 42) gives a smoothed mean square of Σ λ_f(n)χ_{8d}(n); range not visible [A]. Their 2509.20085 needs GRH.
- Pascadi, Forum Math. Pi 14 (2026): savings at the edge only for exceptional-spectrum large sieves with sparse-Fourier weights.
- **No power saving for any thin quadratic-twist family; none for twists of the cubic theta.**
- **Payoff 0.**

### 3.6 Möbius sums over twist families below the diagonal (M_κ type)

- Gao–Zhao 2212.09241 assumes GRH and needs Y ≪ X^{1−ε}, i.e. above the diagonal.
- Toma, arXiv 2405.18778: unconditional asymptotics for k-th moments of Σ_{n≤Y}μ(n)χ_{8d}(n) via de la Bretèche's multivariate Tauberian theorem in a restricted range [A]. Range not visible; by the method it is short Y.
- Banks–Shparlinski (powerful moduli, TAMS 2020) [A].
- Nothing with Y > X (family smaller than length) without GRH, consistent with OBSTRUCTIONS §3 (Prop O1, move calculus).
- **Payoff 0.**

---------------------------------------------------------------------------------------------------------------------

## 4. Map

| input (BRIEF table) | nearest literature | gap to 41/47 [M] |
|---|---|---|
| numerator 6th moment λ ≲ 0.09 | large sieve λ = 0.34; de Faveri λ = 0.45; GRH (Munsch–Toma-type) λ = 0 | needs LoA 0.25 beyond the diagonal, unconditionally |
| DIAG / marked + plain capacities | de Faveri (diagonal to B ≤ A^{1/2}); sharpness by dFDH and DR(GRH) | coefficient-specific; false for general and prime-supported vectors |
| row counts (any form of zero density) | ζ: A ≈ 2.22–2.31 on [0.65, 0.76]; family ≥ 7/3 | needs A ≤ 2.125 on [0.59, 0.80] |
| lossless wall W(0.045), Gram g ≤ 0.75 | DDHL (1.11) at 3/2; DDFW LoA section (GRH paper) | no LoA below the large-sieve exponent for twisted Kubota series |
| M_κ, Möbius below the diagonal | Gao–Zhao (GRH), Toma (short range) | GRH-type |

**Conclusion.** The literature to 10 Oct 2026 contains no provable input for any row of the BRIEF table. The cheapest-looking rows, the lossless wall and the 0.8/0.4 capacities, have no published analogue in any family. For counts, the model makes precise that a density estimate would have to beat ζ's best on [0.59, 0.76].

## Sources (URLs; all access via search snippets unless marked)
- Tracker: https://github.com/GTimG/QuasiRiemannTracker (cloned; PRs #6–#10, #14; issues #11–#13), https://github.com/HaileyStorm/HaileyStorm.github.io
- https://github.com/the-omega-institute/trureturing/pull/14952 (card `Library/Analytic/faveri2026fixedorder.md`, saved in `scratch/omega/`)
- ANTEDB: https://github.com/teorth/expdb (cloned, `scratch/expdb`)
- arXiv IDs cited above: 2610.04045, 2607.07911, 2109.07463, 2507.08296, 2609.27715, 2604.02960, 2206.11729, 2501.16779, 2503.05269, 2604.11941, 2409.01457, 2306.09474, 2301.10979, 2511.08783, 2608.06286, 2502.18366, 2403.13151, 2608.25812, 2506.22256, 2509.20085, 2402.01497, 2603.21739, 2212.09241, 2405.18778, 2609.25885, 2503.18009, 2603.10746
