# Prior work, competing results and literature

This file holds citations and links only. No third-party paper text is vendored in this repo. Titles marked † were filled in from general knowledge rather than from a source we read, so check them before citing.

**Li Xiang.** This project has a standing rule not to copy Li Xiang's work. His repository is listed below as public metadata only: the README and repository metadata were read, and no proof code was fetched or reproduced.

## 1. The starting point: OpenAI's QRH papers and Lean development

| | |
|---|---|
| **P78** | OpenAI, *The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane Re s > 7/8*, preprint, 30 Sep 2026. [PDF @ fd4aeeb](https://github.com/openai/math/blob/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/paper.pdf). Source `build/paper.tex`. Part I gives 11/12; Part II, with prime compensation, gives 7/8. |
| **P2** | OpenAI, *The Quasi-Riemann Hypothesis* (alternate 11/12 proof), preprint, 5 Oct 2026. [PDF @ fd4aeeb](https://github.com/openai/math/blob/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/paper2.pdf) |
| Lean | OpenAI Lean formalization of the 7/8 theorem: 2924 modules, about 486k lines, Lean 4.34.1, in [openai/math `lean/`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/lean). Doc: [`lean/docs/003.md`](https://github.com/openai/math/blob/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb/lean/docs/003.md). Licence Apache-2.0. Commits used: `adc7f124` (initial), `fd4aeeb2e` (Oct 7). |
| lean-base | Yongxi Lin, [CoolRmal/quasi-riemann-seven-eighths](https://github.com/CoolRmal/quasi-riemann-seven-eighths) @ `6ff3bd5`: an isolated extraction of the 7/8 Lean sources with a Linux Comparator workflow |

Other OpenAI manuscripts in the same repository were scanned (R5) and found not usable:
- "Uniform exclusion of Landau–Siegel zeros" (Oct 1);
- "Primitive roots for every admissible integer base" (Oct 4; Re s > 1 − 10⁻⁶);
- "An unconditional first moment for cubic Gauss sums" (Sep 25);
- and others.

## 2. Competing and comparable QRH bounds (state as of 2026-10-10)

The public leaderboard is the [QRH Bounds tracker](https://gtimg.github.io/QuasiRiemannTracker/) ([GTimG/QuasiRiemannTracker](https://github.com/GTimG/QuasiRiemannTracker)). Snapshots: [`data/tracker_catalogue_main_20261009.json`](data/tracker_catalogue_main_20261009.json) and [`data/public_search_final_20261009.json`](data/public_search_final_20261009.json). Outside claims were read, not re-run.

| Claimant | θ | Status / link |
|---|---|---|
| OpenAI | 7/8 = 0.875 | Tracker: framework-verified (baseline) |
| Argonaut Math | 3499999/4000000 | [repo @ 5971383](https://github.com/Argonaut-Math/argonaut-math-quasi-riemann-boundary/tree/5971383edcb0b1cf863927eed3b359d2d3ba346e). Tracker: verification pending ([PR #4](https://github.com/GTimG/QuasiRiemannTracker/pull/4)) |
| Baiying Liu | 34999/40000, and (1507 − 2√921)/1653 ≈ 0.874957069799 (lemma-set optimum) | [arXiv:2610.12234](https://arxiv.org/abs/2610.12234), *Slightly improved zero-free half-planes for the quasi-Riemann hypothesis*; Lean in [liubaiying101/Slightly-improved-…](https://github.com/liubaiying101/Slightly-improved-zero-free-half-planes-for-the-quasi-Riemann-hypothesis). Tracker: pending ([PR #3](https://github.com/GTimG/QuasiRiemannTracker/pull/3)) |
| Zhongpai Gao | 20999/24000 | [Gaozhongpai/seven-eighths-refinement-certificates](https://github.com/Gaozhongpai/seven-eighths-refinement-certificates/blob/309a757ea539e2bf11558d8b5539e1269b792386/README.md) |
| CrossNox / LW56e2 | 4374785349/5·10⁹ | [CrossNox/quasi-riemann-hypothesis-extension](https://github.com/CrossNox/quasi-riemann-hypothesis-extension/blob/f677cc0a820bb2aff86bb4cbe783d076f5c097db/README.md); [LW56e2/quasi-riemann-improved](https://github.com/LW56e2/quasi-riemann-improved/blob/b83dbbd914c5751acf9bc3389cde58970437278f/README.md) |
| GPT-6 Astra & J. M. Shunia | cubic root θ\* ≈ 0.8749570194200989 (weaker rational 162861/186136) | [jshunia/rhtm `compensated_half_plane_refinement.pdf`](https://github.com/jshunia/rhtm/blob/c3abed0cae8fdb0090c12b0e4d20aea68b6e067d/compensated_half_plane_refinement.pdf) |
| Jizhou Guo | cubic root; weaker 43747851/50000000; conditional 328/375 (Hypothesis 11.1) | [aster2024/openai-math-followups `zero-free-beyond-7-8`](https://github.com/aster2024/openai-math-followups/blob/5750cc81e489729a2111d6c9c073f4adecdb9de7/zero-free-beyond-7-8/README.md) |
| Tim Gehrunger / ProofCouncil | 874957019421/10¹² | Tracker: framework-verified |
| Nielstron | 874957019420098946128604623/10²⁷ | Tracker: framework-verified, current best ([PR #1](https://github.com/GTimG/QuasiRiemannTracker/pull/1), merged) |
| Nielstron | cubic root, seven-target three-kernel | [PR #2](https://github.com/GTimG/QuasiRiemannTracker/pull/2) (open) |
| Steven E. Creech | 437479/500000 (draft) | [stevenecreech/math](https://github.com/stevenecreech/math/blob/b1ebd8d8661538992cbc2f27f97f0aa37a15525c/README.md) |
| Khamzat Chakhkiev | 10937/12500 (conditional) | [Zenodo 23221313](https://zenodo.org/records/23221313) |
| Li Xiang (lixiang90) | cubic σ\* = 11/12 − e\*/4; the Lean theorem is conditional on `ArithmeticProbeObligation` | [lixiang90/RH-Zero-Free-Formalization](https://github.com/lixiang90/RH-Zero-Free-Formalization) (metadata only) |
| Shiva Kintali | 47/48 (short proof) | [*A Short Proof of the Quasi-Riemann Hypothesis*](https://shivakintali.github.io/papers/QRH.pdf), 7 Oct 2026 |
| tomoto0 | Lean rebuild of the 7/8 proof; heuristic figures 0.87496, 0.87451, 0.86979, 13/15 | [tomoto0/quasi-riemann-hypothesis-7-8-verification](https://github.com/tomoto0/quasi-riemann-hypothesis-7-8-verification) |
| **This project (Akash Levy)** | **10499/12000 ≈ 0.8749166** | [PR #6](https://github.com/GTimG/QuasiRiemannTracker/pull/6) (open). If merged, it would be the best listed entry by about 4·10⁻⁵. |

Recorded as non-comparable (applications, conditional work or different statements):
- Puglisi, arXiv:2210.03121;
- arXiv:2606.16097, 2605.24363, 2610.08965;
- Rubinstein, `mrubinst/L-function_zero_free_strips`;
- Zenodo 23228131, 23248611, 23269461, 23269135;
- several GitHub issues.

## 3. Literature by topic

### Large sieves for characters of fixed order
- D. R. Heath-Brown, *A mean value estimate for real character sums*, Acta Arith. 72 (1995) 235–275 (the quadratic large sieve; the Lean development uses this style).
- D. R. Heath-Brown, *Kummer's conjecture for cubic Gauss sums*, Israel J. Math. 120 (2000) 97–124 (the cubic large sieve Q + N + (QN)^{2/3}).
- L. Goldmakher, B. Louvel, *A quadratic large sieve inequality over number fields*, Math. Proc. Camb. Phil. Soc. 154 (2013) 193–212, [arXiv:1112.1642](https://arxiv.org/abs/1112.1642). Cited by P2; never used in Lean.
- V. Blomer, L. Goldmakher, B. Louvel, *L-functions with n-th-order twists*, IMRN 2014(7) 1925–1955, [arXiv:1112.1650](https://arxiv.org/abs/1112.1650). n-th order large sieve; zero density (Cor. 1.6).
- S. Baier, M. P. Young, *Mean values with cubic characters*, J. Number Theory 130 (2010) 874–890, [arXiv:0804.2233](https://arxiv.org/abs/0804.2233).
- A. de Faveri, *Optimal large sieve for fixed order characters*, [arXiv:2610.04045](https://arxiv.org/abs/2610.04045) (2 Oct 2026). Improves P78 Part I's sextic row count; changes θ in neither paper.
- A. de Faveri, A. Dunn, J. Hoffstein, *Non-orthogonality of the cubic and quartic large sieves via Rankin–Selberg*, [arXiv:2607.07911](https://arxiv.org/abs/2607.07911). The cubic (AB)^{2/3} term is sharp.
- P. Gao, L. Zhao: quartic large sieve, [arXiv:1101.3646](https://arxiv.org/abs/1101.3646); *Mean squares of quadratic twists of the Möbius function*, J. Number Theory 247 (2023), [arXiv:2212.09241](https://arxiv.org/abs/2212.09241) (GRH-conditional, above the diagonal only); sextic Hecke L-functions of Q(ω), [arXiv:2201.01885](https://arxiv.org/abs/2201.01885); first moments via double Dirichlet series, [arXiv:2306.10726](https://arxiv.org/abs/2306.10726).
- L. Guth, J. Maynard, *New large value estimates for Dirichlet polynomials*, Ann. of Math. 203 (2026), [arXiv:2405.20552](https://arxiv.org/abs/2405.20552). Covered by the coefficient-blind obstruction.
- M. N. Huxley, *The large sieve inequality for algebraic number fields*, Mathematika 15 (1968).

### Cubic, sextic and quartic theta functions; Gauss sums
- T. Kubota, *On automorphic functions and the reciprocity law in a number field*, Lectures in Math. Kyoto Univ. 2 (1969), [hdl:2433/84907](http://hdl.handle.net/2433/84907). Also *Ein arithmetischer Satz über eine Matrizengruppe*, J. reine angew. Math. 222 (1966).
- S. J. Patterson, *A cubic analogue of the theta series*, J. reine angew. Math. 296 (1977) 125–161 (and part II). *The cubic Shimura correspondence*, Asian J. Math. 2 (1998). Arithmetic-progression and prime equidistribution results: Proc. LMS (3) 54 (1987); [arXiv:1411.7594](https://arxiv.org/abs/1411.7594).
- D. R. Heath-Brown, S. J. Patterson, *The distribution of Kummer sums at prime arguments*†, J. reine angew. Math. 310 (1979).
- A. Dunn, M. Radziwiłł, *Bias in cubic Gauss sums: Patterson's conjecture*, Ann. of Math. 200 (2024) 967–1057, [arXiv:2109.07463](https://arxiv.org/abs/2109.07463).
- A. Dunn, *Metaplectic cusp forms and the large sieve*, Algebra & Number Theory 19 (2025), [arXiv:2403.13151](https://arxiv.org/abs/2403.13151).
- C. David, A. de Faveri, A. Dunn, J. Stucky, *Non-vanishing for cubic Hecke L-functions*, [arXiv:2410.03048](https://arxiv.org/abs/2410.03048).
- C. David, A. Dunn, A. Hamieh, H. Lin, *Quartic Gauss sums over primes and metaplectic theta functions*, [arXiv:2306.11875](https://arxiv.org/abs/2306.11875).
- C. Castillo, A. de Faveri, A. Dunn, *Non-vanishing for quartic Hecke L-functions and ranks of elliptic curves*, [arXiv:2604.01316](https://arxiv.org/abs/2604.01316).
- R. Bröker, J. Hoffstein, *Fourier coefficients of sextic theta series*, Math. Comp. 85 (2016), [arXiv:1312.0568](https://arxiv.org/abs/1312.0568).
- D. Kazhdan, S. J. Patterson, *Metaplectic forms*, Publ. Math. IHÉS 59 (1984) 35–142.
- F. Gao, *Distinguished theta representations for certain covering groups*, Pacific J. Math. 290 (2017), [arXiv:1602.01880](https://arxiv.org/abs/1602.01880).
- G. Wellhausen, *Fourierkoeffizienten von Thetafunktionen sechster Ordnung*, PhD thesis, Göttingen 1996.
- T. Suzuki, *Some results on the coefficients of the biquadratic theta series*, J. reine angew. Math. 340 (1983) 70–117.
- C. Eckhardt, S. J. Patterson, *On the Fourier coefficients of biquadratic theta series*†, Proc. LMS (3) 64 (1992).
- J. Hoffstein, *Theta functions on the n-fold metaplectic cover of SL(2) — the function field case*, Invent. Math. 107 (1992).
- S. Friedberg, D. Ginzburg, *Descent and theta functions for metaplectic groups*, JEMS 20 (2018), [arXiv:1403.3930](https://arxiv.org/abs/1403.3930); *Metaplectic theta functions and global integrals*, JNT 146 (2015), [arXiv:1403.3929](https://arxiv.org/abs/1403.3929).
- P. Deligne, *Sommes de Gauss cubiques et revêtements de SL(2)*†, Sém. Bourbaki 1979/80, exp. 539.
- N. Proskurin, *Cubic metaplectic forms and theta functions*, LNM 1677 (1998).
- Y. Flicker, *Automorphic forms on covering groups of GL(2)*†, Invent. Math. 57 (1980).

### Kubota series, Eisenstein series, functional equations
- B. Brubaker, D. Bump, *On Kubota's Dirichlet series*, J. reine angew. Math. 598 (2006). It assumes μ_{2n} ⊂ F, so it does not apply to n = 6 over Q(ω).
- J. Elstrodt, F. Grunewald, J. Mennicke, *Groups acting on hyperbolic space*, Springer 1998.
- W. Müller, *The trace class conjecture in the theory of automorphic forms*†, Ann. of Math. 130 (1989).
- J. Bernstein, E. Lapid, *On the meromorphic continuation of Eisenstein series*, [arXiv:1911.02342](https://arxiv.org/abs/1911.02342).
- C. Mœglin, J.-L. Waldspurger, *Spectral decomposition and Eisenstein series*, CUP 1995.
- D. Zagier, *The Rankin–Selberg method for automorphic functions which are not of rapid decay*†, J. Fac. Sci. Univ. Tokyo 28 (1981).
- A. Diaconu, *Mean square values of Hecke L-series formed with r-th order characters*, Invent. Math. 157 (2004). This is the "Kubota–Diaconu" functional equation used in the quartic refutation.

### Multiple Dirichlet series and Weyl groupoids
- B. Brubaker, D. Bump, S. Friedberg, J. Hoffstein, *Weyl group multiple Dirichlet series III: Eisenstein series and twisted unstable A_r*†, Ann. of Math. 166 (2007); also *Coefficients of the n-fold theta function and Weyl group multiple Dirichlet series* (2012).
- G. Chinta, S. Friedberg, J. Hoffstein, *Multiple Dirichlet series and automorphic forms*, PSPM 75 (2006); *Double Dirichlet series and theta functions*, Patterson Festschrift (2012).
- G. Chinta, P. Gunnells, *Constructing Weyl group multiple Dirichlet series*, JAMS 23 (2010), [arXiv:0803.0691](https://arxiv.org/abs/0803.0691); *Weyl group multiple Dirichlet series of type A₂*, [arXiv:math/0703040](https://arxiv.org/abs/math/0703040); *Weyl group multiple Dirichlet series constructed from quadratic characters*, [arXiv:math/0604449](https://arxiv.org/abs/math/0604449).
- W. Sawin, A. Whitehead, *Functional equations of axiomatic multiple Dirichlet series, Weyl groupoids, and quantum algebra*, [arXiv:2507.08662](https://arxiv.org/abs/2507.08662).
- A. Diaconu, B. Ion, V. Paşol, A. Popa, twisted moments of r-th order Hecke L-functions, [arXiv:2607.27131](https://arxiv.org/abs/2607.27131).
- S. Friedberg, J. Hoffstein, Ann. of Math. 142 (1995). B. Brubaker, S. Friedberg, J. Hoffstein, *Cubic twists of GL(2) automorphic L-functions*, Invent. Math. 160 (2005). S. Friedberg, J. Hoffstein, D. Lieman, Math. Ann. 327 (2003). A. Diaconu, D. Goldfeld, J. Hoffstein, Compositio 139 (2003).
- A. Bucur, A. Diaconu, *Moments of quadratic Dirichlet L-functions over rational function fields*, Moscow Math. J. 10 (2010) (an affine D₄^{(1)} example).

### Spectral and bilinear methods
- B. Louvel, *On the distribution of cubic exponential sums*, Forum Math. 26 (2014), [arXiv:1107.0197](https://arxiv.org/abs/1107.0197).
- R. Livné, S. J. Patterson, *The first moment of cubic exponential sums*, Invent. Math. 148 (2002).
- P. E. Herman, [arXiv:1305.1243](https://arxiv.org/abs/1305.1243), [arXiv:1208.1705](https://arxiv.org/abs/1208.1705).
- N. Watt, weighted spectral large sieve for SL(2, ℤ[i]), Funct. Approx. 48 (2013), [arXiv:1302.3127](https://arxiv.org/abs/1302.3127).
- S. Friedberg, O. Offen, metaplectic Kuznetsov distributions on the cubic cover of SL₃, [arXiv:2507.06963](https://arxiv.org/abs/2507.06963).
- A. Pascadi, *Large sieve inequalities for exceptional Maass forms and applications*, Forum Math. Pi 14 (2026), [arXiv:2404.04239](https://arxiv.org/abs/2404.04239).
- W. Duke, H. Iwaniec, *Bilinear forms in the Fourier coefficients of half-integral weight cusp forms and sums over primes*†, Math. Ann. 286 (1990).
- R. Hamdar, *Hecke L-functions away from the central line*, [arXiv:2509.06152](https://arxiv.org/abs/2509.06152).

### Background on quasi-Riemann hypotheses and zero-free regions
- H. Iwaniec, E. Kowalski, *Analytic Number Theory*, AMS 2004. H. Davenport, *Multiplicative Number Theory*. K. Ireland, M. Rosen, *A Classical Introduction to Modern Number Theory*.
- Cited by the exposition paper: Bettin–Gonek (2017), Bhowmik–Ruzsa (2018), Montgomery–Vaughan, Ford (2002), Vinogradov and Korobov (1958), Landau (1903; the prime ideal theorem, avoided via b⁶ rows).
