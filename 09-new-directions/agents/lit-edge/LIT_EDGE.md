# lit-edge: targeted literature scan for beyond-the-edge inputs (10 Oct 2026, ~11:15 UTC)

Scope: the five targeted questions in the brief. I did not repeat anything in `agents/lit-combos/LIT.md`
(de Faveri 2610.04045, dFDH 2607.07911, DdFDS, DDHL, Castillo–dF–D, Dunn ANT 2025, DIPP, Hamdar, Pascadi–Thorner,
Friedberg–Offen, BGL, Bröker–Hoffstein); I cite them only where a mapping needs them.

**Access.** arxiv.org, doi.org, Springer, AMS, Semantic Scholar, ems.press and researchgate are blocked. WebFetch also fails
(DNS) for every non-GitHub host. So the sources are WebSearch snippets plus GitHub. I found three useful GitHub-hosted sources:
- `teorth/expdb` (Tao–Trudgian–Yang ANTEDB, commit c8eda5e, 26 Sep 2026): full text of the zero-density chapter.
- `AlexenderSokolov/cross-theory-formal-database`: about 15k CC-BY arXiv LaTeX sources, readable via raw.githubusercontent.
  It is only a random sample. From it I read Conrey–Rodgers in full.
- The local `openai/math` preprints (749 dirs).

**Labels.**
- [F]: I read the statement in full text this session.
- [A]: abstract or snippet level.
- [K]: standard background that I did not re-verify this session.
- [I]: my own inference or computation.

Excerpts are saved under `scratch/`, and `calib.py` / `calib.out` hold the calibrations.

**Levers.**
- W(w), with w ≥ 0.064.
- (H1)/(H2): Lindelöf on average plus decorrelation for K quadratic twists of the level-L² marked cubic theta. At the 27/31 geometry, K = Z^{0.86} and L = Z^{0.204} = K^{0.237}.
- M_κ with κ ≤ 0.895 (0.964 with DIAG).
- The plain sixth moment.
- The Gram saving (V′).

---------------------------------------------------------------------------------------------------------------------

## 0. Bottom line

**No item in the 2023–Oct 2026 literature supplies any lever.** The scan sharpens three facts.

1. **Every known beyond-the-edge Lindelöf-on-average (LoA) result for a thin character family uses one of two structures that our family lacks.**
   - A full character group mod d on POWERFUL moduli, where χ is p-adically analytic on cosets (Postnikov). Examples: Petrow–Young, Dasgupta, Garcia–Young.
   - A full spectral family with reciprocity: Conrey–Iwaniec, Petrow–Young, Hu–Petrow–Young.

   Our residual twists ψ_k have squarefree k, because the Möbius μ-absorption forces it, and one character per modulus.
   So the coset mechanism is excluded structurally, not just technically.
2. **Asymptotic large sieves stop AT the edge.**
   - CIS covers N ≤ Q^{2−ε}, which is the family size. The CLMR eighth moment sits exactly at the edge.
   - For quadratic families, even under Lindelöf for every L(s,χ_d), the best general asymptotic (Conrey–Rodgers [F], Čech [A]) controls only "0/1-swap" terms, i.e. polynomial length < D^2.
   - W(w) is a first moment of length K^η with η = 2(1+w/M) = 2.149 at w = 0.064. That is the 2-swap regime, where Conrey–Rodgers' error term (with LH) already exceeds the main term by D^{0.037} [I, calib.out].
3. **"Ingham below σ = 7/10" is open for ζ itself.** ANTEDB [F] lists Ingham's 3/(2−σ) as the best known A(σ) for 1/2 ≤ σ ≤ 7/10, unchanged since 1940.
   - M_κ, κ ≤ 0.895, would therefore be the first improvement in this range for any family, ζ included.
   - The plain sixth moment is open even for the quadratic Dirichlet family: unconditional sharp bounds stop at the fourth moment.

**Calibration [I].** (H1) at the 27/31 geometry implies the pointwise bound D_P(½,k) ≪ K^{1/2} = C^{0.202}, with C = (KL)² the conductor of
each twist. At the 19/22 optimum it implies C^{0.190}. This is a subconvex bound, between convexity 0.25 and Burgess 0.1875, for twists of
the cubic-cover theta. No subconvexity of any strength is known on the cubic cover: Dunn [A] notes that nothing beyond Rankin–Selberg is known even for single coefficients.
In the GL2 analogue, such a bound is inside the Blomer–Harcos hybrid range. So (H1) is plausible, but it needs GL2 tools (Euler product, amplification over
a full spectral family) that do not exist on the cover.

---------------------------------------------------------------------------------------------------------------------

## 1. Q1 — Lindelöf on average beyond the edge for thin twist families

### 1.1 Results that DO go beyond the edge (family < AFE length), and why

| item | statement | family vs AFE length | mechanism |
|---|---|---|---|
| Good 1982 (Mathematika); Iwaniec 1980 | [K, Good via Dasgupta's abstract [A]] ∫_T^{T+T^{2/3}} \|L(f,½+it)\|² ≪ T^{2/3+ε} for f of level 1; ∫_T^{T+T^{2/3}}\|ζ\|⁴ ≪ T^{2/3+ε} | T^{2/3} vs T | t-aspect, Kuznetsov/Voronoi (archimedean "coset") |
| Petrow–Young, *Fourth moment of Dirichlet L-functions along a coset and the Weyl bound*, Duke 172 (2023), arXiv 1908.10346 | [A] LoA for Σ_{ψ mod d}\|L(½,χψ)\|⁴ (χ mod q) when q* \| d, q* minimal with q² \| q*³. For q = p^k, d ≈ q^{2/3}. Gives Weyl for all χ. | d ≈ q^{2/3} vs q (4th moment) | coset of a subgroup on powerful moduli |
| Dasgupta, *Short second moment bound for GL(2) L-functions in q-aspect*, arXiv 2309.14593 | [A] second moment of L(½, f⊗χψ), f of level 1, ψ in a coset of characters mod q^{2/3}, q = p³; "q-aspect analogue of Good's T^{2/3} bound" | q^{2/3} vs q | same (Postnikov) |
| **Garcia–Young**, *Asymptotic second moment of Dirichlet L-functions along a thin coset*, Forum Math. Sigma 13 (2025) e83, arXiv 2312.08482 | [A] **power-saving asymptotic** for Σ_{ψ mod d}\|L(½,χψ)\|² when v_p(d) ≥ v_p(q)/2 for all p, and still power-saving for d > q^{2/5}; secondary main term of size ≈ q^{1/2} (or ≈ d), not predicted by CFKRS | d ≥ q^{2/5} vs q^{1/2}: **asymptotic beyond the edge** | same |
| Conrey–Iwaniec, Annals 151 (2000); Petrow–Young cubic moments; **Hu–Petrow–Young**, arXiv 2506.14741 (v2 Mar 2026) | [K]/[A] Σ_{f} L(½,f⊗χ)³ ≪ q^{1+ε} over spectral families of level q (HPY: square-full q, local components specified, including supercuspidal). Weyl bounds. | family q vs cube length ≈ q^3 | nonnegativity (Waldspurger) + Petersson/Kuznetsov + Voronoi; full spectral family |
| Chen–Lin, arXiv 2608.25812 (26 Aug 2026) | [A] power-saving asymptotic for "the second moment of the product" L(½,f⊗χ)·L(½,χ) over primitive even χ mod q (I could not verify the exact form; only the abstract and the reference list were indexed). f has level 1, the moduli q are restricted, and the result is **conditional on Selberg's eigenvalue conjecture**. | if the moment is \|L(f⊗χ)L(χ)\|²: q vs q^{3/2} | full group mod q, i.e. an unbalanced shifted convolution with spectral theory and Kloosterman-product sums |
| Dunn–Zaharescu, *Twisted second moment of modular half-integral weight L-functions*, JEMS 27 (2025) 1–69, arXiv 1903.03416 | [A] power-saving asymptotic for Σ_{χ mod p}\|L(½,f,χ)\|² with f a half-integral-weight Kohnen newform (**no Euler product**), unconditional | at the edge (family p, length p) | full group mod p; shifted convolution, Salié sums |

**Common structure [I].** Each result has either a full character group mod d on powerful moduli, where χ(1+dx) is additive (Postnikov), so
coset families behave like short t-intervals, or a full spectral family with positivity and reciprocity. Our ψ_k are quadratic characters of
squarefree moduli, one per modulus. The squarefreeness is forced by μ-absorption (NOTES §22, step 2). Kazhdan–Patterson cube marks p³, the only
powerful-modulus freedom, amplify nothing (NOTES §9). **Payoff for W(w)/(H1)/(H2): 0.** The Dunn–Zaharescu result is the right template *if* the
residual twist could be made full (NOTES §22 option (b)), but it lives at the edge, not beyond.

### 1.2 Quadratic-twist families: nothing beyond the edge, nothing uniform in the level

- **Li**, *Moments of quadratic twists of modular L-functions*, Invent. Math. (2024), arXiv 2208.07343 [A]. This is the unconditional SY asymptotic, with only log savings (NOTES §10).
  The key input is an "optimal large-sieve-type bound".
- **Kumar–Mallesham–Sharma–Singh** (first main term for L′) and **Jiang–Shen–Tang**, arXiv 2603.21739 (Mar 2026; three main terms for the second moment of L′) [A].
  Both run on Li's large sieve. No power saving is advertised.
- **Shen–Stucky**, *Fourth moment of quadratic Dirichlet L-functions II*, arXiv 2402.01497 [A]. Unconditional, four main terms. The extra input is a
  fourth-moment version of Li's large sieve plus cross-cancellation between diagonal and off-diagonal. Log savings only.
- **OpenAI**, *Mean analytic rank of quadratic twists* (6 Oct 2026), `sections/analytic.tex` [F, `scratch/oai_meanrank_analytic_props.tex.txt`].
  This is the cleanest public write-up of Li's method: induction on the height, prime-square inflation d ↦ dp² (p ≍ P), and compact Fourier support.
  - Prop. "Dirichlet-polynomial second moments": (1/M) Σ_{d~M} |Σ_n λ_F(n)χ_d(n)n^{−½−it}G(n/N)|² ≪ log^B(3M)(1+|t|)^A **(1+N/M)^A**.
  - Beyond the edge (N > M) it loses a power of N/M. Inflation must take P ≥ N/M so that the dual height N²/(MP²) ≤ M, and recovery costs P.
    That is the same loss as the large sieve [I]. **So Li's technique does not cross the edge.**
- **Unconditional power-saving SY-type asymptotic at the edge since 2023: none found.**
- **Hybrid / uniform-in-level LoA for quadratic twists with D < N^{1/2}: none found.**
  - Jackson–Knightly 1204.2253 [A] averages over forms (relative trace formula), not over twists.
  - Ye 1404.2336 [A] is Rankin–Selberg, not quadratic twists.
- Higher moments for the quadratic family:
  - Unconditional sharp upper bounds exist only for moments ≤ 4: Gao, arXiv 2101.08483 [A], for 0 ≤ k ≤ 2 in his normalisation, plus Heath-Brown's large sieve [K].
  - Asymptotics exist for moments 1–4: Jutila, Sono, Young (3rd, power saving), Shen–Stucky (4th, log).
  - The 5th and 6th moments are open unconditionally. **This is the quadratic analogue of our "plain sixth moment" lever.**

### 1.3 Calibration of (H1) [I, calib.py]
- Each twist θ⊗ψ_kχ_P^{−2} has conductor C ≈ (KL)², per wall-spectral §0.
- (H1) gives D_P(½,k) ≪ K^{1/2} = C^{M/(4(M+ℓ))}: C^{0.2021} at (M,ℓ) = (0.860, 0.204) and C^{0.1897} at the 19/22 optimum.
- Reference exponents: convexity 0.25, Burgess-type 0.1875, Weyl 0.1667.
- So (H1) does not imply sub-Weyl. It does imply a subconvex bound in the twist aspect for quadratic twists of the cubic theta, and none of any strength is known.
  - Dunn (ANT 2025) [A] remarks that he knows no bound for cubic metaplectic cusp coefficients ρ_f(·) beyond Rankin–Selberg.
  - Dunn–Zaharescu's LoA for half-integral weight gives only convexity pointwise.

---------------------------------------------------------------------------------------------------------------------

## 2. Q2 — Asymptotic large sieve for thin families

- **Conrey–Iwaniec–Soundararajan**, *Asymptotic large sieve*, Crelle (2012), arXiv 1105.1176 [A].
  - It covers all primitive χ mod q, averaged over q ≤ Q, for N up to Q^{2−ε}, with off-diagonal main terms coming from switching moduli.
  - **N ≤ Q² is the family size, so this is the edge, not beyond.**
  - The sixth moment (CIS) is within the edge (L³ has length Q^{3/2} < Q²).
  - The eighth moment (Chandee–Li under GRH; Chandee–Li–Matomäki–Radziwiłł, Duke 173 (2024) 3453–3493, unconditional [A]) is exactly at the edge.
- **Switching moduli** (q ↔ (m−n)/q) needs the full group mod q, so that orthogonality gives a congruence. For one-character-per-modulus families the only
  substitute is Poisson in the modulus (Soundararajan). That gives Gauss-sum-weighted duals, which is exactly P78's reflection. **No thin-family
  version of CIS exists.** Searches for quadratic or cubic "asymptotic large sieve" found nothing beyond the following:
  - Conrey–Farmer–Soundararajan, *Transition mean values of real characters*, JNT 82 (2000) [A]: an asymptotic for Σ_{m≤X}Σ_{n≤Y}(m/n) valid for all X, Y, with constant coefficients only.
  - Gao–Zhao, *Mean values of some Hecke characters*, Acta Arith. 187 (2019), arXiv 1708.08164 [A]: quadratic, cubic and quartic analogues, comparable lengths, smooth coefficients.
- **Conrey–Rodgers**, *Averages of quadratic twists of long Dirichlet polynomials* [F, Thm 1, `scratch/conrey_rodgers_thm.tex.txt`].
  - Statement: assume LH for L(s,χ_d). Then Σ_d μ²(2d)Ψ(d/D) Σ_n W(n/N) τ_A(n)χ_{8d}(nℓ)n^{−1/2} equals the CFKRS recipe restricted to |U| ≤ 1
    (0- and 1-swap terms), plus O(N^{1/4}D^{1/2+ε}ℓ^{1/4+ε}).
  - This is a power saving iff N = D^η with η < 2. Terms with |U| = j are negligible only for η < j.
  - **Čech**, arXiv 2402.07473 [A]: MDS continuation for the k-th moment (conditional on GLH for k ≥ 5); unconditional power-saving asymptotics for
    the first three moments; detection of the 0/1-swap terms for every k under GLH; recovers Conrey–Rodgers.
- **Mapping to W(w) [I].**
  - Σ_k |S_k|² is the first moment over k of a product polynomial of length N² = K^η, with η = 2(1+w/M).
  - At w = 0.064, η = 2.149, which is the 2-swap regime. There the Conrey–Rodgers error, even under LH, is D^{+0.037} × main term.
  - In MDS language:
    - finite type corresponds to η < 2, with power savings (Young's third moment, Čech);
    - affine corresponds to η = 2, with log savings (Li, Shen–Stucky);
    - η > 2 is untouched.
  - This agrees with NOTES §16/§22 (affine C̃₂). W(w) sits just past the affine point.
- Thin families with extra averaging over moduli, all at or within the edge:
  - Baluyot–Pratt, JEMS (2022) [A]: L(½,χ_p) for prime p, more than 9% non-vanishing, order of the second moment (one quadratic character per prime modulus).
  - Basak(–Zaharescu), arXiv 2409.12474 [A]: q in short intervals or APs.
  - Drappeau–Pratt–Radziwiłł [A]: one-level density with extended support, full family.
- **Payoff: 0.**

---------------------------------------------------------------------------------------------------------------------

## 3. Q3 — Metaplectic spectral theory over Q(ω)

- **Louvel**, *On the distribution of cubic exponential sums*, arXiv 1107.0197 [A]: a Bruggeman–Kuznetsov formula for the cubic cover over Q(ω) at fixed level,
  applied to the distribution of cubic exponential sums.
- **Proskurin**, *Cubic metaplectic forms and theta functions*, LNM 1677 (1998) [A].
- A thesis "Twisted Kloosterman sums and cubic exponential sums" (d-nb.info) refines Livné–Patterson via Kuznetsov [A].
- **Not found:** any level-aspect Kuznetsov formula on the 3-fold cover (needed for twists of level k²); any Kuznetsov or large sieve on the 6-fold cover;
  any spectral large sieve for cubic metaplectic Maass forms; any nontrivial pointwise bound for cubic-cover coefficients (Dunn's remark above).
- **Half-integral weight analogues exist**, and that is the contrast:
  - Iwaniec, Invent. 87 (1987) and Duke, Invent. 92 (1988) [K]: nontrivial coefficient bounds via averaging over levels plus Kloosterman/Salié sums.
  - Blomer (2008) spectral second-moment lemma, improved in arXiv 2301.11901 [A].
  - Lam, JTNB (2014) [A]: large sieve for half-integral weight.
  - Ahlgren–Andersen [K]: Kuznetsov for half-integral weight.
  - Dunn–Zaharescu (above).
- **Why they do not transfer [I].** Half-integral weight gains come from Salié sums being explicitly evaluable (sums over √ mod c), which feed a
  full spectral family. Cubic Kloosterman sums have no such evaluation (Patterson/Livné give only distribution results), and our family is not full.
- Quadratic twists of the cubic theta: only DdFDS Lemma 6.6 (σ > 1) and the DDHL quartic (1.11) at 3/2, both already in LIT.md.
  Nothing else exists, including no subconvexity for Kubota series ψ(h,s) in the h-aspect (search found none).
- **Payoff: 0** for (H1)/(H2)/W/(V′).

---------------------------------------------------------------------------------------------------------------------

## 4. Q4 — Zero density below σ = 7/10; sixth moments

- **ANTEDB** [F, `scratch/antedb_zero_density_table.tex.txt`], "Current best upper bound on A(σ)" for ζ:
  - A(σ) = 3/(2−σ) for **1/2 ≤ σ ≤ 7/10** (Ingham 1940);
  - 15/(3+5σ) for 7/10 ≤ σ < 19/25 (Guth–Maynard);
  - further entries above that.
- The ANTEDB remarks:
  - Conrey's mollified bounds N(σ,T) ≪ T^{1−θ(2σ−1)} hold for θ < 4/7. "An improvement to Ingham's bound would only occur once θ exceeded 2/3."
  - Heuristic: Ingham corresponds to the L² mean-value theorem (τ₀ = 2−σ). Beating it needs LV(σ,τ₀) ≤ 3−3σ at some τ₀ > 2−σ.
- **Consequence for M_κ (κ ≤ 0.895) and the row-count lever [I].** "Ingham below 7/10 for the sextic family" would beat a bound that has not
  been improved for ζ itself. Any gain must come from sextic-family structure (theta/Gauss-sum arithmetic), not from large-values technology.
  This matches R5's (O2): M_κ is beyond the density hypothesis (DH).
- Family zero-density results (all Ingham-type or weaker below 7/10):
  - Heath-Brown 1995, quadratic family: Q^{3(1−σ)/(2−σ)+ε}T^{…} [A, as quoted in 2402.12133]. This is exactly Ingham-type.
  - Bin Chen, arXiv 2507.08296 [A]: Guth–Maynard for χ mod q, (qT)^{7(1−σ)/3+ε}, improving Huxley's 12/5 (the useful range is near 3/4).
  - Corrigan–Zhao, Bull. Aust. MS 108 (2023), arXiv 2211.00260 [A]; Corrigan, arXiv 2310.16518 [A] (via Petrow–Young's Weyl bound; better only far from ½).
  - Cossaboom, arXiv 2609.27715 [A]: GL2, Guth–Maynard, uniform exponent 8/3 → 5/2.
  - Sun–Wang–Yu, arXiv 2505.23573 [A]: Selberg-type, near ½.
- **Sixth moment.** No LoA sixth moment for any cubic or sextic family.
  - Güloğlu (cubic over the Eisenstein field) [A] has a first-moment asymptotic and a second-moment upper bound.
  - For the quadratic Dirichlet family the 5th/6th moments are open unconditionally (§1.2).
- **Payoff: 0.**

---------------------------------------------------------------------------------------------------------------------

## 5. Q5 — SOTA and new QRH work (checked 11:00–11:15 UTC)

- **Tracker** GTimG/QuasiRiemannTracker:
  - `main` = c480675 (12:36 +0200, "Merge PR #4 nielstron/verify-argonaut"). PRs #1–#8 exist (`scratch/tracker_refs_now.txt`).
  - Catalogue `as_of` 2026-10-10. Best framework-verified entry is still Nielstron's 0.874957019420… No record change.
- **PR #8** (HaileyStorm = Hailey Collet, "developed with assistance from GPT-6 Pro"), *A zero-free half-plane from a single compensating prime* [F, abstract, §1, Lemma 4.1]:
  - **29/33 = 0.8788 > 7/8, a non-record.** It is a simplified P78 with scales (x,y,ℓ) = (13,15,5)/33 and "one compensating prime".
  - Two shared witnesses at the same zero.
  - A "signed additive Gram bound" that is the same as P78's: (Q₀/Y)(1 + P_a^{1/6} + P_a²/Y).
  - It avoids the prime-marked moment. **No new mechanism below 7/8**, though it may be a useful simpler template.
  - It cites an unseen manuscript "[C] A zero-free half-plane from shared detector witnesses and signed Gram estimates" (8 Oct 2026).
- **openai/math** HEAD is still fd4aeeb (Oct 7).
- Competitor HEADs:
  - Gao 8c99552 (unchanged).
  - aster2024 11838f1: two new commits since 867b0dc, both non-QRH (hot spots, rough pairs).
  - Argonaut 5971383, Liu 7d10420, CrossNox f677cc0, tomoto0 b83a73c, Shunia c3abed0, Creech b1ebd8d, lixiang90 ebf33a7 (not compared with earlier HEADs).
- New repos:
  - **mrubinst/L-function_zero_free_strips** (Michael Rubinstein, 9 Oct): expository guides to 11/12 and 7/8, plus a CM-curve adaptation (conditional).
    Its §"How tight" re-evaluates P78's exponent and concludes "**7/8 appears to be close to the limit of this particular design**".
  - **LW56e2/quasi-riemann-improved**: Lean re-tuning to 0.8749570698 (above SOTA).
  - GitHub repository search for "quasi-Riemann", "zero-free half-plane" and "QRH" finds only these plus spam ("qrhX" content repos).
- Web: only Liu 2610.12234 and press pieces. **Nothing below 0.8749, and no new mechanism anywhere.**

---------------------------------------------------------------------------------------------------------------------

## 6. Map: item → lever → payoff

| item | (H1)/(H2), W(w) | M_κ / counts | plain 6th | Gram (V′) | payoff |
|---|---|---|---|---|---|
| Petrow–Young coset, Dasgupta, Garcia–Young | beyond edge only for powerful-modulus cosets; excluded by squarefree μ-absorption | — | — | — | 0 |
| Conrey–Iwaniec / Hu–Petrow–Young cubic moments | needs a full spectral family + positivity | — | — | — | 0 |
| Chen–Lin 2608.25812, Dunn–Zaharescu JEMS 2025 | full family; DZ is the template for option (b) "make the residual twist full" | — | — | — | 0 |
| Li / KMSS / JST / Shen–Stucky; OAI mean-rank Props | at the edge, log savings; inflation loses (N/M)^A beyond the edge [F] | — | — | — | 0 |
| CIS, CLMR | at the edge (N ≤ family) | — | — | — | 0 |
| Conrey–Rodgers [F], Čech | 0/1-swaps only (η < 2), even under LH; W needs η = 2.149 | — | — | — | 0 |
| ANTEDB [F] | — | Ingham unbeaten for σ ≤ 7/10 even for ζ | — | — | 0 (calibration) |
| Chen 2507.08296, Corrigan(–Zhao), Cossaboom, Heath-Brown 1995 | — | ≥ Ingham below 7/10 | — | — | 0 |
| Gao 2101.08483 | — | — | quadratic analogue open beyond the 4th moment | — | 0 |
| Louvel, Proskurin, half-integral weight Kuznetsov/large sieve | no level-aspect or 6-fold-cover version; Salié trick has no cubic analogue | — | — | — | 0 |
| PR #8 (Collet, 29/33) | — | — | — | same Gram bound as P78 | 0 |

## 7. What the literature says a successful input would have to look like [I]
1. **Beyond-edge LoA requires a "complete" structure.**
   - Option 1: a full group of characters mod d with d powerful relative to the conductor (Postnikov).
   - Option 2: a full spectral family with positivity and reciprocity (Motohashi / Conrey–Iwaniec / Kwan type).

   For our family both are blocked at the source: μ-absorption forces squarefree moduli, and only cubes k³ appear in the Kubota index.
   A route must either find positivity plus a spectral family containing the k-twists of the cubic theta (no metaplectic level-aspect Kuznetsov exists), or change the detector.
2. **In swap language, W(w) is a 2-swap statement just past the affine point (η − 2 = 0.15).** No paper controls 2-swap terms for any thin
   family, even assuming GLH. Of all the levers, a "lossless 2-swap" for theta coefficients (where reflection is exact) is the one closest to existing technology.
   Edge-theta's Proposition U is already of this type: a power-saving asymptotic for Z ≤ H ≤ Z^{3/2}, unmarked, CONDITIONAL on P78's energy lemma.
3. **Row-count levers (M_κ, "Ingham below 7/10") would be firsts even for ζ.** They should be deprioritised unless a sextic-specific identity is found.
