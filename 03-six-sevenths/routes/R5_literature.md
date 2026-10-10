# R5 — Literature search for an input giving theta <= 6/7

(working notes; started 2026-10-08)

## Target inputs (what would suffice)

* (M_kappa) sum_{0<N(u)<=D^kappa} |A_u(D)|^2 << D^{1+kappa+eps}, kappa = 6/7,
  A_u(D) = sum_n mu(n) nu(n) chi_n(u) W(N n / D)  (sextic symbol over Z[omega]).
  Gives theta = max((6+5k)/12, 1-k/6) [paper2 prime extraction] or (6+5k)/12 [7/8 paper, U~1].
* Equivalent dual target (paper2 Step 3): sum_{N(h) << D^2/H} |B_h(D)|^2 << D^2 with H = D^kappa,
  i.e. dual row range calH = D^{2-kappa} = D^{8/7} > D (columns).  Below-diagonal for the original;
  above-diagonal (rows > columns) for the dual.  Known to be FALSE in the positive form (prior note (3)).

## Local manuscript scan

### Local OpenAI manuscripts (719 dirs, 8007 .tex files) — grep results
* No local manuscript claims a zero-free half-plane better than 7/8. grep for 6/7, frac{6}{7}, 0.857: only
  non-number-theory hits.  Papers mentioning zero-free regions use 7/8 (QRH Sep-30) as an INPUT:
  A-Power-Saving-for-Polynomial-Differences-at-Prime-Arguments (Thm zero-free input = QRH Thm 1.1, 7/8),
  the-additive-indecomposability-of-the-primes (Remark: optional 7/8 refinement),
  The-E-times-CM-Component-of-Zilber-Pink (uses 7/8).
* Primitive-roots-for-every-admissible-integer-base (Oct 4): Thm hecke-zero-free: Re s > 1 - 10^{-6} for Hecke
  L-functions of cyclotomic F containing mu_12.  Same reflection/quadratic-large-sieve mechanism over F;
  exponent far weaker than 7/8.  Not usable.
* Uniform-exclusion-of-Landau-Siegel-zeros (Oct 1): (1-beta) log q >= c via interpolation determinants.
  Only real zeros, log-scale; irrelevant to a fixed half-plane.
* An-unconditional-first-moment-for-cubic-Gauss-sums (Sep 25): Patterson first moment over primes,
  sum_{N pi<=X} G(pi) ~ (6/5)c_* X^{5/6}/log X; angular modes l != 0 cancel.  Inputs: HB cubic LS
  (Q+N+(QN)^{2/3}), HB metaplectic mean square int|f_r(1/2+eps+it)|^2 << (Nr)^{1/2+4eps}T^2,
  DR level Voronoi.  New: Prop 'exceptional character moments' sum_{p,q}|A(p,q)|^2 << Y^{7/3-eta}
  (unspecified eta>0) for PRIME-CONVOLUTION coefficients at two configurations of the CUBIC sieve.
  Not a mean square for Mobius against sextic characters below the diagonal; eta ineffective-small.
  Not usable for M_kappa.

## Web / arXiv
### arXiv 2607.07911 (de Faveri, Dunn, Hoffstein, Jul 2026) "Non-orthogonality of the cubic and quartic large sieves via Rankin-Selberg"
* Thm 1.1 (UNCONDITIONAL): Xi_3(A,B) >>_eps (AB)^{-eps}(A + B + (AB)^{2/3})  [squarefree, Z[omega]].
  => Heath-Brown/BGL cubic large sieve is OPTIMAL; the "cubic large sieve conjecture" (remove (QN)^{2/3})
  is FALSE for general coefficients.  (Under GRH this was DR24 Thm 1.4.)
* Thm 1.2: Xi_4(A,B) >> (AB)^{-eps}(A+B+A^{3/4}B^{1/2}+A^{1/2}B^{3/4}).
* Conj (1.8): Xi_n(A,B) = (AB)^{o(1)}(A + B + A^{1-1/n}B^{2/n} + A^{2/n}B^{1-1/n}).
  For n = 6: A + B + A^{5/6}B^{1/3} + A^{1/3}B^{5/6}  (at A=B=X: X^{7/6}, vs BGL X^{4/3}).
* Mechanism: extremal vector beta_b = normalized cubic Gauss sum g~_3(b); polar term (Patterson residue at 5/6)
  P_B(a) ~ g~_3(a) N(a)^{-1/6} B^{5/6}  [DdFDS24 Lemma 6.7], remainder controlled by Lindelof-on-average
  second moment of Kubota series psi^{(3)}(a,s) in the twist aspect [DdFDS24 Lemma 6.6].

### Baiying Liu, "Slightly improved zero-free half-planes for the quasi-Riemann hypothesis" (arXiv submit/8200716, 8 Oct 2026;
github liubaiying101/Slightly-improved-zero-free-half-planes-for-the-quasi-Riemann-hypothesis)
* Same lemma set as the 7/8 paper; varies b and ell only.  B_r = 34999/40000 (b=1/8, ell=1/6+1/10^4) — identical to
  our QRH-analysis.md section 3; B_new = (1507 - 2 sqrt 921)/1653 = 0.874957069799 (algebraic b, ell).
* Consistent with our exact LP lower bound 33216105/37963129 = 0.8749569879 for the paper's lemma set
  (Liu's B_new exceeds it by ~8e-8).  Liu also "expects B_new optimal for these estimates".
* Lean 'V2Conditional.lean' is conditional on 'PhysicalProbeEstimates' at the new geometry (not 6/7).
* => independent confirmation that the 7/8 lemma set stops at ~0.87496.  No 6/7.

### Blomer-Goldmakher-Louvel, "L-functions with n-th order twists", IMRN 2014 (arXiv 1112.1650)
* Thm 1.3 (n>3, K contains mu_n): sum*_{Na<=M} |sum*_{Nb<=N} lambda_b chi_a(b)|^2 << (MN)^eps (M + N + (MN)^{2/3}) ||lambda||^2.
* Cor 1.6 (zero density, family of n-th order chi_a, squarefree a):
  sum*_{Na<=N} N(sigma,T,a) << N^{g(sigma)} T^{1+d(1-sigma)/(3-2sigma)} (NT)^eps,
  g = 2(10s-7)(1-s)/(24s-12s^2-11) (5/6<s<=1),  8(1-s)/(7-6s) (1/2<s<=5/6).
  Compare density hypothesis 2(1-s): at s=0.9, g=0.4545 vs 0.2.  Counting only (g->0 as s->1, no exclusion).
* Plug-in: only as a row-count input in the 7/8 architecture (Prop detector-counts).  Even PERFECT counting
  (density hypothesis for every row family) gives 0.86979 (sensitivity.py), and GRH for rows gives 13/15.
  => cannot reach 6/7.  Not usable.

### Baier-Young, "Mean values with cubic characters", JNT 2010 (arXiv 0804.2233)
* Thm 1.4/1.5: large sieve for order-3 (and order-6) DIRICHLET characters mod rational q:
  (QM)^eps min{Q^{5/3}+M, Q^{4/3}+Q^{1/2}M, Q^{11/9}+Q^{2/3}M, Q+Q^{1/3}M^{5/3}+M^{12/5}} ||a||^2.
  Different family (rational moduli), general coefficients, always contains the column term >= M.
  Not usable for M_kappa (needs beating the column term).

### Gao-Zhao, "Mean squares of quadratic twists of the Mobius function", JNT 247 (2023) (arXiv 2212.09241)
* Thm 1.1 ASSUMES GRH: S(X,Y) = sum*_d (sum_n mu(n) chi_{8d}(n) Phi(n/Y))^2 W(d/X) = c X Y + O(X^{1/2+eps}Y^{3/2+eps} + X Y^{1/2+eps}).
  Asymptotic only for Y << X^{1-eps}, i.e. rows >= columns (above the diagonal).  GRH used for 1/L bounds.
* Quadratic analogue of M_kappa is the opposite range (Y = columns > X = rows) and is open even
  as an asymptotic under GRH (as an upper bound it is trivial under GRH).  Not usable.

### David-de Faveri-Dunn-Stucky, arXiv 2410.03048 (v2 Mar 2026) "Non-vanishing for cubic Hecke L-functions"
* Mollified 2nd moment asymptotic for cubic family; Lemma 6.6 Lindelof-on-average 2nd moment of Kubota
  series psi^{(3)}(a,1+eps+it) over twists a (sum_{Na<=A}|psi|^2 << A^{1+eps}(1+|t|)^R); Lemma 6.7 polar term.
* Mollifier lengths < family size (rows > columns).  No below-diagonal Mobius mean square.  Not usable.
### Castillo-de Faveri-Dunn arXiv 2604.01316 (quartic analogue, non-vanishing + ranks): same regime.  Not usable.
### David-Dunn-Hamieh-Lin arXiv 2306.11875 (ANT): quartic Gauss sums over primes; Type I/II; Suzuki at squares. Not usable.
### Dunn arXiv 2403.13151 (ANT 19 (2025)) "Metaplectic cusp forms and the large sieve": cusp-form coefficient
  bilinear large sieve, power saving for sum over primes of rho_f; level of distribution past 2/3. Cusp (not theta)
  coefficients; first moments.  Not usable for a mean square below the diagonal.
### Broker-Hoffstein arXiv 1312.0568 "Fourier coefficients of sextic theta series": partially proven/conjectured
  tau(pi)^2, tau(pi^2), tau(pi^4).  Paper2 needs only the CUBIC theta (Patterson, proven).  Not an input.
### Patterson's quartic conjecture: still OPEN (dFDH 2026 p.7: "Determining tau_4(a) ... is a central open problem").

## Structural obstruction for ANY literature input (why the search must come back empty)
(O1) M_kappa (kappa<1) contains the single row u=1: |A_1(D)|^2 <= D^{1+kappa+eps} => zero-free Re s > (1+kappa)/2
     (paper2 Sec. 3 Mellin argument).  So M_kappa is itself a quasi-RH; no pre-Sep-2026 paper can contain it.
(O2) Sharper: M_kappa + a zero detector gives #{u: N(u)<=H, L(s,nu chi_u) has a zero with Re >= sigma}
     << H^{1-(2 sigma-1)/kappa+eps}.  For kappa = 6/7 this is H^{1-7(2sigma-1)/6}, STRICTLY STRONGER than the
     density hypothesis H^{2(1-sigma)} for every sigma>1/2 (at sigma=7/8: H^{1/8} vs H^{1/4}).
     Large-sieve/moment methods give at best DH-type counts (BGL Cor 1.6 is weaker than DH).  So M_{6/7} is a
     "beyond density hypothesis" statement for the sextic family: out of reach of every published technique.
     Consistency check (heuristic, zeros-only model): with DH counts, sum_u |A_u|^2 ~ max_s H^{2(1-s)} D^{2s}
     exceeds D^{1+kappa} by D^{(2s-1)(1-kappa)} for kappa<1, so DH-quality input reproduces only the diagonal kappa=1.
(O3) 7/8 architecture: perfect (DH) row counting -> 0.86979; GRH for rows -> 13/15 > 6/7 (barriers.py);
     low side alone >= 13/15 - w/60 under marked-wall relaxation W(w) (energy_lp.py) -> needs w >= 4/7.
     Every literature item above is a row-count or large-sieve input => cannot reach 6/7 there.

### ** de Faveri, "Optimal large sieve for fixed order characters", arXiv 2610.04045 (v1 2 Oct 2026) ** — NEW, PROVEN
* Thm 1.1: for every n>=3, K containing mu_n,
  Theta_n(A,B) := sup_{||lambda||=1} sum^flat_{Na<=A} |sum^flat_{Nb<=B} lambda_b chi_a(b)|^2
               << (AB)^eps (A + B + A^{1-1/n} B^{2/n} + A^{2/n} B^{1-1/n}),
  sums over n-th-POWER-FREE ideals (stronger than squarefree).  Proves the upper half of dFDH Conj (1.8);
  optimal for n=3,4 by dFDH lower bounds.  n=3: same as Heath-Brown (no gain).  n=6: A+B+A^{5/6}B^{1/3}+A^{1/3}B^{5/6}.
* Cor 1.5 zero density (n-th power free family): sum N(sigma,T,a) << N^{1-delta_n(sigma)+eps} T^{...}.
  n=6 exponents (computed): sigma=0.75: 0.727 | 0.8: 0.641 | 7/8: 0.466 | 0.9: 0.390  (DH: 0.5,0.4,0.25,0.2).
* PLUG-IN, 7/8 paper: the sextic sieve (Lemma sextic-large-sieve, line 4707, K+D+(KD)^{2/3}) is used ONLY in Part I
  (Prop sextic-row-count, line 5181; applied at 5921/6621).  Replacing it: R(delta)=min{1,max(1-delta/2, 4/3-delta)}
  improves to min{1,max(1-delta/2, 7/6-delta)} (derivation: #B U^{delta r} <= U+U^r+U^{5/6+r/3}+U^{1/3+5r/6}, r in [1/2,1]).
  Part I stays pinned at 11/12 by its low side (C_I(s)=s-2/3, l_x=l_y=1/2).  Part II never uses the sextic sieve
  (its counts come from Lemmas marked/plain/inverse-amplification via Prop detector-counts) and its counts are already
  stronger (R_* = 0.665 at delta=0.39,x=1/2 and 0.692 at x=0, vs 0.805 from the de Faveri count).
  Part II low side uses Heath-Brown's CUBIC sieve (Lemma hybrid-energy, eq. old-eq:4.6b, (NL)^{2/3}) — optimal (n=3).
  => theta of the 7/8 paper unchanged.
* PLUG-IN, paper2: none (needs to beat the column term; sixth-power rows u=p^6 are excluded by 'flat' anyway).
* LEAD (not either paper): quartic transplant over Q(i) via the Kubota/Eisenstein route.  Prior project script
  strengthening/quartic_routes.py modelled the "known" quartic sieve as (AB)^{3/4}; that is outdated: BGL Thm 1.3 already
  gives (AB)^{2/3} for n=4, and de Faveri gives A+B+A^{3/4}B^{1/2}+A^{1/2}B^{3/4}.  Re-run (r5/fine_n6.py, R-table 151x26,
  t 201, r 401; model calibrated: n=6 gives 0.874978 vs paper lemma-set optimum 0.874957-0.874975):
      n=4 Eisenstein + BGL (AB)^{2/3}   : 0.85665   (ell ~ 0)
      n=4 Eisenstein + de Faveri n=4    : 0.85391   (ell ~ 0)
      n=4 theta route (Patterson quartic conj) : 0.81022
  (coarse 61x21 grids gave 0.85717 for BGL — the apparent "6/7" was a grid artifact; ell=0 high-res: 0.85654.)
  CAVEATS: no quartic analogue of the 7/8 paper exists; the model transplants every lemma with 1/6->1/4, alpha=3/4,
  same capacities c_M=1/2, c_P=2/9, and ASSUMES the Kubota-FE reflection has rows Z^M, dual columns Z^{2M+ell-1} and a
  quartic twist.  My own check of the quartic dual: Poisson + absorption gives mu(n)gamma_{-1}(n) = Phi(n) conj(alpha(n))
  gamma_1(n) (verified at 538 primes, strengthening/quartic_absorption_check.py), and the dual twist conj(chi_n(h)) is
  fully absorbed into the SHIFTED quartic Gauss sum g~_4(h,n): B_h = partial sums of the Kubota series psi^{(4)}(h,s)
  (angular-twisted).  Its FE conductor in h (N(h) vs N(h)^2) and the polar/residual terms were not derived.  So: a
  research programme, not an input.

## Verdict for R5
* Usable proven inputs giving theta <= 6/7 in either paper: NONE.
* Best literature-backed exponent: 0.874957069799 = (1507-2 sqrt 921)/1653 (Liu, Oct 8 2026; same lemma set; our LP floor
  33216105/37963129 = 0.8749569879 shows nothing better is possible from those lemmas).
* Genuinely new proven input since the 7/8 paper: de Faveri Thm 1.1 / Cor 1.5 — no change to theta of either paper.
* Structural reason (O1)-(O3) above: M_{6/7} is itself a quasi-RH and a beyond-density-hypothesis statement for the
  sextic family; published techniques give at best DH-type counts (BGL, de Faveri Cor 1.5 are weaker than DH), and the
  7/8 architecture with any counting input is capped at 0.86979 (DH) / 13/15 (GRH rows).
* Scripts: r5/quartic_routes_r5.py, r5/probe_model.py, r5/fine_quartic.py, r5/fine_n6.py (all in wf67/r5/).
  PDFs/text: wf67/lit/ (2607.07911, 2610.04045, 1112.1650, 0804.2233, 2212.09241, liu1/liu2).

### Rigor status of the obstruction
* (O1) RIGOROUS: one nonnegative term (u=1) of the M_kappa sum; for each fixed u, M_kappa at all D >= N(u)^{1/kappa}
  gives A_u(D) << D^{(1+kappa)/2+eps} for all large D and all W, hence zero-free Re s > (1+kappa)/2 (paper2 Sec. 3).
* (O2) HEURISTIC: the paper's detector witnesses have length D' = U^r, r in [1/2,1] (<= conductor), whereas M_kappa only
  covers rows N(u) <= D'^kappa < D'; turning M_kappa into a family density bound needs a long-witness detector
  (D' >= U^{1/kappa}).  Stated as a consistency/size statement, not a theorem.
* (O3) relies on the project's exact LP certificates (barriers.py: 13/15 GRH-rows; sensitivity.py: 0.86979 DH counts;
  energy_lp.py: low side >= 13/15 - w/60).
