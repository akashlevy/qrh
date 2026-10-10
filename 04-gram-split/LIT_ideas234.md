# Literature scan for Ideas 2–4 (sextic theta, metaplectic Kuznetsov/large sieve, double Dirichlet series, μ-absorption)

*Compiled 9 Oct 2026. Labels:*
- **[V]** I read the statement in the source (full text or the abstract).
- **[S]** Taken from a secondary source (another paper's summary of it, or a search snippet), not from the original.
- **[I]** My own inference. It is not stated in any source.

---

## 0. Most decision-relevant facts

1. **The Gram term P_a^{1/6} appears to be the true size of the sextic large sieve, so no orthogonality argument can remove it.**
   - [V] de Faveri–Dunn–Hoffstein, arXiv:2607.07911 (July 2026), prove **unconditionally** that Heath-Brown's cubic large sieve is sharp: Ξ_3(A,B) ≫ (AB)^{−ε}(A+B+(AB)^{2/3}). This removes the GRH from DR Thm 1.4. They also show the quartic large sieve is not perfectly orthogonal.
   - [V] They **conjecture**, for every n ≥ 3 and K ⊇ Q(ζ_n), that Ξ_n(A,B) = (AB)^{o(1)}(A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}).
   - [V] They write that the lower bound "should be within reach of current tools for all n ≥ 5". They also write that "any power-saving improvement on the (AB)^{2/3+ε} term for any n ≥ 4 would be a substantial advance".
   - [I] Take n = 6, with A = rows (size Q) and B = columns (size Y). The excess over the diagonal is A^{5/6}B^{1/3}/A = (B²/A)^{1/6} = (Y²/Q)^{1/6} = **P^{1/6}**. That is exactly the 7/8 paper's Gram excess.
   - [I] The obstruction is the same rank-one "theta vector" that IDEA1_NOTES measured numerically. In DR/DFDH it is the polar term P_B(a) ∝ g̃(a) N(a)^{−1/6} B^{5/6}, a single Patterson-theta vector.
   - **Consequence [I]:** P^{1/6} can only be avoided by using the theta structure, as in Idea 1's direct theta pairing. A better large sieve cannot remove it.
   - **Correction to QRH-analysis.md §4 [V/I]:** That section says "removing (QN)^{2/3} is the cubic large sieve conjecture". The cubic large sieve conjecture is **false**: DR showed this under GRH in 2021, and DFDH showed it unconditionally in 2026. For n = 6 the conjectured truth is Q^{5/6}N^{1/3} + Q^{1/3}N^{5/6}. That is smaller than (QN)^{2/3} only when Q^{1/2} < N < Q², and it is unproved.

2. **The sextic theta coefficients at squarefree indices are unknown, and nobody has a conjecture for them.**
   - [V] Kazhdan–Patterson / Gao: the theta representation on GL_r^{(n)} has a unique Whittaker model iff n = r or n = r+1. For GL_2 that means only n = 2, 3.
   - [V] For θ_6, the Hecke relations and KP periodicity fix only τ_6(p⁵) = 0 and τ_6(p⁴) = G_1(1,p). They give two relations:
     - τ_6(p) = G_2(1,p) τ_6(p³);
     - τ_6(mp²) = G_3(m,p) τ_6(mp²), where G_3 is a quadratic Gauss sum, so this is a quadratic-residue vanishing condition.
   - [V] The CFH conjecture covers **squares** only: τ_6(m₁²m₂⁴) is given in terms of cubic Gauss sums. It is proved only over F_q(t).
   - [V] Bröker–Hoffstein have numerics for τ_6(π)². It equals a root of unity × a cubic Gauss sum × 3^{l} × elements of norm 7 or 13, and they found "no pattern".
   - **Consequence [I]:** Idea 1 step (V) cannot use arithmetic facts about ρ(m) for squarefree m. That includes "main terms only for cube μ" and the size of a′(n). Only automorphy-based analysis is available: functional equations, Voronoi, Rankin–Selberg averages. None of it has been written down for n = 6.

3. **The Voronoi machinery in the level aspect exists, but only for the cubic cover.**
   - [V] DR §5 / App. A is uniform in the modulus r. It covers **arbitrary periodic twists ψ mod r** (Lemma 5.2, Prop 5.1), with the dual twisted by a cubic-Kloosterman transform ψ̃. The residue is ∝ ψ̃(0).
   - [V] Dunn (ANT 2025, §8) does the same for cubic cusp forms with additive twists.
   - [V] David–Dunn–Hamieh–Lin (§8–9) do it for quartic Kubota series with level structure.
   - **Nothing exists for the 6-fold cover.**

4. **Mean values of Kubota series, the closest thing to step (D).**
   - [V] DdFDS Lemma 6.6 (cubic): Σ_{N(m)≤M} |ψ^{(3)}(hm, σ+it)|² ≪ N(h)^{1/2} M^{1+ε}(1+|t|)², for 1 < σ ≤ 5/4. This is a Lindelöf-on-average bound.
   - [V] DFDH prove the quartic analogue with Zagier-regularised Rankin–Selberg, and expect the same "for all n ≥ 5".
   - [V] DDHL (1.11) treat a case with **fewer rows than the conductor**: Σ_{N(a)≤A} μ²(a)|ψ̃^{(4)}(½+ε+it, a²)|² ≪ A^{3/2+ε}. Here the conductor is about A². They say A^{1+ε} "seems likely" with more work.

5. **Kuznetsov over Q(ω) on the cubic cover exists.**
   - [V] Louvel (Forum Math 2014, Thm 4.2) states it, built from Bruggeman–Motohashi and Lokvenec-Guleska.
   - [V] Louvel Thm 3.4: the exceptional spectrum of the cubic cover lies in (0, s(d)] ∪ {1/3}, with s(d) ≤ 1/6. The proof uses Flicker's correspondence.
   - [V] Livné–Patterson used the Goldfeld–Sarnak version.
   - **No spectral large sieve for metaplectic Maass forms over Q(ω) was found.** No Kuznetsov-based result for the Gram regime was found either.

6. **Follow-ups to the 7/8 paper.**
   - Nothing rigorous is below 0.874957.
   - Nothing treats the P_a^{1/6} Gram term.
   - Liu and Li Xiang only re-optimise parameters (details in §6).

---

## 1. The sextic theta function θ_6 on the 6-fold cover over Q(ω)

### What is known about its coefficients

- **Kazhdan–Patterson, "Metaplectic forms", Publ. IHÉS 59 (1984) 35–142.** ([numdam](https://www.numdam.org/item/PMIHES_1984__59__35_0))
  - [S] They prove periodicity τ_n(m pⁿ) = N(p)^{1/2} τ_n(m) and the Hecke eigen-relations in general.
  - [V, via Gao's restatement] dim Wh(Θ(GL_r^{(n)})) = 1 iff n = r (any KP cover), or n = r+1 (a special type of KP cover).
  - So for GL_2, the Whittaker model is unique only for n = 2, 3. For n ≥ 4 it is not unique, and Deligne (Bourbaki 1980) identifies this as the reason n ≥ 4 is inaccessible.
- **F. Gao, "Distinguished theta representations for certain covering groups", Pacific J. Math. 290 (2017) 333–379.** ([arXiv:1602.01880](https://arxiv.org/abs/1602.01880))
  - [V] Gives upper and lower bounds for dim Wh(Θ), and recovers the KP criterion.
  - Relevance: the precise representation-theoretic reason the n = 6 coefficients are not "Gauss sums".
- **Hecke relations for n = 6.** [V] From Bröker–Hoffstein §1 and CFH (2.5): for (m,p) = 1,
  - τ_6(m p^i) = G_{i+1}(m,p) τ_6(m p^{4−i}) for 0 ≤ i ≤ 4;
  - τ_6(m p⁵) = 0.

  Hence:
  - τ_6(p⁴) = G_1(1,p), determined;
  - τ_6(p) = G_2(1,p) τ_6(p³), linked but undetermined;
  - τ_6(mp²) = G_3(m,p) τ_6(mp²), with G_3 quadratic, so τ_6(mp²) = 0 unless m is a QR mod p (up to normalisation).
  - [V, DDHL §1.2] In general, n/2 − 1 coefficients among indices 1, π, …, π^{n−1} are undetermined for n even. For n = 6 that is **2**.
- **Bröker–Hoffstein, "Fourier coefficients of sextic theta series", Math. Comp. 85 (2016); [arXiv:1312.0568](https://arxiv.org/abs/1312.0568).**
  - [V] Computations over Q(ζ_6) with Patterson's algorithm.
  - Lemma 5.1, proved for |π| ≡ 1 mod 4: τ(π⁴)/τ(1) = |π|^{−1/2} g_6(1,ε,π).
  - Conj 5.2: for |π| ≡ 3 mod 4 there is an extra factor (−1,π)_6/√3.
  - τ(π²) is about 2G^{(3)}(1,π)/√3 up to a 12th root of unity. The evidence is overwhelming for |π| ≡ 7 mod 12.
  - Conj 5.7: (τ(π)/τ(1))² = ζ_6^{k} · g_3(1,ε²,π)/|π|^{1/2} · 3^{l} · (1−3ζ_6)^{m}(−2+3ζ_6)^{n}, or the analogous norm-13 elements. This rests on numerics for |π| < 900, and the authors found "no pattern" in k, m, n.
  - Relevance: the **only data on θ_6 at primes**. τ_6(π) behaves like a square root of a cubic Gauss sum, with sporadic algebraic factors.
- **G. Wellhausen, "Fourierkoeffizienten von Thetafunktionen sechster Ordnung", PhD thesis, Göttingen 1996.** [S] Cited by Bröker–Hoffstein. Earlier numerics; it gives τ(1,V) ≈ 0.1358547858696091.
- **Chinta–Friedberg–Hoffstein, "Double Dirichlet series and theta functions", Patterson Festschrift, Springer Proc. Math. 9 (2012) 149–170.** ([pdf](https://chinta.ccny.cuny.edu/publ/patterson.pdf))
  - [V] Conj 2.1 (the CFH conjecture) is a relation between the GL_2 theta on the double cover and the cubic objects: ζ(3u−½) Σ τ^{(6)}(m²) N(m)^{−u} = ζ(3u−½) Σ G^{(3)}(1,d) N(d)^{−u} · Σ τ^{(3)}(m) N(m)^{−u}.
  - Consequence: τ^{(6)}(m₁²m₂⁴) = G^{(3)}(1,m₂)² G^{(3)}(1,m₁) (m₂/m₁)_3² Σ_{m₁=d₁d₂} (d₂/d₁)_3.
  - It is proved over the **rational function field** for 2n-fold covers with n odd. That includes the 6-fold case (n = 3).
  - They note its naive general-n form is "certainly not true over a number field for n ≥ 5". Here n ≥ 5 means 10-fold covers and up. The sextic case over number fields is still a conjecture, and Bröker–Hoffstein's data show it needs a √3 correction at the prime 3.
  - Relevance: the only n-fold ↔ n/2-fold relation. It concerns **square** indices only.
- **Friedberg–Ginzburg, "Descent and theta functions for metaplectic groups", JEMS 20 (2018) 1913–1957.** ([arXiv:1403.3930](https://arxiv.org/abs/1403.3930))
  - [V] They run a descent from a small representation on the 3-fold cover of Sp_4, built from the cubic theta on GL_2. The result is an L² function on the **6-fold cover of SL_2**.
  - Its coefficients are arithmetic for infinitely many p. It has nonzero projection to θ_6, and they conjecture it *is* θ_6.
  - The quartic analogue comes from the 4-fold cover of SO_7.
  - Relevance: a structural link between θ_6 and the cubic theta, but only partial.
- **Friedberg–Ginzburg, "Metaplectic theta functions and global integrals", JNT 146 (2015) 134–149.** ([arXiv:1403.3929](https://arxiv.org/abs/1403.3929)) [V] A new Eulerian integral on the 3-fold cover of GL_2. For the quartic case, a Dirichlet series combining half-integral-weight coefficients with quartic Gauss sums.
- **Hoffstein, "Theta functions on the n-fold metaplectic cover of SL(2) — the function field case", Invent. Math. 107 (1992) 61–86.** [S] Explicit formulas over F_q(t). Patterson re-derived them in [arXiv:1411.7594](https://arxiv.org/abs/1411.7594) (2014) and found a hidden Aut(F_q(t))-action.
- **Brubaker–Bump–Friedberg–Hoffstein, "Coefficients of the n-fold theta function and Weyl group multiple Dirichlet series" (2012).** ([UMN](https://experts.umn.edu/en/publications/coefficients-of-the-n-fold-theta-function-and-weyl-group-multiple/))
  - [S] For n large relative to r, the Whittaker coefficients of the GL_r theta in a Hecke orbit at p equal the p-parts of stable WMDS.
  - Relevance: this is the **stable range**, not GL_2 with n = 6.
- **Quartic theta (n = 4), for comparison.**
  - **Suzuki, "Some results on the coefficients of the biquadratic theta series", Crelle 340 (1983) 70–117.**
    - [V, via DDHL (1.10)] ψ^{(4)}(a²) = g̃_4(a)/N(a)^{1/4} for squarefree a. This is the quartic theta at squares.
    - [V, via Friedberg–Ginzburg] Suzuki, Ryukyu Math. J. 25 (2012) 21–106, uses Rankin–Selberg and gets arguments but not moduli.
    - [V, via CFH] Suzuki 1997 proved the Bump–Hoffstein argument conjecture over function fields.
  - **Patterson's quartic conjecture (1984) and Eckhardt–Patterson, PLMS 64 (1992) 225–264.** [V, via Bröker–Hoffstein] D_2 = D_1², i.e. τ_4(p)² = 2G_4(p). Still open over number fields.

### Voronoi summation and twisted sums of theta coefficients (growing modulus)

- **Dunn–Radziwiłł §5 and Appendix A** (details in §2 below). [V]
  - Lemma 5.2: let r ≡ 1 mod 3, and let ψ be any **periodic function mod r** whose ψ̂ is supported on (Z[ω]/r)^×. Then the ψ-twist of the cubic theta at ∞ equals (1/(3N(r))) (F_1 + ωF*_{19} + ω²F*_{10}) at the reflected point, twisted by ψ̃.
  - Here ψ̃(u) = Σ_x ψ(x) S_{1/3}(x,u;r), the "cubic-Kloosterman (K_{1/3}) transform" (5.22–5.23).
  - Prop 5.1 gives the functional equation with N(r)^{2s}. The pole at 5/6 has residue ∝ ψ̃(0).
  - Props 5.2–5.3 give the explicit Voronoi for ψ = (·/r)_3, r squarefree, with dual length about N(r)²/X.
  - **It is fully uniform in r, so a growing modulus is already handled.**
  - Appendix A tabulates Patterson's d_j(μ) for the 27 coset reps and k_j(P,T,E).
  - [I] Take ψ = quadratic character mod r. Then ψ̃(u) = g_2(r) Σ_a χ_2(a)(a/r)_3 ě(ua/r), a **sextic** Gauss sum, and ψ̃(0) = 0, so there is no pole. This is the cubic-cover mirror of Idea 1's (V) heuristic, where a quadratic twist plus the sextic multiplier gives a cubic Gauss sum.
- **A. Dunn, "Metaplectic cusp forms and the large sieve", ANT 19 (2025) 1823–1880.** ([arXiv:2403.13151](https://arxiv.org/abs/2403.13151))
  - [V] §8 gives functional equations for cubic cusp forms twisted mod λ^ℓ q, and a level-aspect Voronoi (Prop 8.17).
  - [V] Dunn: "no 'Waldspurger-type' formula known for the coefficients of cubic metaplectic cusp forms". No bound better than Rankin–Selberg is known for ρ_f.
- **David–Dunn–Hamieh–Lin, "Quartic Gauss sums over primes and metaplectic theta functions", to appear in ANT.** ([arXiv:2306.11875](https://arxiv.org/abs/2306.11875)) [V] §8–9 treat the quartic Kubota Dirichlet series with level structure, a Voronoi formula, and average Type-I sums.
- **Yoshimoto (1987); Patterson (1977, 1978).** [S, via DR] Functional equations for Σ g̃(c)ψ(c)N(c)^{−s} with ψ³ nonprincipal.
- **Hamdar, "Hecke L-functions away from the central line", [arXiv:2509.06152](https://arxiv.org/abs/2509.06152).** [V, abstract] Uses infinitely many residues of twisted cubic Kubota series, i.e. cubic theta coefficients, in the dual sum. It notes this "is only known for cubic Gauss sums".
- **Sextic: nothing found.** [I] No paper gives functional equations or Voronoi for θ_6 twisted by characters. The ingredients for the cubic case are Proskurin's cusp expansions for the cubic Kubota multiplier and Patterson's table. Their 6-fold analogues are the scattering data of the sextic Kubota Eisenstein series at all cusps of Γ_1(λ^k). Those analogues exist in principle, but nobody has written them out.

---

## 2. Dunn–Radziwiłł, "Bias in cubic Gauss sums: Patterson's conjecture", Ann. of Math. 200 (2024) 967–1057 ([arXiv:2109.07463](https://arxiv.org/abs/2109.07463))

### Conditional results (GRH for Hecke L-functions over Q(ω))

- [V] Thm 1.1: Patterson's asymptotic Σ_p g̃(p) W(p/X) ~ c X^{5/6}/log X.
- [V] Thm 1.2: the bias in equidistribution.
- [V] Thm 1.3: Σ g̃(ϖ)^k = o(X^{5/6}/log X) uniformly for 1 < |k| ≤ X^{1/100}.
- [V] The **lower bound** in Thm 1.4 (sharpness of the cubic large sieve). DFDH (2026) have since made this lower bound unconditional.

### Unconditional parts

- [V] Section 4 (Poisson) and Section 5 (level-aspect Voronoi).
- [V] The upper bounds in Thm 1.4, which are Heath-Brown's.
- [V] Heath-Brown's bound Σ g̃(p) ≪ X^{5/6+ε} (Israel J. Math. 120 (2000)) is unconditional.

### Where GRH enters

- [V] **Type-I:** cancelling the contribution of cubes, pointwise for C > N(r)^{1+ε} and on average for C > R^{2−ε}.
- [V] **Type-II:** the corrected cubic large sieve, a "dispersion estimate" (Props 9.1–9.2). It subtracts the main term (2π)^{2/3}/(3Γ(2/3)) · g̃(a)N(a)^{−1/6} Σ β_b N(b)^{−1/6}. It saves (AB)^{−ρ} below (AB)^{2/3} **only for sequences in the class C_η(B,w)**: squarefree, w-rough, and square-root cancelling against every non-principal cubic character twist (Def 3.1, axiom (3.2)). GRH is used to verify (3.2) for rough-number indicators.
- [V] **Type-III**, at A = B = C = X^{1/3}: bounded by ε × the main term via a GRH-tight large sieve (§7).

### Method

- [V] Heath-Brown/Vaughan combinatorial decomposition into Types I, II and III.
- [V] Voronoi for Type I.
- [V] Cauchy–Schwarz plus the corrected cubic large sieve for Type II.
- [V] GRH-tight large sieve for Type III.

### Sextic analogue

- Nothing found.
- [V] Patterson, PLMS 54 (1987) 193–215, proves equidistribution with a power saving for n-th order Gauss sums at primes over K ⊇ Q(ζ_n), for any n ≥ 3. [S] DDHL cite a quartic bound of X^{19/20}.
- Quartic analogue [V]: DDHL prove, unconditionally, Σ_c g̃_4(c)(c/|c|)^ℓ Λ(c) R(N(c)/X) ≪ X^{5/6+ε} + X^{3/4+ε}|ℓ|^{3/2+ε}.
  - They conjecture a main term b_β X^{3/4} for ℓ = 0.
  - Type-II uses the quadratic large sieve over Q(i) (Onodera, Goldmakher–Louvel).
  - Type-I uses Suzuki's τ_4 at squares.
- [I] By the pattern X^{1/2+1/n} (cubic 5/6, quartic 3/4), the expected sextic bias would be X^{2/3}.
- Cuspidal analogue [V]: Dunn (ANT 2025), Thm 1.2 gives Σ ρ_f(ν)Λ(λ³ν)W ≪ X^{1−1/34}. The ingredients are:
  - a new bilinear large sieve: B_f ≪ ((AB)^{1/2} + A^{3/2}B^{1/4}) ‖μ²α‖_∞ ‖β‖_2, proved with the Browning–Vishe δ-method plus Voronoi, **not Kuznetsov**;
  - a linear-form estimate through Heath-Brown's cubic large sieve, which beats level 2/3.

---

## 3. Idea 2: Kuznetsov, spectral theory and large sieves on metaplectic covers over imaginary quadratic fields

### Kuznetsov and spectral theory

- **B. Louvel, "On the distribution of cubic exponential sums", Forum Math. 26 (2014) 987–1028.** ([arXiv:1107.0197](https://arxiv.org/abs/1107.0197))
  - [V] §4, Thm 4.2: a **Bruggeman–Kuznetsov formula on the cubic metaplectic cover over Q(ω)**. It covers general unitary multipliers and essential cusps, with a Bessel transform B_h built from J_{s,p}, K_{s,p}. It is derived from Bruggeman–Motohashi (SL_2(Z[i])) and Lokvenec-Guleska (any imaginary quadratic field).
  - [V] Thm 3.4: exceptional parameters on the cubic cover lie in (0, s(d)] ∪ {1/3}, with s(d) ≤ 1/6. The proof combines Flicker's cubic Shimura correspondence (s̃ = s/3) with Gelbart–Jacquet.
  - [V] Thm 1.3 is an asymptotic for Σ_{c≡0 (d)} S(1,c)/√N(c) with explicit level dependence. The "theta term" is the pole at 4/3.
  - [I] The formula itself does not depend on the multiplier. For the sextic multiplier one would need Flicker's correspondence for 6-fold covers; [I] the same argument should put the exceptional spectrum in (0, small] ∪ {theta point}.
- **Livné–Patterson, "The first moment of cubic exponential sums", Invent. Math. 148 (2002) 79–116.** [V, via Louvel] They use the Goldfeld–Sarnak version of the sum formula, not a full Bruggeman–Kuznetsov formula. The first moment is ~ c_θ X^{5/6}, with error X^{1/4+ε}-type and constants for squarefree d.
- **P. E. Herman, "On Patterson's conjecture: sums of quartic exponential sums", [arXiv:1305.1243](https://arxiv.org/abs/1305.1243)**, and **"On the analytic properties of a cubic Dirichlet series associated to a cubic metaplectic form", [arXiv:1208.1705](https://arxiv.org/abs/1208.1705).**
  - [V, abstract] A Kuznetsov-like trace formula on the quartic cover over Q(i). A secondary main term appears in the quartic asymptotic.
  - [V, abstract] A cubic "beyond endoscopy" spectral average.
- **N. Proskurin, "Cubic metaplectic forms and theta functions", LNM 1677 (1998).** [V, via Dunn] The standard reference for cusp expansions (Thm 0.3.1), the Hecke algebra C[T_{ν³}] (§0.3.12), and Kloosterman-type evaluations.
- **Patterson, "The cubic Shimura correspondence", Asian J. Math. 2 (1998); Flicker, Invent. Math. 57 (1980).** [S] Spectral parameters of cubic cusp forms lie on 1 + iR.
- **N. Watt, "Weighted spectral large sieve inequalities for Hecke congruence subgroups of SL(2,Z[i])", Funct. Approx. 48 (2013).** ([arXiv:1302.3127](https://arxiv.org/abs/1302.3127)) [V, abstract] Covers only the **non-metaplectic** case. It is the closest spectral large sieve over an imaginary quadratic field.
- **Friedberg–Offen, [arXiv:2507.06963](https://arxiv.org/abs/2507.06963) (2025).** [V, abstract] Metaplectic Kuznetsov distributions on the cubic cover of SL_3, toward a relative trace formula. Not directly usable here.
- "Venkatesh–?": nothing relevant was found.

### Second moments and large sieves for metaplectic coefficients, and the regime with fewer rows than columns

- **Heath-Brown's cubic large sieve, and BGL for all n ≥ 3 over any K ⊇ μ_n.** [V] Ξ_n ≪ (AB)^ε(A + B + (AB)^{2/3}). Reference: Blomer–Goldmakher–Louvel, IMRN 2014, [arXiv:1112.1650](https://arxiv.org/abs/1112.1650). The 7/8 paper's sextic large sieve is this bound.
- **DFDH 2026, [arXiv:2607.07911](https://arxiv.org/abs/2607.07911).**
  - [V] Thm 1.1 (unconditional): Ξ_3 ≫ (AB)^{−ε}(A+B+(AB)^{2/3}).
  - [V] Thm 1.2 (unconditional): Ξ_4 ≫ (AB)^{−ε}(A+B+A^{3/4}B^{1/2}+A^{1/2}B^{3/4}).
  - [V] Conjecture (1.8) for general n is stated in §0.
  - [V] The extremal sequence is β_b = g̃(b)H(N(b)/B). By the polar term, the inner sum is P_B(a) ∝ g̃(a)N(a)^{−1/6}B^{5/6}, which is the theta vector.
  - [V] Their two proof ingredients are Lindelöf on average for Kubota series, and a Rankin–Selberg lower bound Σ_{N(a)≤A} |τ_n(a)|² over **squarefree** a. The latter needs a Zagier-type regularisation.
- **DR Props 9.1–9.2.** [V] Conditional power-saving *asymptotics* in the large-sieve range: Σ_a μ²(a)V |Σ_b β_b g̃(b)(a/b)_3 − main|² ≪ A^{2/3−ρ}B^{5/3−ρ} + …. This requires β ∈ C_η, which means GRH-type cancellation.
- **Duke–Iwaniec 1990, half-integral weight.** [V, via Dunn] Σ_{a≤A,b≤B} μ²(a)α_a β_b c_g(ab) ≪ (AB)^ε(B^{1/2} + AB^{1/4})‖α‖‖β‖. This uses Kloosterman sums, which are available for holomorphic forms but not Maass forms.
- **Pascadi, "Large sieve inequalities for exceptional Maass forms…", Forum Math. Pi 14 (2026).** ([arXiv:2404.04239](https://arxiv.org/abs/2404.04239)) [V, abstract] Non-metaplectic. A model for sequences with sparse Fourier transform.
- **Not found:** any Kuznetsov-based second moment for θ_3/θ_6 coefficients in the regime with fewer rows than columns, or any power-saving asymptotic there that is unconditional.

---

## 4. Idea 3: double Dirichlet series and mean values

### Double and multiple Dirichlet series

- **Friedberg–Hoffstein, Ann. of Math. 142 (1995) 385–423.** [S] Quadratic twists of GL_2 over **any** number field, via the metaplectic Eisenstein series on the double cover. Used for non-vanishing.
- **Brubaker–Friedberg–Hoffstein, "Cubic twists of GL(2) automorphic L-functions", Invent. Math. 160 (2005) 31–58.**
  - [V, via CFH survey §4.2] Over Q(√−3) there are three linked series:
    - Z_1, the cubic twists of L(s,π);
    - Z_6, in which "the square of a cubic Gauss sum introduces, via Hasse–Davenport, a conjugate 6-th order Gauss sum". Its functional equation uses the **6-fold Kubota series** Σ G^{(6)}(m,d)|d|^{−w};
    - Z_3, made of cubic twists of the Rankin–Selberg product **π × cubic theta**.
  - Remark: the group of functional equations has order 384. A G_2 cubic-cover interpretation is suggested.
  - Relevance: **the closest published precedent for μ-absorption in reverse**, i.e. cubic² → sextic.
- **Friedberg–Hoffstein–Lieman, Math. Ann. 327 (2003) 315–338.** [S] n-th order twists of Hecke L-series, paired with a series of n-th order Gauss sums (Kubota).
- **Diaconu–Goldfeld–Hoffstein, Compositio 139 (2003)**; **Diaconu, Invent. Math. 157 (2004) 635–684.** [V, abstract] Mean square of r-th order twists, via Rankin–Selberg of the r-fold metaplectic Eisenstein series with itself.
- **Chinta–Friedberg–Hoffstein, "Multiple Dirichlet series and automorphic forms" (PSPM 75, 2006).** ([pdf](https://chinta.ccny.cuny.edu/publ/cfh2.pdf))
  - [V] The survey records Σ_{|d|<x} |L(½, π⊗χ_d)| ≪ x^{5/4} for GL_3 (Chinta–Diaconu), valid only over Q because it needs Heath-Brown's large sieve.
  - [I] Goldmakher–Louvel later supplied the number-field large sieve.
- **Brubaker–Bump–Friedberg (WMDS series); Chinta–Gunnells, JAMS 2010.** [S] Weyl group MDS. The theta link is in BBFH 2012 (§1).
- **Sawin–Whitehead, [arXiv:2507.08662](https://arxiv.org/abs/2507.08662) (2025/26).** [V, abstract] A classification of axiomatic MDS whose functional equations include the Kubota/Gauss-sum type. It is a checklist for which "mixed" moments, e.g. quadratic twists of a cubic Kubota series, admit full continuation.
- **Diaconu–Ion–Paşol–Popa, [arXiv:2607.27131](https://arxiv.org/abs/2607.27131) (July 2026).**
  - [V] Gives asymptotics for the first and second twisted moments of r-th order Hecke L-functions over K ⊇ μ_{2r}, for r ≥ 3, by the MDS method. Error terms: X^{(r+A_r)/(r+1)} for the second moment and X^{(r+1+A_r)/(r+3)} for the first.
  - [V] Their first-moment secondary-term constant is "a series involving the Whittaker–Fourier coefficients of theta functions on the r-fold cover". They captured that term only for r = 3.
  - [V] For r = 3 over Q(√−3) the second moment "falls just short" of the conjectured X^{5/6} secondary term. The bottleneck is Heath-Brown's cubic large sieve, which DFDH show is sharp.
- **Gao–Zhao.** Moments and one-level density of **sextic** Hecke L-functions of Q(ω), Funct. Approx. 70 (2024), [arXiv:2201.01885](https://arxiv.org/abs/2201.01885). Also first moments for quadratic, cubic, quartic and sextic families over imaginary quadratic fields via double Dirichlet series, [arXiv:2306.10726](https://arxiv.org/abs/2306.10726). [V, abstracts] The sextic paper is already cited by the 7/8 paper.

### Lindelöf on average, and families smaller than the conductor

- [V] DdFDS Lemma 6.6 (cubic Kubota, twist aspect): Σ_{N(m)≤M} |ψ(hm, σ+it)|² ≪ N(h)^{1/2} M^{1+ε}(1+|t|)². The proof goes through Heath-Brown's cubic large sieve after the functional equation. Reference: [arXiv:2410.03048](https://arxiv.org/abs/2410.03048).
- [V] DFDH: the quartic analogue via Rankin–Selberg, "expected for all n ≥ 5".
- [V] DDHL (1.11), with family A and conductor A²: Σ μ²(a)|ψ̃^{(4)}(½+ε+it, a²)|² ≪ A^{3/2+ε} via the approximate functional equation plus the quadratic large sieve. Remark 9.7 says A^{1+ε} "seems likely".
- **Not found:** any mean-value theorem for L(½, θ_3 ⊗ χ_d), i.e. quadratic twists of the cubic theta's "L-function" Σ τ_3(m)χ_d(m)N(m)^{−s}, or for its sextic analogue. [I] By DR Prop 5.1 with ψ = χ_d, each twist individually continues with a functional equation of conductor about N(d)², and the dual is twisted by ψ̃ (sextic Gauss sums). A large-sieve second moment over d ≤ D would be Lindelöf-on-average only when D² ≲ the number of d. That is the same balanced/unbalanced issue as Goldmakher–Louvel.

---

## 5. Idea 4: other μ-absorbing families

- [V, CFH / BFH] The relation is G^{(3)}(1,m)² = G^{(6)}(1,m) by Davenport–Hasse, up to normalisation and bad primes.
- [I] **General structure.** Davenport–Hasse with p = 2 gives g(χ)g(χ^{m+1}) = χ̄(4) g(χ²) g(χ^m) for χ of order 2m.
  - Suppose m is odd. Then χ^{m+1} is a power of ψ = χ², and Stickelberger/Weil makes g(ψ)^a/g(ψ^a) a Jacobi-sum **Hecke character**. So g̃_{2m} = (Hecke character) × g̃_m(ψ^j) × g̃_2, with j ≡ (1−m)/2 mod m. That exponent is forced by matching the reciprocity cross-terms.
  - m = 3 gives j = 2, a conjugate cubic Gauss sum, which matches the 7/8 identity.
  - m = 5 (n = 10) gives a quintic Gauss sum g̃_5(ψ³). Quintic theta coefficients are unknown.
  - m even (n = 4, 8, 12) does not linearise. For n = 4 one only gets g̃_4² = (Hecke) × g̃_2. That is "square-absorption", the identity DDHL use for k ≡ 2 mod 4.
  - **So n = 6 is the unique order where Möbius × Gauss sum collapses linearly to a Gauss sum of an order whose theta is fully known.** Cubic is known: Patterson, and Kazhdan–Patterson for n = r+1 over any K ⊇ μ_3. The quadratic case is classical. The quartic case is known only at squares (Suzuki), which suits square-absorption rather than linear absorption.
- **Higher fields ⊇ μ_6, e.g. Q(ζ_12), Q(ζ_9).** [I]
  - The identity and all imported tools exist with the same exponents:
    - KP cubic theta over any K ⊇ μ_3;
    - the BGL n-th order large sieve, A+B+(AB)^{2/3};
    - the Goldmakher–Louvel quadratic large sieve;
    - the DFDH conjecture, stated for any K ⊇ Q(ζ_n).
  - The obstruction exponents (P^{1/6}, (QN)^{2/3}, Q²/Z) are expressed in norms and do not change.
  - A higher degree adds archimedean places and Hecke characters with infinity type, plus higher-dimensional Poisson steps.
  - I see **no mechanism by which a larger field improves θ**.
  - The transfer to ζ is unaffected, since ζ divides ζ_K.

---

## 6. Follow-ups to the 7/8 paper (30 Sep–9 Oct 2026)

- **Baiying Liu, "Slightly improved zero-free half-planes for the quasi-Riemann hypothesis", [arXiv:2610.12234](https://arxiv.org/abs/2610.12234) (8 Oct 2026).** ([GitHub](https://github.com/liubaiying101/Slightly-improved-zero-free-half-planes-for-the-quasi-Riemann-hypothesis))
  - [V] Varies only b and ℓ.
  - First result: 34999/40000 = 0.874975, at b = 1/8, ℓ = 1/6 + 10⁻⁴.
  - Second result: **(1507 − 2√921)/1653 = 0.874957069799…**, at b = −4/29 + 230√921/26709, ℓ = (33+8√921)/1653.
  - Lean-formalised.
  - [I] This sits 8·10⁻⁸ above our LP barrier of 33216105/37963129 = 0.8749569879 (QRH-analysis.md).
- **Li Xiang, [lixiang90/RH-Zero-Free-Formalization](https://github.com/lixiang90/RH-Zero-Free-Formalization) (created 9 Oct 2026).**
  - [V, README] Claims σ* = 11/12 − e*/4 = **0.874957019420…**, where 657e³ − 954e² + 21e + 20 = 0. That is between Liu and our barrier. It comes from a "kappa-feedback cubic boundary" paper.
  - [V, README] The proof is **conditional**: "the complete physical probe and its low/raw-high estimates remain open". Lean checks only the algebraic certificate and the reductions.
- **S. Kintali, "A Short Proof of the Quasi-Riemann Hypothesis" (7 Oct 2026, [pdf](https://shivakintali.github.io/papers/QRH.pdf)).**
  - [V] Proves the weaker **47/48**.
  - Inputs: the probe framework of the 7/8 paper, DR §5 and App. A (unconditional cubic theta formulas), Goldmakher–Louvel Thm 1.1/Cor 1.2, and a Hecke zero-density lemma.
  - It explicitly does **not** use DR's conditional prime-sum or dispersion estimates.
  - No Gram-term innovation.
- **[tomoto0/quasi-riemann-hypothesis-7-8-verification](https://github.com/tomoto0/quasi-riemann-hypothesis-7-8-verification) (7 Oct 2026).**
  - [V, README and reports] Rebuilds the Lean proof; it finds 0 errors.
  - Heuristic model figures (all "optimistic"):
    - a geometry-only optimum of 0.87496;
    - a "bootstrap" with α = 3/4, i.e. 7/8 fed back into the row count, giving 0.87451, and iteration converging to about 0.8745;
    - a density-hypothesis-strength row count giving 0.86979;
    - a floor of 13/15.
  - The bootstrap needs `lem:bin-contour-accounting`, `lem:local-euler` (currently only Re x ≥ 7/8) and `lem:plain` (currently only κ ≥ 3/4) to be re-proved below 7/8. **Unproved.**
- **Others.** The OpenAI repo was updated 8 Oct ("Update manuscripts and Lean formalizations"); I did not diff it. arXiv searches for "quasi-Riemann", "zero-free half-plane" and "sextic Gauss/theta" since 30 Sep return only Liu.
- **Nothing found below 0.8749 rigorously, and nothing touching the P_a^{1/6} Gram term.**

---

## Reference URLs (consolidated)

- DR (Annals 2024): https://arxiv.org/abs/2109.07463
- de Faveri–Dunn–Hoffstein (2026): https://arxiv.org/abs/2607.07911
- David–de Faveri–Dunn–Stucky: https://arxiv.org/abs/2410.03048
- Castillo–de Faveri–Dunn (quartic): https://arxiv.org/abs/2604.01316
- David–Dunn–Hamieh–Lin (quartic): https://arxiv.org/abs/2306.11875
- Dunn (ANT 2025): https://arxiv.org/abs/2403.13151
- Bröker–Hoffstein (sextic theta): https://arxiv.org/abs/1312.0568
- CFH, double Dirichlet series and theta: https://chinta.ccny.cuny.edu/publ/patterson.pdf
- CFH survey: https://chinta.ccny.cuny.edu/publ/cfh2.pdf
- Friedberg–Ginzburg descent: https://arxiv.org/abs/1403.3930
- Friedberg–Ginzburg global integrals: https://arxiv.org/abs/1403.3929
- Gao (distinguished theta): https://arxiv.org/abs/1602.01880
- Patterson, F_q(t) theta: https://arxiv.org/abs/1411.7594
- Louvel (Forum Math 2014): https://arxiv.org/abs/1107.0197
- Herman: https://arxiv.org/abs/1305.1243, https://arxiv.org/abs/1208.1705
- Watt: https://arxiv.org/abs/1302.3127
- BGL: https://arxiv.org/abs/1112.1650
- Goldmakher–Louvel: https://arxiv.org/abs/1112.1642
- Baier–Young (cubic and sextic large sieve over Q): https://arxiv.org/abs/0804.2233
- Gao–Zhao (quartic large sieve): https://arxiv.org/abs/1101.3646
- Gao–Zhao sextic: https://arxiv.org/abs/2201.01885
- Gao–Zhao first moments: https://arxiv.org/abs/2306.10726
- Diaconu–Ion–Paşol–Popa: https://arxiv.org/abs/2607.27131
- Sawin–Whitehead: https://arxiv.org/abs/2507.08662
- Hamdar: https://arxiv.org/abs/2509.06152
- Liu: https://arxiv.org/abs/2610.12234
- Kintali: https://shivakintali.github.io/papers/QRH.pdf
- Li Xiang: https://github.com/lixiang90/RH-Zero-Free-Formalization
- tomoto0: https://github.com/tomoto0/quasi-riemann-hypothesis-7-8-verification
