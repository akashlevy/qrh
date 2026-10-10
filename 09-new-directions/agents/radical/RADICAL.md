# RADICAL: triage of radically different routes to a zero-free half-plane (agent `radical`, 10 Oct 2026)

Files in this directory:
- `RADICAL.md` (this file);
- `dictionary.py` → `runs/dictionary.out`: exact numbers for the equivalent reformulations at 7/8, the record, 27/31 and 6/7, plus the bootstrap threshold;
- `model_radical.py` (a copy of `agents/hls/model_hls.py`), `run_ceiling.sh` → `runs/ceiling*.jsonl`: the "truth ceiling" of the P78 architecture;
- `ext/qrh-debruijn-newman/`: a clone of the de Bruijn–Newman spinoff (read-only, for §2.8).

Status labels:
- **PROVED**: a short argument written here or standard textbook;
- **STD**: known in the literature, with the reference given;
- **HEUR**: heuristic or exponent-level;
- **MODEL**: the r4/hls exponent model.

The record to beat is 0.874887 (tracker PR #9, pending). The bare target is 27/31 = 0.870968.

---

## 0. Bottom line

1. **No radically different strategy gives any half-plane θ < 1, even 0.999.** Of the 24 strategies triaged (§1), every one falls into one of five classes.

   | class | strategies | why it fails |
   |---|---|---|
   | (U) unitary-axis-locked | Eisenstein / Maass–Selberg / Sarnak, Langlands–Shahidi, Kuznetsov, transfer operators | Positivity or self-adjointness lives on the unitary axis, which maps to Re s = 1 for 1/L. Off the axis the inequalities point the wrong way (§2.2). |
   | (B) Beurling-blind | pretentious, Turán, elementary and sieve methods, Landau positivity | They use only the Euler product plus integer counts. Diamond–Montgomery–Vorhauer Beurling systems satisfy all of that with zeros on σ = 1 − a/log t. |
   | (H) bounded-height | restricted Weil positivity, Li coefficients, de Bruijn–Newman | The criteria grade zeros by (β − ½)/γ², not by β, so any partial version sees only bounded height. |
   | (E) exact reformulations | horocycle rates, Nyman–Beurling L^p, resonance gap, Lax–Phillips decay, de Branges/Hermite–Biehler | No partial result exists beyond the trivial endpoint except those *deduced from* P78. |
   | (M) measure-theoretic or model-only | probabilistic/RMT, function fields | Exceptional sets of density zero are allowed, or the statement is true only in a model that has Frobenius. |

2. **Families bootstrap.** The bootstrap over the family *is* P78's logical form (`lem:continuation-criterion`, the sup β* over all Hecke characters). The extra information it supplies is self-Lindelöf on Re s > β*, and it is provably inert below (1+β*)/2 ≈ 0.937 (§3.A).

3. **Truth ceiling of the only viable mechanism.** P78-type μ-absorption with every random-matrix-true input (H1, H2, GRH rows, refined case A) is capped at **6/7 exactly** (MODEL, §3.C). This follows from the law θ ≥ 1 − h/6, where 1/6 is the copy density of the sextic family (copies u = a⁶). The cap is 19/22 with (H1)+(H2) alone, and 0.85955 with DH counts. A radical input that only supplies P78's open mean values can therefore reach 6/7 at best. Beating 0.874887 still needs one of the known frontier levers.

4. **Necessary features of any half-plane proof.** These can be read off the counterexamples:
   - the Euler product (else Davenport–Heilbronn);
   - structure beyond multiplicativity plus counting (else Beurling/DMV);
   - an object with positivity or analytic control at a point *off* the unitary axis. On covers this is the residual theta at s′ = 1 + 1/n; nothing else is known;
   - a replication (copy) mechanism, because density bounds never exclude a single zero. Copies come only from imprimitive characters, at density 1/n.

   P78 is, as far as we can tell, the unique known object with all four.

---

## 1. Triage table

"Half-plane?" means: can the method alone give *any* fixed θ < 1 for all Dirichlet / Hecke L-functions?

| # | strategy | half-plane? | best plausible θ | first precise obstruction |
|---|---|---|---|---|
| 1 | Eisenstein scattering φ(s) = ξ(2s−1)/ξ(2s); Maass–Selberg; Sarnak's spectral non-vanishing (also Gelbart–Lapid–Sarnak lower bounds) | **No** | 1 (region 1 − c/log^A t) | Positivity of ‖Λ^T E(s)‖² forbids zeros of the denominator only on Re s = ½, i.e. Re 2s = 1. For Re s < ½ the Maass–Selberg identity gives a *lower* bound on \|φ(s)\| (§2.2), which allows poles. Zeros of ζ(2s) *are* the poles of E there. |
| 2 | Langlands–Shahidi, including residual/Speh inducing data and Kubota scattering on n-fold covers | **No** | 1 | Gindikin–Karpelevich denominators are ξ(⟨λ,α∨⟩+1). The unitary axis ⟨λ,α∨⟩ ∈ iℝ gives Re = 1. On covers the constant term ∝ ζ_F(ns′−n)/ζ_F(ns′−n+1) is again Re 1 on Re s′ = 1 (cB-theory §1). Telescoping residual data leaves the top denominator at Re ≥ 1 (HEUR). |
| 3 | Resonance-free strip / essential spectral gap / FUP (Bourgain–Dyatlov), Lax–Phillips local-energy decay | Equivalent | — | A gap g on SL₂(ℤ)\ℍ is equivalent to θ = 1 − 2g. The record needs g ≥ 0.06256. Every gap theorem needs a porous trapped set (δ < 1, infinite area) or pressure P(½) < 0. In finite area the trapped set is everything and P(½) = ½ > 0. Under Phillips–Sarnak/Wolpert deformations, dissolved cusp forms give resonances at distance ∝ ε² → 0, so no geometric gap exists and the gap must be arithmetic. |
| 4 | Kuznetsov / Kloosterman zeta; Selberg's 3/16 mechanism | **No** | 1 | Z(m,0;s) = σ_{1−2s}(m)/ζ(2s). Selberg's 3/16 comes from Weil's bound making Z(m,n;s) converge absolutely beyond the spectral line. For n = 0 the sums are Ramanujan sums, \|c_c(m)\| ≤ (c,m), with no Weil-type saving. Absolute convergence then gives only Re 2s > 1, and further cancellation in c *is* Möbius cancellation. |
| 5 | Fourier / constant-term coefficients of metaplectic Eisenstein series; scattering on covers with 1/L(Hecke) | Only via the theta residue, i.e. P78 | P78's 7/8; truth ceiling 6/7 | On an n-fold cover the zeros of ζ_F sit as resonances at distance (1−β)/n left of the axis. The only non-tempered object off the axis is the theta pole at distance 1/n to the right; its Gauss-sum coefficients are what μ-absorption uses. Other splits are ≥ 10/11 (rounds 3–4). |
| 6 | Weil positivity: restricted support (Yoshida; Connes–Consani), trace-formula versions | **No** | 1 | Support length L resolves zeros only up to height ≈ e^{cL}. An off-line zero shows up only when L ≳ log γ, so any fixed L means bounded height and all L means RH. A θ-graded positivity exists (ξ(s+θ−½) is Hermite–Biehler ⇔ QRH(θ), Lagarias-type) but positivity does not survive analytic continuation, so it has no partial form. |
| 7 | Li / Keiper coefficients | **No** | 1 | The partner zero 1−ρ̄ of an off-line zero contributes \|1−1/(1−ρ̄)\|^n ≈ exp(n(2β−1)/(2γ²)). Zeros are graded by (2β−1)/γ², not β. QRH(θ) does not even bound λ_n polynomially. |
| 8 | de Bruijn–Newman, including the repo `DrMurphyIsIn/qrh-debruijn-newman` | **No** (one-way) | — | The repo is a Lean composition of OAI's 7/8 theorem with de Bruijn's strip theorem: θ ⇒ Λ ≤ 2(θ−½)², giving 9/32, and 0.28108 at the record. Nothing flows back. Ki–Kim–Lee: for every t > 0, all but finitely many zeros of H_t are real, so Λ-bounds see bounded height only. An off-line zero (β, γ) is pushed to ℝ in time ≈ 4(β−½)/log γ (HEUR). |
| 9 | Nyman–Beurling / Báez-Duarte; Beurling's L^p theorem (ζ ≠ 0 on Re s > 1/p ⇔ B_p = L^p(0,1)) | Equivalent | — | QRH(θ) ⇔ density in L^{1/θ}. The record needs p ≥ 1.14301. Hardy-space point evaluation gives d_N(θ)² ≳ (β−θ)/\|ρ\|², so the obstruction lives at bounded height. On a fixed window, d_N → 0 is equivalent to the Möbius series converging on Re s = θ. Approximants are mollifiers; the known rates (BBLS) are 1/log N. |
| 10 | Landau's theorem (nonnegative Dirichlet coefficients) | **No** | 1 | Landau sees only *real* singularities, which suits Siegel zeros. A nonnegative series encoding μ has positive density of support, so its abscissa is 1. |
| 11 | Turán power sums | **No** | 1 | A zero gives a large prime power sum. To contradict it one needs Σ Λ(n)χ(n)n^{−it} ≪ x^θ uniformly, which is the conclusion. In families it gives density near σ = 1, worse than Ingham at σ ≈ 0.65. It also falls to the Beurling counterexamples. |
| 12 | Pretentious methods (Halász, Granville–Soundararajan, Koukoulopoulos) | **No** | 1 | Savings are only (log x)^{−c}, or at best dlVP. DMV Beurling systems with N(x) = cx + O(x^θ), θ ∈ (½,1), have infinitely many zeros on σ = 1 − a/log t. |
| 13 | Sieve, Vaughan/Heath-Brown identities, Type I/II for one character | **No** | 1 | For a single character a Type II sum factorises, so no cancellation is available. Family averaging (BV, dispersion) gives only GRH-on-average. Parity also blocks. |
| 14 | Exponential sums (Vinogradov–Korobov); Lindelöf | **No** | 1 | Even LH gives only N(σ,T) ≪ T^ε for σ > ¾ (Halász–Turán), and never excludes a zero. |
| 15 | Zero density (Ingham, Huxley, Guth–Maynard) for individual L-functions | **No** | 1 | Counts cannot exclude a zero. The coefficient-blind / copy obstruction is OBSTRUCTIONS §3. |
| 16 | Families bootstrap: β* = sup, plus density near β* and rigidity | Already P78's logic | no gain | See §3.A. Self-Lindelöf on Re s > β* acts only above (1+β*)/2. The R2 consistency barrier applies. Deuring–Heilbronn repulsion needs the pole at s = 1. Replication comes only from imprimitive copies (density 1/6, already used). |
| 17 | Frobenius / interpolation determinant (OpenAI Oct 1 Landau–Siegel preprint: θ_n^p ≡ σ(θ_n) mod p) | **No** | 1 | It needs log-mass bias Σ_{χ(p)=1} log p/p ≪ ℓ + δ(log X)²/ℓ with δ = (1−β)log q = O(1). A zero at fixed β < 1 has δ ≍ log q, so there is no bias. Complex zeros have no algebraic Frobenius (n^{iγ} is not Galois). |
| 18 | Horocycle equidistribution with rates (Zagier, Sarnak, Strömbergsson); Farey (Franel–Landau); Möbius disjointness with rates | Equivalent | — | Error O(y^a) for all f ⇔ θ = 2(1−a). The record needs a ≥ 0.56256. Dynamics and mixing see only the L² spectrum, so a = ½ is stuck; beyond ½ means continuing E(z,s) into Re s < ½, which is #3. P78 now *gives* a = 9/16. |
| 19 | Transfer operators / dynamical zeta (Mayer: Z_{PSL₂ℤ}(s) = det(1−L_s²)); Dolgopyat strips | **No** | 1 | L_s is defined only for Re s > ½, while the ζ-resonances sit at Re s = β/2 ≤ ½. Dolgopyat strips near σ = 1 control the Laplace spectral gap, not ζ. |
| 20 | Function-field transfer (Weil/Deligne, Stepanov, tensor-power trick, BDPW homological stability, Sawin–Whitehead MDS) | Model only | 19/22 → 6/7 if P78 is transplanted (§3.C) | Over ℚ there is no Frobenius or cohomology. The tensor-power trick acts on Satake parameters, not zeros. Stepanov needs q ≫ g², which fails in the large-conductor regime. Horizontal (varying-modulus) families are not geometric. |
| 21 | Probabilistic / RMT (random multiplicative functions, Bohr–Jessen, CFKRS / ratios) | **No** | 1 | Measure-theoretic: density-zero t-sets are invisible. Derandomisation is the statement itself. |
| 22 | Operator theory: Hilbert–Pólya with a bounded skew part (zeros of ξ(½+iz) = spec(S+iB) ⇒ \|Im\| ≤ ‖B‖) | Graded, but no construction | — | No operator is known. Lax–Phillips realises the resonances as the spectrum of a contraction semigroup generator, but that is #3 (decay rate (1−θ)/2 needed). |
| 23 | Axiomatic Selberg class (Euler product + FE) | Unknown | — | No method is known. The Davenport–Heilbronn and Beurling counterexamples show that both axioms are needed, and degree 1 is just Dirichlet L. |
| 24 | Combinations with the metaplectic probe: spectral theory on the cubic cover (Louvel, Dunn) for W/(H1)/(H2)/(HLS); a β* feed-in; Turán in the Mellin step; Heath-Brown–Patterson / Livné–Patterson prime equidistribution | Inside P78 only | ≤ 19/22 (spectral), 6/7 (all truth) | There is no level-aspect metaplectic Kuznetsov. The family is a thin quadratic-twist family, and (H2) is a freeness statement (h2, hls). β* is inert (§3.A). Mellin continuation is already optimal, so Turán buys nothing. Prime equidistribution of Gauss sums is first-moment, while the levers are second-moment and thin. |

---

## 2. Notes behind the table (key computations)

### 2.1 Beurling and Davenport–Heilbronn meta-obstructions (STD)

**Euler product without FE.** Diamond–Montgomery–Vorhauer (Math. Ann. 334, 2006) build, for each θ ∈ (½,1), a Beurling system with N(x) = cx + O(x^θ) whose zeta has infinitely many zeros on σ = 1 − a/log t and none to its right.
- So any argument using only the Euler product, |coefficients| ≤ 1 and a power-saving integer count gives at best dlVP.
- This kills pretentious, Turán, elementary/sieve, Landau-positivity and Hadamard–dlVP-kernel methods *as standalone half-plane tools*.

**FE without Euler product.** Davenport–Heilbronn-type functions have FEs and zeros in Re s > 1.

So a proof must couple multiplicativity with additive or automorphic structure. P78 does this with Gauss sums plus reciprocity plus theta.

### 2.2 Why Maass–Selberg / Sarnak cannot leave the line (PROVED, elementary)

For SL₂(ℤ), with s = σ + it, s₁ = s, s₂ = s̄, the Maass–Selberg relation gives

  ‖Λ^T E(·,s)‖² = (T^{2σ−1} − |φ(s)|² T^{1−2σ})/(2σ−1) + (φ(s̄)T^{2it} − φ(s)T^{−2it})/(2it) ≥ 0.

- **At σ = ½** the first term becomes 2 log T − φ′/φ(½+it). Holomorphy of E on the unitary axis, from the spectral theorem, then controls the non-constant coefficients ∝ 1/ξ(1+2it). That gives ζ(1+2it) ≠ 0 and 1/ζ(1+it) ≪ log^A t (Sarnak 2004).
- **For σ < ½** the first term is (|φ|²T^{1−2σ} − T^{2σ−1})/(1−2σ). Positivity therefore yields only a **lower** bound |φ(s)|² ≳ T^{4σ−2}·(1 − O(T^{1−2σ}/|t|)). Hence |ξ(2s)| ≲ |ξ(2s−1)|·T^{1−2σ}, which is harmless.
- A zero of ξ(2s) is a **pole** of φ and of E. That is perfectly compatible with ‖Λ^T E‖² ≥ 0.

So the positivity is one-sided. On the unitary axis it is supplemented by holomorphy (self-adjointness). Off the axis, holomorphy is precisely the unknown. Langlands–Shahidi (#2) inherits this through Gindikin–Karpelevich denominators L(1+⟨λ,α∨⟩).

On the n-fold cover (constant term ∝ ζ_F(ns′−n)/ζ_F(ns′−n+1), unitary axis Re s′ = 1):
- the zeros of ζ_F become resonances at Re s′ = 1 − (1−β)/n;
- the theta pole sits at s′ = 1 + 1/n.

P78 is exactly "Sarnak at the theta point". Coefficients there are Gauss sums, and Jacobi/Stickelberger turns μ·Hecke into Gauss-sum triples.

### 2.3 Selberg-3/16 analogue (PROVED)

Selberg's λ₁ ≥ 3/16 comes from S(m,n;c) ≪ c^{½+ε} (Weil), which makes Z(m,n;s) = Σ S(m,n;c)c^{−2s} absolutely convergent for Re s > ¾.

For the frequency n = 0, which carries ζ, S(m,0;c) = c_c(m) and Z(m,0;s) = σ_{1−2s}(m)/ζ(2s). Since Σ_{c≤X}|c_c(m)| ≍_m X, absolute convergence gives only Re 2s > 1.

Any improvement is cancellation of Σ_c c_c(m)c^{−2s}, i.e. of Σ μ. There is no "Weil bound for Ramanujan sums" to exploit. On the cubic cover the analogous sums are Gauss sums of modulus N(c)^{½}, *with* a residual pole (Patterson's 5/6). That is the extra room P78 uses.

### 2.4 Bounded-height criteria (HEUR, standard mechanisms)

- **Li:** |1 − 1/ρ′|² = 1 + (2β−1)/|ρ′|² ≈ 1 + (2β−1)/γ² for ρ′ = 1 − ρ̄. A zero at height γ is invisible until n ≳ γ²/(2β−1).
- **Restricted Weil:** frequency resolution of support L is 1/L, against zero spacing 2π/log γ.
- **de Bruijn–Newman:** the backward heat flow moves an off-line zero to ℝ at speed ≈ 2πρ_z ≈ ½ log γ, so Λ ≳ 4(β−½)/log γ only. Ki–Kim–Lee: for t > 0, H_t has finitely many non-real zeros.

All three are lossy in the height direction. A half-plane is uniform in height.

### 2.5 Dictionary (`runs/dictionary.out`, exact)

| statement equivalent to QRH(θ) | 7/8 (now unconditional via P78) | record 0.874887 | 27/31 | 6/7 |
|---|---|---|---|---|
| horocycle rate a, error O(y^{a−ε}) (Zagier; a = 1 − θ/2) | 0.5625 | 0.56256 | 0.56452 | 0.57143 |
| Beurling density in L^p, p = 1/θ | 1.14286 | 1.14301 | 1.14815 | 1.16667 |
| SL₂(ℤ) resonance gap / Lax–Phillips rate (1−θ)/2 | 0.0625 | 0.06256 | 0.06452 | 0.07143 |
| ⇒ Λ ≤ 2(θ−½)² (one-way) | 0.28125 | 0.28108 | 0.27523 | 0.25510 |

- Pre-P78 endpoints for each reformulation were the trivial ones: a = ½, p = 1, gap 0.
- Beating Λ ≤ 0.2 (Platt–Trudgian) through dBN needs θ ≤ 0.81623; beating 0.22 needs θ ≤ 0.83166.

---

## 3. Deeper analysis of the three most promising candidates

### 3.A Families bootstrap with self-Lindelöf: compatible with P78, zero payoff (PROVED / HEUR)

**Setup.** P78's `lem:continuation-criterion` (paper.tex l.400) is already a family bootstrap: β* = sup Re ρ over *all* finite-order Hecke L over F. One shows the signal 1/L^S(s,η) has a uniform saving Z^{C(β*)−ε*}, which contradicts β* > σ₀. The question is whether H(β*) ("all such L are zero-free on Re s > β*") gives *extra* input.

**Lemma A1 (standard; Borel–Carathéodory plus three circles).** H(β*) ⇒ for σ ≥ β* + ε,
  |log L(σ+it,ψ)| ≪ (log q(|t|+2))^{1−c(ε)},
hence L^{±1}(σ+it,ψ) ≪ (q(|t|+2))^{ε′}, uniformly over the family.

**Key inequality A2 (Halász–Montgomery under H(β*)).** Take distinct ψ₁,…,ψ_R of conductor ≤ Q and D(ψ) = Σ_{N n∼N} b_n ψ(n) with |b_n| ≤ 1 and |D(ψ_r)| ≥ V. Duality gives
  R V² ≤ N (N + R·max_{r≠s}|Σ_n w(n/N)ψ_rψ̄_s(n)|).
Under H(β*) the off-diagonal is ≪ N^{β*+ε}Q^ε (shift the Mellin contour to Re s = β*+ε). Zero detection at σ gives V ≍ N^σ, so
  R ≪ N^{2−2σ} provided N^{1+β*−2σ} ≤ ½, i.e. **σ > (1+β*)/2**.

- At β* = 0.874887 the threshold is σ > 0.93744.
- P78's binding row zeros are at a = (1+δ)/2 ∈ [0.63, 0.70] (NOTES §26).
- Without H(β*), the off-diagonal bound N^{½}Q^{¼} needs σ > ¾.

Either way the binding range is untouched.

**Other entry points.**
- The ingredients of the counts are the 4th/6th moments and the large sieve on Re s = ½, where H(β*) says nothing.
- Under H(β*) the witnesses (inverse Möbius polynomials M_r) are short approximants of 1/L only for σ > β*.
- The low side (theta twists D_P) has no Euler product, so A1 does not apply. h2 found that GRH on the dual still misses by Z^{0.67}.
- Deuring–Heilbronn repulsion uses positivity anchored at the pole s = 1.
- Replication of the target zero across the family happens only through imprimitive copies u = a⁶, already used.

**Black-box impossibility (HEUR).** A Beurling/DMV-type family whose rightmost zeros lie exactly on Re s = β* satisfies A1 and every multiplicative input. So any self-improvement H(β*) ⇒ H(β*−ε) must use additive structure, which is P78's probe.

**Exponent payoff: 0.**

### 3.B Spectral gap via arithmetic: the resonance route, made quantitative (PROVED reductions, HEUR conclusion)

**Target.** A resonance-free strip of width g ≥ 0.06256 for SL₂(ℤ)\ℍ, or (1−θ)/n on the n-fold cover. Equivalently: exponential decay of the Lax–Phillips semigroup; equivalently a horocycle rate a ≥ 0.56256.

**Geometric tools give g = 0.**
- Pressure P(½) = ½ > 0 (the trapped set has full dimension δ = 1).
- FUP needs δ < 1.
- Phillips–Sarnak/Wolpert: for character or ε-deformations of Γ₀(N), dissolved Maass forms become resonances with Re s − ½ ≍ −ε²|L(½+it_j, f_j × …)|², which tends to 0.

So no finite-area gap holds uniformly, and the input must be Hecke-arithmetic.

**Arithmetic tools.**
- Hecke amplification of Maass–Selberg: E is a Hecke eigenfunction with λ_n(s) = Σ_{ad=n}(a/d)^{s−½}. Amplifying ‖Λ^T Σ_n x_n T_n E(s)‖² ≥ 0 reproduces §2.2 with the same one-sided sign for σ < ½. PROVED: the identity is Maass–Selberg multiplied by |Σ x_nλ_n(s)|².
- Selberg-3/16 transfer: trivial (§2.3).

The only arithmetic object giving control off the axis is the residual theta on covers (§2.2). The best plausible outcome is P78 itself: 7/8 proved, 6/7 as the truth ceiling (§3.C).

**Exponent payoff: none outside P78.**

### 3.C Function-field oracle and the truth ceiling of μ-absorption (MODEL, `runs/ceiling*.jsonl`)

**Transfer heuristic.**
- Over F_q(T) with q ≡ 1 (6), the cubic theta coefficients are explicit (Hoffstein 1992).
- P78's open levers are moments of a thin quadratic-twist family: W(w), (H1), (H2), (HLS). Their MDS has an *affine* group C̃₂ (wall-spectral).
- Over function fields such MDS are rational in q^{−x}, and all moments of quadratic families follow from homological stability for large q (Bergström–Diaconu–Petersen–Westerland, with Miller–Patzt–Petersen–Randal-Williams; Sawin–Whitehead axiomatic MDS). So the levers are plausibly theorems there.
- Over ℚ, affine MDS are expected to stop at the imaginary-root hyperplane (DGH 2003 conjecture). wall-spectral needs I = x₁+x₂+y ≤ 0.263, below the initial hull I ≥ ½.

The first obstruction is therefore **rationality/finite-dimensional cohomology of horizontal families**, which ℚ lacks.

**What the transplanted architecture would give** (r4/hls model; Gram term and Patterson bias kept, since they are genuine main terms):

| inputs assumed (all numerically true per edge-theta / h2 / hls) | θ |
|---|---|
| paper lemmas (sanity) | 0.874957 |
| (H1)+(H2) | 19/22 = 0.863636 |
| (H1)+(H2) + GRH rows, or + DH counts (paper case A, which binds) | 19/22 |
| (H1)+(H2) + refined case A + DH counts | 0.859551 |
| (H1)+(H2) + refined (or no) case A + GRH rows | **6/7 = 0.857143** (h = 6/7, the bound θ ≥ 1 − h/6 is tight) |
| GRH rows only (refined case A) | 13/15 |
| W(0.064) | 0.870955 |
| unphysical: Gram term also deleted | 16/19 = 0.8421 |

So a radical input feeding P78 can reach **6/7 at most**. The 1/6 in θ ≥ 1 − h/6 is the copy density of order-6 twist families. A larger copy density needs order < 6, which is incompatible with μ-absorption (alt-detect, twovar). To beat 0.874887 one still needs a frontier lever: W(≥0.0638), (H2) with saving ≥ L^{−0.93}, or Ingham below 7/10.

---

## 4. What remains genuinely open (honest residue)

1. **A θ-graded positivity that survives continuation.** Hermite–Biehler for ξ(s+θ−½) is equivalent to QRH(θ). Any operator realisation (#22) would be graded. Nothing is known, but the half-plane problem is "the right shape" here, unlike RH-style Li or Weil criteria.
2. **An arithmetic resonance gap on SL₂(ℤ)\ℍ or the Bianchi orbifold.** This is equivalent to the problem. The only known non-tempered off-axis objects are covering thetas, which lead back to P78.
3. **Residual spectrum of other groups.** Exceptional or minimal representations, and odd-orthogonal double-cover thetas, have divisor-type (not Gauss-sum) coefficients, so no μ-absorption (HEUR). Higher-rank Kazhdan–Patterson thetas are already worse (wall (r−1)(M+ℓ) ≤ 1).

None of these yields a concrete inequality today.

## 5. References (verified by search unless marked)

- **Beurling systems.**
  - Diamond, Montgomery, Vorhauer, *Beurling primes with large oscillation*, Math. Ann. 334 (2006).
  - Broucke, Debruyne, Vindas, *Beurling integers with RH and large oscillation* (arXiv 2004.11501).
- **Nyman–Beurling.** Beurling, PNAS 41 (1955) 312–314 (L^p version); Bercovici–Foias (p = 1).
- **Horocycles.** Zagier (horocycle ⇔ RH at y^{3/4}); Sarnak (1981); Strömbergsson (Duke 2004).
- **de Bruijn–Newman.** Ki, Kim, Lee (2009): finitely many non-real zeros for t > 0. DrMurphyIsIn/qrh-debruijn-newman (Λ ≤ 9/32, Lean; read here).
- **Function fields.** Bergström, Diaconu, Petersen, Westerland, and Miller, Patzt, Petersen, Randal-Williams (all moments over F_q(T), large q).
- **OpenAI preprints.**
  - *Uniform exclusion of Landau–Siegel zeros*, Oct 1 2026: determinant method, Lemma `lem:primebias`.
  - *Primitive roots…*, Oct 4 2026: the same metaplectic mechanism over F ⊇ μ₁₂, giving 1 − 10⁻⁶.
- **From memory, check before citing:**
  - Sarnak, *Nonvanishing of L-functions on Re s = 1* (2004);
  - Gelbart–Lapid (2006);
  - Lagarias (2006), de Branges spaces for ξ(s+a);
  - Phillips–Sarnak (1985);
  - Mayer (1991);
  - Bombieri–Lagarias (1999);
  - Yoshida (1992);
  - Connes–Consani (2021);
  - Diaconu–Goldfeld–Hoffstein (2003).
