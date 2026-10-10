# theta-pairing: the Gram exceptional term P_a^{1/6} and the pairing T = ⟨A^sig, B⟩

Agent `theta-pairing`, 10 Oct 2026. Direction: get a power saving on
T(Q,Z) = Σ_m K(N m/Q) ρ(m) B_m(Z) (P78 Prop. `probe-gram`, l.8340; theta split of `04-gram-split`).
Labels: **[P]** proved here (given the cited lemmas), **[R-]** routine/derivation not fully written,
**[N]** numerical, **[M]** model computation (my copy `model_ext.py` of r4 `model.py`), **[H]** heuristic.

**Verdict: no saving proved. θ stays at 0.874957.** The main new finding is a structural one (§1):
at every LP optimum the theta pairing sits exactly on the energy wall, where the marked reflected
column count equals the row count Q. That kills every mean-square (large-sieve/Bessel/sparsity/restricted-energy)
route at the paper's geometry, and it quantifies how strong an off-wall input must be (§2).

Files: `model_ext.py` (knobs `gram_s`, `gram_sQ`, `gram_ls`), `model_summary.txt` (all runs),
`ext_*.log`, `model_gram_b.log`, `kp_forced_zero_density.py`.

---

## 0. Checks of the three observations in the task

(i) **[R-, from the structure of the exceptional-frequency coefficient in P78's proof of Prop. probe-gram]** For exceptional k = a⁶e the joint coefficient in the two residual
columns is a fixed finite sum of products ψ(m₁)ψ'(m₂)·(zero masks) over a fixed ray group. Hence, for *any* bounded
s-weights c (ray indicators, slot-marked s, primes, …), the exceptional excess is a fixed finite-rank PSD Hermitian form
H(c) = Σ_C L_C Σ_{ψ,ψ'} a_{ψψ'}⟨c(C·),ψ⟩conj⟨c(C·),ψ'⟩, while the signal is a linear functional of c. So
|signal|²/excess ≤ const(ray group) for every c. Ratio "invariant up to a constant", never a power. The constant can
be improved only by projecting c onto the principal isotypic part; exponents do not move.

(ii) **[P, = Lemma C2 + Prop. C3 of vls/C]** With ŴT = Y^{1/3}f + T_E: removing T′ := T_E/Ŵ (by β* or anything else)
leaves the comparison "f against f". The theta-dominated low constraint carries information only through the
signal-blind bound |S + T′| ≤ ‖ρ_Q‖‖B‖; T′ is the Y^{-1/3} part of the high-side error E, which the paper already controls
with its row counts. The a-priori half-plane β* enters only through T′, i.e. it is already used. LP gain: 0.

(iii) **[P]** Theta part of the low bound = Q^{-1/2}Y^{-1/3}‖ρ_Q‖‖B‖ = Y^{-1/3}Q^{5/12}; signal = Q^{1/3}Y^{-1/3}Z^{β-5/6}P^{1/6}.
Constraint σ ≥ 5/6 + M/12 − ℓ/6. On the wall M+ℓ = 1 this is **σ ≥ 2/3 + M/4**; at the paper's optimum
M = 0.83317 gives 0.874959 = θ. 27/31 would need M ≤ 0.8172 on the wall.

## 1. The wall observation (new) [P for the identity, M for the optimum]

Marked reflected row (P78 reflection, κ = 1 marks): dual length L^{mk} = Q²P/Z = Q·Z^{M+ℓ−1}. The paper's
energy wall e = M+ℓ−1 is exactly the column term of the quadratic large sieve (Q + L^{mk}).

**Every LP optimum has M+ℓ = 1.0000** (paper; Gram halved; Gram removed; any uniform saving s ≤ 0.0127 — see
`model_summary.txt`). So at the binding geometry **L^{mk} = Q**: T^{mk} = Σ_{col≍Q} β_col F′(col) is a square Q×Q form.
Consequences:

1. **Mean-square routes are void at the wall.** Any bound of the form CS-over-columns + Σ_col|F′(col)|² is ≥
   (L^{mk}‖π‖²)^{1/2} = Q^{11/12} if Σ|F′|² is at its true (random) size, which `vls/E` measured (R_rand = 1.00±0.05 up to L = Q) [N].
   (V′_LS) at L = Q is *false*, not just unproved.
2. **Modulus-blind barrier [P, given flatness of |Kβ′|, which is the generic case]**. For the square kernel, take β′ with |β′| = |β| and π′(m) = |π(m)|·phase((Kβ′)(m)):
   |⟨π′,Kβ′⟩| = Σ|π||Kβ′| ≍ Q^{11/12} when |Kβ′| is flat. So no inequality that sees only |π|, |β| and the kernel beats CS
   at the wall (this extends `vls/B`'s "β-structure used up" to both sides).
3. **Sparsity / weighted CS / restricted energy cannot help at the wall [P].** For any row set S with |S| ≤ L,
   ‖K_S‖²_op ≥ ‖K_S‖²_F/|S| = L = Q, so the restricted large sieve gives no gain however sparse supp ρ is.
   Independently, ρ is not sparse: the KP-forced zeros (v_p ≡ 5, or v_p ≡ 2 with G₃ ≠ 1) cover
   ≤ 5.95% of rows (`kp_forced_zero_density.py`), i.e. a CS gain ≤ 3%; Bröker–Hoffstein's prime values
   τ(π)² = ζ₆^k·(cubic Gauss sum)·F(π), F in a finite nonzero set, are nonvanishing [S]. Under |τ(m)| ≪ N(m)^ε,
   ‖ρ‖₁ ≥ ‖ρ‖₂²/‖ρ‖_∞ ≫ Q^{11/12−ε} = Q^{1/2}‖ρ‖₂·Q^{−ε}: the ℓ¹ question has answer **no**.

So a saving at the paper's geometry must use the *phases* of π = ρζ (automorphy of θ₆ and θ₃): a first-moment /
bilinear (V′_β) input. The only first-moment structure known is the A₂/4-node series, whose four FE images
(vls/A) are: CS-level, worse, worse, and one (m-Voronoi) containing the signal — circular (vls/C, C3).

## 2. What an off-wall mean-square input would have to deliver [M]

Knob `gram_ls = η`: saving η(1−M′−ℓ′)/2 (Z-units) = a fraction η of the random rate (L^{mk}/Q)^{1/2}; it vanishes on
the wall, so the LP must leave the wall to use it.

| η (fraction of random-rate (V′_LS)) | 0.10 | 0.15 | 0.16 | 0.20 | 0.25 | 0.277 | 0.30 | 0.40 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|
| θ | 0.874957 | 0.874957 | 0.874823 | 0.873583 | 0.872159 | ≈27/31 | 0.870182 | 0.867776 | 0.863101 |
| M+ℓ at optimum | 1 | 1 | 0.935 | 0.920 | 0.862 | ~0.85 | 0.870 | 0.886 | 0.953 |

**No gain below η₀ ≈ 0.155; 27/31 needs η ≈ 0.277.** (P²/Y kept; dropping it changes only η = 1: 0.8581.)
So "any s₀ > 0 beats Liu" (04 REPORT) holds only for a saving that survives at L = Q; a large-sieve-type saving must
reach ~16% of square-root cancellation over a beyond-the-edge family before it does anything.

Uniform saving `gram_s = s` (valid at the wall, i.e. a genuine (V′_β)/first-moment saving), reproduced:
s = 0.002 → 0.873332; 0.004 → 0.871707; **0.005 → 0.870895 (< 27/31; threshold s ≈ 0.0049)**; 0.01 → 0.866832;
0.0127 → 0.864692. Gram ×0.5 → 0.869495, removed → 0.859796 (`model_gram_b.log`).

## 3. Angles tried and the precise obstruction for each

1. **Period / Petersson / spectral projection.** Pairing a single residual vector θ₆ against anything is Bessel for one
   vector = CS. Writing T = ⟨θ₆, P_B⟩ (Poincaré series at height Q^{−1/2}) gives |T| ≤ ‖θ₆‖‖P_B‖ and ‖P_B‖² is the full
   spectral large sieve, which is no better. The "projection of the θ₃-built object onto θ₆" is the z = 5/6 polar
   coefficient of the double series Ξ(z,u) = Σ_m ρ(m)N(m)^{−z}L_m(u), whose residue is H(u)/L(u+5/6) (vls/C Lemma C1):
   **the projection is the signal itself**. B is not the coefficient vector of an automorphic form on the 6-fold cover:
   there is no z = 11/12 pole, so there is no "orthogonal-projection" main term to remove. [P/R-]
2. **β* on the non-principal part.** Vacuous, see §0(ii). A genuinely non-circular use of β* would be through the
   *row* Möbius factor of π = (ray)·μ(m)ᾱ(m)N^{1/6}ρ(m³) (07/independent_low §7). Two obstructions: ρ(m³) has a
   patternless ζ₆^{k(p)} factor (BH), so no L-function; and ᾱ is infinite order, outside P78's β*.
   *Side computation [H]*: in the fictitious model π = μ·(finite-order Hecke)·N^{−1/12}, mean square over columns plus a
   log-free density N(σ) ≪ L^{c(1−σ)} near σ = β* needs cλ < 2. That replaces 07's λ < 2(1−B). Heath-Brown's quadratic
   density c = 3/(2−σ) = 2.67 fails narrowly at λ = 0.79 and works for λ < 0.75. It is moot at the slotted wall (λ = 1).
3. **Exploiting KP / s-side amplification.** For c = c₀·1_{p|s}, twisted multiplicativity and reciprocity give the
   theta part g̃(p,−m)N(p)^{−2/3}ρ(mp⁴); with KP τ(mp⁴) = χ_p(m)γ₅(p)τ(m) (BH orientation) this is (const_p)·ρ·1_{p∤m}.
   It is collinear with ρ, with the signal's phase. So s-side marks cannot amplify (cf. §0(i), which does not need the orientation).
   B-side marks amplify because the B^{(p)} ∝ χ_p(m)B_{mp^j} are orthogonal while ρ is fixed. Their gain is
   P^{1−5/6} = P^{1/6}, i.e. the theta pole 5/6 against energy exponent 1, and it is capped by the wall M+ℓ ≤ 1. [R-]
4. **Two geometries.** T depends only on (Q,Z), but signal and theta-CS both scale as Y^{−1/3}, so Y-combinations kill
   both (known). A combination over Q-scales has signal/CS ≤ max_i Q_i^{−1/12}·O(1) (orthogonal ρ_{Q_i}). Bounding T′ at one Y and the
   theta-CS at another reproduces the paper's single-geometry condition. Z-window families (orthogonal B, coherent signal)
   gain only log Z, because signal ∝ energy per window. [P]
5. **Standalone theta detector** (bound T′ by evaluation, use only CS on T). T′ is the Voronoi dual of T at
   Q^∨ = (ZP)²/Q > ZP. Its only bound is trivial (≫ Q^{11/12}) or the high side. At M < 1/2, T ≈ 0 forces T′ ≈ −S. No gain.

## 4. Precise obstruction (summary)

At the binding geometry of P78 Part II (and of every improved LP tested), the marked reflected family has
L^{mk} = Q. The theta pairing is therefore a square bilinear form ⟨π, Kβ⟩ at the edge of the quadratic large sieve.
There the following are all capped at Q^{11/12}: (a) any modulus-only inequality (§1.2); (b) any mean-square/large-sieve
input including the true random rate (§1.1); (c) any sparsity or restricted-energy argument (§1.3).
A uniform saving s ≥ 0.0049 (Z-units), which is what 27/31 needs, must come from the joint phases of
π = ρ·conj(γ₂)⁴ and of Patterson's cubic theta columns, i.e. from a first moment of the 4-node A₂-type series.
Evaluating that series exposes the Möbius signal (circular), and its FE images give nothing better than CS. Leaving the
wall to use a mean-square input requires ≥ 15.5% of square-root cancellation in (V′_LS) for any gain and 27.7% for 27/31.
That is a beyond-the-edge large sieve for a family with conductor ~L³ (vls/C), for which no tool exists.
