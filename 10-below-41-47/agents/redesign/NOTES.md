# Agent `redesign`: probe redesigns that interact with the weighted numerator (10 Oct 2026)

Scripts in this directory (copies of `model/gmodel.py` and `model/wmodel.py`, plus):
- `rmodel.py`: weighted high side with a pluggable numerator bound, one function per design;
- `ctable.py`: per-bin dimension count for direction (c);
- `c_optimistic.py`: (c) with no Euler overhead and a wider b range;
- `frac.py`: the (b) upper envelope (diagonal moments of every real order) and θ reached;
- `oracle.py`: the per-bin geometry oracle for (d);
- `thetas.py`: θ reached by each design.

Labels:
- **PROVED**: a proof is written here, relative to the stated P78 lemmas;
- **MODEL**: an exact statement inside the exponent model of BRIEF.md (an adversary consistent with the proved lemmas);
- **CONDITIONAL / HEURISTIC**: as usual.

Target: ℓ = 0.1773 on the wall, B = 41/47. Worst bin (0.45, 1/2). Baseline worst E = +0.0117 (reproduced).

## 0. Verdict

**None of (a)–(d) gives any gain.** The precise obstruction is the same in all four.
- The numerator on row u is a *linear* image of the s-weight (§1).
- Provable inputs exist only for two classes of weights:
  - the fourth moment (lem:plain), for weights whose numerator is a plain smooth sum times short prime polynomials;
  - the Gram bound (prop:probe-gram), for weights in a fixed finite span of ray characters.
- Within that class the numerator is exactly the paper's: reflected length n₀ = 1 − l_y/d, moment degree 4.
- Every other weight is controlled only by the large sieve, which is worse on the high side. At b > 0 it is also worse on the low side, by up to Z^{5b/12}.
- The target-dependent and coupled versions either fail a dimension count at the worst bin, or need a spectral input that, taken alone, already closes the numerator deficit.

| design | worst E at 41/47 | status |
|---|---|---|
| paper (weighted) | +0.0117 | reproduced |
| (a1) s-mollifier / Selberg c = 1*λ, length Z^{0.005} / Z^{0.02} | +0.0145 / +0.0228 | MODEL, worse |
| (a2) c = Λ (prime s) | +0.0146 | MODEL, worse |
| (a3) squares/cubes, APs, rough, angular or short-interval supports | — | violates (H1), (H4), or has a longer dual |
| (b) numerator = S(n₁)S(n₂), n₁+n₂ = n₀ (lem:plain capacity 1−n₀) | +0.094 (cP=2/9), +0.074 (cP=1) | MODEL, much worse |
| (b) upper envelope: diagonal moments of every real order k ≤ 1/n₀ | +0.0047 (cP=2/9), +0.0019 (cP=1); θ ≥ 0.87273 | HYP, still fails |
| (c) target-dependent annihilation, optimistic (all side costs ignored) | +0.0118 (no Euler overhead: +0.0114) | MODEL, the worst bin is not annihilable |
| (d) per-bin oracle over geometries (ℓ ≥ 25/141, b) | +0.0113 | MODEL, upper bound for any coupling |
| HYP (sLS): Σ_{u∈C}\|S̃_u\|² ≪ U^ε(#C + U^{l_y/d}) | −0.029 | the input (c) would implicitly need |

**θ reached on the wall** (`thetas.py`). The grid is coarse (8 values of b, 15 values of δ, x ∈ {0, 1/4, 1/2}), so the paper's design shows 0.87505 here, against the certified 0.874917. Compare rows with each other, not with the certificate.

| design | ℓ_max | θ |
|---|---|---|
| paper (weighted) | 0.1665 | 0.87505 |
| (a1) mollifier, θ = 0.01 | 0.1626 | 0.87602 |
| (a2) Λ | 0.1663 | 0.87508 |
| (b) product, cP = 2/9 or cP = 1 | 0.0774 | 0.89731 |
| HYP: fractional diagonal moments, cP = 1 | 0.1757 | 0.87274 |
| HYP: (sLS) on the bad rows | 0.2016 | 0.86626 |

## 1. Exact high identity for an arbitrary s-weight (PROVED)

**Setup.** Replace the ray coefficient c_σ(s) = χ_s(b_*)/(τξ(s)) in eq. `general-probe` by c_σ(s)λ(s), where λ is any function with finite support in s (non-multiplicative and target-dependent weights are allowed). Write I_η[λ] for the result.

**Local factors graded by the s-valuation.** In the proof of lem:local-euler, the full coefficient K(c,n,s,a;u) of eq. `local-full-coefficient` factors over primes. This uses only the ray structure of c_σ, which is kept here. Write P_p(k) for the part of eq. `local-P` with v_p(s) = k, including its Q^{-kw}, so that P_p = Σ_k P_p(k). From eq. `local-scalar-values` and the table `local-family-table`:
- for k ≥ 2 only t = 0 survives, so P_p(k) = 1_{j=0}·ρ^k·Q^{-kw};
- P_p(1) is ρQ^{-w}·1_{j=0} plus the k = 1 rows of the table, each divided by (1−R);
- P_p(0) = 1/(1−V) plus the k = 0, t > 0 rows divided by (1−R).

Since λ multiplies the coefficient only through s, in absolute convergence

  F_{η,u}[λ](x,w,z) = Σ_s λ(s) Π_p P_p(v_p(s)) = [Π_p P_p(0)] · Σ_s λ(s) π_u(s;x,z) q_s^{-w},

where π_u is multiplicative with π_p(k) := Q^{kw}P_p(k)/P_p(0). This is **independent of w**.

**The two factors.** At j = 0, P_p(0)(1−V) = (1 − D − R/Q + VD)/(1−R). Hence

  Π_p P_p(0) = ζ^S(6z) L^S(x, ηχ̄_u)^{-1} H⁰_u(x,z),   H⁰_p = [1 + (VD − R/Q)/(1−D)]/(1−R) = 1 + O(Q^{4−6x−6z} + Q^{−x−6z}),

which is holomorphic on x_r ≥ 51/100, z_r ≥ 17/50. The multiplicative weight satisfies:
- p ∤ u: π_p(1) = (ρ + O(Q^{1−x−6z} + Q^{4−6x−6z}))/P_p(0) = χ_p(u) + η(p)p^{−x} + (smaller);
- p | u: π_p(k ≥ 2) = 0 and π_p(1) = O(Q^{1−x}) (the strict ramified term).

**The identity.** Doing the w-integral by Mellin inversion of W₁, with no contour shift since the s-sum is finite:

  I_η[λ] = (2πi)^{-2} ∫∫ X^{1/2−z} Z^{x+z−1} Φ(x+z−1) M(z) Σ_u^{(6)} q_u^{-z} ξ̄(u) ζ^S(6z) H⁰_u(x,z) L^S(x,ηχ̄_u)^{-1} N_u[λ](x,z) dz dx,

  N_u[λ](x,z) := Y^{-1} Σ_s λ(s) π_u(s;x,z) W₁(q_s/Y).

For λ ≡ 1 this is eq. `probe-poisson-identity` after its w-integral. Expanding π_u(s) = Σ_{k|s}A_u(k)·χ_{s/k}(u) recovers the 08 reconstructed numerator Y^{-1/2}Σ_k A(k) q_k^{-1/2} S_{χ_u}(Y/q_k). Here A(p) ≈ η(p)p^{-x} for p ∤ u, and the strict ramified terms are those of the 08 lemma.

**The principal row** is the finite sum N_1[λ](x,z). For λ = 1 it is Ŵ₁(1)·Res + O(Y^{-1/2+ε}). For a general λ, "keeping the pole" means |N_1[λ](ρ,1/6)| ≫ 1 at the target zero ρ.

**Consequence.** The numerator is linear in λ: N_u[λ] = ⟨λ, v_u⟩ with v_u(s) = π_u(s)W₁(q_s/Y)/Y. Everything in (a)–(d) is a statement about which λ admit provable bounds for ⟨λ, v_u⟩ on the bad rows C.

## 2. Low side for an arbitrary s-weight (PROVED, relative to P78)

1. **Exact separation.** eq. `low-separated` holds verbatim with A_{m,σ,v}[λ]: its derivation is Mellin inversion in W₀ only.
2. **Any bounded λ.** For any |λ| ≤ 1, lem:planar-additive-sieve and the proof of lem:balanced-additive-norm give Σ_m |A[λ]|² ≪ (Q+Y²)/Y = (Q/Y)(1+P_a). That proof uses only |c(s)| ≤ 1 and the distinctness of the reduced fractions d/s.
3. **The Gram bound** prop:probe-gram, (Q/Y)(1 + P_a^{1/6} + P_a²/Y), is proved only for "fixed finite linear combinations of the ray coefficients". Its proof needs the column coefficient a_{C,d,k} to be periodic modulo a fixed 𝔩_fixed·rad(C)·rad(k).
4. **Corollary (free s-weights at b = 0).** If P_a = Y²/(q_{b_*}Q) ≍ 1, i.e. b = 0, the compensated low estimate holds unchanged for every |λ| ≤ 1, including target-dependent λ. This is because the B-side (lem:probe-row-norm) does not involve s, and item 2 equals the Gram bound up to constants. On the wall, B = 11/12 − ℓ/4.
5. **Cost for b > 0.** A non-admissible λ costs up to (P_a^{5/6})^{1/2} = Z^{5b/12} on the low side unless a new Gram lemma is proved; at b = 0.06 that is 0.025 in B.
6. **The admissible class gives the paper's numerator.** For λ = Σ_ψ a_ψ ψ with fixed ray ψ, N_u is a fixed combination of plain sums for χ_uψ, of conductor ≍ q_u. The reflected length is n₀ = 1 − l_y/d up to constants, and lem:plain applies with the same capacity 1 − 2n₀.

So **at b > 0, changing the numerator requires leaving the class that the Gram lemma covers**.

At b = 0 the high side is worse: worst E = +0.0132 at (0.55, 1/2). The freedom there does not help either (§3).

## 3. (a) s-weights

### (a1) Mollifier or Selberg weights, λ = 1*λ′ with λ′ supported on q_d ≤ D = Z^θ, |λ′_d| ≤ 1

**Identity.** By §1,

  N_u = Y^{-1/2} Σ_d λ′_d π_u(d) q_d^{-1/2} S̃_u(Y/q_d) + (non-coprime corrections of the same shape).

Equivalently, in Mellin form, the numerator is L(w,χ_u)·M_{λ′}(w;u)·(Euler correction).

**What the lemma set gives.** lem:plain accepts only plain sums times prime polynomials with fixed ray coefficients. M_{λ′} is neither, so it can only be removed by the coefficient norm, exactly as the 08 lemma removes A(k). This costs Σ_d |λ′_d| q_d^{-1/2} ≍ D^{1/2}. Each plain sum also has reflected length n₀ + θ/d, so the capacity drops to 1 − 2n₀ − 2θ/d. Hence E_wt and E_old both get at least +θ/2 (MODEL; numbers in §0).

**A twisted fourth moment does not help either.** Expanding |M S̃|⁴ gives Σ_u χ_u(e) × (four plain sums) with e = d₁d₂/(d₃d₄). lem:plain gives no cancellation in χ_u(e), so the cost (Σ|λ′|q^{-1/2})⁴ is the same.

**Low side.** Not Gram-admissible (each piece is periodic mod d, with d growing); free at b = 0 (§2.4).

**No moment gain is possible in principle (HEURISTIC).** For the full row family, the fourth moment of any numerator with a fixed signal is ≍ U (orthogonality: Σ_u |⟨λ, v_u⟩|² ≈ U‖λ‖² ≥ U|Σλ|²/Y). A mollifier lowers only log factors.

### (a2) λ = Λ (s prime, weighted by log)

- **High side.** The principal row keeps its w = 1 pole: −ζ′/ζ gives Ŵ₁(1). On a row in bin a, the explicit formula and the bin's zero-free rectangle give only a pointwise bound, |Σ_{p~Y}χ_u(p)log p·W| ≪ Y^{a+ε}. In U-units that is σ_Λ = δ l_y/(2d). Then E_Λ − E_old = δ(l_y − d/2) > 0, since l_y > d/2, so it is worse than the paper's pointwise branch. There is no fourth moment at length Y.
- **Low side.** Not Gram-admissible.
- **Model:** +0.0146.

### (a3) Restricted supports

- **s ≡ b mod r, or a ray class.** With r fixed, this is a combination of ray characters and changes nothing. With r = Z^θ, the conductor is U·r and the dual length U r/Y is longer.
- **s = t² or t³.** The s-Gauss sum at p² has character χ_p² (a = 2), so (H1) fails: every cube row u = b³ is principal, with denominator 1/L(x, ηψ_b) where ψ_b is quadratic. These rows carry q_u^{-1/6} = q_b^{-1/2} with no further cancellation, so the extra signals total Σ_{q_b ≤ Z^{h/3}} q_b^{-1/2} ≈ Z^{h/6}. That is fatal.
- **Rough s (sieve weights).** Reduces to (a1) through the fundamental lemma.
- **Angular localisation.** α^k weights kill the w = 1 pole (D_probe_redesign C6). A smooth sector of angle 1/T expands into Σ_{|k|≤T} L(w, α^kχ_u), so there is no gain.
- **Short norm intervals |q_s − Y| < H₀.** The dual length U/H₀ > U/Y is longer.

**The plain smooth W₁ minimises the reflected length.** Among s-weights smooth at scale ≥ Y with fixed-conductor twists, every numerator has dual length ≥ U/Y·O(1). Any finer structure in λ lengthens the dual.

### (a4) Two-scale radial weights, W₁(q_s/Y₁) + βW₁(q_s/Y₂)

This is a sum of probes. The numerator is governed by the smaller scale, and the low side by the superposition lemma (07 independent_low §3). The signal and theta coefficients stay proportional, as in 04 IDEA1 §2.

### Remark on the Gram term with arithmetic λ (HEURISTIC, outside the mission)

For radial weights the Gram and signal coefficients are both Ŵ₁(1)Ŵ₀(5/6). For λ = 1*λ′ the theta part of A_m becomes ρ_λ = Σ_d λ′_d ρ^{(d)}, so the Gram term becomes a quadratic form in λ′ while the signal stays a linear form. This is not ruled out by 07 §2, which covers only radial kernels.

However, since the signal is the 5/6-pole of Σ ρ_λ(m)B_m, the Cauchy–Schwarz slack is the angle between ρ_λ and B. A D-dimensional λ′ space can improve that angle by at most about D^{1/2}. This is the same factor the high side pays in (a1), so the trade is break-even at best.

Making it rigorous is a (V′)-type problem, which is signal-circular.

## 4. (b) Splitting Y = Y₁Y₂, or s-slot primes

**Exact identity.** Take c = W₁(·/Y₁) ⋆ W₁(·/Y₂), a Dirichlet convolution. By §1, with the local factor Σ_{k₁,k₂} π_p(k₁+k₂) Q^{−k₁w₁−k₂w₂}:

  N_u = (2πi)^{-2} ∫∫ Y₁^{w₁−1} Y₂^{w₂−1} Ŵ₁(w₁)Ŵ₁(w₂) L(w₁,χ_u)L(w₂,χ_u) H₂(w₁,w₂;u,x,z) dw₁dw₂,

that is, Y^{-1/2}·Σ A(k₁,k₂)·S̃(Y₁/q_{k₁})·S̃(Y₂/q_{k₂}). This is a product of two plain sums.

**Lengths (PROVED).** Each factor has length n_i = min(y_i, 1−y_i), where y₁ + y₂ = l_y/d. So n₁ + n₂ ≥ n₀, with equality only for the degenerate split (the paper's). Splitting never shortens the total length.

**Moment degree (MODEL; the Hölder computation is exact).**
- lem:plain bounds Σ_u |S(n₁)S(n₂)Q|² ≪ U whenever n₁ + n₂ + 6κz ≤ 1. For a product numerator this is only its **second** moment. Hölder then gives σ_b = (1−R)/2 − q·min(L, (1−n₁−n₂)cP), against σ_wt = (1−R)/4 − q z_N/2.
- The extra capacity is worth at most qL ≤ δL/2 ≈ 0.11δ. The lost Hölder degree is worth (1−R)/4 ≈ 0.25(1−R), which is larger at every binding bin.
- **Mixed Hölder for the product.** The other lem:plain moments, Σ|S_i|⁴|Q|², do not help. Matching the S₁- and S₂-degrees, 2α′ + 4β = 2α′ + 4γ = 1, forces the Hölder weights to sum to Σθ = 1/2 in every combination. The capacity average is again (1 − n₁ − n₂)cP, so σ_b is unchanged.
- **Mixed Hölder for the paper's single sum.** Mix the fourth moment (weight α, capacity z₁ = (1−2n₀)cP) with the second moment, i.e. lem:plain with n₂ = 0 (capacity z₂ = (1−n₀)cP). Then dσ/dα = −(1−R) + 2q(2z₂ − z₁) = −0.41 + 0.10 < 0 at the worst bin, so the pure fourth moment, α = 1/4, is optimal.
- **Model:** +0.094; +0.074 with DIAG capacity; +0.106 for an unreflected Y₁Y₂.

**Low side.** c is a divisor-type weight, so it is not Gram-admissible. Its theta part is a double (A₂-type) Kubota series, which 04 wf2 shows is not controlled.

**Upper envelope (HYP).** A product with n₁ ≠ n₂ could beat the single sum only through higher diagonal moments Σ|S₁^a S₂^b|² ≪ U with a n₁ + b n₂ ≤ 1. These are not available: lem:plain has degree 4, and the sextic large sieve is diagonal only up to length U^{1/2}.

Even granting every real order k ≤ 1/n₀ (the ideal σ = n₀(1−R)/2 when capacity is exhausted), `frac.py` gives:
- worst E = +0.0047 with cP = 2/9, i.e. θ ≥ 0.87326;
- worst E = +0.0019 with cP = 1, i.e. θ ≥ 0.87273.

**No diagonal-type numerator moment reaches 41/47.** One needs an off-diagonal input, such as the sixth moment Σ|S̃|⁶ ≪ U^{1.10} at length 1.32, or (sLS) of §5.

**s-slot primes, c(s) = Σ_{p|s, p∈𝒫′} W′(q_p/P′)** (MODEL):
- **High side.** The numerator becomes Q′(u)·S̃(n₀ + ℓ′/d), where Q′ = Σ_p χ_u(p)p^{-1/2}·… is in the numerator orientation.
  - lem:plain accepts Q′ only squared, so Σ_C |Q′S̃| ≤ (Σ|S̃|⁴|Q′|²)^{1/4}(Σ_C|Q′|^{2/3})^{3/4}, which gives σ = (1−R)/2.
  - Alternatively, |Q′| bounded pointwise costs ℓ′/2.
  - Unlike the A-marks, Q′ has no lower bound on C, so it cannot buy capacity.
- **Low side.** The signal factor is Σ_p W′/q_p ≍ 1/log P′, with no P^{1/6} amplification (the A-marks' Σ W q_p^{-5/6} has it). There is no low-side gain.

So s-slot primes are a pure loss.

## 5. (c) Target-dependent (non-effective) s-weights

**Formulation.** Write c₀ for the paper's weight, and v_{u,k} for the vectors s ↦ 1_{k|s}χ_u(s/k)W₁(q_s/Y), one for each Euler-correction label k (§1). Take

  λ = 1 − (Π_W c₀)/c₀,   W = span{v_{u,k} : u ∈ C, q_k ≤ U^{θ_D}}.

**Condition 1: dimension (PROVED count; MODEL thresholds).**
- In the reconstructed numerator, A(k) ≍ q_k^{-a}. Truncating at q_k ≤ U^{θ_D} leaves an untouched remainder of at most U^{θ_D(1/2−a)}·U^{σ_wt}.
- Reaching σ* therefore needs θ_D = (σ_wt − σ*)/(δ/2), and #C·U^{θ_D+ε} < Y, i.e. R + θ_D < l_y/d.
- At the worst bin (0.45, 1/2), b = 0.06: #C = U^{0.610} = Z^{0.486} against Y = Z^{0.441}, and θ_D = 0.066.
- Under Gram domination (b < 0.1755), l_y/h ≤ 0.585 < R = 0.610. So **the worst bin is never annihilable, even with θ_D = 0.** Only the bins with δ ≥ 0.6 pass (`ctable.py`).
- Model: optimistic (c), with every side cost ignored, gives +0.0118; with no Euler overhead, +0.0114.

**Condition 2: the signal survives.** One needs ‖Π_W c₀‖² ≤ ‖c₀‖²/4.
- The available bounds are ‖Π_W c₀‖² = S*G_C^{-1}S ∈ [‖S‖²/λ_max(G_C), ‖S‖²/λ_min(G_C)]. Here S_u = ⟨c₀, v_u⟩ is the numerator vector and G_C is the Gram matrix of the v's.
- The only sufficient condition is therefore Σ_{u∈C}|S̃_u|² ≤ λ_min(G_C)/4 ≤ Y/4: a reverse large sieve plus the bound (ii) Σ_C|S̃|² ≲ U^{l_y/d}. Neither is provable.
- **Self-domination (PROVED inequality, MODEL check).** Hypothesis (ii) alone, with Cauchy–Schwarz, gives σ ≤ (l_y/d − R)/2. At every annihilable bin this is ≤ σ*: at δ = 0.6 it is 0.050 against 0.122, and at δ = 0.65 it is 0.078 against 0.141. So annihilation adds nothing beyond its own admissibility hypothesis.
- In sparse-LS form, (sLS) Σ_C|S̃|² ≪ U^ε(#C + U^{l_y/d}), the hypothesis closes the whole target (−0.029).
- The minimal version needed at the worst bin is Σ_C|S̃|² ≤ U^{R+2σ*} = U^{0.770}, against the proved U^{(1+R)/2} = U^{0.805}.

**Condition 3: the other bins.** The correction c₁ = −Π_W c₀ is controlled on all other rows only by the large sieve. The sextic large sieve constant at length Y is ≍ (UY)^{2/3} = U^{1.04}, giving σ_LS = (1.04 − R′)/2 ≫ σ* at the δ ∈ [0.25, 0.45] bins. The fix would be ‖c₁‖²/Y ≤ U^{R′+2σ*−1.04}, which contradicts the size ‖c₁‖² ≈ Σ_C|S̃|² that annihilation needs in the adversary scenario. On the low side c₁ costs Z^{5b/12}·(‖c₁‖²/Y)^{1/2}; this cost vanishes at b = 0 (§2.4).

**Partial projections (MODEL, exact).** Annihilating a subset C′ with #C′ ≤ Y·U^{-θ_D} removes a fraction ≤ U^{-0.09} of the worst bin. The adversary, who has uniform numerators U^{σ_wt} on C, keeps (1 − o(1))#C rows with the same numerators. The Hölder bound on C∖C′ is unchanged, and the c₁ part only adds a large-sieve term. **No gain**, at any level of projection.

**Other target-dependent choices.**
- Slot windows (choosing among K′ windows): the adversary uses amplitude 1/2 in all of them, which is consistent because the count lemma's prime capacity (≤ 7/37) is already saturated.
- The s-scale Y_j: large values of L(1/2+it,χ_u) make every smooth scale large at once.
- ±1 discrepancy colourings: these reduce to the same kernel and dimension condition.

## 6. (d) Coupled probes sharing rows

**Linearity (PROVED).** For fixed Z and slots, I_η is linear in the joint weight Ψ(s,m) = c_σλ(s)W₁(q_s/Y)W₀(q_m/(q_{b*}q_sX))/(Y√(q_{b*}q_sX)). So Σ_j a_j I_η(X_j,Y_j,Z) is one probe with Ψ = Σ_j a_jΨ_j:
- **Same X:** it is a single probe with the multi-scale or arithmetic s-weight Σ_j a_j c_j, which is §3.
- **Different X_j:** the Mellin weight is no longer a product. The row-u numerator at fixed z is Σ_s χ_u(s)π_u(s)ψ_z(q_s), where ψ_z = Σ_j a_j X_j^{1/2−z} W₁(·/Y_j)/Y_j. This is radial at scales ≥ min Y_j, so the numerator cannot be shorter.
- **Signal and theta part:** both are proportional to Σ_j a_j Y_j^{-1/3}(…) (07 §2). Cancelling the theta part cancels the signal.

**Joint low bound.** A joint bound better than the triangle inequality needs cancellation between the low forms Σ_m A^j_m B^j_m. With a common Q and a common B this is just the linear case. Otherwise it needs a cross-correlation estimate that nothing supplies. Products and positive powers give weighted averages (07 §3).

**Oracle (MODEL).** Even allowing a joint low bound as good as each probe's and a per-bin choice of geometry (ℓ ∈ [0.1773, 0.2], b ∈ [0, 0.2]) on the high side, the worst bin gives +0.0113 at (0.45, 1/2), with the best b* = 0. No coupling of these probes reaches 41/47.

## 7. What the numerator side would actually need

Every redesign above either needs, or is dominated by, an off-diagonal large-values input for short sextic character sums over the bad rows. Any one of these suffices in the model:
- (sLS);
- or Σ_C|S̃_u|² ≤ U^{R+2σ*}, i.e. U^{0.770} at the worst bin;
- or the sixth moment of NOTES §2.4.

Diagonal inputs, including all fractional diagonal moments (§4), stop at θ ≥ 0.8727.

The probe's remaining free parts that touch the numerator are:
- the s-weight (this note);
- w(m), which must be ray by (H4);
- W₀, Φ, ξ, b_* and T, which change only constants;
- the geometry, which is already optimised.

The slots (a, j, κ) are rigid by D_probe_redesign Theorem 6.1.
