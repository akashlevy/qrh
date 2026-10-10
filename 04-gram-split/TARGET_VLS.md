# Target: the theta-pairing estimate (V′) needed to beat Liu's 0.874957

*9 Oct 2026. This is the precise target for the "prove the large sieve estimate" effort. Background:*
- *`REPORT.md`, including the CORRECTION banner;*
- *`wf/AUDIT_chain.md` and `wf/V_gaps.md` (Lemma Z, F1);*
- *`wf2/A2_theory.md` and `wf2/A2_numerics.md`;*
- *`wf/FE_LoA6_proof.md` for (FE)/(LoA6), which hold.*

## 1. Objects

Notation:
- K = Q(ω), O = Z[ω], λ = 1+2ω.
- q_x = N(x).
- e(z) = exp(4πi Im z/√3).
- Primary means ≡ 1 mod 3.
- χ_n(u) = (u/n)_6, zero on non-units.
- g(s,h) = Σ_{d mod s} χ_s(d) e(hd/s), and g̃ = g/√q_s.
- γ_j(c) is the normalised Gauss sum of χ_c^j.
- α(a) = a/|a|.

**Sextic theta coefficients.**
- D_c(w;m) = Σ_{s primary, (s,6)=1} c(s) g̃(s,−m) q_s^{−w}, where c is a fixed finite-order ray-class twist (the 7/8 paper's c(s) = 1_{s∈σ} χ_s(b_*)/(τξ(s))).
- D_c continues meromorphically; its only pole in Re w > 1/2 is at w = 2/3 (`wf/FE_LoA6_proof.md`, Theorem A).
- ρ_c(m) := Res_{w=2/3} D_c(w;m).
- Σ_{q_m≍Q} |ρ_c(m)|² ≍ K₀ Q^{5/6}. This mean square is proved (D_derivation §6.3); numerically ρ is "unitary" after the factor q_m^{1/12}.

**Reflected completed row.** Paper78 eq. reflection, ~l.1763; row phase eq. reflection-row-phase, l.1832. For the base probe (no slots) and row m:

  B_m = ζ_m · B⁰_m,  B⁰_m = Σ_{μ ∈ λ^{-4}O∖0} d(μ) α(μ) ϑ(λ⁴μ) q_μ^{−1/2} ∏_{p|m} B_p(λ⁴μ) · V♯(q_μ Z / q_{c(m)}²),

where:
- c(m) = c_F · rad(m), so the dual length is L_m ≍ q_{rad m}²/Z ≤ L := Q²/Z;
- for squarefree m, ∏_{p|m} B_p(x) = χ_m(x)^{-3}, i.e. the quadratic symbol;
- d(μ) is Patterson's cubic theta coefficient: 3^{5/2}|b| g̃₃(n) at μ = λ^{-3}nb³ with n squarefree;
- ζ_m is the reflection scalar. **Lemma Z:** for squarefree m, ζ_m = κ₀ ψ₁(m) conj(γ₂(m))⁴ = (ray)·μ(m)ᾱ(m)·conj(γ₂(m)), a cubic Gauss sum times an angle character.

**Theta pairing.** With Ω the row weight, ξ the fixed residue character and v ∈ R:

  T(Q,Z;v) := Σ_{m≠0} Ω(q_m/Q) ξ(m) (q_m/Q)^{−iv} ρ_c(m) B_m = Σ_m π(m) B⁰_m (weighted), where π(m) := ρ_c(m) ζ_m.

## 2. What is needed

Geometry: Q = Z^M with M ∈ [0.80, 0.86]; at the slot-free optimum M ≈ 0.826. Then:
- L = Q²/Z = Q^{2−1/M} ≈ Q^{0.79} < Q;
- the cubic-theta column weights have Σ_μ |β_μ|² ≍ 1 per dyadic block (β_μ = d(μ)α(μ)ϑ q_μ^{−1/2}).

Trivial (Cauchy–Schwarz) bound: |T| ≤ ‖ρ‖_Ω ‖B‖_Ω, with ‖B‖² ≪ Z^{M+ε} = Q^{1+ε} (Part I's unmarked completed moment), so |T| ≪ Q^{5/12}·Q^{1/2} = Q^{11/12}.

**(V′) Target.** For some fixed s₀ > 0, uniformly in |v| ≤ Z^ε (polynomial dependence allowed),

  |T(Q,Z;v)| ≪ Q^{11/12} Z^{−s₀+ε}.

Any s₀ > 0 beats Liu. With the slot-free low side of `wf/PARTII_slotfree.md` and (LoA6), the Gram excess term becomes g = max(0, b/6 − 2s₀):

| s₀ | θ |
|---|---|
| 0.002 | 0.8733 |
| 0.005 | 0.8709 |
| 0.01 | 0.8668 |
| 0.02 | 0.8615 |
| 0.04 | 0.8532 |

(`wf/lp_const_saving.py`.) A saving that holds only for M in a sub-range is also useful; the LP can use any range.

**Sufficient large-sieve form (V′_LS).** For some η > 0,

  Σ_{q_μ ≍ L} |Σ_m Ω(q_m/Q) π(m) ψ_μ(m)|² ≪ Q^{1−η} Σ_m Ω|π(m)|²,  for L = Q^{0.79±0.02},

where ψ_μ(m) = (λ⁴μ/m)_2 · (fixed ray factor). Compare the general large sieve (Goldmakher–Louvel), (L + Q)‖π‖², whose Q-term is exactly Cauchy–Schwarz. Since L < Q this is a "fewer moduli than length" estimate. It is false for general π and must use the structure of π.

**Weaker, also sufficient (V′_β).** It is enough to bound the specific bilinear form Σ_μ β_μ Σ_m π(m) ψ_μ(m), using the structure of β_μ, i.e. Patterson's cubic theta coefficients, which are Gauss sums at squarefree μ.

## 3. What is known or ruled out (do not redo)

- **The ζ-free version is fine.** Σ_m ρ(m)(μ/m)₂ has a Voronoi formula with dual length ≍ q_μ²/Q and a cubic-character dual (V_derivation.md), and Heath-Brown's cubic large sieve then saves. **But the true row carries ζ_m.**
- **The A₂ functional-equation route does not save.** π is the theta residue of the sextic twisted A₂ WMDS. The quadratic twist is a twisting parameter μ³, and the functional equation z ↔ 5/6−z has factor N(μ)^{3(1/2−z)}, so the dual has length about L³/Q > Q. Its dual coefficients ρ(sμ^k) couple s and μ, so the trivial bound on the dual gives no saving for any k ≥ 2 at M ≈ 0.826 (`wf2/A2_theory.md` and its referee report).
- **Numerics.** π is at random level against quadratic twists: Σ_μ|F′(μ)|² ≈ L‖π‖²·O(1), so (V′_LS) is plausibly true. π is not "simpler" than ρ, and the cubic-type degeneracy is excluded (`wf2/A2_numerics.md`).
- **The unreflected route.** In the unreflected form, B_m = Σ_{A≍Z} γ₂(A)ᾱ(A)Ψ(A) χ_A(m). The j=1 twist is coherent and its dual is principal (Ramanujan sums). By eq. signal, the pairing becomes S + T′, where S is the Möbius (zero-detecting) signal. S cannot be bounded unconditionally on its own, so that route can only help if T′ is bounded and S is handled together with the high side.
- **Larger families.** Ideas 2 and 3 of the original goal (Kuznetsov, energy wall) need second moments "beyond the edge" of the same type.

## 4. Tools available

Throughout, for (s,m) = 1: g̃(s,−m) = χ̄_s(−m) g̃(s); the Kubota multiplier on Γ(72) is κ = (c/d)_6; and γ₁γ₂ = μαG (Hasse–Davenport).

- Heath-Brown cubic large sieve and de Faveri's Θ_n(A,B) ≪ A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}. This is sharp (DFDH).
- The Goldmakher–Louvel quadratic large sieve over K.
- Kubota/Patterson functional equations, and the cubic theta coefficients via Patterson / Dunn–Radziwiłł (Lemma 5.2 Voronoi, uniform in the modulus).
- The paper's planar additive large sieve.
- Kazhdan–Patterson relations for θ₆: τ(mp^n) = N(p)^{1/2}τ(m), τ(mp^{n−1}) = 0, τ(mp^j) = G_{j+1}(m,p) τ(mp^{n−2−j}).
- The Louvel Bruggeman–Kuznetsov formula on the cubic cover.
- Numerics library `zomega.py`; fast Gauss sums in `D_numerics/fastgauss.py`; ζ_m in `wf2/a2lib.py`; data `A_Q*.npy` and `D_numerics/`.
