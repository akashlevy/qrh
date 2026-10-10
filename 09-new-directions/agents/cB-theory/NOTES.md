# cB-theory: the conductor exponent c_B of the cubic architecture's dual (10 Oct 2026)

Labels: **PROVED** (finite algebra), **R-** (rigorous modulo standard rank-one Kubota theory, not written out),
**HEUR**, **NUM**. Files are listed in sec 8.

## 0. Verdict

**c_B = 4 (R-), and it is sharp.** The cubic θ₃⊗θ₃ architecture therefore caps at 10/11 (alt-detect's table).

- The coprime part of an n-th order Kubota series at index p^l has a vector-valued local FE at p.
  - Its diagonal factor is exactly N(p)^{(l mod n)(1/2−s)}, with weight 1.
  - The only other entries are off-diagonal, of weight N(p)^{−1/2}.
  - Unlike σ_{1−2s}(n) for GL2 Eisenstein series, nothing is carried away by the p-part polynomial when l < n.
- B_h is the coprime sextic series at index h⁴, so l = 4 < 6, which gives N(h)⁴.
- The rule reproduces every case we can check rigorously:
  - P78's reflection at all six twist exponents j_p;
  - the A₂ referee's n = 2 drop (index μ³ → conductor |μ|).
- **Calibration premise corrected.** P78's dual is the *Möbius-twisted* sextic series. By Davenport–Hasse (DH) that is a quadratic twist of the cubic theta series, with conductor N(k)², and it is not a sextic Kubota index shift.
  - The plain sextic coprime part at index k³ has conductor N(k)³.
  - B_h is the reverse case. It is Möbius-twisted cubic, which by DH is plain sextic at index h⁴, giving N(h)⁴.
- **Orientation (task item 4): alt-detect's claim is confirmed (PROVED).**
  - Index h⁴ is forced in both shapes: in P2 shape by Poisson, in P78 shape by μ-absorption.
  - The index-h² object (c = 2) is exactly the orientation in which γ²·γ̄ = γ, so it detects nothing.
- **No c_B = 2 variant with a single Poisson step** (sec 6).
  - The triple can be split as (Poisson, Poisson, θ₃) with a two-variable probe. That gives c_B = 2 rigorously, with a full-index cubic dual.
  - But the second Poisson step multiplies the frequency count by Z.
  - In P2 shape the split variant is provably no better (θ ≥ 7/6). P78 shape is untested (HEUR: likely no gain).

## 1. Exact functional equations (task item 1)

**Normalisation.** ψ_n(r,s) = Σ_c g̃_n(r,c) N(c)^{−s}, with g̃ = g/N(c)^{1/2}. It converges absolutely for Re s > 1, has centre 1/2, and n | #μ(F) (n = 3 or 6 over Q(ω)).

**Archimedean input (R-).**
- The r-th Fourier coefficient of the n-fold metaplectic Eisenstein series E(w,s′), with s′ = s + 1/2, is ∝ |r|^{s′−1} ψ_n(r,s) K_{s′−1}(4π|r|y)/Γ(s′).
- The cover splits at the complex place, so this is the ordinary K-Bessel function.
- E(s′) = Φ(s′) E(2−s′), where Φ is the finite scattering matrix over the S-classes (S = primes above 6).
- Its constant term is ∝ Γ(s′−1)/Γ(s′) · ζ_F(ns′−n)/ζ_F(ns′−n+1).
- Gauss multiplication then gives:

  Ψ_n(r,s) := (2π/n)^{−ns} Π_{j=1}^{n−1} Γ(s − 1/2 + j/n) · ζ_F(ns − n/2 + 1) · ψ_n(r,s),

  **Ψ_n(r,s) = N(r)^{1/2−s} · Σ_{Ψ′} τ_{Ψ,Ψ′}(s) Ψ_{n,Ψ′}(r, 1−s).**

  Here τ is a finite matrix of exponentials and roots of unity. It depends on r only through r mod F_S^{×n} (Hilbert symbols); cf. Brubaker–Bump 2006, Thm 1, "D(s,m) = |m|^{1−s} Σ T_ij D(2−s,m;i)".
- The factor N(r)^{1/2−s} comes *only* from the Bessel argument |r|y. It is therefore the complete r-dependence of the full series.

**Gamma factors (spherical case).**
- n = 2: Γ(s) (Hecke L-function) ✓.
- n = 3: Γ(s∓1/6).
- n = 6: Γ(s−1/3)Γ(s−1/6)Γ(s)Γ(s+1/6)Γ(s+1/3). Up to the ζ normaliser the object has "degree n−1".
- Grössencharacter twists change the K-type and shift these Gamma factors. P78's ᾱ-twisted cubic has R(t) = Π_±Γ(1+t±1/6)/Γ(1−t±1/6), i.e. Γ(s+1/2∓1/6).
- **None of this affects the r-aspect.** Label: R- (Gamma shifts HEUR).

**Cubic case made fully explicit (rigorous).**
- P78 Prop. `completed-reflection` (l.1676) *is* this FE for the coprime cubic series twisted by Ψ = Π χ_p^{j_p}. It gives:
  - the scale q_c² (c = c_F · Π_{active} p);
  - the local factors B_p;
  - the scalar ζ with |ζ| ≤ 1/81 and the explicit ω_{p,j};
  - three fixed cusp-coefficient functions.
- By reciprocity, j_p = 4 ↔ cubic numerator p, and j_p = 2 ↔ numerator p².
- j_p odd is a quadratic or sextic twist, which is not an index shift.

## 2. p-part polynomials and the local FE of coprime parts (task item 2)

Fix one prime p (q = N(p), X = q^{−s}). Let C_j = Σ_{(c,p)=1} (p/c)^{−j} g̃(m′,c) N(c)^{−s} (j mod n) be the coprime part at index p^j m′.

**Gauss-sum p-parts (PROVED; NUM check (A) in `orient_check.out`).** For g̃(p^l, p^k):
- k ≤ l with n | k: q^{k/2}(1−1/q);
- k = l+1: q^{l/2}·G̃_{l+1}(p), or −q^{(l−1)/2} if n | l+1;
- all other k: 0.

Twisted multiplicativity then gives, for 0 ≤ l ≤ n−1,

  ψ(p^l m′) = C_l + a_l X^{l+1} C_{l′},   l′ = −l−2 mod n,   |a_l| = q^{l/2}.

| n | v_p(r) = l | p-part relation |
|---|---|---|
| 3 | 1 | ψ = C₁ + q^{1/2}G̃ X² C₀  (pair {0,1}) |
| 3 | 2 | ψ = (1 − q^{1/2}X³) C₂  (singleton) |
| 3 | 3 | ψ = (1+q^{3/2}(1−q^{−1})X³) C₀ + q^{3/2}G̃ X⁴ C₁ |
| 3 | 4 | ψ = (1+q^{3/2}(1−q^{−1})X³) C₁ + q²G̃ X⁵ C₀ |
| 3 | 5 | ψ = (1+q^{3/2}(1−q^{−1})X³ − q²X⁶) C₂ = (1+q^{3/2}X³)·ψ(p²) (KP) |
| 6 | 1 | ψ = C₁ + q^{1/2}G̃₂ X² C₃  (pair {1,3}) |
| 6 | 2 | ψ = (1 + qG̃₃ X³) C₂  (singleton; G̃₃ is a quadratic Gauss sum: the L(3s−1)-type normaliser) |
| 6 | 3 | ψ = C₃ + q^{3/2}G̃₄ X⁴ C₁ |
| 6 | 4 | **ψ = C₄ + q²G̃₅ X⁵ C₀  (pair {0,4})** |
| 6 | 5 | ψ = (1 − q²X⁶) C₅  (singleton) |

**Local FE (PROVED given sec 1; `local_fe.py` / `local_fe.out`).** Insert the full-series FE into the 2×2 block and solve. With T_j = ζ^{(p)} C_j:

  T_l(s) = G(s)·[ q^{l(1/2−s)}(1−q^{−1}) T_l(1−s) + q^{−1/2}(e·q^{−(1/2−s)} − e′·q^{(n−1)(1/2−s)}) T_{l′}(1−s) ] / (1 − q^{−1}q^{−n(1/2−s)}).

- The diagonal exponent is l, with weight 1.
- The off-diagonal terms have weight q^{−1/2}. They are just the shifted conductor-q^{l′} series a_l X^{l+1}ψ_{l′}.
- Singletons are pure: q^{l(1/2−s)} times an Euler-factor ratio.
- Coprime parts depend only on l mod n. So:

**Rule (PROVED given sec 1): the coprime part at index r has effective conductor Π_p q_p^{(v_p(r) mod n)}.**

## 3. c_B = 4

**Exact identity (PROVED).** Invert the tensor product of the local 2×2 blocks (pair {4,0}). For squarefree h:

  B_h(s) = Σ_{d|h} Π_{p|d}(−a_{4,p} N(p)^{−5s}) · ψ̃₆((h/d)⁴; s) / Π_{p|h}(1 − N(p)^{2−6s}).

Here ψ̃₆ is the full (Hecke-, ray- and Ψ-twisted) sextic Kubota series, whose FE has conductor N(h/d)⁴.

**Upper bound.** The AFE length is N(h)², i.e. c_B ≤ 4, as alt-detect stated.

**Lower bound (R-).** The d = 1 term contributes a dual sum Σ t₄*(m) K(N(m)N(h)⁴/x) with coefficients of full size: Gauss sums at squarefree m, times the S-matrix, which is invertible.
- No other term has a dual at that scale with those coefficients. The characters (p/m)^{−j} are independent, and the exponentials q^{e(1/2−s)} are independent.
- So the smoothed sums Σ b_h(n)V(N(n)/x) stop oscillating only at x ≍ N(h)⁴.

**Robust to completions (PROVED).** Consider any completion F = T₄ + Σ_{k≥1} X^k Σ_j c_{kj} T_j.
- Its ψ₄-component has the coefficient (P(X) − a₀X·Q(X))/D with constant term 1, so it cannot vanish.
- The q^{4(1/2−s)}-dual of T₄ could only be cancelled with |c_{50}| ≍ q³. That completion is larger than the row itself on the critical line.

**General h.** h = h₁h₂²h₃³ gives conductor N(h₁)⁴N(h₂)², since 4e mod 6 ∈ {4, 2, 0}.
- Cubes have conductor 1. These are exactly the copy rows u = a³, a consistency check.

**Heuristic picture (HEUR).** ψ_n behaves like a degree-(n−1) object, and the index p^l "ramifies l of its n−1 components".
- For the sextic family, the theta object is the cubic theta: degree 2, and every twist has conductor ≤ q².
- For the cubic family, the dual lives on the degree-5 sextic Kubota series, where index h⁴ costs h⁴.
- The GL2×GL2 Rankin–Selberg heuristic ("3") does not apply to metaplectic thetas.

## 4. Calibration on P78 (task item 3)

μ moves a series between orders (P78 Lemma `prime-gauss-identities`; DH):
- μ·γ̄₁ = (ray·Hecke)·γ₂;
- μ·γ₂ = (ray·Hecke)·γ̄₁ = (ray·Hecke)·conj G̃(sextic).

| object | Möbius-free form | c (rule) | rigorous check |
|---|---|---|---|
| Möbius-sextic, numerator k³ (P78 dual) | cubic theta ⊗ (k/·)₂ (not an index shift) | 2 | P78 j_p = 3: B_p = χ_p^{−5}, scale q_p² ✓ |
| Möbius-sextic, numerator P⁴ (P78 marks χ_P^{−2}) | cubic, numerator P | 1 (+ q^{−1/2}·q² piece) | P78 j_p = 4: B_p = q^{−1/2}(−1+q·1_{p\|x}) ✓ (both pieces match sec 2) |
| Möbius-sextic, numerator P² | cubic, numerator P² | 2, pure | P78 j_p = 2: B_p = χ_p^{−4} ✓ |
| plain sextic, numerator k³ | — | 3 | (not P78's object) |
| **Möbius-cubic, numerator h (= B_h)** | **plain sextic, index h⁴** | **4** | — |
| Möbius-cubic, numerator h² | plain sextic, index h² | 2 | — |
| quadratic Gauss sums, index μ³ (wf2 referee, n = 2) | — | 1 (KP) | referee's n2check ✓ |

## 5. Orientation (task item 4) — alt-detect confirmed (PROVED)

**P2 shape.** Poisson of χ_n(u) gives g(χ_n,h) = χ̄_n(h)g(χ_n). So the dual is Σ μν g̃(χ_n)χ̄_n(h): numerator h, hence sextic index h⁴. The conjugate family gives the conjugate dual.

**P78 shape.** The row coefficient is c(A)·ψ_A(m)^a. Poisson in m supplies γ(ψ^a), and absorption needs c·γ(ψ^a) = μ·Hecke.
- The only N = 3 solution is c = γ(ψ^a)², i.e. the twist has the same orientation as the Gauss sums.
- γ₁(A)²κ_A(m) = Hecke·conj(μγ₁(A; m)), and by DH this is sextic index m⁴.
- The index-m² coefficient is γ₁²κ̄(m). Absorption would then need Poisson to supply γ̄₁, and γ₁²γ̄₁ = γ₁ contains no μ. **So c = 2 ⇔ no detection.**

**NUM** (`orient_check.out`, sextic Gauss sums computed by brute force over O/n):
- μγ₂(n;h)/conj γ₁(n;h⁴) is independent of h, for n = p, pq, pqr.
- The same ratio with h² in place of h⁴ varies with h.

## 6. Alternatives (task item 5)

- **(a) Twisting one theta factor.** γ₁(A)γ₁(A; m²) = γ₁(A; m)² on coprime A. The Dirichlet series is identical, so c = 4 (PROVED).
- **(b) Other pairings.** For N = 3 the only options are (1,1,1) and its conjugate (2,2,2). Both give c = 4.
  - For other N, γ_b² is not a single Gauss sum, since the DH duplication needs a sextic σ with σ² = ψ^b. So there is no automorphic row (alt-detect sec 2).
- **(c) θ₃ on one side, Gauss-sum polynomials on the other.**
  - One Poisson step yields only one γ. Putting the second γ₁(A) into the C–S partner kills absorption off the diagonal A = A′.
  - The only working split is a **two-variable probe**: I = Σ_{m₁,m₂} A_{m₁}A′_{m₂}B_{m₁m₂}.
  - Here B_n = Σ_A τ(A)κ_A(n) is the single cubic theta, and the two Poisson steps give γ₂(A)²·τ(A) = μα.
  - Its row is P78's j_p = 2 case: pure conductor N(n)², rigorous, and full-index with a cubic dual (no residual quadratic twist). **This is the only c_B = 2 cubic variant.**
  - Cost: the modulus A enters both Poisson steps.
    - In P2 shape: ρ_O ≥ (l−f₁)+(l−f₂) = 2−h > 0 for h < 2, so θ ≥ 7/6 (PROVED within the move calculus).
    - In P78 shape the high-side frequency count becomes (Z/X₁)(Z/X₂). That is an extra factor Z, a GRH-rows loss of Z^{1/2} (HEUR: probably cancels the gain).
    - It needs its own exponent model; not done here.
- **(d) Square rows h = k².** B at index k⁸ ≡ k² has conductor N(k)², but no probe produces square-indexed duals. The family over squares is the conjugate family, whose dual is conj(B), so c = 4.

## 7. Function-field by-product (NUM, not used for conclusions)

- Over F_q[T], for squarefree monic c: **g_n(1,c) = (−1)^{deg c − ω(c)} G_q(χ₀)^{deg c} χ₀(Res(c,c′))**, by Hasse–Davenport plus the product formula for Res.
  - Verified by brute force: 21/21 cases (`ff_validate.py`).
  - With Pellet's formula, (−1)^{d−ω} = (disc c / q). So every Kubota coefficient is a character of the discriminant.
- `ff_kubota.c` (7^9 polynomials in 16 s) shows the predicted p-part pattern: ψ(T^l) has its non-polar terms at degrees 0 and l+1 only.
- At q = 7, −1 is not a sixth power, and the place at ∞ adds a pole near s = 5/6 (coefficient growth ≈ q^{0.86 d}). So FE checks there would need the ∞-scattering, which was not done.

## 8. Files

| file | content |
|---|---|
| `local_fe.py`, `local_fe.out` | symbolic local FE of coprime parts for n = 2, 3, 6, all l (sec 2) |
| `orient_check.py`, `orient_check.out` | Z[ω]: sextic p-part table at p = 7, 13; orientation test R4 vs R2 (sec 5) |
| `ff_validate.py`, `ff_kubota.c`, `ff_analyze.py`, `ff_D9.txt` | function-field by-product (sec 7) |
