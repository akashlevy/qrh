# Attack D: redesigning the 7/8 probe to remove the reflection scalar ζ_m

*9 Oct 2026. This is Attack D of the (V′) effort in `TARGET_VLS.md`. It builds on:*
- *`paper78.tex`: eq. `general-probe` (l.3397), `low-separated` (l.3470), `probe-poisson-identity` (l.3831), Lemma `local-euler`, Prop. `completed-reflection` with eqs. `reflection`, `reflection-local-fourier`, `reflection-multiplier`, `reflection-row-phase`, `reflection-cross-prime-phase`, and eq. `signal`;*
- *`wf/V_gaps.md` (Lemma Z, §1.6);*
- *`wf2/A2_theory.md` with its referee report.*

*Scripts in this directory (run with the scratchpad venv):*
- *`D_class_table.py` (N1): the Gauss-class / Möbius-sign table;*
- *`D_cross_term.py` (N2): the s–A cross term of the Poisson coefficient;*
- *`D_lemmaZ_general.py` (N3): Lemma Z for every twist exponent j.*

*LP: `wf/lp_const_saving.py`, rerun at s₀ = 0.*

Labels:
- **[R]** rigorous, given the paper's displayed formulas and cited results;
- **[R-]** routine, not fully written;
- **[S]** taken from a secondary statement;
- **[H]** heuristic;
- **[G]** gap;
- **[N]** numerical evidence.

---

## 0. Verdict

**Attack D fails, and the failure is structural.**

Consider the probe's "slots" that can be changed while keeping the high side:
- the power a of the character in the s-Gauss sum;
- the twist exponent j in χ_A(m)^j;
- the column coefficient κ(A);
- the weights in s and m.

Within this space the paper's probe is **rigid**, up to complex conjugation and ray/Hecke twists (Theorem 6.1). In every admissible design the reflected theta pairing is again

  Σ_μ β_μ Σ_m ρ(m) · ᾱ(m)² γ₁(m) γ₃(m) · (ray) · (λ⁴μ/m)₂,

or its complex conjugate. So F1 comes back unchanged. No candidate gives a saving s₀ > 0, and the LP stays at **0.874957** (rerun today).

The new results, in order of usefulness:

1. **The trade-off theorem (§5) [R, given the paper's formulas; N].** For the paper's completed row, a uniform twist exponent j on the free row primes produces:
   - the column kernel χ_p^{−(j+2)};
   - the row scalar ζ_R ∝ ᾱ(R)² γ_j(R) γ_{j+2}(R), which has Gauss class 2j+2 mod 3 (Lemma Z is the case j = 1).

   Consequences:
   - **The signal (Möbius absorption) holds if and only if j = 1.**
   - **The row scalar is ζ-free (ray × Hecke) if and only if j ∈ {2, 5}.**

   The two conditions are incompatible. The twists for which the Voronoi + Heath-Brown argument would work are exactly the ones that detect no zeros. The signal also forces the column kernel to be quadratic, so the short-dual rows of the A₂ conductor table (j_col = 4, 5) are unreachable.
2. **The Euler-product constraint j ≡ a (mod 6) (§3.3) [R-, N].** The s–A cross term of the Poisson coefficient is χ̄_s(A)^a χ_A(s)^j = (ray)·χ_A(s)^{j−a}. The paper's R(A,s) insertion cancels it only when j ≡ a. Otherwise the w = 1 residue loses the Möbius A-sum, and with it the signal.
3. **Gauss-class calculus (§2) [R, N].** γ_j γ_i is Hecke × ray if and only if i + j ≡ 0 (mod 3). It carries the Möbius sign if and only if {i, j} = {1, 2} or {4, 5}. This is checked on 1422 primes to 2·10⁻¹⁵.
4. **Rigidity (Theorem 6.1) [R-].** Constraints (H2) (scalar Euler product), (H3) (Möbius signal) and (L1) (shortening theta reflection) together force (a, j, i) ∈ {(1,1,2), (5,5,4)}. These are the paper's probe and its conjugate.
5. **Candidate-by-candidate analysis (§7).** This covers every modification named in the brief, plus several more. The closest near miss is a cubic s-slot (C5). It makes the reflected row sequence *Möbius-explicit* (π = μ·ψ), so (V′) would become a zero-density statement. But C5 violates (H2) and (H1).
6. **A structural reading (§8) [H].** Within the admissible designs, the Möbius information of the signal has to sit in the reflected row sequence. It is either explicit (π = μψ: a zero-density problem) or hidden in a Gauss sum of nontrivial class (π = ρ·ᾱ²γ₁γ₃: F1). The class arithmetic places every admissible design on the hidden side.

**Status: failed** (no saving; θ = 0.874957). Attack D is now closed in the following precise sense. A probe that removes F1 must give up one of (H2), (H3) or (L1), and §3–§4 show what each loss costs:
- dropping (H2) loses the signal;
- dropping (H3) loses the signal;
- dropping (L1) loses the ‖B‖² ≪ Q bound, i.e. a power Z^{(1−M)/2} on the low side.

---

## 1. The design space

Strip the normalisations from the probe (eq. `general-probe`). It is a triple sum

  I = Σ_s c(s) W₁(q_s/Y) Σ_m w(m) g̃_a(s,−m) W₀(q_m/(q_s X)) Σ_{A = c n³} κ(A) χ_A(m)^j · (V_G-weight in q_A/Z) / √q_A,

where g̃_a(s,h) = q_s^{−1/2} Σ_{d mod s} χ_s(d)^a e(hd/s).

**The paper's choice:**
- a = 1, j = 1;
- κ(A) = γ₂(c) ᾱ(A) Ḡ(A) ν(A), which has Gauss class 2 (§2);
- c(s) = 1_σ χ_s(b_*)/(τξ(s)), a ray function;
- w(m) = ξ(m), a ray function.

**Three views of the same I:**
- **High side.** Poisson in m (eq. `probe-poisson-identity`) gives the scalar Euler product of Lemma `local-euler`. The principal row u = 1 is ζ_F(6z) ζ_F(w)/L(x,η) · 𝓗. Its residue at w = 1, z = 1/6 is the signal f_η.
- **Low side.** I = Q^{−1/2} Σ_m A_m B_m. Here A_m is the s-sum, with theta part Y^{−1/3} ρ(m), and B_m is the A-sum, reflected by eq. `reflection`.
- **Theta pairing.** T = Σ_m ρ(m) B_m = Σ_μ β_μ Σ_m ρ(m) ζ_m (λ⁴μ/m)₂ (`TARGET_VLS.md` §1).

**Requirements for a redesigned probe.** Items (a)–(b) of the brief are made precise as follows.

| | requirement | why |
|---|---|---|
| (H1) | only u = 1 is principal in w | otherwise further rows carry ζ(w)-poles, and with them further "signals" |
| (H2) | the Poisson coefficient factors into a scalar Euler product | Lemma `local-euler`; needed to isolate 1/L^S(x,η) |
| (H3) | the principal row contains 1/L(x, η′), η′ in the target family | the zero detector |
| (H4) | Poisson in m is available and keeps the ζ_F(6z) pole | the z = 1/6 residue is half of the signal |
| (L1) | B_m has a reflection to length ≤ L = Q²/Z < Q, so that Σ_m \|B_m\|² ≪ Q^{1+ε} | Part I's completed moment, the Cauchy–Schwarz scale Q^{11/12} |

**How the brief's candidates map onto the parameters:**
- C1: a weight in m that cancels ζ_m, i.e. the parameter w;
- C2: conjugate or other orientation, i.e. the parameters (a, j);
- C3: cubic instead of quadratic structure in the reflection, i.e. the parameter j through the column kernel;
- C4: a different μ-absorbing combination, i.e. the parameters (j, κ).

Section 7 adds C5–C9.

---

## 2. Gauss-class calculus

**The group 𝓗.** Let 𝓗 be the multiplicative group of unimodular functions on squarefree good ideals generated by:
- Hecke characters, including the angle α;
- functions on a fixed ray group, including the quadratic refinements G of `lem:fixed-gauss-phase`;
- the Möbius sign μ.

All twisted-multiplicativity cross terms are taken as in the paper. A function f has **Gauss class** cl(f) = k ∈ Z/3 if f·γ₁^{−k} ∈ 𝓗.

**Lemma 2.1 [R, given Heath-Brown–Patterson].** cl(γ_j) = j mod 3. Moreover γ_j ∉ 𝓗 for j ∈ {1, 2, 4, 5}.

*Proof.*
- At primes, γ_a γ_b = (J(χ^a, χ^b)/√q) γ_{a+b} whenever χ^{a+b} ≠ 1, and J is a Hecke character up to sign (Weil).
- γ₃ is a ray function (eq. `fixed-gauss-phases`, γ₃ = χ(4)G).
- γ_j γ_{−j} = χ(−1)^j.
- Cross terms: γ_j(ab) = (χ_a(b)χ_b(a))^j γ_j(a)γ_j(b), and the cube of (χ_a(b)χ_b(a)) is R(a,b), which lies in 𝓗.
- γ₂ ∉ 𝓗 by the Heath-Brown–Patterson equidistribution of cubic Gauss sums at primes [S]. ∎

**Lemma 2.2 (Möbius sign) [R].** For i, j ∈ Z/6, γ_j γ_i ∈ 𝓗 if and only if i + j ≡ 0 (mod 3). In that case:
- γ_j γ_i = μ·(Hecke)·(ray) **if and only if {i, j} = {1, 2} or {4, 5}**: eq. `signal` and its conjugate; the sign is the Jacobi sign J ≡ −1 mod (1−ζ)²;
- γ_j γ_i is Hecke × ray with no μ in the other cases (i + j ≡ 0 mod 6, or {i, j} = {0, 3}).

**Numerics (N1, `D_class_table.py`).** Setup:
- 1422 split primes of norm ≤ 12000, classes mod 72;
- products with i + j ≡ 3 (mod 6) are divided by the quadratic refinement G before the multiplicativity test, since G is not a character.

| j | i | i+j mod 6 | best α-power | within-class spread | c(x)c(y)/c(xy) | verdict |
|---|---|---|---|---|---|---|
| 0 | 3 | 3 | 0 | 2·10⁻¹⁶ | +1 | Hecke × ray (G), no μ |
| **1** | **2** | 3 | +1 | 2·10⁻¹⁵ | **−1** | **μ × Hecke × ray** |
| 1 | 5 | 0 | 0 | 2·10⁻¹⁵ | +1 | Hecke × ray |
| 2 | 4 | 0 | 0 | 2·10⁻¹⁵ | +1 | Hecke × ray |
| 3 | 3 | 0 | 0 | 3·10⁻¹⁵ | +1 | Hecke × ray |
| **4** | **5** | 3 | −1 | 2·10⁻¹⁵ | **−1** | **μ × Hecke × ray** |
| all others, incl. single γ₁, γ₂, γ₄, γ₅ | | ≢ 0 (mod 3) | none | 5.81 rad (random) | — | nontrivial Gauss class |

---

## 3. High-side constraints

### 3.1 (H4): the m-weight must be a ray function [R-]

The Poisson identity needs m ↦ w(m)·g̃_a(s,−m)·χ_A(m)^j to be periodic modulo s·b_*·A times a fixed modulus. Two kinds of weight break this.

**A Gauss sum with modulus m** (for example γ̄₁(m)γ̄₃(m), the Gauss part of ζ̄_m):
- it is not periodic, so Poisson is unavailable;
- the m-node becomes a Kubota node and the coefficient is no longer an Euler product.

**An angle weight α(m)^k with k ≠ 0** (for example the ᾱ(m)² inside ζ_m):
- Poisson goes through, but the planar Fourier transform of α(y)^k W(|y|²/R) is α(ξ)^{∓k} times a Hankel transform;
- so the frequency H = u a⁶ carries α(H)^{∓k}, and the a-sum becomes L^S(6z, α^{∓6k}) in place of ζ_F^S(6z);
- this is entire, because α^{6k} is a nontrivial Grössencharacter;
- so there is no z = 1/6 residue and no signal.

### 3.2 (H1): the s-slot must be sextic, a ∈ {1, 5} [R-]

After Poisson, the s-node of row u is L^S(w, c·χ_•(u)^a) (the W-factor of Lemma `local-euler`, with χ_p(u) replaced by χ_p(u)^a). It is principal exactly when χ_•(u)^a is trivial:
- a ∈ {1, 5}: only u = 1;
- a ∈ {2, 4}: u = b³ (cube rows);
- a = 3: u = b² (square rows).

**The s-weight c(s) must also be ray.**
- A Gauss-sum weight in s turns every row's s-node into a Kubota series. Each row then has a theta pole, so every row is "principal".
- An angle weight α(s)^k with k ≠ 0 removes the ζ(w) pole, i.e. the signal.
- By the Heath-Brown–Patterson mechanism [S] it plausibly also removes the theta pole of A_m, i.e. the Gram excess. This is an instance of the signal–excess link of `REPORT.md` §0: excess and signal disappear together.

### 3.3 (H2): the scalar Euler product forces j ≡ a (mod 6) [R-; N]

**The cross term.** Take s = r and A = p prime, and H a unit at both. Eq. `local-F`, with the s-character χ_r^a and the A-twist χ_p^j, gives

  F(r, p, H) = χ_r(H)^a χ̄_p(H)^j G_j(p) · X(r,p),   X(r,p) = χ̄_r(p)^a χ_p(r)^j = R(p,r)^{−a} χ_p(r)^{j−a}.

The paper's inserted R(A,s) cancels R^{−a}. The factor χ_p(r)^{j−a} is a genuine sextic symbol coupling A and s, so it can be cancelled only if j ≡ a.

**What happens when j ≢ a.** In the principal row, for squarefree A ≠ 1 the s-sum is L(w, c·χ_A^{j−a}), which is entire. So the w = 1 residue comes only from A without squarefree part:
- A = 1 always;
- A = n³ as well when j − a is even.

Neither carries the Möbius factor μ(c). So 1/L(x,η) disappears from the residue and the signal is lost [R-].

**Numerics (N2, `D_cross_term.py`).**
- The closed form above matches a brute-force finite Fourier transform in 402 cases, with relative error ≤ 6.2·10⁻¹⁵.
- The cross term X(r,p) over 150 × 150 prime pairs has within-class spread (mod 12) as follows:

| (a, j) | j − a mod 6 | spread | |
|---|---|---|---|
| (1,1), (5,5), (2,2) | 0 | ≤ 6·10⁻¹⁵ | ray, cancellable |
| (1,5), (5,1), (2,1), (1,2), (4,1) | ≠ 0 | 4.2–6.3 rad | not cancellable |

**Repair by insertion fails [R-].** One could insert χ_A(s)^{a−j} into the completed row to cancel the cross term. Then the A-sum becomes the completed row at the index m^j s^{a−j}. The low side no longer factors as Σ_m A_m B_m, and the s-theta structure is lost.

Example: a = 2, j = 1 gives B_{ms}.
- The s-sum collapses to a divisor sum.
- B is evaluated at rows of norm ≍ QY.
- Cauchy–Schwarz loses Y^{1/2} relative to the paper.
- Moreover the cube rows u = b³ become principal with denominators 1/L(x, η·ψ_b), ψ_b quadratic. These are uncontrolled extra signals of size ≫ f_η, against eq. `high-probe-contract`.

### 3.4 (H3): Möbius absorption [R]

The signal term of the local table (row (1, 2r), k = 0, j = 0 of eq. `local-family-table`) has coefficient γ_j(p)κ(p)Q^{−x}. It must equal −η′(p)Q^{−x}.

**Gauss-type κ.** With κ = γ_i·(Hecke)·(ray), Lemma 2.2 gives

  **(j, i) ∈ {(1,2), (2,1), (4,5), (5,4)}.**

**κ with an explicit Möbius factor.** With κ = μ·γ_{−j}·(Hecke), the completed row Σ_A μ(A)γ_{−j}(A)χ_A(m)^j is a Möbius × Gauss sum. It is not a theta L-series [R-], so (L1) fails:
- Patterson's coefficients d(nb³) = 3^{5/2}|b|g̃₃(n) carry no Möbius factor;
- μ·γ_i is not in γ_i·(Hecke·ray), since μ is not itself a Hecke × ray function;
- for example μγ₂ ∝ γ₁²·(Hecke) and μγ₄ ∝ γ₂²·(Hecke) are squares of Gauss sums.

---

## 4. (L1): the completed row must be a cubic-theta L-series [R- / H]

**What the low side needs.** The paper's low side needs Σ_{q_m ≲ Q} |B_m|² ≪ Q^{1+ε} with Z = Q^{1/M} > Q. It obtains this from the reflection, which has length L = Q²/Z < Q, together with the quadratic large sieve: (L + Q)Σ|β|² ≍ Q.

**What this requires of κ.** κ must be the coefficient sequence of an automorphic theta function whose L-series has a Hecke-type functional equation with Voronoi twists. With known tools this means Patterson's cubic theta: θ̄ for class i = 2, or θ for class i = 4.

**The sextic alternative fails.** For i ∈ {1, 5} (κ = γ₁ or γ₅), write γ₁(A)χ_A(m)² = g̃₆(A, m⁴)·1_{(A,m)=1}. Then B_m is a length-Z partial sum of the sextic Kubota series D(·; m⁴).
- Its approximate functional equation has length q_m² ≍ Q² > Z, so there is no shortening [R-].
- The best available mean-value bound is the cubic large sieve, Σ_m |B_m|² ≪ (Q + Z + (QZ)^{2/3})Z^ε. For M < 2 this is ≥ (QZ)^{2/3} > Q.
- Even the random model gives max(Q, Z) > Q for M < 1.
- So the low side loses at least Z^{(1−M)/2} [R- for the bound; H for "no better"].

---

## 5. The trade-off theorem: generalised Lemma Z

Take the paper's completed row (θ̄, class i = 2) and let every free row prime be active with the same exponent j ∉ {0, 4}. Prop. `completed-reflection` gives three ingredients [R]:
- **Column kernel.** B_p(x) = χ_p(x)^{−(j+2)}. The shift +2 comes from the multiplier, eq. `reflection-multiplier`.
- **Row scalar.** ζ = (fixed)·ᾱ(c)² ∏_p χ_p(σ_p)^{−2} ω_{p,j}, with ω_{p,j} = τ⁻_{p,j} τ⁺_{p,j+2} χ_p(ε_p)^{−j−2} and τ⁻_{p,m} = χ_p(−1)^m γ_m(p).
- **Cross-prime exponent.** χ_p(σ_p)^{−2} ω_{p,j} = ξ_{p,j} χ_p(c/p)^{2j+2} (eq. `reflection-cross-prime-phase`).

**Theorem 5.1 [R, given those formulas; N].** For squarefree R with c_F fixed,

  ζ_R = κ_j · ψ_j(R) · ᾱ(R)² · γ_j(R) · γ_{j+2}(R),

with κ_j a constant and ψ_j a ray function. Hence cl(ζ) = 2j + 2 (mod 3).

*Proof.*
- Write ∏_{p|R} γ_j(p) = γ_j(R) ∏_{p<p′} (χ_p(p′)χ_{p′}(p))^{−j}, and likewise for γ_{j+2}.
- The moving cross factor is ∏_{p≠p′} χ_p(p′)^{2j+2}.
- The total pair exponent is 2j + 2 − j − (j + 2) = 0.
- The remaining factors χ_p(−1), χ_p(λ), χ_p(c_F) are ray functions. ∎

For j = 1 this is Lemma Z: ᾱ²γ₁γ₃ = (ray)·γ̄₂⁴.

**Numerics (N3, `D_lemmaZ_general.py`).** On 1212 products R = p₁p₂ with N(R) < 20000, with γ_j(R) computed by brute force over O/R, the ratio of the two sides is a function of the classes of (p₁, p₂) mod 12. The spread is ≤ 5·10⁻¹⁵ for j = 1, 2, 3, 5.

**The full table for the paper's row (θ̄, i = 2).**

| j | column kernel χ^{−(j+2)} | γ_j γ₂: signal? | ζ ∝ | cl(ζ) | ζ-free? |
|---|---|---|---|---|---|
| 0 (active) | χ^{−2}, cubic (·q^{−1/2}) | no coupling | γ₂ | 2 | no |
| **1** | **χ³, quadratic** | **μαG: yes** | **γ₁γ₃** | **1** | **no** |
| 2 | χ², cubic | μαγ₄: no (class 1) | γ₂γ₄ = 1 | 0 | **yes** |
| 3 | χ, sextic | class 2: no | γ₃γ₅ | 2 | no |
| 4 | Ramanujan | γ₄γ₂ = 1: no μ | γ₄ | 1 | no |
| 5 | χ^{−1}, sextic | class 1: no | γ₅γ₁ = χ(−1) | 0 | **yes** |

**Corollary 5.2 (incompatibility) [R].**
- The probe detects zeros (H3) if and only if j = 1, i.e. j ≡ 1 (mod 3).
- The theta pairing is ζ-free, so that Voronoi + Heath-Brown would apply with s(M) = −½max(M−1, 7M/3−2, −M/6), if and only if j ≡ 2 (mod 3).

These never hold together.

The conjugate row (θ, i = 4) behaves the same way: the signal needs j = 5 and ζ-freeness needs j ∈ {1, 4}.

**In general.** Suppose the completed row is the L-series of a theta function whose coefficients and multiplier both have class i. (This is the "same symbol" principle of Patterson and Kazhdan–Patterson. The paper's own computation gives shift 2 for the class-2 row θ̄.) Then:
- cl(ζ) = 2j + i;
- the signal needs i + j ≡ 0 (mod 3) with j ≢ 0, which forces cl(ζ) = j ≢ 0;
- ζ-freeness together with the signal would force j ≡ i ≡ 0, and then no Gauss sum is present and there is no Jacobi sign, hence no μ.

**Corollary 5.3 (the quadratic kernel is forced) [R].** The signal forces the column kernel to be quadratic: j = 1 gives χ³. In the A₂ conductor table (`A2_theory.md` §2.4), the rows with small conductor (j_col = 4: N(μ)²; j_col = 5: N(μ)) need j ∈ {0, 5} here. That holds up to the j_col ↔ 6 − j_col orientation caveat of that table; neither member of the swapped pair, {0, 5} or {2, 3}, contains j = 1. None of these is Möbius-absorbing.

In particular j = 5 gives a **ζ-free row with a sextic kernel**, which is the best low side one could wish for. It belongs to a probe with **no signal**.

---

## 6. Rigidity

**Theorem 6.1 [R-].** Let P = P(c, a, w, j, κ) be a probe of the form in §1, where:
- c and w are ray functions (Hecke twists allowed);
- a, j ∈ Z/6;
- κ ∈ {γ_i·(Hecke)·(ray), μ·γ_i·(Hecke)·(ray)}.

If P satisfies (H2), (H3), (H4) and (L1), then (a, j, i) ∈ {(1,1,2), (5,5,4)}, i.e. the paper's probe or its complex conjugate. On squarefree rows the theta pairing is then

  Σ_μ β_μ Σ_m ρ_a(m) ζ_m (λ⁴μ/m)₂ W,   with ζ_m = (ray)·ᾱ(m)²γ₁(m)γ₃(m) (resp. its conjugate),

and ρ_a is the sextic theta coefficient (resp. its conjugate). So F1 holds for every admissible P.

*Proof.*
1. (H4) forces w to be a ray function (§3.1).
2. (H3) with Lemma 2.2 excludes Gauss-type κ with i ∈ {0, 3}: a Hecke-type completed row gives no Jacobi sign. (Such a row would satisfy (L1) trivially: Σ_A φ(A)χ_A(m)^j over A ≍ Z > Q is negligible by Poisson in A.) For the other Gauss-type κ, (H3) leaves (j, i) ∈ {(1,2), (2,1), (4,5), (5,4)}.
3. (L1) excludes the explicit-μ type (§3.4) and i ∈ {1, 5} (§4). So i = 2 with j = 1, or i = 4 with j = 5.
4. (H2) gives a = j (§3.3).
5. Theorem 5.1 (and its conjugate) gives ζ. ∎

**Remarks.**
- (H1) is not needed; it holds automatically.
- Ray twists and Hecke twists (on c, κ, w and ν) change ζ only by Hecke × ray. They cannot change its Gauss class.

---

## 7. The candidates, one by one

| # | modification | Möbius absorption / high side | low side | verdict |
|---|---|---|---|---|
| C1 | weight w(m) cancelling ζ_m (Gauss part γ̄₁γ̄₃, angle α²) | Gauss part: no Poisson; angle: L(6z, α^{∓12}) entire, no z = 1/6 pole (§3.1) | would give the ζ-free sum | ✗ (H4) [R-] |
| C2a | global conjugation (a, j, κ) = (5, 5, γ₄) | ✓ | π̄-type pairing, ζ ∝ α²γ₅γ₃ | admissible, **same obstruction** |
| C2b | partial conjugation: twist j = 5 with a = 1, κ = γ₄ | local μ ✓, but cross term χ_A(s)⁴ (N2) | — | ✗ (H2) [R-] |
| C2c | j = 5 with the paper's θ̄-row | γ₅γ₂ has class 1: no μ | ζ-free, sextic kernel (the best low side) | ✗ (H3) [R] |
| C3 | cubic column kernel (j ∈ {0, 2}) | j = 2: γ₂² = μαγ₄ has class 1; j = 0: no coupling | j = 2 is even ζ-free | ✗ (H3) [R] |
| C4 | other μ-absorbing pairs (j, i) = (2,1), (4,5) | needs a = j ∈ {2,4} (H2), and then cube rows; κ is sextic | B is a sextic Kubota sum, no shortening | ✗ (L1) [R-/H]; see C7 |
| C5 | cubic s-slot a ∈ {2,4}, j = 1, κ = γ₂ | cross term χ_A(s)^{1−a} (N2); cube rows u = b³ principal with 1/L(x, ηψ_b) | **near miss**: for one a ∈ {2,4}, cl(ρ_a) = −1 cancels cl(ζ) = 1, so π = ρ_a ζ = μ·ᾱ·(ray)·q^{−1/6} is Möbius-explicit | ✗ (H2), (H1) [R-] |
| C6 | angle twists: α^k on s / on A / on m | s: kills the ζ(w) pole (and plausibly the excess); A: changes the target to a Grössencharacter; m: as C1 | A-twist changes ζ by a Hecke character only | no gain |
| C7 | "swapped" probe (a, j, i) = (2, 2, 1) | (H2) ✓; cube rows give coherent 1/L(x,η), producing a new pole at z = 1/3 with signal X^{1/6}Z^{x−2/3} [R-] | ‖A^sig‖² excess P^{1/3}; B a sextic Kubota sum; Heath-Brown cubic LS ⇒ the low side needs β > 1; with an ideal LS at best β > 5/6 + … before the high side | ✗ [H] |
| C8 | regrouping the low side (Cauchy–Schwarz over s or over μ), or combining probes | — | rows s: dual of the additive sieve, same bound; rows μ: that *is* (V′_LS); combinations: signal–excess link (`REPORT.md`) | no gain [H/R-] |
| C9 | handle small μ separately (A₂ FE: F′(μ) negligible for N(μ) ≤ Q^{1/3−ε}) | — | each dyadic μ-block carries ≍ the full CS mass (Σ_{μ∈block} \|β_μ\|² ≍ 1), so the blocks near L = Q^{0.79} still give Q^{11/12} | no power saving [R- given the A₂ FE (H)] |

**Notes.**

**C1.** ζ_m contains the angle ᾱ(m)². Even this Hecke part cannot be moved into the m-weight. Its Gauss part γ₁γ₃ is a Gauss sum of modulus m.

**C5.** This is the only design found in which the reflected row sequence is "simple".
- The class arithmetic needs the s-theta to be cubic with class −j = 2.
- (H2) then fails because a ≢ j.
- The repair of §3.3 costs Y^{1/2} and creates uncontrolled signals 1/L(x, ηψ_b) from the cube rows.

Which of a = 2 or a = 4 gives the cancellation depends on the orientation of the cubic theta coefficients. Bröker–Hoffstein's two displayed statements point in opposite directions: τ₃(d) = g₃(1,d)/√N on one hand, but τ₃(mp) = G₂(m,p)τ₃(m) from their Hecke relation on the other (`scratchpad/bh.txt` l.209 vs l.254). The statement above is orientation-free.

**C7.** This is a genuinely different probe. It is cubic where the paper is sextic and sextic where the paper is cubic. Two obstacles:
- **The B-norm.** It needs a beyond-the-edge cubic large sieve, without the (QZ)^{2/3} term, for sextic Kubota partial sums. De Faveri–Dunn–Hoffstein show that term is sharp for general coefficients.
- **Circularity for M > 1.** Its unreflected theta pairing is the signal plus a cubic dual of length Z/Q < 1 (the KP collapse τ₃(mp²) = 0, `A2_theory.md` §3.2). So for M > 1 the theta part *is* the signal, which makes the low side circular.

The rows u = u₀b³ also need a full redo of the zero detector, which has not been done [G].

---

## 8. Why: where the Möbius information goes [H, with R ingredients]

**The signal is always inside T.** For every admissible design the unreflected theta pairing is T = 𝒮 + T′ (`V_gaps.md` §1.6).
- 𝒮 is the Möbius sum over A coming from the coherent main term.
- Its contribution to I is Q^{−1/2}Y^{−1/3}Q^{5/6}·Z^{x−5/6}·(…) = X^{1/3}Z^{x−5/6}·(…).
- This is exactly the scaling X^{1/3}Y⁰Z^{x−5/6} of the high-side principal residue (`W(X,Y,Z; x, 1, 1/6)`) [R- as exponent bookkeeping].

So 𝒮 is the signal and cannot be split off: the continuation criterion needs |I| small, and I − 𝒮 is signal-free on both sides.

**Two ways to hide it.** In a reflected form that does not display 𝒮, the Möbius information must be carried by the row sequence π:
- **explicitly**, π = μ·ψ (only in C5-type designs). Then (V′_LS) becomes Σ_{μ≍L}|Σ_{m≍Q} μ(m)ψχ_μ(m)W|² ≪ Q^{2−δ}, a zero-density statement for the family {L(s, ψχ_μ)}_{N(μ)≍L};
- **implicitly**, through a Gauss sum of nontrivial class, π = ρ·ᾱ²γ₁γ₃ (F1).

**What the explicit form would need.** Inside the contradiction argument one may use that every zero has real part ≤ β_*. Then the explicit form needs a density exponent c(σ) < 2/ℓ on σ ∈ [½, β_*], where ℓ = log L/log Q = 2 − 1/M ≈ 0.79, i.e. c < 2.53.
- The density hypothesis c = 2 would suffice.
- Montgomery-type large-sieve exponents c(σ) = 3/(2 − σ) suffice only for σ < 2 − 3ℓ/2 ≈ 0.82, short of β_* ≈ 0.87 [H/S].

The A₂/Sawin–Whitehead picture agrees: the reflected row node acquires M_mm = 4, i.e. γ₁γ₃ (`A2_theory.md` §4.2). The class calculus of §5 shows that every admissible design lands on the implicit side.

---

## 9. LP

No candidate supplies a uniform saving s₀ > 0, so the LP is unchanged. Rerun of `wf/lp_const_saving.py` at s₀ = 0 (`th("const", 0.0)`):

  θ = **0.874957** at (l_x, l_y, ℓ) = (0.3549, 0.4783, 0.1668), M = 0.8332, b = 0.1234.

This is Liu's value. The reference table s₀ ↦ θ of `TARGET_VLS.md` §2 stands: 0.002 → 0.8733, 0.01 → 0.8668, 0.04 → 0.8532.

---

## 10. Conclusions and what is left

**Attack D is closed within the probe family of §1.** F1 is not an artefact of the paper's choices. It is forced by the same Davenport–Hasse alignment (γ₁γ₂ = μαG at the Poisson node) that produces the signal. Two facts make this precise:
- Theorem 5.1 / Corollary 5.2: the signal needs j ≡ 1 (mod 3) and the ζ-free pairing needs j ≡ 2 (mod 3);
- Theorem 6.1: no other slot (a, w, κ, c) can repair this without violating (H2), (H3), (H4) or (L1).

**Gaps and caveats.**
- **[G]** Probes outside the family are not classified: for example several Poisson variables, several completed rows, or a signal not produced by a Davenport–Hasse collapse at a Poisson node. The partial checks are §7 C8 and the coupled-row remark in §3.3; the latter loses its factorisation.
- **[R-/H]** (L1) for sextic completed rows rests on "no known shortening plus large-sieve bounds". A beyond-the-edge cubic large sieve for sextic Kubota sums would reopen C7. That is an input of the same difficulty as (V′_LS).
- **[S]** That an angle twist on s also removes the theta pole of A (§3.2) is quoted from the Heath-Brown–Patterson mechanism. It is not needed, because the loss of the w = 1 pole already excludes the candidate.
- **Orientation.** The orientation of the cubic theta coefficients is convention-sensitive (Bröker–Hoffstein, §7 C5). All conclusions here are orientation-free.

**What remains for the (V′) program.** The only routes left are Attacks A–C on the actual pairing:
- (V′_LS) or (V′_β) for π = ρ·ᾱ²γ₁γ₃, i.e. a beyond-the-edge large sieve with few moduli;
- or a probe of an essentially different shape.

The most informative by-product of this attack is Corollary 5.2's twist table, which should be cited whenever a "change the twist" repair of F1 is proposed.

---

## Referee report

*Independent skeptical review, 9 Oct 2026. Scope:*
- *this file;*
- *`paper78.tex`: Prop. `completed-reflection` (l.1677–1860), the multiplier and cross-prime phase (l.2170–2290), eq. `signal` (l.965), eq. `general-probe`/`low-separated` (l.3383–3480), eq. `local-F` (l.3742), Lemma `local-euler` and its table (l.4010–4130);*
- *`wf/V_gaps.md` §1.2–1.7 and `wf2/A2_theory.md` §2.4.*

*Re-runs (scratchpad venv):*
- *`D_class_table.py`, `D_cross_term.py`, `D_lemmaZ_general.py`: all outputs reproduced to the digit (1422 primes / 402 cases / 1212 pairs, spreads ≤ 6·10⁻¹⁵).*
- *`wf/lp_const_saving.py` at s₀ = 0: θ = 0.874957 at (l_x, l_y, ℓ) = (0.3549, 0.4783, 0.1668), M = 0.8332, b = 0.1234. Reproduced.*

*New referee scripts:*
- *`ref_D_n3_power.py`: tests the discriminating power of N3;*
- *`ref_D_j0_node.py`: the j = 0 Poisson node.*

### Verdict

The file claims no saving, and nothing here produces one. **The negative result stands; it is not refuted.** The core mathematics checks out:
- **Theorem 5.1 (generalised Lemma Z)** is correct relative to the paper.
- **The incompatibility in Corollary 5.2** is correct.
- **The reduction of rigidity to (H2) + (H3) + (L1)** is correct.

There are defects in one lemma, in one proof step, in the numerics claims and in several labels. None of them changes the verdict. Details follow.

### 1. What I verified independently

**(a) Theorem 5.1 [R relative to the paper].** Prop. `completed-reflection` explicitly allows arbitrary exponents j_p ∈ {0,…,5}, so the input is not an extrapolation from j = 1. I re-derived the cross exponent by hand:
- ε_p = −λ⁻⁵(c/p)⁻² mod p, so χ_p(ε_p)^{−j−2} contributes χ_p(c/p)^{2j+4};
- χ_p(σ_p)^{−2} contributes χ_p(c/p)^{−2};
- the total is 2j + 2, matching the paper's eq. (`reflection-cross-prime-phase`).

Then τ⁻_{p,j}τ⁺_{p,j+2} = χ_p(−1)^j γ_j(p)γ_{j+2}(p), and twisted multiplicativity γ_k(pp′) = (χ_p(p′)χ_{p′}(p))^k γ_k(p)γ_k(p′) has total cross exponent j + (j+2) = 2j + 2. These cancel exactly, so ζ_R = (ray)·ᾱ(R)²γ_j(R)γ_{j+2}(R).

The i = 4 row is the complex conjugate of the i = 2 row with exponent 6 − j. This gives cl(ζ) = 2j + 4 rigorously for i = 4, not just as the "same symbol" heuristic.

**(b) Lemma 2.2 at i + j ≢ 0 (mod 6), i, j ≠ 0.** The Möbius sign for {1,2} and {4,5} and its absence for {1,5}, {2,4}, {3,3} are consistent with eq. `signal` and the table in Lemma `local-euler`. The paper uses exactly γ₁γ₂ = −α χ̄(4)γ₃ and γ₄/γ₁ = −ᾱ/G there.

**(c) (H2) cross term.** The closed form X(r,p) = χ̄_r(p)^a χ_p(r)^j follows from the CRT lift in eq. `local-F`. With R(A,s) inserted, the residual coupling χ_p(r)^{j−a} is a genuine sextic symbol unless j ≡ a. Reproduced.

**(d) Rigidity, reformulated.** This is a cleaner statement of what the file proves:
- the signal needs {j, i} ∈ {{1,2},{4,5}};
- (H2) needs a = j;
- together these give (a, i) ∈ {(1,2), (2,1), (4,5), (5,4)}.

So **exactly one of the two theta nodes (the s-slot or the completed row) must be sextic and the other cubic.**
- The paper puts the sextic theta in the s-slot. There its coefficients are unknown, but its mean square is known (LoA6).
- The swapped choice puts it in the completed row, which has no known reflection.

F1 is then forced by Theorem 5.1. I find this argument convincing.

**(e) Bookkeeping.** All of these are correct:
- Corollary 5.3's map to the A₂ table: kernel exponent −(j+2) ≡ 4 ⇔ j ≡ 0, and ≡ 5 ⇔ j ≡ 5; the swapped orientation gives {2, 3}; neither set contains 1.
- s(M) at M = 0.826 is ½·0.073 > 0.
- ℓ = 2 − 1/M ≈ 0.79, so 2/ℓ ≈ 2.53.
- The Montgomery threshold is σ < 2 − 3ℓ/2 ≈ 0.82.
- C9: each dyadic μ-block carries Σ|β_μ|² ≍ 1, so dropping small μ cannot save a power.

### 2. Errors and gaps found

**E1. The j = 0 node is misdescribed (Lemma 2.2, Theorem 6.1 step 2). [Error, harmless to the verdict]**
- The paper's convention is χ_p(x)⁰ = 1_{p∤x}. So the j = 0 Poisson node is the Ramanujan sum c_p(H) = −1 at unit H, not "γ₀ = 1".
- `ref_D_j0_node.py` checks this: T = q_r χ_r(H)^a χ̄_r(p)^a·(−1), on 84 cases, rel. err. 4·10⁻¹⁵.
- So j = 0 **does** carry a Möbius sign, with a q_p^{−1/2} deficit. With a unimodular κ the detector becomes 1/L(x+½, ·); recovering 1/L(x, η) needs |κ| ≍ √q_A, which is not a theta row.
- So the stated reason in step 2 ("Hecke-type rows give no Jacobi sign") is wrong for j = 0.
- The conclusion survives anyway, because (H2) with j = 0 forces a ≡ 0.

**E2. Theorem 6.1 is missing a hypothesis. [Gap, harmless]**
- As stated, a ∈ Z/6 is arbitrary and (H1) is dropped. Then (a, j, i) = (0, 0, 0) with κ = η′√q is not excluded by (H2), (H3) or (H4), and (L1) is trivial because B_m is a long smooth Hecke sum.
- That design has no theta pole in s, hence no ρ, no Gram excess and no probe in the paper's sense.
- (a, j) = (3, 3) is excluded only via (L1) for explicit-μ rows, which is [H].
- **Fix:** add the hypothesis (H0): the s-slot carries a sextic or cubic theta pole, i.e. a ∈ {1, 2, 4, 5}.

**E3. N3 is not independent evidence for Theorem 5.1. [Overclaim in numerics]**
- `ref_D_n3_power.py` shows that the N3 ratio is **identically 1** (a single global value), not merely a class function. The exponent 2j + 2 is taken from the paper, so N3 only checks the sign convention of CRT twisted multiplicativity of Gauss sums. It does discriminate that convention: exponent −4 fails for j = 1.
- For j = 2 and j = 5, 2j + 2 ≡ 0 (mod 6). The test is then vacuous: exponents 0 and −6 also pass.
- "Checked numerically for j = 1, 2, 3, 5" should read "the Gauss-sum multiplicativity step is checked; the reflection formula itself is the paper's".
- Likewise the (0,3) row of N1 is trivial by construction, since one divides by G = χ̄(4)γ₃. Only the *relative* sign between (1,2) and (0,3) is informative.

**E4. Label inflation in the (L1) parts.**
- Theorem 6.1 step 3 uses two inputs:
  - "sextic completed rows admit no shortening";
  - "explicit-μ rows are not theta L-series".

  Both are statements about *known tools* ([R-/H] in §4 and [H] respectively). So Theorem 6.1 should be labelled **[R-/H]**: rigidity relative to known reflection formulas. It is not [R-].
- In §4, "even the random model gives max(Q, Z) > Q" is unjustified. A random-phase model for γ₁(A)χ_A(m)² gives ≈ Q. The real obstructions are:
  - no provable shortening or large sieve;
  - possibly a theta main term, since B_m is a partial sum of the sextic Kubota series D(·; m⁴), which has its w = 2/3 pole and a Kazhdan–Patterson-degenerate index m⁴.

  The conclusion "(L1) unavailable" stands; that sentence should be removed.

**E5. "ζ-free" for j ∈ {2, 5} is overstated. [Minor]**
- ζ = ᾱ²·(ray) there still carries the angle ᾱ(m)², and the column kernel is cubic (j = 2) or sextic (j = 5), not quadratic.
- The quoted saving s(M) = −½max(M−1, 7M/3−2, −M/6) was derived for a ray weight and the quadratic kernel. So "Voronoi + Heath-Brown would apply" is [H] for those j.
- This is moot, since neither j has a signal.

**E6. Coverage. [Acknowledged gap; I agree with its scope]**
- I looked for designs outside the slot family. A twist involving both m and m̄, χ_A(m)^j χ_A(m̄)^{j′}, is allowed by the Proposition (arbitrary prime set 𝒫). But it doubles the conductor: c ∝ rad(m m̄), so L_m ≍ Q⁴/Z ≫ Q. (L1) fails, so this is not a loophole.
- The (H2) argument examines only the w = 1 pole. Away from j ≡ a, the principal-row coefficient is Σ_s c(s)q_s^{−w}/L(x, η′χ_•(s)^{j−a}). It has no isolated 1/L(x, η′) anywhere, so this is not a loophole either.
- Several completed rows, Rankin–Selberg-type rows, and signals not produced by a Davenport–Hasse collapse remain unclassified, as the file says.

**E7. Checks requested in the brief.**
- **μ₁₂ ∉ Q(ω):** no 12th-order characters or large sieves are used anywhere in the file, so the issue does not arise. Angle characters α^k are Grössencharacters only for 6 | k, which is how the file uses them (L(6z, α^{∓6k})).
- **Estimator bias:** all numerics test exact algebraic identities by within-class phase spread. Spread ~10⁻¹⁵ versus ~4–6 rad leaves no estimator bias, and the populations (≥ 2–3 primes per class) suffice for identities.

### 3. Recommended edits
1. Lemma 2.2 and Theorem 6.1 step 2: replace "γ₀ := 1" with the Ramanujan node −q^{−1/2}, and exclude j = 0 via (H2) (a ≡ 0), not via (H3).
2. Theorem 6.1: add (H0), a ∈ {1, 2, 4, 5}. Relabel the theorem [R-/H] (rigidity relative to known reflections). State the cleaner form: "one of the s-node and the row node must be sextic".
3. N3: describe it as a check of Gauss-sum multiplicativity. Note that it is vacuous for j = 2, 5 and that the j-general exponent rests on the paper's formula (re-derived by hand above).
4. §4: delete the "random model gives max(Q, Z)" sentence.
5. Corollary 5.2: label "Voronoi + HB would apply" for j ∈ {2, 5} as [H].

### Surviving claim

Attack D fails: no saving, and θ stays at Liu's 0.874957 (LP reproduced).

Rigorous relative to the paper's Prop. `completed-reflection`: for the paper's cubic-theta completed row with a uniform active exponent j ∉ {0, 4}, the reflection scalar is ζ_R = (ray)·ᾱ(R)²γ_j(R)γ_{j+2}(R), of Gauss class 2j + 2. Within the slot family:
- Möbius absorption needs j = 1 (j = 5 for the conjugate row);
- ζ-freeness needs j ∈ {2, 5} (resp. {1, 4});
- so the two never coincide.

Combined with (H2) (a ≡ j), and with (L1) taken as the requirement of a known shortening reflection for the completed row ([R-/H]), this forces one of the two theta nodes to be sextic. With the extra hypothesis a ∈ {1, 2, 4, 5}, the only admissible designs are the paper's probe and its conjugate, so F1 is intrinsic to this family.
