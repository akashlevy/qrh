# wall4 — can the energy wall be relaxed losslessly, W(w) for some w > 0? (10 Oct 2026)

Agent `wall4`. P78 = `openai/.../September-30-2026/build/paper.tex`; line numbers refer to it.
Labels: **[P]** proved (given P78's lemmas and identities), **[B]** proved as exponent bookkeeping,
**[C]** conditional, **[R-]** rigorous modulo routine unwritten lemmas (04's labels), **[H]** heuristic, **[M]** model.

## 0. Verdict

* **No proof of W(w) for any w > 0.** All four suggested routes fail, each at a step identified below
  (§3–§6). Exponent payoff of everything proved here: **0** (worst high exponent at 41/47 stays +0.0114 on the fine grid).
* **New reduction [C on 04's Corollary D, status R-] (§7).** Beyond the wall, the low-side loss of the
  separated probe is confined to **one direction**, the sextic theta vector ρ(m). Write A = A^sig + A^rest
  (04, `gram-split/wf/FE_LoA6_proof.md` §5). Then the A^rest part tolerates a wall excess w ≤ b/6 with *no new input*.
  The only input still missing is a **rank-one wall W_θ(w)**: the single theta pairing
  T = Σ_m φ(q_m/Q) ξ(m) ρ(m) B_m must be ≪ Q^{11/12+ε} past the wall.
  W(w) ⇒ W_θ(w), and with W_θ the low side is lossless for w ≤ b/6.
  W_θ fails at the same step as W (quadratic large sieve for (βX, π) at N = Q·Z^w), and T carries the zero-detecting
  signal, so it is not a cheaper target in practice (§7.4). The model payoff of W_θ is in §7.3.
* **Localization [B] (§2).** Past the wall the only positive branch of lem:reflected-energy is
  (u = z_a, max = dual), with loss at most
  w − O/2 − 2Δ_H − 2(ℓ′ − z_a) − 3ℓ_b − (3/2)(T_d − y).
  So W(w) is equivalent to a lossless bound for one "core" sum (§2) with all slots active, no dual cube part,
  rows in the top Z^{w/2} of the row range, and dual length within Z^{2w/3} of T_d.
* **Model [M] (§8).** A lossless W(w_W) closes 41/47 only for w_W ≈ 0.0415 or more (fine grid; W(0.040) gives +0.0006,
  W(0.043) gives −0.0004). An excess that grows past the wall at slope f gives *nothing* for f ≥ 1/4, and closes only for
  f ≈ 0.1–0.15, at a far vertex (ℓ ≈ 0.24, w ≈ 0.12). Every partial mechanism checked here either has f ≥ 1/2 or does not act at the binding dyad: kernel saving
  (none at v = T_d), Mellin-height averaging (loses T₀^{1/2}), slot transfer (break-even), cube part (f = 1/2), Hölder.
* **Rank-one wall payoff [M].** W_θ alone, for all w, gives +0.0047, worth about W(0.03). W_θ together with a full W(0.02)
  closes (−0.0020). So the full-W requirement drops from 0.0415 to ≈ 0.014.

## 1. The object and the target

Probe, J = ∅ summand (all slots marked; the rescaled summands J ≠ ∅ have excess w − 3d and an extra Z^{−d},
prop:probe-low l.8620). Low separation (eq. low-separated, l.3470):

  I = (Q^{−1/2}/2π) Σ_σ ∫ Ŵ₀(iv) Σ_m Ω(q_m/Q) ξ(m)(q_m/Q)^{−iv} A_{m,σ,v}(Y) B_{m,σ}(Z) dv,   Q = q_{b*}XY ≍ Z^M.

Cauchy–Schwarz in m with ‖A‖² ≪ (Q/Y)(1 + P_a^{1/6} + P_a²/Y) (prop:probe-gram) and
‖B‖² ≪ Z^{M + w + ε}, w := (M + ℓ − 1)_+ (lem:probe-row-norm via lem:reflected-energy + lem:hybrid-energy),
gives |I| ≪ X^{1/2} P_a^{1/12} Z^{w/2+ε}, i.e. the low side

  B_low = 5/6 + M/12 − ℓ/6 + w/2.

After reflection, at the binding dyad (§2), the energy is

  𝒩 = Σ_{k sf, q_k ≍ K} | Σ_P c_P 1_{(P,k)=1} Σ_{n sf ~ N} β(n) ψ_k(n) χ_P(n)^{−2} |²,
  K = Z^{M}, L = q_P = Z^{ℓ}, N = K²L/Z = K·Z^{w}, c_P = a(P) q_P^{−1/2}, |a|, |β| ≤ 1, ψ_k = χ_k³.

lem:hybrid-energy gives 𝒩 ≪ (K + N)·N·Z^ε. **W(w)** is 𝒩 ≪ K·N·Z^ε for N ≤ K·Z^w.
Weighted geometry: ℓ ≈ 0.18–0.20, b ≈ 0.06–0.15, M ≈ 0.82–0.86, so N/K = Z^w with w ≈ 0.04 and L = Z^{0.19} ≫ N/K.

**Model thresholds** (`wm4.py`, `scanpar.py`, `runall.py`, `runW.py`; fine grid δ step 0.005, t 101 pts, r 201 pts):
see table in §8.

## 2. Localization of the excess [B]

Let E_ref be P78's exponent (eq. old-eq:4.6, l.7787) at the probe (f = 1, ρ = 1, N_* = 1 + ℓ′,
2A₀ ≤ O, N₀ ≤ A₀, z_a ≤ ℓ′; l.8292). 03's exact LP certificates (`03-six-sevenths/routes/r4/energy_dual_cert.py`,
R4_partII_map §6b) bound the six branches (u ∈ {v, z_a, (v+z_a)/3}) × (max ∈ {H, dual}) by

  0, (1+2ℓ′−2M′)/3 [refined count; P78's own (1+3ℓ′−2M′)/4], (1−2M′+ℓ′)/3, M′+ℓ′−1, (M′+4ℓ′−2)/3, −M′/2+ℓ′.

At (M′, ℓ′) = (0.85, 0.195) these are 0, −0.103 [−0.029], −0.168, **+0.045**, −0.123, −0.230. Only the branch
(u = z_a, max = v + ℓ_b) is positive. In that branch, with y = v + 3ℓ_b + e_λ ≤ T_d and Δ_H = M′ − O − H ≥ 0,
substituting T_d = 2H + 2A₀ + 2z_a − 1 − ℓ′ − N₀ − 3B₀ and 2A₀ ≤ O gives

  **E_ref − M′ ≤ w − O/2 − 2Δ_H − 2(ℓ′ − z_a) − 3ℓ_b − (3/2)(T_d − y) − S₀ − B₀ + O(ε).**   (2.1)

Consequence. W(w) for the probe is equivalent to the "core" estimate 𝒩_core ≪ K N Z^ε. Here 𝒩_core is 𝒩 restricted to:
- powerful row part O < 2w;
- residual rows H > M′ − w/2;
- active slot length z_a > ℓ′ − w/2;
- dual cube part ℓ_b < w/3;
- dual length v > T_d − 2w/3.

Inactive slots (coefficient q_p^{−1}, absent from the conductor) cut the excess at rate 2. A dual cube part cuts it at
rate 3, which is why the cube-part idea (§6) has nothing to act on. The J ≠ ∅ summands have low exponent
−d + (w−3d)_+/2 < w/2.

## 3. Route (1): a Lemma-D split of the probe energy — fails at the target

Split 𝒩 = 𝒩_diag + 𝒩_off by tuples, P = P′ versus P ≠ P′:
  𝒩_diag = Σ_P |c_P|² E_P,  E_P := Σ_k |S_{k,P}|²,  S_{k,P} = Σ_n β(n)ψ_k(n)χ_P(n)^{−2}.

**3.1 The tuple-diagonal is NOT lossless with available tools [P].**
- (a) *Lemma D's mechanism gives exactly the wall.* Primal form: τ(PA′)χ_{PA′}(m) = τ(P)χ_P(m)·τ(A′)χ_{A′}(mP⁴)
  (twisted multiplicativity; |τ(P)χ_P(m)| ≤ 1 dropped by positivity, as in Lemma D step 2). So E_P is the
  *unmarked* theta energy over the rows mP⁴, i.e. the tuple becomes a j = 4 label. In lem:reflected-energy a j = 4
  label prime contributes 2A₀ − N₀ = +1 per log to T_d, exactly as an active mark does: +2z_a in T_d against
  +z_a in N_*. Hence T_d = 2M + ℓ − 1 is unchanged. In lem:canonical-moment the same conversion is free because the
  target there is the column length F = N′ + V, which grows with the label length V. Here the target is the row
  count Z^{M′}, which does not grow (eq. completed-moment, l.8098: max{M − O, 2M − O + Q + V − N_*}).
  **This is the exact failing step**, as the marked-cap agent anticipated.
- (b) *Every single-mark energy has the wall.* lem:hybrid-energy with one P gives E_P ≪ (K+N)N·Z^ε, so
  Σ_P|c_P|²E_P ≪ (K + N)N, the full column term. A lossless 𝒩_diag needs the theta-blind hybrid sieve
  (HLS) Σ_{k,P}|Σ_n a_n ψ_kχ_P^{−2}(n)|² ≪ (KL + N)^{1+ε}‖a‖² at N ∈ (K, K Z^w], or (H1) for each P.
  Both are open. The best proved form (hls agent) is min{L(K+N), K(L+N+(LN)^{2/3}), N + KL² + KL√N}, which here
  is L(K+N): the loss is again N/K = Z^w.
- (c) *Möbius-dual check.* §6's identity writes the off-diagonal of the marked energy as
  Σ_H|Π_H|²|M_H|² with Π_H = Σ_P c_P(H/P)₆. One might hope that its P-diagonal is the unmarked Möbius family,
  which has a power-saving asymptotic (edge-theta Prop. U). It is not. For P = P′ the two copies share P, so
  χ_{A₁}χ̄_{A₂} is principal at P. Poisson in m then runs modulo A₁′A₂′ (dual length Z²/Q, not (ZL)²/Q), and μ-absorption
  gives the χ_P⁴-twisted family Σ_{h ≍ Z²/Q}|Σ_n μ(n)χ_P(n)⁴(h/n)₆|². Its Poisson dual in h is E_P again.
  So there is no mark removal: the single-mark energy is self-dual, and its wall is real for every method in the move calculus.

**3.2 Cross-mark terms and "extra averaging over marks" [P].** For independent random signs ε_P,
E_ε[𝒩(εc)] = 𝒩_diag exactly. So averaging over marks removes 𝒩_off completely but leaves the wall in 𝒩_diag.
Even "W on average over mark signs" is equivalent to (H1)-average, and is open. 𝒩_off itself is (H2).
The h2 agent proved there is no (off-diagonal)^{1/2} inequality, since beyond the wall 𝒩_off = −𝒩_diag on a subspace of dimension ≥ N − K.
Route (1) therefore needs (HLS) or (H1), together with (H2), which is the 19/22 package. It cannot give a partial w.

## 4. Route (2): large-values dichotomy for B_m + Hölder on A_m — fails at the flat configuration

- **A is flat [R-/N].** By Corollary D, A = A^sig + A^rest with A^sig_m = Y^{−1/3}Ŵ₁(1/6+iv)ρ(m)
  and |ρ(m)| ≈ q_m^{−1/12} for typical m (04 numerics: ρ is unitary after the factor q^{1/12}).
  So |A_m| ≍ Y^{−1/2}P_a^{1/12} for most m, and Σ_S|A_m|² ≍ (|S|/Q)‖A‖² for every set S that is not
  adapted to ρ. Hölder with Σ|A|^{2k} gains exactly (|S|/Q)^{1−1/k}.
- **Dichotomy.** Put S = {m : |B_m| > Z^ε}. The rows off S cost ‖A‖√Q Z^ε, which is lossless. The rows in S need
  |S|·Σ_S|B_m|² ≤ Q², i.e. |S| ≤ Q·Z^{−w}.
- **Failing step [P].** The flat configuration b := βX ∈ span{ψ_k}, b = Σ_k d_k ψ̄_k with generic d, gives
  |B_m|² ≍ N/K = Z^w for *all* m. It is consistent with every input that sees β only through |β|: the hybrid sieve, pointwise
  convexity, and Halász–Montgomery with Pólya–Vinogradov off-diagonals. Halász–Montgomery gives
  Σ_S|B|² ≤ (N + |S|K)‖b‖², which beats N only for |S| < N/K = Z^w ≪ Q. Then |S| = Q and nothing is gained.
  - Excluding the flat configuration needs #{m : |B_m| > Z^{w/2−η}} ≤ Q Z^{−η′}, which by Chebyshev is ⇔ W(w − 2η).
  - Excluding |B_m| > Z^ε on all but Q Z^{−w} rows needs a lossless 2k-th moment of B. B^k has dual length N^k ≫ Q, which is hopeless.
  - Hölder in m with any exponent pair also gains nothing when both factors are flat: (Σ|A|^p)^{1/p}(Σ|B|^q)^{1/q} = Q^{−1/2}‖A‖‖B‖·Q^{1/2}.

## 5. Route (3): slot groups, slot powers, many tiny slots — the wall sees only q_P

- **Slot group moved to the rows [B].** Write X = X₁X₂ (lengths ℓ₁ + ℓ₂ = ℓ). Cauchy–Schwarz over P₂ (Σ|c_{P₂}|² ≍ 1)
  gives 𝒩 ≤ Σ_{k,P₂}|Σ_n βX₁ ψ_kχ_{P₂}^{−2}|². Even granting (HLS) for the (k, P₂) family this is ≤ (KL₂ + N)N.
  The loss is max(0, w − ℓ₂) + ℓ₂ ≥ w: break-even at best (O2 exchange rate). It is strictly worse with proved tools.
- **Duality form [P].** 𝒩 = sup_{‖d‖=1}|Σ_n b(n)D(n)X(n)|², with D = Σ_k d_kψ_k. Any Cauchy–Schwarz or Hölder in n is phase-blind
  (O1, energy-wall) and gives at best N·(Σ_n|D|⁴)^{1/2}(Σ|X|⁴)^{1/2} ≍ KN^{3/2}, which is far worse.
- **Slot powers / polarisation [P].** E_f Π_c Q^{(f,c)} = (k!/k^k)e_k replaces disjoint windows by k-subsets of one
  window, but the product of active primes still has norm Z^{z_a}. In (2.1) only z_a enters, never the number or
  shape of the slots.
- **Many tiny slots [B].** The J-sum has low exponent max_J [−d_J + (w − 3d_J)_+/2] = w/2, attained at J = ∅,
  independently of the slot mesh. Dropping J = ∅ destroys the compensation identity (the scalar prime terms on the
  high side return).
- **Inactive branches** reduce the excess at rate 2, but they are fixed by the inclusion–exclusion
  1_{q|A} = 1 − χ_q(A)⁰ and cannot be chosen.

## 6. Route (4): the cube part of cn³ — nothing at the binding dyad

- **Dual cube part [B].** By (2.1) a dual cube part of length ℓ_b cuts the excess by 3ℓ_b, so every dyad with
  ℓ_b ≥ w/3 is already inside the wall. At the binding dyad b = 1, and χ_P(b³)^{−2} = 1 is then empty.
  The refined case-A count (03) is (1+2ℓ′−2M′)/3 = −0.10. That branch is not active, and its improvement does not touch the dual branch.
- **Primal cube part [P].** Write B_m = Σ_n q_n^{−1}ᾱ(n)³(m/n)₂·B^{sf}_m(ZL/q_n³). The n-th piece is a squarefree-index
  theta of length ZL/n³ twisted by (·/n)₂. Its reflection conductor is (QLn)², so its dual length is Q²Ln⁵/Z, worse
  for every n > 1. Splitting the cube part off the completed sum breaks the completion that makes the reflection
  exact, and every piece has a worse wall. Marks on the cube part give no amplification (09 §9).
- **Cube part as an averaging variable [B].** It is decoupled from the marks (energy-wall O2), so the excess is w/2, a slope of 1/2. That is useless (§1).

## 7. New: the loss is rank one (theta-direction reduction)

### 7.1 Statement [C on Corollary D of 04, R-]

Fix a rescaled subset J with (d, M′, ℓ′, Q, Y′, X′, P_a) as in prop:probe-low, and put w_J = (M′+ℓ′−1)_+.
For a ray class σ let c_σ(s) = 1_{s∈σ}χ_s(b*)/(τξ(s)), D_σ(w;m) = Σ_s c_σ(s)g̃(s,−m)q_s^{−w}, and
ρ_σ(m) = Res_{w=2/3}D_σ(w;m). (04 Theorem A gives continuation and a single pole at w = 2/3 in Re w > 1/2.) Put

  φ(r) := Ω(r)·(2π)^{−1}∫ Ŵ₀(iv)Ŵ₁(1/6+iv) r^{−iv} dv ∈ C_c^∞(0,∞),
  T^J_σ := Σ_{m≠0} φ(q_m/Q) ξ(m) ρ_σ(m) B^J_{m,σ}(Z).

**Proposition T.**

  |I_J| ≪ Z^ε [ X′^{1/2} Z^{w_J/2} + Q^{−1/2}Y′^{−1/3} Σ_σ |T^J_σ| ].

In particular:
- (i) Cauchy–Schwarz on T together with Σ_mφ|ρ_σ|² ≍ Q^{5/6} gives back P78's X′^{1/2}P_a^{1/12}Z^{w_J/2}, now without
  P78's Gram-domination hypothesis P_a²/Y′ ≪ P_a^{1/6} (i.e. l_y − d − 11b/6 ≥ 0).
- (ii) Under **W_θ(w): |T^J_σ| ≪ Q^{11/12}Z^ε whenever M′+ℓ′ ≤ 1 + w**, the bound is

  |I_J| ≪ X′^{1/2}(Z^{w_J/2} + P_a^{1/12})Z^ε.

  So the low side is lossless as long as w ≤ b/6.

*Proof.*
1. Corollary D (04, `wf/FE_LoA6_proof.md` §5) splits A_{m,σ,v}(Y′) = Y′^{−1/3}Ŵ₁(1/6+iv)ρ_σ(m) + A^rest_{m,σ,v}, with
   Σ_mΩ|A^rest|² ≪ (1+|v|)^K(Q/Y′)Z^ε.
2. Insert this into eq. low-separated.
3. *Rest part.* Cauchy–Schwarz in m, Lemma probe-row-norm's energy bound Z^{M′+w_J+ε} (§2), and integrability of
   (1+|v|)^{K/2}|Ŵ₀(iv)| give Q^{−1/2}(Q/Y′)^{1/2}Z^{(M′+w_J)/2+ε} = X′^{1/2}Z^{w_J/2+ε}.
4. *Theta part.* The m-sum is finite, Ŵ₀(i·) is Schwartz and Ŵ₁(1/6+i·) is bounded, so the v-integral can be taken inside.
   B^J does not depend on v, and the integral produces φ. ∎

Remarks.
- (a) The theta part depends on (Q, Z, L) only, and Y′ enters it only through the scalar Y′^{−1/3}. In the
  B_low formula the theta part alone gives 5/6 + M/12 − ℓ/6 + w_θ/2, with no b and no l_x.
- (b) Since the rest part has spare room P_a^{1/12} = Z^{b/12} at the wall, this is the first place where the Gram
  excess and the wall excess can be traded.
- (c) [H] The zero-detecting signal sits in the theta part: 04 found that the P_a^{1/6} Gram term is the signal's
shadow, and vls/C that T carries the full signal. Consequently, for b > 0, P78's whole low side is in effect the single
inequality |⟨ρ,B⟩| ≤ ‖ρ‖‖B‖. The A^rest part has spare room Z^{b/12}, which is wasted unless the wall is relaxed.
For b < 0 the A^rest part dominates: it carries no signal and wastes room, which is why b ≥ 0 is optimal.

### 7.2 W_θ is strictly weaker than W, and what it is

- W(w) ⇒ W_θ(w) by Cauchy–Schwarz.
- In reflected form, T = Σ_{μ~N}b(μ)·Σ_m π(m)(μ/m)₂·(kernel) with π = ρζ (Lemma Z of 04), which is the 04 theta pairing in the
  marked geometry.
- W_θ ⟺ |⟨Ψ*π, Π_V b⟩| ≲ √K‖π‖‖b‖, where Π_V is the projection onto the K-dimensional row space V of Ψ. This is a
  **rank-one** version of W: it asks that Π_V b have no excess along the single vector Ψ*π.
- **Random model:** |T| ≈ Q^{5/12}, against the needed Q^{11/12}.

### 7.3 Payoff [M]

The low side becomes B_low = 5/6 + M/12 − ℓ/6 + max(w − w_θ, w − b/6, 0)/2, with no Gram constraint.
With a full W(w_W) added for the rest part, the excess is max(w − w_θ, w − w_W − b/6, 0) (mode `thetaW`).
(`runall2.py thetaInf thW01 …`; table in §8.)

### 7.4 Failing step for W_θ [P/B]

- (a) Cauchy–Schwarz over μ: Σ_{μ~N}|Σ_m π(m)(μ/m)₂|² ≍ N‖π‖² is the *truth* (rows N > Q), so this loses N/Q = Z^w
  even in the random model.
- (b) Cauchy–Schwarz over m is W itself.
- (c) The four functional-equation images R0–R3 of the 3-node model (04 vls/A, an exhaustive [RC] search) are, past
  the wall:
  - R0 (ZL, Q): sextic sieve, loss ZL/Q;
  - R1 (QZ^w, Q): quadratic sieve, loss Z^{w/2};
  - R2 (coupled, trivial bound);
  - R3: contains the Möbius signal.
  None reaches Cauchy–Schwarz quality at M + ℓ = 1 + w.
- (d) Cauchy–Schwarz over the marks, Σ_P|⟨ρ,B_P⟩|², needs ‖[B_P]_P‖²_op ≪ Q, i.e. single-mark W for every P plus (H2). That is worse.
- (e) [H] T contains the zero-detecting signal (Remark c). An unconditional proof of W_θ at a geometry with
  B_low < 0.8749 is therefore itself a zero-density statement for the family M_h (the copies h = a⁶), and pointwise
  zero-free inputs (β*) give only Z^{2β*−1} there (edge-theta O2).

## 8. Numbers [M]

Target B = 41/47, weighted high side (E_old/E_wt, paper counts, κ = 3/4). Each row minimises the worst high exponent
over ℓ, b ∈ [0, 0.45] and M ≤ M_max(ℓ, b), where M_max is the largest M with B_low ≤ 41/47.
Coarse search, then a fine grid: δ ∈ [0.05, 0.75] in 0.01 steps (0.005 in the first three rows), t with 101 points, r with 201 points.
Results are in `runall.out`.

| low-side input | worst exponent | at (ℓ, b, M, w = M+ℓ−1) | binding (δ, x, R) |
|---|---|---|---|
| paper | +0.01139 | (0.1775, 0.065, 0.8226, 0.000) | (0.495, 0.5, 0.566) |
| lossless W(0.040) | +0.00060 | (0.191, 0.080, 0.8492, 0.040) | (0.46, 0.5, 0.600) |
| lossless W(0.043) | −0.00036 | (0.1915, 0.070, 0.8511, 0.043) | (0.47, 0.5, 0.590) |
| lossless W(0.045) | −0.00083 | (0.1925, 0.090, 0.8526, 0.045) | (0.45, 0.5, 0.609) |
| Corollary D + W_θ(∞) (rank-one wall only) | **+0.00471** | (0.190, 0.225, 0.8476, 0.038) | (0.25, 0.45, 0.794) |
| Cor. D + W_θ(∞) + full W(0.01) | +0.00131 | (0.1925, 0.215, 0.8531, 0.046) | (0.29, 0.5, 0.756) |
| Cor. D + W_θ(∞) + full W(0.02) | **−0.00204** | (0.196, 0.215, 0.8599, 0.056) | (0.29, 0.5, 0.756) |
| Cor. D + W_θ(∞) + full W(0.025) | −0.00369 | (0.1975, 0.215, 0.8631, 0.061) | (0.29, 0.5, 0.756) |
| Cor. D + W_θ(∞) + full W(0.03) | −0.00535 | (0.200, 0.225, 0.8676, 0.068) | (0.26, 0.5, 0.782) |
| slope f = 0.5 / 0.25 (excess f·w past the wall) | +0.01133 / +0.01128 | optimum stays on the wall (w ≈ 0) | (0.47, 0.5, 0.59) |
| slope f = 0.15 | +0.00910 | (0.243, 0.055, 0.8607, 0.104), far vertex | (0.47, 0.5, 0.59) |
| slope f = 0.10 | −0.00164 | (0.243, 0.070, 0.8802, 0.123), far vertex | (0.47, 0.5, 0.59) |

Reading.
- A lossless W closes 41/47 at w_W ≈ 0.0415 (between 0.040 and 0.043), the shared wall alone, without the lem:marked transfer. This agrees with marked-cap's
  0.041.
- The rank-one wall alone is worth about W(0.03). With W_θ the full-W requirement drops from ≈ 0.042 to ≈ 0.014.
  The large b (≈ 0.22) needed to make the A^rest part lossless moves the binding bin to small δ ≈ 0.25–0.3, where
  the high side gets worse. That is why W_θ alone falls short.
- W_θ, like W, is an open beyond-the-edge statement (§7.4).

## 9. Files

| file | content |
|---|---|
| `wm4.py` | weighted high side (vectorized gmodel) with a general low side: paper / W / slope / theta modes |
| `scanpar.py`, `runall.py`, `runall2.py` | parallel coarse-to-fine scans; all results in `runall.out` |
| `runW.py`, `runT.py` | first (serial) scans |
| `gmodel.py`, `wall.py`, `wmodel.py` | copies of the coordinator's models |
