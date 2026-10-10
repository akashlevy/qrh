# single-mark — the single-mark wall SM(w) (energy-wall program, step 2)

Agent `single-mark`, 10 Oct 2026. P78 = `openai/math/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex`;
line numbers refer to it. Labels as in `../BRIEF.md`: **[P]** proved (given P78's lemmas), **[B]** exponent bookkeeping,
**[C]** conditional, **[N]** numerical, **[H]** heuristic, **[L]** literature (from memory; arXiv was not consulted).
Every exponent quoted below is reproduced by `exponents.py` (output in `exponents.out`).

## 0. Verdict

1. **No proof of SM(w) for any w > 0.** No partial statement found here gives a lossless w > 0 for the coherent
   sum either. Exponent payoff of everything proved here: **0**.
2. **Exact location of the wall [P + B, §3].** SM(w) is equivalent [B] to extending the single-twist version of
   Proposition U, which I prove [P, §3.2]:

       F_P(H) := Σ_{h~H} | Σ_{n~Z} a(n) χ_P(n)² (h/n)₆ W(q_n/Z) |²  =  Diag_P(H) + O(Z^{2+ε} (1 + ZL/H)),

   from its natural range H ≥ ZL down to H = ZL·Z^{−w} (a = μḠ; H = Z²/K). The factor ZL/H = KL/Z = N/K comes from
   exactly one step: Goldmakher–Louvel's quadratic large sieve applied to the reflected dual of length K²L/Z > K
   (both pieces of the Ramanujan factor at P). All earlier steps are exact identities.
3. **The current bound is exactly (K+N)N [P, §2].** Its column term N² is the V-term of P78 eq. (quadratic-large-sieve)
   (l.2644) = the NG/C² term of lem:hybrid-energy (l.7598). Reflection plus the sextic sieve is lossless only for
   N ≥ K^{3/2}L (§2.2); the probe would need M ≥ 2 for that.
4. **Routes (a)–(e) each fail at an identified step** (§4–§8). New facts on the way:
   - (a) The needed short-range (H1) is SM itself on N′ ∈ [KZ^{w/3}, KZ^w] with allowance (N/N′)^{1/2} (kernel decay).
     No mean-value statement for D(s) on a vertical line captures it: (H1) at the central point is equivalent to
     closing the whole gap (K, K^{3/2}L), and its proved form loses K^{1/5}L^{2/5}. For every Mellin kernel at least as
     wide in height as P78's Gaussian, the hybrid (k,t) sieve gains exactly what the high-side signal loses: an SNR
     invariant (uncertainty principle). That is the T₀^{1/2}.
     Interpolation (Phragmén–Lindelöf) between the two lossless regimes fails by Paley–Wiener.
     Individual Lindelöf would give SM(∞), but subconvexity is useless unless δ > 0.2394 (Lindelöf: 0.25).
   - (b) **The mark is invisible to the Hecke structure:** the p³-shift eigenvalues of θ⊗ψ_kχ_P^{−2} and of θ⊗ψ_k agree
     for every p ≠ P. Hecke amplifiers for the k-family are quadratic-character polynomials: they lengthen the sum.
   - (c) Prop U breaks at the reflection step (above). There is a **j = 4 / j = 2 dichotomy**: the two cubic twists of the
     Möbius family have walls at H = ZL and H = ZL². W needs only the P-*average* of E_P, not SM for each P
     (a correction to the brief's wording). A **new duality (Poisson in P**, via β·conj γ₂ = ᾱ) turns that average into
     a cubic family of GL(1) Hecke L-functions, a finite-FE-group family. But it needs a power saving δ ≥ 0.44,
     i.e. nearly square-root, uniformly in a twist of conductor H′^{0.53}.
   - (d) The cubic-first-moment mechanism beats *theta terms*, i.e. the (QN)^{2/3} term of the cubic sieve. SM's
     obstruction is the column term of the quadratic sieve, which has none. The flat configuration satisfies every
     moment, BHM and near-orthogonality constraint.
   - (e) No Kloosterman sums arise under any available exact transformation. θ⊗χ_P^{−2} is residual spectrum, so
     positivity gives the trivial bound. The family spans K different levels.
5. **What a proof must supply (minimal form):** the projection bound ‖Π_V b_P‖² ≪ K Z^ε for the *single* theta vector
   b_P = (βχ_P^{−2}W(q_n/N))_n against the K-dimensional quadratic space V = span{ψ_k : k ~ K}. It is true numerically
   with margin ≈ K^{1/2} (energy-wall §4). It is the theta analogue of an open GL(2) problem: Lindelöf on average for
   K quadratic twists of one form of level L² (or L), at lengths between K and the AFE length.

## 1. The statement

### 1.1 Objects

Notation as in P78 and wall4 §1:
- Z is the probe scale; K = Z^M; L = q_P = Z^ℓ for one good prime P.
- β = γ₂ᾱ is the cubic-theta coefficient on squarefree n.
- ψ_k = χ_k³ is the quadratic row character.
- W is a fixed smooth bump.

    S_{k,P}(N) := Σ_{n sf} β(n) ψ_k(n) χ_P(n)^{−2} W(q_n/N),      E_P(K,N) := Σ_{k sf, q_k~K} |S_{k,P}(N)|².

The diagonal is Σ_k Σ_n |βW|² 1_{(n,k)=1} ≍ KN.

**SM(w).** E_P(K,N) ≪_ε K N Z^ε whenever N = K²L/Z and K < N ≤ K Z^w. The bound must be uniform in P with q_P = L,
for L ≥ Z^w (the brief's "L ≫ N/K").

**Parameter relations [B].** N/K = KL/Z = Z^w, hence:
- M = 1 + w − ℓ, N = Z^{1+2w−ℓ}, KL = Z^{1+w};
- the Möbius-dual family size is H = Z²/K = Z^{1+ℓ−w} = ZL·Z^{−w}.

At the weighted geometry (ℓ, M, w) = (0.19, 0.85, 0.04):

| quantity | value |
|---|---|
| K, N, L | Z^{0.85}, Z^{0.89}, Z^{0.19} |
| KL | Z^{1.04} |
| K√L | Z^{0.945} |
| H | Z^{1.15} (against ZL = Z^{1.19}) |

**Conductors [P, from Prop. completed-reflection, l.1676ff].** For Ψ = ψ_kχ_P^{−2} = χ_k³χ_P⁴:
- every p | k has j = 3, so B_p = χ_p^{−5} = χ_p;
- P has j = 4, so it is active, with the Ramanujan factor B_P(x) = q_P^{−1/2}(−1 + q_P 1_{P|x}) (l.1791);
- the kernel scale is q_c² ≍ (q_k q_P)².

The Ramanujan factor splits every dual into two pieces:

| piece | length | normalised coefficient |
|---|---|---|
| 1_{P|μ} | K²L/X | O(1), after twisted multiplicativity at P |
| −1 | (KL)²/X | L^{−1/2} |

So the effective analytic conductor is K²L (AFE point K√L), the "U_P picture" of edge-theta §3.2. The family has K
members against an AFE length K√L. In the probe the sums have length N ∈ (K, K√L): **beyond the edge by N/K.**

### 1.2 Equivalent forms [B unless stated]

- **(E1) Reflected.** SM(w) itself.
- **(E2) Primal (P78's marked primal).** Let 𝒫_P(K) := Σ_k |𝒯_P(ZL;k)|², where 𝒯_P is P78's completed sum with the
  whole-index mark 1_{P|nb³} (eq. marked-completion, l.7318). By Cor. marked-reflection (l.7361), 𝒯_P is the sum of:
  - an active branch: q_P^{−1/2} × the E_P-sum;
  - an inactive branch: q_P^{−1} × an unmarked quadratic-row sum of length K²/(ZL) < K, which is lossless by the
    quadratic sieve.

  Hence E_P ≪ KN Z^ε ⟺ 𝒫_P ≪ KZ·Z^ε. Twisted multiplicativity τ(PA′) = τ(P)τ(A′)χ_{A′}(P)⁴ (wall4 §3.1a) makes 𝒫_P
  the unmarked theta energy over the thin row set {kP⁴} (j = 4 at P), up to cube-part terms.
- **(E3) Möbius dual.** Apply Poisson in k to 𝒫_P. The modulus is P·A₁′A₂′; at P the character is principal, and its
  Fourier transform is the Ramanujan sum c_P(h) = −1 + L·1_{P|h}. Then apply μ-absorption τγ₁ = μG (P78 eq. signal,
  l.965). This gives

      off-diag(𝒫_P) = (K/Z)·[F_P(H) − Diag_P] − (K/(LZ))·[F^{(1)}_P(LH) − Diag^{(1)}],   H = Z²/K.

  Here F^{(1)}_P is the family twisted by the sextic χ_P^{±1} at range LH. The proof of Prop U_P (§3.2) gives
  F^{(1)}_P − Diag ≪ Z^ε(Z² + HN), which is ≤ LZ² because K ≤ Z, so the companion family is harmless. Therefore

  **SM(w) ⟺ F_P(H) = Diag_P(H) + O(Z^{2+ε}) at H = Z²/K ∈ [ZL·Z^{−w}, ZL).**

  *Convention warning [P].* With a = μḠ and rows (h/n)₆ (edge-theta's normalisation), the twist that is Poisson-dual
  to the j = 4 theta is χ_P(n)². wall4 writes χ_P(n)⁴ = conj χ_P(n)² with conjugated rows. The *other* cubic twist
  (χ_P⁴ in edge-theta's normalisation) is dual to the j = 2 theta. It has a different wall (§3.4), so a numerics run
  implementing "the χ_P⁴ family" in edge-theta's code would test the wrong family.
- **(E4) Short (H1).** See §4.2.
- **(E5) Projection form [P, one direction].** E_P = ‖Ψ b_P‖² = ‖Ψ Π_V b_P‖² ≤ ‖Ψ‖²_op ‖Π_V b_P‖² ≤ (K+N) Z^ε ‖Π_V b_P‖².
  Hence ‖Π_V b_P‖² ≪ K Z^ε ⟹ SM. The converse needs the lower Riesz bound σ_min(Ψ|_V)² ≫ N Z^{−ε}, which is true
  numerically (hls, h2) and not needed here.
- **(E6) Large values [B].** SM ⟺ #{k : |S_{k,P}|² ≥ N Z^{2η}} ≪ K Z^{−2η+ε} for all η > 0 (dyadic Chebyshev, up to ε).

## 2. The current bound and where N² comes from

### 2.1 Proposition 1 [P]

For all K, N ≥ 1, every good prime P and all coefficients |b(n)| ≤ 1 on squarefree n ~ N,

    Σ_{k sf, q_k ≤ K} |Σ_n b(n) ψ_k(n)|²  ≪_ε  (KN)^ε (K + N) Σ_n |b(n)|².

In particular E_P(K,N) ≪ (K+N) N Z^ε, i.e. **excess max(0, w) over the diagonal.**

*Proof.* This is P78 eq. (quadratic-large-sieve), l.2644 (Goldmakher–Louvel), with U = K, V = N and
c_n = b(n)χ_P(n)^{−2}W, which is independent of k. P78 obtains the same bound, with the coherent mark sum inside the
square, as lem:hybrid-energy (l.7510), where B = 1 and {N + L + (NL)^{2/3}} ≍ N. ∎

**Which step produces N² [P].** It is the V-term of the quadratic sieve, the column term: ‖Ψ‖²_op ≥ max_k ‖ψ_k‖² = N.
In P78's proof of lem:hybrid-energy it is the term NG/C² of eq. (hybrid-after-quadratic-reduction), l.7598, at G = C = 1.
The term is attained: b = conj ψ_{k₀}W gives |S_{k₀}|² ≍ N². The *flat* configuration b ∈ span{ψ_k} gives
|S_k|² ≍ N²/K for every k (wall4 §4). So no argument that treats b = βχ_P^{−2} as a general bounded sequence can beat
it; by energy-wall O1, neither can any argument that sees only |β|.

### 2.2 The full provable set [P]

**Lemma 2.** With the available tools (quadratic sieve, P78's sextic sieve lem:sextic-large-sieve
K + D + (KD)^{2/3} (l.4706), completed reflection), E_P(K,N) ≪ KN Z^ε holds for N ≤ K and for N ≥ K^{3/2}L.
In the gap K < N < K^{3/2}L the best bound is min{(K+N)N, reflection branch}, which equals (K+N)N throughout the SM range.

*Proof.* Reflect E_P using §1.1 (j = 3 at k gives sextic dual rows χ_k; j = 4 at P gives the Ramanujan factor):

    N^{−1/2} S_{k,P}(N) = ζ [ N₂^{−1/2} S₂(N₂) − L^{−1/2} N₁^{−1/2} S₁(N₁) ],    N₂ = K²L/N,  N₁ = L N₂.

Both S_i have bounded coefficients and sextic rows. Apply the sextic sieve to each (fixed sectors absorb the
reciprocity sign 𝔯; cube parts as in lem:reflected-energy):

    E_P ≲ N [ (K + N₂ + (KN₂)^{2/3}) + L^{−1}(K + N₁ + (KN₁)^{2/3}) ].

This is lossless iff N₂ ≤ K^{1/2}, i.e. N ≥ K^{3/2}L.

In the SM range N₂ = Z > K. The excess of this branch is (2−M)/3 = 0.383, against w = 0.04 for the quadratic sieve.
The maximal excess over the whole gap is K^{1/5}L^{2/5} (= Z^{0.246}), attained at N = K^{6/5}L^{2/5}. ∎

*Remark [P].* With L = 1 the same lemma says the **unmarked** E₀(K,N) is lossless only for N ≤ K or N ≥ K^{3/2}. In the
unmarked probe N = K²/Z ≤ K always (M ≤ 1). **The mark matters only in that it pushes the probe's N = K²L/Z past K.**

### 2.3 What the probe actually needs at each dual dyad [B]

At the core (wall4 §2), lem:reflected-energy (l.7783ff, branch u = z_a, max = v) gives exponent v − ½(T_d − v).
So a dual dyad of length N′ = N Z^{−δ} has excess w − 3δ/2: the column excess w − δ, minus δ/2 from the kernel factor
m_ker² = Z^{−D_ker/2} (l.7884). The probe therefore needs exactly

    (SM*)_w :  E_P(K, N′) ≪ K N′ (N/N′)^{1/2} Z^ε    for all N′ ∈ [N Z^{−2w/3}, N],  N = K Z^w.

At N′ = N this is SM(w). For N′ ≤ N Z^{−2w/3} Proposition 1 already suffices.

The other core restrictions of wall4 (2.1) carry over, since every one of them cuts the excess:
- row dyads q_k ~ KZ^{−Δ} have excess w − 2Δ, so only Δ < w/2 matters;
- the powerful row part must satisfy O < 2w;
- the dual cube part must satisfy ℓ_b < w/3, since it cuts the excess at rate 3.

So the exact dyadic box for the single-mark problem is
q_k ∈ [KZ^{−w/2}, K], O < 2w, ℓ_b < w/3, N′ ∈ [NZ^{−2w/3}, N], with allowance (N/N′)^{1/2}.
The worst corner is the full dyads q_k ~ K, N′ ~ N with O = ℓ_b = 0, where the needed statement is exactly SM(w).

## 3. Proposition U and its single-mark extension

### 3.1 Anatomy of Proposition U (edge-theta §3.2, re-derived)

F(H) = Σ_h W |M_h|², with M_h = Σ_{n~Z} a(n)(h/n)₆ W and a = μḠ.

1. **Poisson in h, modulus n₁n₂ (n₁ ≠ n₂) [identity].**
   - The Gauss sums give g(χ_n;k) = conj χ_n(k)√q_n γ₁(n).
   - With aγ₁ = conj β this gives F − Diag = (H/Z)·[Θ(K) − diag], where K = Z²/H and
     Θ(K) = Σ_{q_k ≲ K} |Σ_{n~Z} β(n) χ_k(n)|² (sextic rows, j = 1).
2. **Completed reflection [identity].** j = 1 gives B_p = χ_p^{−3}, quadratic. So Θ(K) ≍ (Z/N₀)·E₀(K, N₀) with
   N₀ = K²/Z.
3. **Quadratic sieve.** E₀ ≪ (K + N₀)N₀. For H ≥ Z we have N₀ ≤ K, so Θ ≪ ZK and F − Diag ≪ Z².
4. **Range.** H ≤ Z^{3/2} (M ≥ ½) keeps P78's powerful-row term (1−2M)/4 out (completed-moment, l.8053).

### 3.2 Proposition U_P [P, same inputs and status as Prop U]

For Z ≤ H ≤ Z^{3/2}, every good prime P with q_P = L, and F_P as in §0:

    F_P(H) = Diag_P(H) + O_ε( Z^{2+ε} (1 + ZL/H) ).

*Proof.* We follow §3.1 step by step.

1. **Step 1 is unchanged.** The twist rides on the coefficient. By the convention check in §1.2 (the Poisson dual of
   ξ is conj ξ on the theta side), the dual energy is the j = 4 theta energy
   Θ_P(K) = Σ_k |Σ_{n~Z} β(n) χ_P(n)^{−2} χ_k(n)|².
2. **Step 2.** P is active, since j = 4 ≠ 0 forces activity (l.1777). With N := K²L/Z the Ramanujan factor gives:
   - a piece of length N with O(1) coefficient: the twisted multiplicativity of the dual Gauss-sum coefficients
     turns 1_{P|μ} into the column mark χ_P(μ′)^{−2}, exactly as in Cor. marked-reflection;
   - a piece of length NL with coefficient L^{−1/2};
   - both pieces have quadratic rows.
3. **Step 3.** Apply the quadratic sieve to each piece:

       Θ_P ≲ Z [ (K + N) + L^{−1}(K + NL) ] ≍ Z (K + N).

   Hence F_P − Diag_P ≲ (H/Z)·Z(K + N) = Z² + HN = Z² + Z³L/H. ∎

**Consequences.**
- (i) For H ≥ ZL, i.e. KL ≤ Z, i.e. M + ℓ ≤ 1, the twisted family has Prop U's error O(Z²): it is inside the wall.
- (ii) For Z√L < H < ZL it is still a power-saving asymptotic, with relative error (Z/H)(ZL/H).
- (iii) By (E3), SM(w) is exactly the statement that (ii) improves to (i) on [ZL·Z^{−w}, ZL).

At the weighted geometry the needed error is Z^{2.00} and the proved one is Z^{2.04}. The diagonal is Z^{2.15}.

### 3.3 The exact failing step [P]

The failing step is **Step 3: the quadratic sieve at reflected length N = K²L/Z > K.**
- Both Ramanujan pieces contribute the same column term ZN. For the −1 piece: L^{−1}·(NL)·Z = ZN.
- Steps 1–2 are exact identities.
- In P78's own chain this is the NG/C² term of lem:hybrid-energy.
- In the Möbius language, the relative error ZL/H = N/K is precisely the reflected family's beyond-edge ratio.

### 3.4 The j = 4 / j = 2 dichotomy [P]

Let F′_P(H) := Σ_h |Σ_n a(n) χ_P(n)⁴ (h/n)₆|², the other cubic twist.
- Its Poisson dual is the j = 2 theta Σ_k |Σ_n β(n)χ_P(n)²χ_k(n)|².
- For j = 2, B_P = χ_P^{−4} = χ_P² is a character: there is **no** Ramanujan splitting, and the dual length is
  (KL)²/Z = NL.
- The proof of §3.2 then gives F′_P(H) = Diag + O(Z^{2+ε}(1 + ZL²/H)). **Its wall is at H = ZL².**

Only the twist dual to the Ramanujan ("U_P", level-L) theta has its wall at ZL. P78's j = 0 primal mark becomes
j = 4 on the reflected column, and U_Pθ ∝ θ⊗χ_P⁴ (twisted multiplicativity) has level P, not P².

### 3.5 Why SM is "Prop U-type": the mark moves the dual into the diagonal regime [B]

Fix (K, N) with N > K.
- **Unmarked.** E₀(K,N) is dual (reflection, then Poisson in k) to the Möbius family at H₀ = K³/N² with length
  K²/N > H₀. That is *beyond* the Möbius edge: a GRH or zero-density type statement (edge-theta O2).
- **Marked.** The dual is F_P at H = K³L²/N² with length Z = K²L/N and H/Z = KL/N > 1: the *diagonal* regime. What is
  needed is an asymptotic with relative error N/(KL).

So the mark trades a beyond-the-edge problem for an asymptotic in the diagonal regime, the class where Prop U lives.
This is why the brief calls SM the least hopeless lever. The required precision, however, is exactly what the theta
side cannot reach (§3.3).

## 4. Route (a): Perron, short (H1), hybrid (k,t) sieve, Mellin kernels

### 4.1 Pointwise inputs [P/B]

- **(i) Individual Lindelöf would suffice [P].** Assume D(½+it; ψ_kχ_P^{−2}) ≪ (KL)^ε for |t| ≤ Z^ε.
  - Mellin: S_{k,P}(N) = (1/2πi)∫_{(1/2)} D(s) N^s Ŵ(s) ds, with D = T/L₀; 1/L₀ is harmless on Re s = ½ + 1/log Z.
  - This gives |S|² ≪ N Z^ε, hence E_P ≪ KN Z^ε for all N: SM(∞) for that P.
  - For a single mark, individual bounds lose nothing. This differs from the coherent sum, where they lose L
    (wall-spectral Route 3).
- **(ii) Subconvexity is useless [B].** The AFE with Ramanujan splitting gives the convexity bound |D(½)| ≪ C^{1/4+ε}
  with C = K²L. A bound C^{1/4−δ} gives E_P ≪ KN·C^{1/2−2δ}. This beats (K+N)N only if (½ − 2δ)(2M+ℓ) < w, i.e.
  **δ > ¼ − w/(2(2M+ℓ)) = 0.2394** at the weighted geometry. Any amplification or Kuznetsov argument that outputs
  individual bounds must land within 0.011 of Lindelöf.

### 4.2 The short-range (H1) [P/B]

- **(H1).** Σ_k |D(½+it)|² ≪ K^{1+ε}L^ε.
  - By Perron, (H1) implies the diagonal bound for E_P at *every* length.
  - Conversely, the AFE splits D(½) into quadratic-row pieces of length ≤ X and sextic-row dual pieces of length
    ≲ K²L/X (plus the L^{−1/2}·(KL)²/X piece). Reflection identifies sextic-row energies at length N″ with E_P at
    length K²L/N″.
  - So (H1) is equivalent [B] to closing the **whole** gap of Lemma 2: E_P ≪ KN′Z^ε for all N′ ∈ (K, K^{3/2}L).
  - **Proved [B]:** Σ_k |D(½+it)|² ≪ K^{6/5}L^{2/5}·Z^ε, from the optimal split X = K^{6/5}L^{2/5} with the quadratic
    sieve below X and the sextic sieve K + D + (KD)^{2/3} above. The excess K^{1/5}L^{2/5} = Z^{0.246} is the gap
    maximum of Lemma 2.
- **What the probe needs** is only (SM*)_w of §2.3: lengths in [KZ^{w/3}, KZ^w], with allowance (N/N′)^{1/2}.
- **No L-function form [P].** No mean value of D on a single vertical line captures (SM*)_w:
  - From Re s = σ > ½, Perron loses N^{2σ−1}, because the diagonal Σ_k |D(σ+it)|² ≍ K/(2σ−1).
  - On Re s = ½ the mean square is (H1), which involves every length in the gap (K, K^{3/2}L).
  - So the short-range (H1) has no L-function formulation: it is a mean value of Dirichlet polynomials of length
    ≤ KZ^w, i.e. SM itself.

### 4.3 Hybrid (k,t) sieve and the Mellin kernel: an SNR invariant [B]

Grant even the optimal hybrid inequality

    Σ_k ∫_{|t|≤T} |Σ_n c_n ψ_k(n) q_n^{−it}|² dt ≪ (KT + N) Z^ε ‖c‖².

(Gallagher's lemma reduces it to the quadratic sieve on intervals of relative length 1/T.)

- **(i) Averages over heights only.** It is lossless on average over heights iff T ≥ N/K = Z^w. At one height,
  Gallagher's |f(0)|² ≤ ½∫_{−1}^{1}(|f|² + |ff′|) returns only (K+N)N.
- **(ii) SNR invariant.** In P78 the zero-detecting signal is the Mellin integral f_η(Z) of
  Z^{C(s)}e^{(s−5/6)²}H_η/L (eq. common-residue, l.~416). It is attached to *one* target and *one* O(1) Mellin window
  (the Gaussian). The high-side contract eq. (high-probe-contract) and the low-side contract (low-probe-contract) must
  hold per target.

  Suppose the kernel is changed to Ŵ, with a = ‖Ŵ‖_∞, b = ‖Ŵ‖₁, T_eff = b/a, and support radius T ≥ T_eff/2. Then:
  - the signal is ≤ S·a;
  - the noise² certified by Cauchy–Schwarz is ‖A‖²·b·N·min{a(KT+N), b(K+N)}, using the hybrid sieve or the pointwise
    sieve, whichever is better.

  Hence the certified ratio is

      SNR² ≤ (S²/‖A‖²) · max{ 1/(T_eff(KT_eff/2+N)N),  1/(T_eff²(K+N)N) } ≤ (S²/‖A‖²)·O(1)/((K+N)N)

  for every kernel with T_eff ≥ 1. Those are the only kernels on which the hybrid sieve can act. Narrower kernels
  (T_eff < 1, i.e. Z-independent Gaussian widths λ > 1) change only Z-independent constants, through the factors
  e^{λ²(·)} of the shifted Gaussian. So the best ratio is attained at T_eff ≍ 1, which is P78's Gaussian.

  **This is where the T₀^{1/2} sits.** The signal's amplitude is ≤ ‖Ŵ‖₁/T_eff (uncertainty), and the hybrid sieve
  returns exactly T_eff^{1/2} in amplitude. A gain would
  require the hypothetical zeros of each bad target to fill a positive proportion of the T₀ heights. The adversary
  places one zero per target, consistent with all density inputs; the Part II counts grow only like T₀^R,
  R ≈ 0.66 < 1 (09 NOTES §9, §11).
- **(iii) Rough kernels also lengthen the reflected sums.** V^♯ obeys eq. (reflection-kernel) with the seminorm
  𝔐_{A,⌈4A⌉+j+2}(V), and the gamma ratio is R(t) ≍ |t|^{4σ} (l.2414). So a primal kernel of Mellin height T₀
  multiplies the effective dual length by T₀^{O(1)}: the wall moves the wrong way.

### 4.4 Interpolation between the lossless regimes fails (outside the move calculus) [B]

Lemma 2 gives two lossless regimes, N ≤ K and N ≥ K^{3/2}L.
- For a weight with an entire log-extension, z ↦ (S_{k,P}(e^z))_k is an entire ℓ²-valued function. So log‖S(e^z)‖ is
  subharmonic, and Hadamard three-lines would interpolate log E_P linearly across the strip of width
  Δ = log(K^{1/2}L) (= 0.615 log Z).
- **The trade-off is fatal.**
  - The endpoint bounds hold only if the weight samples the gap with log-tails |φ(x)| ≲ e^{−|x|} out to |x| ≍ Δ: the
    dyad at length Ke^x carries excess up to e^x.
  - Such decay forces exponential type σ ≫ 1. By Poisson–Jensen in the upper half-plane,
    log|φ(iy)| ≤ σ|y| + (y/π)∫log|φ(t)|dt/(t²+y²). With |φ(i/2σ)| ≥ ½ (Bernstein) this gives
    a ≤ πσ(½ + log 2)/log(2σΔ) for decay rate a, so σ ≫ log Δ is needed.
  - The boundary growth e^{σ|v|} costs, after harmonic interpolation, exp(c σ Δ) in energy with
    c = 8G/π² = 0.742, where G is Catalan's constant (the harmonic measure of the strip at its midline).
  - This gives (K^{1/2}L)^{0.742σ} ≥ Z^{0.457σ}, which exceeds the maximal gap loss K^{1/5}L^{2/5} = Z^{0.246}.
  - Conversely, small σ makes the weight sample the gap with only polynomial-in-log tails.
- This is the same uncertainty principle as in 4.3(ii).

## 5. Route (b): Hecke structure [P unless stated]

- **5.1 Periodicity is exactly the completion.**
  - P78's completed coefficient at nb³ is β(n)Ψ(n)·ᾱ(b)³Ψ(b)³√q_b (normalised by √q). The b-factor is independent of
    n mod b and of shared primes (eq. reflection, l.1728–1736: Kazhdan–Patterson periodicity in P78's form).
  - Hence T(s,Ψ) = L₀(s,Ψ)·D(s,Ψ), where L₀ = L(3s−½, ᾱ³Ψ³) is a GL(1) Hecke L-function.
  - The p³-shift acts on θ⊗Ψ by the scalar ᾱ(p)³Ψ(p)³√q_p. [L/H] The Kazhdan–Patterson T_{p³} eigenvalue is a fixed
    function of this.
- **5.2 The mark is invisible to the Hecke structure.**
  - For Ψ = ψ_kχ_P^{−2}, Ψ³ = ψ_k·1_{(·,P)=1}. So for every p ≠ P the p³-eigenvalue of θ⊗ψ_kχ_P^{−2} equals that of
    θ⊗ψ_k, and L₀ is the same apart from the Euler factor at P.
  - The unmarked family at the probe's length N₀ = K²/Z ≤ K is inside the wall. Any argument that uses only Hecke data
    (amplifiers, the relations between β(n) and β(np³), periodicity) treats the marked and unmarked families alike.
  - **The obstruction lives entirely in the non-Hecke squarefree part χ_P(n)^{−2}.** [P for the eigenvalue identity;
    H for the moral.]
- **5.3 Hecke amplification lengthens the sum.** Since λ_{p³}(k) = λ_p·ψ_k(p), an amplifier is
  A(k) = |Σ_{p≤Λ} x_p ψ_k(p)λ_p|². Then

      Σ_k A(k)|S_k(N)|² = Σ_k |Σ_m (Σ_{pn=m} x_pλ_p b(n)) ψ_k(m)|²,

  a quadratic-family energy at length NΛ, beyond the edge by NΛ/K. For mean values this is strictly worse. For
  individual twists it gives at best subconvexity, useless by 4.1(ii).
- **5.4 No Hecke relation at P.**
  - The only Hecke-active variable is the cube part. It is decoupled from the mark (χ_P(b³)^{−2} = 1), cuts the excess
    at rate 3, and is absent at the binding dyad (wall4 §2, §6).
  - At P itself: θ⊗χ_P⁴ ≈ U_Pθ/τ(P), but U_P is not a Hecke operator of the cubic cover, and U_P(θ⊗χ_P⁴) = 0.
  - The Fricke/Atkin–Lehner involution at P *is* the reflection, i.e. the Ramanujan factor. So the Hecke algebra at P
    yields no new relation. [P/H]

## 6. Route (c): Prop U's proof, weaker statements, the coherent sum

### 6.1 Which step of Prop U breaks

The failing step is identified in §3.3: the reflected dual length K²L/Z exceeds the row count K, and the quadratic sieve
pays the column term. It is not Poisson, not μ-absorption, and not the gcd or powerful-row bookkeeping.

### 6.2 Weaker statements

- **(i) "Error Z²L^a".** Prop U_P gives exactly F_P − Diag ≪ Z^{2+ε}·(ZL/H).
  - At the probe's H = Z²/K this is Z²·L^{w/ℓ}, with a = w/ℓ = 0.21 at the weighted geometry. That is the current
    excess Z^w, nothing more.
  - A uniform a < 1 is not available: at H = Z the error is Z²L (a = 1).
  - In E_P language a uniform a would be an *offset* excess aℓ. It is worthless for SM alone. For the coherent sum it
    helps only together with (H2), and only if a ≲ 0.07 (energy-wall §2 offset/slope table). [B]
- **(ii) "P in a sparse set".** Every available estimate depends on P only through q_P: the quadratic and sextic
  sieves, the completed reflection (P's sector data are fixed phases), Prop U_P, and Lemma 2.
  - The only P-specific structures found permute the family and do not shorten it:
    - the conjugation symmetry S_{k̄,P̄} = conj S_{k,P} (h2 §2.1d);
    - inert P, for which χ_P is trivial on Z.
  - Nothing better than the uniform bound is provable on any sparse set of P. [B]
- **(iii) "On average over P, incoherently".**

  Define **AVG-SM**: L^{−1}Σ_{P~L} E_P ≪ KN Z^ε.
  - It is the hybrid large sieve (HLS) for the single vector β at N ∈ (K, KZ^w]. It is open. The best proved bound
    (hls Prop 1) has excess N/K again. [P/L]
  - **New duality: Poisson in P [P modulo squarefree and gcd bookkeeping].**
    - χ_P(n)^{−2} = conj (P/n)₃ (P78 l.1906 and cubic reciprocity), so Poisson in squarefree P ~ L (positivity first
      extends primes to squarefree P) produces the Gauss sum conj γ₂(n).
    - **β·conj γ₂ = ᾱ exactly**, since γ₂³ = μα (P78 eq. signal) and |γ₂| = 1.
    - Hence, for each k,

          Σ_{P~L} |S_{k,P}|² = L N·c + (L/N)·[ 𝒞_k(H′) − Diag′ ],
          𝒞_k(H′) = Σ_{h~H′} |Σ_{n~N} μ²(n) ᾱ(n) ψ_k(n) (h/n)₃ W|²,   H′ = N²/L.

      This is a **cubic family of GL(1) Hecke L-functions**, of the finite-functional-equation-group type where
      power savings exist (David–de Faveri–Dunn–Stucky; wall-spectral §1.2).
    - AVG-SM needs 𝒞_k − Diag′ ≪ N² = H′·L, i.e. a power saving
      **δ ≥ ½ − ℓ/(2·log_Z H′) = 0.44 in H′-units**. That is nearly square-root, against a diagonal H′N.
    - It must also be uniform in the twist ξ_k = ᾱψ_k of conductor K = H′^{0.53}, and the cubic sieve's own theta term
      (H′N)^{2/3} exceeds even the diagonal by L^{1/3}.
    - Available MDS savings are small and for fixed twists. **Dead, but it is the one reformulation with Euler
      products on both sides.** [B/L]
  - **Necessity (a correction to the brief).** wall4 §3.2 gives E_ε𝒩(εc) = Σ_P|c_P|²E_P. So W(w), required for all sign
    patterns of the marks, implies only **AVG-SM**, not SM for each P: a single P may carry E_P up to L·KN without
    violating W. Conversely, AVG-SM ⟸ SM for every P outside an exceptional set 𝔈 with |𝔈| ≤ L Z^{−w}, using E_P ≪ N²
    on 𝔈. [B]

### 6.3 The coherent sum [B]

𝒩 = Σ_P|c_P|²E_P + 𝒩_off.
- Any SM-type statement, full or partial, for each P or on average, controls only the first term.
- With 𝒩_off bounded by Cauchy–Schwarz, 𝒩 ≤ (Σ_P|c_P| E_P^{1/2})² ≲ L·KN. That is worse than (K+N)N whenever L > N/K.
- **No SM-type statement gives a lossless w > 0 for the coherent sum without (H2).** With (H2) at saving L^{−θ₂} the
  payoff is the h2 table: 19/22 at full strength, nothing for θ₂ ≤ 0.895.

## 7. Route (d): the cubic first-moment technique [P/B]

**7.1 Mechanism** (`An-unconditional-first-moment…/build/sections/dual.tex`, Prop. exceptional character moments):
1. The coefficients are bounded convolutions of prime sequences.
2. Λ = (2m − m*m*1)*log splits them into short Möbius pieces and full smooth sums.
3. A full-length factor is shortened by the GL(1) functional equation.
4. Moments of products with multiplicities k_i force each factor of length Y^{a_i} to have the theta profile
   |·| = Y^{5a_i/6} on Y^{2/3} rows. That profile is the extremal of the (QN)^{2/3} term of Heath-Brown's cubic sieve.
5. Near-orthogonality Σ_{p′≠p}|T_{p,p′}|² ≪ Z(P + (P³/Z)^{2/3}), plus Bombieri–Halász–Montgomery, excludes theta
   profiles at product lengths 4/3 < z < 2.

The output is a tiny saving η below the theta term.

**7.2 Why nothing transfers.**
- **(i) No theta term to beat.** SM's wall is the column term V of the *quadratic* sieve (K+N), which has no theta term.
  Its extremal is the flat configuration |S_k|² = N²/K for all k.
  - It satisfies every power-moment constraint Σ_k|S_k|^{2j} ≤ (K + N^j)N^j, since K(N²/K)^j ≤ N^{2j}.
  - BHM with the quadratic Gram gives N + max_k Σ_{k′}|G_{kk′}| ≥ N: the diagonal ‖ψ_k‖² = N *is* the column term.
  - Near-orthogonality controls only off-diagonal Gram entries, and Σ_{k′}|G_{kk′}|² ≈ KN (random), giving
    BHM ≥ N + K√N.
  - So steps 4–5 have nothing to exclude. [P for the constraints, B for the conclusion]
- **(ii) The E_P side has no convolution structure.** The coefficient βχ_P^{−2} is not a prime convolution. By identity
  (∗) it is μG𝔯·conj γ₁(·;k³P⁴): the convolution part (μG) and the character part cannot be separated, because the
  Gauss sum γ₁ is twisted-multiplicative, not a character. energy-wall O1 kills phase-blind use, and diag
  (09 NOTES §19) found theta profiles consistent with all power and BHM constraints for n = 6.
- **(iii) The Möbius side has the structure but needs an asymptotic.** The coefficients μḠχ_P² *are* convolutions
  (Heath-Brown's identity), but the target is F_P = Diag + O(Z²): relative precision Z^{−(ℓ−w)} = Z^{−0.15}.
  - Large-values arguments bound positive moments from above and cannot give second-order precision.
  - The sextic sieve's theta term (HZ)^{2/3}Z = Z^{2.43} exceeds even the diagonal Z^{2.15}. Prop U_P already gives
    the order of magnitude, via theta. [B]

## 8. Route (e): spectral (Louvel, Proskurin, Livné–Patterson) [P/B/H]

- **8.1 No Kloosterman sums arise [P].** Every exact transformation available for E_P produces characters only:
  - Poisson in k: quadratic Gauss sums g₂(n;h) = (h/n)₂g₂(n) are characters times signs. The dual is the same family at
    row length N²/K (hls Prop 2b).
  - Additive expansion of ψ_k, then Voronoi/reflection at the cusps x/k: Σ_x (x/k)₂·κ·e(−x̄m/k) is a sextic Gauss sum,
    a character. That is the reflection.
  - Poisson in h on the Möbius side, or in P: sextic or cubic Gauss sums, then Davenport–Hasse, giving theta or GL(1).

  The Kuznetsov formulas on the cubic cover (Livné–Patterson, Louvel; Proskurin for the double cover) need sums
  Σ_c S(m,n;c)φ of Kloosterman sums over moduli. **None appear.**
- **8.2 Theta is residual spectrum [B].**
  - In a spectral large sieve, θ⊗χ_P^{−2} is the residual (exceptional) form. Dropping the cuspidal and Eisenstein terms
    by positivity bounds |Σ a_nτ(n)|² by (spectral volume + N)‖a‖²: the trivial bound.
  - A gain needs the cuspidal side *evaluated*, which is a GL(2)-type moment beyond its edge.
  - θ is orthogonal to the cusp forms, but S_k = ⟨θ_P, Φ_k⟩ with Φ_k the ψ_k-twisted Poincaré vectors. Unfolding
    returns S_k: circular.
- **8.3 Level aspect [B].**
  - The members θ⊗ψ_kχ_P^{−2} live at K different levels (kP)², or k²P in the U_P picture.
  - A spectral family of a single level containing them has level ≍ lcm, of size ≫ e^{K}.
  - Inside the level-(KL)² spectral family the K twists are a vanishing fraction, and positivity loses the size ratio.
  - A spectral input would have to be level-aspect Lindelöf on average for K quadratic twists of one fixed level-P
    form on the cubic cover. That is (H1) itself. Subconvexity does not help (4.1ii).
- **8.4 The GL(2) analogue is open [L/H].** For a fixed newform g of level q, e.g. g = f⊗χ with f of level 1 and χ of
  prime conductor P:

      Σ_{d~D} |Σ_{n~N} λ_g(n) χ_d(n)|² ≪ DN(Dq)^ε    for D < N < D√q

  is, to my knowledge, open. The known Lindelöf-on-average results for families shorter than their conductor use one of:
  - a continuous parameter with exponential-sum structure (∫_T^{T+H}|ζ|², H ≥ T^{1/3}; p-adic analytic twists,
    Blomer–Milićević);
  - the full character group (BFKMM).

  Quadratic families have neither. The theta analogue additionally lacks an Euler product on the reflected side.

## 9. The minimal missing input

One prime P, one theta vector, and either of the following:

    (MI-1)  ‖Π_V b_P‖² ≪ K Z^ε,      b_P = (β χ_P^{−2} W(q_n/N))_{n sf},  V = span{ψ_k : q_k ~ K},  K < N ≤ K Z^w;
    (MI-2)  Prop U_P with error O(Z^{2+ε}) on H ∈ [ZL·Z^{−w}, ZL).

(MI-1) implies (SM*)_w, hence SM(w) and AVG-SM. (MI-2) is equivalent to (SM*)_w.

- **Truth [N].** energy-wall §4 measured θ·χ_P^{−2} for P = 7, 31, 97, 211, up to N/K = 256 with L ≫ N/K: per-row
  energy 0.74–1.05 of the diagonal, Gaussian row statistics. hls and h2 found random-matrix behaviour. The margin at the
  weighted geometry is about K^{1/2}.
- **What a proof must use.** The joint phases of β against span{ψ_k} (O1), at lengths where reflection lengthens the
  sums. No tool in §4–§8 does this.

## 10. Payoff

Nothing proved here moves θ.
- Prop U_P, Lemma 2, the Poisson-in-P duality and the j = 4 / j = 2 dichotomy are structural.
- The energy-wall threshold is unchanged: a lossless W(≈ 0.0415) closes 41/47 (wall4 §8).
- For W, step 2 now reads: **AVG-SM (needed) ⟸ SM off LZ^{−w} primes; SM ⟺ Prop U_P below ZL; neither suffices
  without (H2).**

## 11. Files

| file | content |
|---|---|
| `NOTES.md` | this file |
| `exponents.py`, `exponents.out` | every exponent and threshold quoted above, at the weighted geometry and three nearby points |
