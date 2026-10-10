# DIAG package — diag agent working notes (10 Oct 2026)

All paper line numbers refer to P78 `paper.tex` (The-Quasi-Riemann-Hypothesis-September-30-2026).
Model runs use my copy `model.py` (= r4 model + one knob `wz`, see §4), baseline config
`{"caseA":"refined","delta0":0.001}` as in lit-combos.

## 0. Verdict

* DIAG is **not provable with current tools**. Neither are the useful partial versions.
* The two losses have **different mechanisms**.
  * In `lem:marked` the cost 2z is "length + conductor". The coordinator's hypothesis is confirmed, with one refinement: the mark is paid **once**, at the terminal theta reflection. The Poisson steps do not pay it again.
  * In `lem:plain` the cost 6κz is "sup-norm via β* × sixth-power exchange rate 6". There is no conductor effect, so the hypothesis is refuted there.
* There is a **rigorous negative result**. Lemma marked in the paper's coefficient-uniform form (any bounded slot coefficients) is **false** for z > m/2, whatever r is: Gauss-phase slots are extremisers of the sextic large sieve. So any DIAG proof must use the Hecke structure of the slot coefficients. This is not yet fatal, because the model only uses z ≤ S ≈ 0.23 < m/2.
* The OpenAI technique (b) **cannot work for n = 6**. Its large-values LP admits the theta profile, and near-orthogonality (BHM) excludes theta profiles only when the excess slope 1/n is above 1/4, i.e. only for n = 3. This holds even with square-root Gram sums.
* **Payoff ceiling.** At the DIAG optimum every binding bin has mark amplitude x = 0. So θ_DIAG = θ(x ≡ 0) = 0.87255, and **no slot-capacity lever of any strength can beat 0.87255**. This is checked in §4.

## 1. Re-verification (all reproduced exactly)

| config (+ caseA refined, δ0 = 1e-3) | θ |
|---|---|
| baseline | 0.8749568 |
| DIAG (mk=1, w2scale=10, cP=1) | **0.8725498** |
| mk=0.75 alone (2r+8z cap kept) | 0.8734852 |
| cP=1 alone | 0.8734526 |
| DIAG + wE=0.0196 | 0.8709686 (27/31 = 0.8709677, so W(0.0196) is right up to ±3e-4) |

New partial points:

| config | θ |
|---|---|
| mk = 0.6 | 0.873988 |
| mk = 0.7 | 0.873486 (the cap 2r+8z ≤ 3 binds from 0.7 on) |
| mk=1, w2scale=10 (marked part of DIAG only) | 0.873159 |
| wz=0.1 (z ≤ (1+wz−r)/2, i.e. r+2z ≤ 1+wz), w2scale=10 | 0.873436 |
| wz=0.2, w2scale=10 | 0.873159 (already saturated) |
| cP = 0.3 / 0.4 / 0.5 | 0.874789 / 0.874523 / 0.874177 |
| mk = cP = 5, w2scale=10 | 0.8725498113 (= DIAG, bit-identical) |
| mark amplitude capped x ≡ 0 (`xcap0.py`) | 0.8725498113 (= DIAG, bit-identical) |

Binding anatomy (`diag_anatomy.py`, `anat_*.txt`):

* **Baseline.** The binding bin is δ = 0.391, x = 1/2. It has two co-binding witnesses:
  * marked inverse with r = 0.7165 and z_I = 0.1417 = (1−r)/2, against plain m = 0.4085 with z_P = 0.0407 = (2/9)(1−2m);
  * the pure long inverse witness r = t = 1.123.
  * Supply S = 0.205.
* **DIAG.** Every binding bin has **x = 0**, near δ ≈ 0.33. The witnesses are r ≈ 0.77, m ≈ 0.38, and the long r = t ≈ 1.15. Slot lengths are irrelevant there, because q = xδ = 0.

## 2. Where the losses are

### 2.1 Lemma marked: r + 2z ≤ m − c₁ and 2r + 8z ≤ 3m − c₂ (l.9254)

**Step 1 — initialization: one Poisson in u (proof of lem:marked, l.11606ff).**
* The slots are merged into the Möbius column. With j = (n,P), write n = jn₀ and P = jP₀; then c = n₀P₀ carries μ(c) and the mark 1_{P₀|c}. The definitions are G = Σ_assigned z_i, z₀ = z − G and D' = r + z − 2G (l.11624–11626). So the column has length D' = r+z (G = 0): **the slot pays its length once here**. That is unavoidable, since it is the length of the product.
* Poisson in u modulo the conductor z₁z₂ gives the dual row scale M_c^init = 2D' − m − 2P₁ (eq. `inverse-initial-row`, l.11751). The canonical lengths are F_c^init = D' − P₁ (eq. `initial-canonical-lengths`, l.11909).
* μ(c) times the Gauss sum is absorbed into the theta-type coefficient ᾱγ₂ (eq. `initial-common-extraction`).
* Then (l.12207–12209):
  `F − M − (P₁+G) − z₀ = m − r − 2z + 2G`, `4F − 3M − 6z₀ = 3m − 2r − 8z + 10G + 2P₁`.
  These are exactly the two hypotheses of lem:marked once they are fed into the canonical invariant.

**Step 2 — the canonical invariant `eq:invariant` (l.9317)** is `3M + 6z₀ + c ≤ 4F` together with the puncture bound `q_r ≤ Z^{F−M−z₀−c}`. At each node it reads `F−M−Q−z₀ ≥ c_node` and `4F−3M−6z₀ ≥ c_node` (eq. `inverse-fixed-invariant`, l.9807).

**Step 3 — the recursion does not pay for marks again.**
* eq. `canonical-transition` (l.11051) gives `F_c − M_c = F − M + D_c` with D_c ≥ 0, and `4F_c−3M_c = 4F−3M + 2(A_i+t+g−θ+V)+j`.
* Also (l.11100): "The nominal cap z₀ is unchanged because the surviving mark is a subcollection."
* So the margin F−M only grows down the recursion. The binding place is the terminal.

**Step 4 — the terminal (reflected energy), where the second z is paid.**
* Cor. `marked-reflection` (eq. `marked-reflection`, l.7374–7383) gives the dual test `V^♯(q_μ X / q_c²)`, with c = c_F(h₀)·r and r ⊇ the *active marked primes*. Each dual term also carries `q_q^{-1/2} χ_q(λ⁴μ)^{-2}`. **An active mark p enters the conductor squared.** For the theta at length X = N (full length, mark included), the dual length is (KP)²/N rather than K²/N.
* In lem:reflected-energy, eq. `old-eq:4.4` (l.7759) reads `T_d = 2H + 2A₀ + 2z_a − N_* − N₀ − 3B₀`.
* The reflected sum is then bounded by the hybrid quadratic–cubic sieve `old-eq:4.6b` (l.7523): `(K+NB)B{N + L + (NL)^{2/3}}`. Its quadratic factor (K + NB) has a **column term**.
* In the canonical proof, eq. `inverse-terminal-width` (l.9884) feeds into `old-eq:4.9` (l.9892): `T_d − S₀ − B₀ ≤ M − O + 2z_a − z₀ − c_node + …`. The column branch u = z_a, eq. `old-eq:4.10` (l.9910), needs the **first invariant**. The branch u ≠ z_a (marks dominate, or the cubic (NL)^{2/3} term), eq. `old-eq:4.11` (l.9921), needs the **second invariant** `3M + 6z₀ ≤ 4F`.

**Clean exponent form.** Canonical rows K = Z^M, column N, marks P = Z^{z₀}:
* reflected length T = (KP)²/N;
* quadratic sieve bound (K + T), needed to be ≤ N;
* hence **KP ≤ N**, i.e. M + z₀ ≤ F.
At the init, K = N²/U with N = Z^{r+z} and U = Z^m, so this is r + 2z ≤ m.

**Equivalent picture.** Twisted multiplicativity τ(pn') = τ(p)τ(n')·χ_{n'}(p)⁴ (eq. `initial-common-extraction`) shows that a mark turns the row k into kp⁴. That is a Kubota shift whose conductor contributes p, and the length drops from N to N/P. **So a mark is an extra row variable, and the canonical (theta) problem is diagonal only when rows ≤ columns.**

### 2.2 Lemma plain: n₁+n₂+6κz ≤ M (old-eq:3.9, l.12562)

* **Pointwise slot bound `old-eq:3.5` (l.12876):** |Q_ω|² ≤ C Z^{κz+ε}, with κ = 2β*−1 = 3/4. It comes from L'/L at Re s = (1+κ)/2 ≥ β* (lem:logarithmic-control), applied to Hecke-twisted prime sums.
* **Slot removal `old-eq:3.16` (l.14256–14276):** "A removed slot of length d decreases A_act + (6κ−1)z_act by 6κd and contributes κd to the squared Z-exponent". The removal is needed when a transformed child violates the affine region by F_act ≤ 6(w+ℓ) + … (l.14236). The parent's slack is only (w+ℓ).
  * The 6 is the sixth-power exchange rate: the extraction ledger has ℓ_p ≥ σ/6 and g ≤ d + (K₀−K) + 6ℓ (l.14190–14205), and the A ≤ 5M/6 threshold of the Θ-rows (l.12594–12598).
* **Exceptional rows are not the problem.** The centered exceptional deficit is ≤ (A−M)₊ (old-eq:2.19, l.14760). Slots there are summed absolutely and need only A ≤ M, so cP = 1 would be compatible with them.
* **Within this architecture**, cP = 1/(6κ) = 1/(6(2β*−1)). Hence **cP = 1 ⟺ β* = 7/12** (a pointwise prime-sum bound P^{7/12}, i.e. a QRH at 7/12), and cP = 1/2 ⟺ β* = 2/3. Mean-square removal instead of sup removal costs rate 1 per unit (the column term), which is worse than κ.
* **Verdict.** The plain loss is not a conductor effect.

### 2.3 Answer to the hypothesis
* **lem:marked: confirmed**, but the mark is paid once, at the terminal reflection (conductor (kp)², length N). It is not paid in each Poisson step: the init Poisson charges only the length, and the recursion charges nothing.
* **lem:plain: refuted**, as explained in §2.2.

## 3. Is DIAG within reach?

### 3.1 Necessary condition: coefficient-uniform DIAG is false for z > m/2 (PROVED, standard Poisson)

**Example.**
* Restrict the primes to one class mod 36, where the units and the reciprocity cross phase are trivial (checked numerically, `recip_check.py`).
* Take r = 0, one slot of length z, and a(p) = conj(g̃₆(p)), the normalised sextic Gauss sum with additive character e(Tr(x/(π√−3))).
* Then Q_u = conj(P^{-1/2} Σ_p g̃₆(p,u)).
* Poisson in u (with a smooth minorant Φ ≤ 1_{q_u≤U} that has Φ̂ ≥ 0) gives
  `Σ_u|Q_u|² = U·(diag) + (U/P²)·Σ_{0<q_h≲P²/U} |Σ_p χ_p(h)|² Φ̂(·) − O(P)`.
  All terms are ≥ 0. The h = w⁶ terms give |Σ_p 1|² = P², and there are ≍ (P²/U)^{1/6} of them.
* Hence `Σ_{q_u≤U}|Q_u|² ≫ U^{5/6}P^{1/3}`, which is ≫ U^{1+δ} once z > m/2. This is de Faveri's theta term, attained, consistent with de Faveri–Dunn–Hoffstein.

**Numerics (`theta_slot_test.py`).**
* Cubic (n = 3, U = 1000, cube-free rows). The excess of the Gauss-phase coefficients over the diagonal is 0.28 / 0.58 / 1.16 / 2.87 at P = 1k / 4k / 16k / 64k, growing like P^{2/3}. Random and trivial coefficients stay at 1.
* Sextic (n = 6, U = 1000). The excess is 0.10 / 0.02 / 0.14 / 0.35 / 0.39 / 0.50 at P = 1k … 256k: slow P^{1/3} growth, and random coefficients stay at 1. See `theta_slot_sextic.txt`.

**Consequences.**
* P78's lem:marked is stated uniformly in |a_i(p)| ≤ 1. Its 8z ≤ 3m cap at r = 0 (z ≤ 3m/8) is therefore not far from the truth (m/2).
* A heuristic alignment construction (a(p) chosen against F(p) = Σ_n μ(n)g̃(n)(p/n)₃) gives a failure threshold that depends on the size of F(p):
  * z > m/2 + 2r if F(p) ~ D^{1/2} (square-root size);
  * z > m/2 for every r if F(p) has the "signal" size D^{5/6}.
* Either way this is outside the model's range (z ≤ S ≈ 0.23), so in the needed range the uniform version may still be true.

### 3.2 A rigorous but useless improvement
* de Faveri Thm 1.1 (n = 6, sixth-power-free; LIT.md row 1), applied by duality, gives lem:marked for all bounded slot coefficients whenever r + z ≤ m/2 − c.
  * All rows are allowed: non-sixth-power-free rows u = w⁶u₀ contribute ≪ U^{1/6}·N ≤ U.
  * At r = 0 this improves z < 3m/8 to z < m/2, and §3.1 shows that is sharp for uniform coefficients.
* **Model payoff 0**, because the inverse witnesses have r ≥ t − 1/2 ≥ 1/2.
* Applying the sextic sieve to the Poisson dual instead reaches only r + z ≤ 2m/3, provided the near-principal dual rows are controlled; the payoff is still 0.

### 3.3 (a) Recursion variants — all reproduce or worsen KP ≤ N
* **Slots as an averaged label.** A mark p extracted from the column is exactly the canonical label twist χ_n(f)⁴ with f = p (eq. `initial-common-extraction`). But p sits inside the square with coefficient a(p)τ(p)χ_p(k), so Cauchy–Schwarz over p is needed. That gives E ≤ P·Z^{N} (with labelled F' = N−z+z), a loss of P. Worse.
* **Hölder before reflecting.**
  * Σ|MQ|² ≤ (Σ|MQ^j|²)^{1/j}(Σ|M|²)^{1−1/j} uses j copies at total length jz, which is the same capacity.
  * The (Σ|M|^{2p})^{1/p} variants need Möbius powers μ^{*p}, which are not μ. μ-absorption into theta is lost, and the Möbius×Möbius cross term τ(d₁)τ(d₂)(d₁/d₂)₃ is again a mark of cost 2.
* **Splitting slots between the two Poisson steps.** Every slot depends on u, so it must be present at the init Poisson. The recursion never touches marks (z₀ is invariant).
* **Vaughan/Heath-Brown on the slot primes.** This produces Möbius pieces whose total length can reach z, so the worst term is unchanged. Smooth pieces ride multiplicatively and gain nothing.
* **Twisted multiplicativity n = pn'.** The theta row shifts to kp⁴: Kubota shift p plus length N/P gives T = K²P²/N, identical to the paper.

### 3.4 (b) OpenAI "moments of different lengths + near-orthogonality" — the LP does not close (PROVED as an LP statement; `lp_bhm.py`)

**Setup** (log_U units). Factors from the short identity have lengths a_j (Σa_j = a = r+z) and normalised excesses W_j on U^ρ rows. DIAG failure means ρ + 2ΣW_j ≥ 1 + ε₀.

**Constraints.**
* **(LS)** for every multiplicity vector k: ρ + 2Σk_jW_j ≤ L₆(Σk_ja_j), where L₆(z) = max(1, z, 5/6 + z/3, 1/3 + 5z/6).
* **(BHM)** with the true Gram size g = 1+z: ρ + 2W ≤ max(z, (ρ+1+z)/2).

**The theta profile is feasible.** The profile ρ = 5/6, W_j = a_j/6 satisfies **all** LS constraints for every factorisation, because it is linear in the length. It also satisfies BHM for every z. Yet it violates the diagonal for a > 1/2.

**General n.** BHM excludes the theta profile (ρ = 1−1/n, W = z/n) iff z(2/n − 1/2) > 1/(2n). That works only for n = 3 (z ∈ (1,2); OpenAI use (4/3,2) with their Gram bound). It never works for n ≥ 4, even with an R-restricted square-root Gram.

**Why.** The sextic theta excess per unit length is 1/6, half the cubic 1/3. That is below the square-root noise floor of any Gram/BHM argument.

**There is also no rigidity.** OpenAI argue at the extremal point (moment = sharp sieve bound), where equality forces every factor onto the theta profile. DIAG needs the whole range between the diagonal and the sieve bound, and there the LP has a continuum of feasible profiles.

**Consistency with the signal.** Möbius correlates with the sextic theta at relative level N^{-1/6}. Its theta-direction contribution to the moment is N^{-1/3}·U^{5/6}N^{1/3} = U^{5/6} < U. So DIAG does not contradict the probe's signal, but an argument that "kills all theta correlation" of the Möbius pieces would. A correct argument must tolerate W_M ≈ 0 with correlation N^{-1/6}, and that is invisible to BHM anyway.

### 3.5 (c) The β* bootstrap gives nothing
* Any pointwise use of β* (on slots, or on M_u) costs κ = 2β*−1 = 3/4 per unit length in the moment exponent.
* A slot's count gain is 2q = 2xδ ≤ δ ≈ 0.39 per unit at the binding bins. Since 0.39 < 0.75, that is a net loss.
* β* says nothing about Gauss-sum correlations (the primal or dual obstruction) or about the quadratic column term (the reflected obstruction).
* κ-feedback with θ_DIAG gives cP = 1/(6·0.745) ≈ 0.224, which is negligible.

### 3.6 (d) Partial relaxations
* Nothing beyond the paper is provable except §3.2, whose payoff is 0.
* What the model would reward (from §1): mk 0.6 → 0.87399; mk 0.7 → 0.87349; marked part of DIAG alone → 0.87316; wz 0.1 → 0.87344; cP 0.5 → 0.87418; full DIAG → 0.87255.
* **Saturation.** At the DIAG optimum all binding bins have x = 0, so θ_DIAG = θ(x ≡ 0) = 0.872550. Checked bit-identically: `xcap0.out` and `sweep_extreme.jsonl` (mk = cP = 5) both give 0.8725498113201074. Adding the x > 0 constraints can only raise θ, so θ ≥ θ(x≡0) for every slot capacity. **No lever acting on slot capacities (mk, cP, wz, w2) can go below 0.87255.** 27/31 needs the wall W(0.0196) or another non-slot lever on top.

## 4. The precise obstruction

**Marked part.** It would suffice to prove the canonical estimate with the invariant **without z₀**:

(★) Σ_{q_k≪Z^M} |Z^{-N/2} Σ_{n sf} ᾱγ₂(n)νρ(n) χ_n(k) 𝔡(n) W(q_n/Z^N)|² ≪ Z^{N+ε}, for M < N < M + z₀ (with labels/punctures as in lem:canonical-moment).

* At the init this gives F − M − Q = m − r − z + G and 4F − 3M = 3m − 2(r+z), which is exactly DIAG-marked.
* In the unreflected picture, (★) is at the edge: rows K < columns N, and the needed bound is the column term. But it needs a coefficient-sensitive sextic sieve, because the theta terms K^{1/3}N^{5/6} and K^{5/6}N^{1/3} exceed N for K > N^{1/2} and are sharp for general coefficients.
* After the (unavoidable) marked theta reflection it becomes a quadratic-character sum with K rows and T = (KP)²/N > N columns, and needs a bound below its column term T. **That is the same shape as the energy wall W(w).** The needed "wall excess" is 2(r + 2z − m) ≤ 2z ≈ 0.46 in family units, much more than the W(0.064) needed for the probe.
* μ-absorption has already used up the Möbius structure, so no arithmetic of μ remains.

**Plain part.** One of two things is needed:
1. a treatment of the sixth-power-full transformed rows that does not shrink the child family (exchange rate below 6); or
2. a pointwise slot bound with κ ≤ 1/6 (β* ≤ 7/12).

Both are far out of reach; the second is circular.

## 5. Files
* `model.py`: r4 copy plus the knob `wz` (marked capacity r+2z ≤ 1+wz, keeping r < 1+w1).
* `run_cfgs.sh`: batch runner.
* Model outputs: `verify.jsonl`, `sweep_*.jsonl`, `xcap0.py`/`xcap0.out`.
* `diag_anatomy.py`, `anat_base.txt`, `anat_diag.txt`: binding witnesses and slot lengths.
* `theta_slot_test.py`, `recip_check.py`, `theta_slot_sextic.txt`: Gauss-phase slot counterexample numerics (n = 3, 6).
* `lp_bhm.py`, `lp_bhm.out`: the (b) LP / BHM barrier.
