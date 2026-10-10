# Agent `amplitude` — exploiting large-prime-amplitude rows (10 Oct 2026)

Scripts in this directory: `amp_model.py` (gmodel extended: supply, auxiliary amplitude, mean-treatment
of marks, ZD hypothesis, prime-only witnesses, vectorised count LP with 201 t-steps × 161 r-steps),
`lpanat.py`, `contour_order.py`, `runs.py`, `runs2.py`, `thresholds.py`, `capcurve*.py`,
`supplyneed.py`, `mkstar.py`; outputs in `results.txt`. The baseline at ℓ = 0.1773 on the finer grid is
**+0.01119** (δ = 0.475, x = 1/2, b = 0.075); gmodel's coarse grid gives +0.0117.

## 0. Verdict

No rigorous improvement. All four directions give **zero exponent payoff** in the weighted model,
and the reasons are proved at the level of the exponent model (§2–§5). One new rigorous lemma came out
of this: **polarised slot powers** in lem:marked (§3.3). It shows that slot *supply* never limits a
coefficient-blind marked moment. It pays nothing at the paper's capacities, but it changes what a
future capacity improvement must deliver (§3.4).

## 1. What binds at x = 1/2 (`lpanat.py`, `runs.py`)

At ℓ = 0.1773 and b = 0.075 we have L = ℓ/h = 0.2207. The count LP at x = 1/2 is a three-way balance:

| δ | R | t* | adversary's split (r, m) | needed R |
|---|---|---|---|---|
| 0.30 | 0.7467 | 1.087 | (0.687, 0.400) short | 0.7359 |
| 0.40 | 0.6563 | 1.130 | (1.130, 0) pure long | 0.6391 |
| 0.47 | 0.5895 | 1.163 | (0.747, 0.416) short ≡ long | 0.5713 |
| 0.55 | 0.5100 | 1.210 | (0.782, 0.428) | 0.4938 |
| 0.60 | 0.4578 | 1.245 | (0.808, 0.438) | 0.4454 |

Three facts about this table:

- **The long branch has no slot capacity.** When r ≥ 1 there is no room for slots. Slots inside the
  sixth-power amplification lose (5/3 − 2q)z.
- **The short branch has small capacities.** At δ = 0.47 the binding split uses marked capacity
  (1−r)/2 = 0.127 and plain capacity (1−2m)/(6κ) = 0.037.
- **Supply is not the constraint.** Both capacities are below L, and below 7/37 for every inverse
  split the LP uses.

## 2. Direction (3): marks entering the high side "in mean" — REFUTED (PROVED)

**Prop 2.1 (the amplitude is a two-sided per-row size).** Fix a retained tuple and an amplitude class
C with slot amplitudes g_i, and let I₊ = {i : g_i > 0}. lem:high-bin defines the class so that
P_i^{g_i} ≤ |Q_i(u)| ≤ P_i^{g_i+ϑ} for i ∈ I₊ and u ∈ C. Then

  Σ_{u∈C} |N_u| ∏_i |Q_i(u)| ≥ Z^{qℓ} Σ_{u∈C} |N_u| ∏_{i∉I₊} |Q_i(u)|.

Hence any bound for the class obtained by the triangle inequality over rows contains the factor
Z^{qℓ}. Such a bound can improve only through #C or through the remaining factors (numerator,
zero-amplitude and error slots).

The LP is attained at uniform amplitude g_i = δ/2, where I₊ is all slots. There the upper and lower
bounds coincide up to Z^{ϑℓ}. *Proof:* this is the definition of the amplitude bins
(P78 §"Prime amplitudes") and multiplication over i. ∎

**Prop 2.2 (a moment in Q instead of pointwise loses).** Replace the pointwise factor by
Cauchy–Schwarz against Σ_u |Q_z(u)|² ≪ U^{1+ε} (lem:marked with r = 0, z ≤ 3/8; or the sextic large
sieve, z ≤ 1/2). Relative to the class bound with #C = U^R, this changes the exponent by
(1 − R − 2qz)/2. That is a loss whenever R + 2qz < 1.

Over all bins δ ∈ [0.05, 0.75], x ∈ [0, 1/2] at the optimal geometry we get
max (R + 2q·min(L, 3/8)) = 0.970 and max (R + q) = 0.984, both attained at δ = 0.05. So the change is
a loss in every bin. *Proof:* Hölder with the stated moment, compared with #C·U^{2qz} on C. ∎

The one mean treatment that does pay is the selected-slot fourth moment of 08, which is already in
E_wt. Every other Hölder or level-set combination of the available moments is log-linear in the same
exponents. The extreme points are the 4th-moment and 2nd-moment numerator bounds, and the
numerator × plain-witness bound (lem:plain applied to S(n₀)·S_m·Q). The last is exactly tight at the
plain/marked crossing: its gain needs z > (1 − m − n₀)/(6κ), which is the maximum allowed. Model check:
`mean_marks=True` gives an unchanged worst exponent, +0.01119.

## 3. Direction (1): explicit-formula supply and rigidity

**Prop 3.1 (capacity, not supply; PROVED in the model).** At ℓ = 0.1773, L = 0.2207 exceeds every
selected capacity used by the LP. Therefore:

- (a) R(δ, 1/2) is unchanged with unlimited supply (all δ; `runs.py`).
- (b) Suppose any auxiliary windows (explicit-formula windows, slot powers, other scales) gave
  *every* row amplitude-1/2 counts with unlimited length, while the high side still pays the physical
  amplitude. The worst exponent is still +0.01119 (`xcount=0.5, Lc=∞`).
- (c) Scanning b ∈ [0, 0.45]: supply matters only for b ≥ 0.3, where L < 7/37. There the worst
  exponent is ≥ +0.023 (`runs2.py`).

*Proof of (b):* E(δ, x) depends on the count only through R. The amplitude-1/2 count at the worst bin
x = 1/2 is unchanged by (a), and E is unchanged there. ∎

**Prop 3.2 (no exponent-level rigidity from the explicit formula).** Fix a row with zero ρ₀ = a + iγ₀
and the prime sum normalised at ρ₀ on the log-scale y = log P:

  F(y) = Σ_ρ Ŵ(ρ − ρ₀) e^{(ρ−ρ₀)y} + small,

where the ρ₀ term is the constant Ŵ(0).

(i) **Nazarov–Turán.** Write F_B for the sum over the box B = {β ≥ a − η, |γ − γ₀| ≤ H} and let E ⊂ I be
intervals. Then sup_I |F_B| ≤ e^{2η|I|} (A|I|/|E|)^{n−1} sup_E |F_B|, where n = #B. Suppose
|F(y₀)| ≍ 1 at one physical slot scale y₀ ∈ I. Then |F| ≳ (2A)^{−n} e^{−2η|I|} on at least half of I.

(ii) **Jensen** (disks centred at a + 6e + iγ, with |L| ≥ U^{−ε} at the centre and the bin's
Phragmén–Lindelöf bound on the disk) gives only n ≲ H·log U / log 2. No better single-row bound on the
number of zeros near a in a unit height window is known for any σ < 1.

(iii) **Tail control.** The zeros outside B must contribute ≤ U^{−c}. This forces H ≳ 1 with
Gaussian log-windows of width √(log U), and H ≥ U^{c'} with compactly supported windows.

Hence the transfer of saturation from one scale to a positive proportion of other scales costs
U^{−κ₀}, with κ₀ ≍ H log(2A)/log 2 = O(1). That is useless at exponent level.

Montgomery–Vaughan ("large on average over log P") needs zero gaps ≫ 1/|I|. With ≍ log U zeros per
unit height it gives only O(log U) information. Zero multiplicities are ≤ O(log U) = U^{o(1)}, so a
dichotomy "few zeros / many zeros" cannot give a count gain on the many-zeros side.

Status: PROVED as a statement about what (Nazarov–Turán, Jensen, MV) can yield. That no other
single-row tool exists is HEURISTIC.

Two remarks, the second HEURISTIC:

- A natural row (one dominant off-line zero) *does* have x = 1/2, with inverse sums saturated at all
  lengths. That gives DH counts.
- The LP's adversary at x = 1/2 must therefore be zero-dominated at the tiny slot scales but saturated
  in M_r and S_m only at the detector's lengths. Nothing above excludes this. It needs ≍ log U zeros
  near a per unit height, which Jensen allows.

**Prop 3.3 (polarised slot powers in lem:marked — PROVED, new).** Let 𝒫 be a slot (bounded
coefficients a(p), scale Z^w) disjoint from the other slots Q_rest, and let k ≥ 1 be fixed. Put

  e_k(u) = Z^{−kw/2} Σ_{p₁<…<p_k ∈ 𝒫} ∏_j a(p_j) ψ_u(p_j) W(q_{p_j}/Z^w).

Under the hypotheses of lem:marked, with z replaced by z_rest + kw,

  Σ_u |M_u(Z^r;W) e_k(u) Q_rest(u)|² ≤ (k^k/k!)² sup_f Σ_u |M_u ∏_{c≤k} Q^{(f,c)} Q_rest|² ≪ Z^{m+ε}.

Here f runs over colourings 𝒫 → {1,…,k}, and Q^{(f,c)} is the slot with coefficient a(p)·1_{f(p)=c}.

Moreover |Q^k − k! e_k| ≤ C_k max(1, |Q|)^{k−2}. So on an amplitude class with g > 0,
|e_k| ≥ |Q|^k/(2·k!).

*Proof.* For a uniform random colouring, E_f ∏_c Q^{(f,c)} = k^{−k} Σ_{distinct ordered tuples} =
k^{−k} k! e_k. Then |E X|² ≤ E|X|² and Fubini give the first inequality. Each colour class is a slot
with disjoint support at scale Z^w and coefficients bounded by 1. lem:marked is uniform in prime
supports and bounded coefficients ("uniform in the prime supports and their bounded coefficients",
P78 l.9263), so the sup is ≪ Z^{m+ε}. The number of colourings is irrelevant because only an average
is used.

For the lower bound, expand the tuples with a repeated prime through set partitions and power sums
p_j = Σ c_p^j. A partition with s singletons and t blocks of size ≥ 2 contributes
≤ |Q|^s Z^{w(t − (k−s)/2)} ≤ |Q|^s, with s ≤ k − 2. ∎

**Consequence.** For the counting witness M_r·∏ e_{k_i}(Q_i), supply in the marked count is unlimited
at the physical amplitude. At the paper's capacity this pays nothing, by Prop 3.1.

Limitation (PROVED): the colour masks are not Θ-combinations, so lem:plain does not allow this.
HEURISTIC: Θ-structured (ray-class) colourings fail, because the adversary can concentrate Q in one
class.

**Prop 3.4 (what a capacity improvement must now deliver; CONDITIONAL, model).** Let mk denote
inverse capacity z ≤ mk·(1 − r) with the plain capacity left at 2/9 (`mkstar.py`, `supplyneed.py`):

| marked hypothesis | worst |
|---|---|
| mk = 1, physical supply L, cubic branch 2r + 8z ≤ 3 dropped | +0.0017 |
| mk ∈ [0.75, 1] needs supply L ≥ 0.2458 (physical 0.2207) | — |
| unlimited supply (= polarised powers) and cubic branch dropped | **mk ≥ 0.744 alone closes 41/47** |
| unlimited supply with cubic branch kept, even mk = 1 | +0.0043 |

So a marked improvement that is uniform in slot coefficients, as lem:marked is now, would need only
mk ≈ 0.744 together with removal of the cubic branch, and no plain improvement. The BRIEF's
mk = 0.8 + cP = 0.4 is replaced by this.

But coefficient-blind DIAG is false at r = 0 (the 09 Gauss-sum example), and the same theta term
plausibly makes (1−r)/2 sharp for all r (HEURISTIC). A Hecke-only improvement cannot use the masked
polarisation. It must itself be stated for products of slots *with repetition* (Q_i^k), which is
natural for any proof that uses the L-function structure of Σ_p λ(p) χ_p(u). Otherwise the supply
L = 0.2207 binds, and mk = 0.8 still needs cP ≥ 0.422 (`capcurve.py`).

## 4. Direction (2): the count LP is complete for the existing lemmas (PROVED, exponent level)

Each item below gives no gain over the LP.

- (a) **Witnesses from all t ∈ [1, 3/2] at once.** The row's ρ and ψ are the same for every t. Without
  a joint moment the count is min_t max_split, which equals the single-t LP.
- (b) **lem:plain with two different plain witnesses S_{m₁} S_{m₂}.** The count
  1 − δA − 2q(1 − A)/(6κ), with A = m₁ + m₂, is decreasing in A because 2q/(6κ) < δ. Two copies of the
  longer witness are better. One copy plus the larger capacity (1 − m)/(6κ) is worse by
  m(δ − 2q/(6κ)) ≥ 0.
- (c) **Slots beyond capacity.** Bounded pointwise by the bin they are neutral at x = 1/2 and lose for
  x < 1/2. Bounded by β* inside the lemmas they lose κ − 2q per unit.
- (d) **Slots in the amplification.** These lose (5/3 − 2q)z. At the capacity boundary they equal the
  marked count.
- (e) **Prime-only witnesses** (marked with r = 0, sextic large sieve, polarised powers, Turán powers).
  Capacity ≤ 1/2 gives count ≥ 1 − δ/2, never binding (`cQ` runs). With auxiliary amplitude-1/2
  supply, a prime-only capacity ≥ 0.904 would be needed.
- (f) **Detector at t < 1.** R_short increases as t decreases.
- (g) **Choosing the w-contour per amplitude class instead of per δ.** Probably legitimate, since the
  w-move is per row at fixed (s, z). Its payoff at ℓ = 0.1773 is zero (`contour_order.py`).

## 5. Direction (4): slot redesign — invariant at exponent level (PROVED, with one sharpness caveat)

E_old, E_wt and the selection lemma depend on the slot system only through four quantities:

- ℓ = Σℓ_i;
- the mean amplitude q;
- the supply L;
- the mesh max_i ℓ_i, which only enters as a vanishing loss.

Hence scales and number are irrelevant. On weights and signs:

- **lem:plain needs Θ-span coefficients.** For these, each character component contributes
  −Σ_ρ P^{ρ−1/2} Ŵ(ρ + z − 1).
- **Amplification needs a nonzero principal component.** The principal amplification
  Σ_p W p^{−5/6} ≍ P^{1/6} requires the principal coefficient c₁ ≠ 0.
- **So the δ/2 bound is reached.** A zero of the c₁-component at a, at a height where Ŵ ≠ 0 (such
  heights exist for any W ≠ 0), gives |Q| ≍ P^{δ/2}.
- **Non-Hecke signs are worse.** They lose the zero-free bound altogether (trivial P^{1/2}).

So the bin information never yields a worst amplitude below δ/2 per unit length. The caveat: this
sharpness is relative to the bin hypotheses, not a theorem about actual L-functions.

Hypothetical caps for reference (`thresholds.py`):

- An amplitude cap x ≤ 0.255 would close 41/47 (x ≤ 0.25 gives −0.0002; x ≤ 0.3 gives +0.0018).
- Equivalently, ZD(x₀): DH counts for every row with x ≥ x₀ closes iff x₀ ≤ 0.255.

Both say that zero-dominated rows obey DH. That is a zero-statistics statement, not a design choice.

## 6. The precise obstruction

At the binding bins (δ ≈ 0.4–0.55, x = 1/2) the count is the balance of two branches.

- **Long branch.** The pure long inverse witness, slope α = 5/6, set by the sparsity of the copies
  u·a⁶.
- **Short branch.** Marked capacity (1−r)/2, from the mark conductor p² plus the cubic branch, and
  plain capacity (1−2m)/(6κ), from β*-pointwise slot removal in child rows.

Amplitude information is per row and single-scale. The explicit formula converts it into other
witnesses only with a U^{−O(1)} loss (Prop 3.2). Supply was never the constraint, and is now provably
unlimited for coefficient-blind marked moments (Prop 3.3).

Any route through large-amplitude rows therefore needs one of three inputs:

- (i) a Hecke-structured marked capacity ≥ 0.744(1−r) that holds for slot powers and removes the cubic
  branch;
- (ii) a single-row local zero bound n(a − η, γ₀ ± 1) = o(log U) for bin rows, plus Möbius log-control,
  which would give rigidity;
- (iii) the DH-for-zero-dominated-rows statement ZD(0.255).

None of these is available.
