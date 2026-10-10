# plain-cap — can lem:plain's prime capacity 1/(6κ) be improved? (10 Oct 2026)

Agent `plain-cap`. Sources: P78 `paper.tex` §"Fourth moments with short prime factors" (l.12477–14985),
Hailey Collet `public_proofs_paper.md` §B3.4–B3.10 (same proof for κ ∈ [7/10,3/4]), model `gmodel.py`.
Scripts here: `anatomy.py`, `scan_cap.py`, `quantify.py`, `game.py`, `game_asym.py`, `model_wr.py`, `thr_wr.py`.

## 0. Verdict

**Negative.** Within lem:plain's proof (reflections, two-transform recursion with sixth-power amplification,
pointwise slot removal) the region n₁+n₂+6κz ≤ M is *exactly optimal* for the symmetric plain products that the
count and the numerator use (Theorem G below, PROVED at the level of the move calculus, confirmed by DP).
Every route in the mission brief either is a re-ordering covered by Theorem G, or needs one of three inputs:
a cheaper slot removal (a twisted fourth moment of self-dual sums "at the edge"), a denser amplifier
(impossible: sixth powers are the only multipliers preserving all Gauss sums), or an asymmetric witness (not
supplied by the detector). Small proved items (§4) have zero payoff. Conditional quantification in §5.

## 1. Where 6κ comes from (line-precise)

* **κ** — eq. (old-eq:3.5), l.12876–12880: on eligible rows k ∈ R_z, |Q_ω|² ≤ C Z^{κz+ε₁}, from
  L'/L(s_κ+e+it,ψ) ≪ log on Re s = s_κ = (1+κ)/2 ≥ β* (old-eq:3.4, l.12838–12848); κ = 2s_κ−1.
* **Amplifier** — l.13448–13490: pool of primes of norm ≍ Z^{ℓ*}, ℓ* = σ/3; Jensen
  |H(h)|² ≤ avg_p{|H(hp⁶)|² + errors}; G(u,hp⁶) = G(u,h); new rows in a ball of length K+6ℓ*; density factor
  Z^{−ℓ*}; "ℓ = ℓ*, g = J₊+2σ" (2σ = 6ℓ*).
* **Errors** — (old-eq:2.9), l.13515: (w,w_o) = (iℓ_p, e_iℓ_p), i = 1,6,7, ℓ = 0, g = J₊+σ.
* **Child allowance** — (old-eq:2.14), l.13906: Δ_child = b₂−p₂+w+B_c+ℓ ≥ w+ℓ.
* **The "6"** — (old-eq:3.14), l.14236–14259: (A'−M')−(A−M) = g+w−d−c₂−(K₀−K)−w_o+g₂−t₂ ≤ 6(w+ℓ)+δ_fr,1.
  Main norm: g ≤ d+(K₀−K)+6ℓ+δ. Errors: w−w_o+σ ≤ 6w because ℓ_p ≥ σ/6.
* **Exchange** — (old-eq:3.16), l.14276–14300: removing a slot of length d lowers A+(6κ−1)z by 6κd and costs κd
  in the squared exponent; κd_z ≤ F_act/6+κη ≤ Δ_child+….
* Hence the induction invariant is Φ := n₁+n₂+6κz−M = A−M+(6κ−1)z ≤ 0: a transform raises Φ by ≤ 6Δ_child,
  pointwise removal lowers Φ by 6 per unit of Δ spent. Terminal (old-eq:3.6): z ≤ M/(6κ).

There are two different sixes. In the errors the 6 is an artefact of ℓ_p ∈ [σ/6, σ/3] (taking ℓ* ≫ σ makes the
i = 1 error ratio σ/ℓ_p → 0, i = 6,7 ratios → 1, 6/7). In the **main amplified norm** it is intrinsic: in the
generic child (c = d = c₂ = 0, K = K₀, so J = 0) one has exactly A' = A, M' = M−6ℓ*, Δ_child = ℓ*. The ratio
is the sixth-power copy ratio: row volume Z^{6ℓ} per density gain Z^{ℓ}. This generic child binds.

## 2. The move calculus of the induction and its exact solution

Only the generic branch is modelled (it must be bounded; other branches are easier). State: represented plain
lengths (c₁,c₂) — the lengths of the plain polynomials actually transformed; children inherit them verbatim —
width M, live slot length z, allowance V (target Z^{M+V}). Moves:

* (R) reflect plain i at the current width: c_i ↦ M−c_i (eq:centered-reflection-length; free).
* (T) transform with ball extension g and amplifier ℓ ≤ g/6: admissible iff the diagonal (old-eq:2.12)
  c₁+c₂+z+g−ℓ ≤ M+V; produces the child (c₁,c₂,z; M−g; V+ℓ) (old-eq:2.13–2.14). The exceptional rows impose
  A ≤ M+V (old-eq:2.19), the same constraint. **The diagonal does not consume V** (diagonal and children are
  separate terms of one sum).
* (X) remove slot length d: z ↦ z−d, V ↦ V−κd (old-eq:3.5/3.16).
* (E) terminal: z = 0 (unrestricted zero-slot lemma), or M → 0.

Invariants: S := c₁+c₂−M, D := c₁−c₂. (T): S ↦ S+g, D fixed. (R₁): (S,D) ↦ (−D,−S); (R₂): (S,D) ↦ (D,S);
both: (S,D) ↦ (−S,−D). Diagonal: S+z+g−ℓ ≤ V.

**Theorem G (PROVED, calculus level).** Start (n,n; M₀; z₀; V₀). Up to O(g) step losses, the slot can be
carried to a terminal iff κz₀ ≤ V₀+(M₀−2n)/6. With V₀ = 0 this is exactly P78's region 2n+6κz₀ ≤ M₀.

*Proof.* (⇐) is P78. (⇒) Initially D = 0. While D = 0, the representations reachable by (R) have
S' ∈ {S, −S, 0}; S' = 0 needs V ≥ z, and S' = −S > 0 needs more. So until the first time τ_c at which
V(τ) ≥ z(τ), no reflection is usable, D stays 0, and S(τ) = 2n−M₀+τ increases by exactly the descent τ; the
diagonal forces S(τ)+z(τ) ≤ V(τ) and V(τ) ≤ V₀+τ/6−κ(z₀−z(τ)). Either z reaches 0 at some τ₀, or V ≥ z
holds at some τ_c ≤ M₀−2n (it must hold at the crossing S = 0). In both cases
V₀+(M₀−2n)/6 ≥ κz₀+(1−κ)z(τ) ≥ κz₀. ∎

Numerical confirmation (`game.py`, exact DP over (S,D,z) with step h = 0.01): for κ = 1 the DP value equals
κz−(1−2m)/6 to 4 decimals; for κ = 3/4 within the step loss h(1−1/6)·κ; ratio a in place of 6 gives
threshold (1−2m)/(aκ) on the grid.

**Asymmetric products (DP, `game_asym.py`).** For (n₁,n₂) the capacity is (M−2n_min*)/(6κ),
n_min* = min_i min(n_i, M−n_i): e.g. (0.2,0.4) → 0.133 (P78's affine region gives 0.089); (0,0.4) → 0.222
(vs 0.133). The crossing that costs is the self-dual crossing of the **shorter** factor. For the two witnesses
used (S_m·S_m and S_u(n₀)·S_u(n₀)) D = 0 and the formula is P78's. A written lemma for asymmetric inputs needs an
invariant in (S,D) instead of A; I did not write it because nothing in the application is asymmetric.

**Corollary.** The capacity (1−2m)/(6κ) (plain witness) and (1−2n₀)/(6κ) (numerator) cannot be improved by
re-ordering the induction, re-tuning (σ, ℓ*, g), delaying removals, reflecting in children, or trading
allowance on diagonals. An improvement needs a move outside {R,T,X,E}.

## 3. The routes in the brief

**(a) Slots in mean over the transformed rows.** By Theorem G all removal can be postponed to the crossing
child: width W = 2n, plain factors self-dual. The needed estimate is
Σ_{h'}|S_{W/2}(ψ_{h'})²Q_y(ψ_{h'})|² ≪ Z^{W+λy}, a fourth moment of self-dual sums twisted by a prime
polynomial *at the edge* (product length W+y > family W). Available tools give λ ≥ κ:
pointwise (old-eq:3.5): λ = κ. Large sieve (BGL, de Faveri) and the double-Poisson diagonal: the column term
Z^{W+y}, so λ = 1, plus theta terms. Hölder (Σ|S|^{4p})^{1/p}(Σ|Q|^{2p'})^{1/p'} with sup|S_{W/2}|⁴ ≤ Z^{κW}
(zero-free region) and p' = W/(2y): λ = 2κ. Large-values split E_τ = {|Q|² ≥ Z^τ}: it needs Σ_{E_τ}|S|⁴, i.e.
an 8th moment or pointwise S; with the BGL bound |E_τ| ≤ Z^{W−τW/(2y)} and pointwise S one needs τ > κy
(computed for y = 0.1, 0.2), so it loses. This is a thin-family twisted fourth moment beyond the large-sieve
edge (the "thin family" barrier of 09 §20/§22); no unconditional analogue with power-length prime twists is
known even for quadratic Dirichlet families.

**(b) Slots outside the reflection.** They already are: slots are passengers in both transforms. The binding
constraint is the second-transform **diagonal** (j = 0 term, size ≍ Z^{A} relative to Z^{M}) at every width.
Reflection shortens only plain factors.

**(c) Structure of the transformed rows.** (c1) Amplifier: H(hb) = H(h) for every column polynomial iff
χ_u(b) = 1 for all u in the column support. By reciprocity this means χ_b(·) is trivial on all good primes,
i.e. b is a sixth power up to a finite group. So the density gain is Z^{λ} per Z^{6λ} of row volume: ratio 6
is forced. Multipliers p, p², p³ twist the columns by χ_u(p)^j. The twists then differ for different p, so the
p-average cannot be merged into one norm. Carried as marks, p adds ≥ ℓ_p of conductor per side: no net gain
(this is the lem:marked mechanism, factor 2, and the energy-wall obstruction). (c2) Gauss-polynomial stage:
G(pl,h) = R(p,l)χ_p(l)²G(p,h)G(l,h) couples slot and plain variables, so there is no slot-only estimate there.
Choosing the pool inside a slot window freezes a slot only in the valuation-1 *error* (at no cost). The main
norm keeps every slot. (c3) Sixth-power copies inside a child family are subdominant (copies of length 6α
contribute Z^{m'−5α}).

**(d) Induction order.** Theorem G quantifies over all strategies of the calculus.

**(e) Hölder with the inverse moment.** Each available moment gives a count U^{E_i−thr_i}. Any Hölder product
∏(Σ|W_i|²)^{θ_i} gives the θ-weighted geometric mean ≥ min_i; this is the log-linearity already noted in 10 §3.
A gain needs a mixed moment Σ|M_rS_mQ|². Its plain/Möbius split needs the sextic sieve without its theta term
(03 R2), which is unavailable. Pointwise |S_m|² ≤ U^{δm} on bin rows amounts to one copy: it loses δm and gains
≤ 2qz ≤ δz ≤ δm.

**(f) Asymmetric witnesses.** Since conj(S_m(ψ)) = ε̄_ψ S*_{1−m}(ψ) (reflection in the same character),
|S_m|⁴ = |S_mS*_{1−m}|². This has D = 1−2m but S = 0 at width 1, so the top diagonal fails by z.
- A lift to width 1+6λ (sixth-power row amplification, cost 5λ) needs λ ≥ z/6, i.e. cost 5z/6 > κz−(1−2m)/6.
  Moving-support lifts cost ≥ z.
- Two values of t₀ in Lemma 1W: the adversary can keep m(t₀) constant (r = t₀−m), and replicate the one-t game.
- Numerator × witness pairs: (1−n₀, m) has D ≈ 0.16 but S = −0.04. Its capacity ≈ 0.2/(6κ+2) ≈ 0.031 is below
  the average 0.036 of the symmetric capacities, which is the break-even for that Hölder. No gain.

## 4. Small PROVED items (all zero payoff)

* **(P1)** Theorem G and the asymmetric capacity (M−2n_min*)/(6κ) (calculus level, DP-verified).
* **(P2) Short witnesses by the paper's own sextic large sieve.** For sixth-power-free q_u ≍ U, arbitrary
  bounded slot coefficients and 2m+z ≤ 1/2: Σ_u|S_m(u)²Q(u)|² ≪ U^{1+ε}. *Proof.* Write u = u₁w (u₁
  squarefree, w powerful, coprime) and n = n₁n₂ likewise. χ_{n₂}(u₁) factors out of the n₁-sum, and χ_n(w) is a
  fixed column twist. Lemma sextic-large-sieve (l.4707) on (u₁,n₁), Minkowski over n₂
  (Σ_{n₂ powerful} q^{−1/2} ≪ log), and summing over w give
  Σ_w(U/q_w+N+(UN/q_w)^{2/3}) ≪ U^{1+ε}+NU^{1/2}+(UN)^{2/3} ≪ U^{1+ε} for N = U^{2m+z} ≤ U^{1/2}.
  It beats (1−2m)/(6κ) only for 2m < 5/14. Binding splits have 2m ≈ 0.77–0.87 (`anatomy.py`). Payoff 0.000.
* **(P3)** κ = 2β*−1 = 35/47 at the target (Collet's variable κ): cP = cN = 47/210 gives worst +0.01166 vs
  +0.01168, a payoff of −0.00002.

## 5. What would be needed (CONDITIONAL / HEURISTIC)

A better lem:plain acts on **both** the plain witness (cP) and the numerator (cN). The brief's
"cP = 1 alone +0.0044" understates the lever (`scan_cap.py`):

| joint c = cP = cN | 2/9 | 0.3 | 0.4 | 0.5 | 0.75 | 1 |
|---|---|---|---|---|---|---|
| worst, mk = 1/2 | +0.0117 | +0.0099 | +0.0076 | +0.0052 | −0.0011 | −0.0028 |
| worst, mk = 0.75 | +0.0030 | +0.0012 | −0.0010 | −0.0029 | — | −0.0078 |

Thresholds (`quantify.py`, fine grid). Alone: c* ≈ 0.73, i.e. an exchange rate 6κ → 1.37. With mk = 0.6, 0.7,
0.75: c* ≈ 0.59, 0.41, 0.37. In exchange-rate terms, the amplifier ratio would have to fall from 6 to
a* = 1/(κc*): 1.8 alone, 3.6 with mk ≥ 0.75. That is impossible by (c1). Alternatively the progressive removal
rate would have to fall from κ = 0.75 to 1/(6c*): 0.23 alone, 0.45 with mk ≥ 0.75. That is a zero-free
half-plane for all transformed families, which is circular.

**Conditional input WR(λ)** (weighted removal / twisted fourth moment): for every lem:plain-type family of
width M, all n₁+n₂ ≤ M and slot products of length z with Θ-coefficients,
Σ_{k∈R_z}|S(n₁)S(n₂)Q|² ≪ Z^{M+λz+ε}. λ = κ is PROVED (zero-slot lemma + old-eq:3.5); λ = 0 is
Lindelöf-on-average (HEURISTIC truth). Used as a terminal (at the top, or after descending to width
2n+(1−λ)z), it gives the moment excess
E(n,z) = min{(κz−(1−2n)/6)₊, ((1+5λ)z/6−(1−2n)/6)₊, λz}. It feeds both the plain count and the numerator
(`model_wr.py`):

| λ | 0.75 (proved) | 0.5 | 0.3 | 0.2 | 0.1 | 0 |
|---|---|---|---|---|---|---|
| worst, alone | +0.0117 | +0.0106 | +0.0031 | −0.0019 | −0.0058 | −0.0078 |
| worst, mk = 0.75 | — | +0.0010 | −0.0034 | −0.0047 | −0.0078 | −0.0078 |

Thresholds (`thr_wr.py`, bisection on the coarse grid): λ* ≈ 0.238 alone; λ* ≈ 0.296 with mk = 0.6; λ* ≈ 0.40
with mk = 0.75. Compare the proved λ = κ = 0.75: a factor ≈ 3 (alone) or ≈ 1.9 (with mk = 0.75) is missing,
and mk = 0.75 is itself unproved (lem:marked gives mk = 1/2 with the 2r+8z cap).

WR(λ) with λ < κ is a twisted fourth moment of self-dual sextic sums with a power-length prime twist, at the
edge. It is of the same type as the thin-family beyond-the-edge inputs in OBSTRUCTIONS §3–4, so it is not
provable by the Poisson/reflection/large-sieve calculus.

## 6. Files
- `anatomy.py`: binding witness splits (2m ≈ 0.77–0.87, slot capacity used 0.03–0.05).
- `scan_cap.py`, `quantify.py`: cP/cN/mk scans, κ = 35/47, large-sieve option, thresholds.
- `game.py`, `game_asym.py`: DP of the move calculus (Theorem G, asymmetric capacity).
- `model_wr.py`, `thr_wr.py`: conditional WR(λ) model and thresholds.
