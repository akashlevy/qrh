# energy-wall: the marked reflected energy beyond M + ℓ = 1

Agent `energy-wall`, 10 Oct 2026. Verdict: **no proof of W(w) for any w > 0. Instead there is a
sharp quantitative obstruction (§2), a mechanism-by-mechanism no-go (§3), and numerics (§4)
showing W is true numerically far past the wall.** Files in this directory:
`dual_energy.py` (+ `dual_energy_N40000.json`), `exchange_rate.py` (+ `.out`),
`model_wall.py` (a copy of r4 `model.py` with a slope knob `wall_factor`), `mw_*.json` and `model_w*.json`.

## 0. The object

At the binding dyad of Lemma `reflected-energy` (O = 0, ℓ_b = 0, v = T_d), with everything in log_Z:

    E = Σ_{k sf, q_k ≤ K} | Σ_{n sf ~ N} β(n) X(n) ψ_k(n) |²,   ψ_k = χ_k^3 (quadratic),
    X(n) = Σ_P ã(P) q_P^{-1/2} χ_P(n)^{-2} = Π_slots X_i(n),
    K = Z^M,  N = Z^{2M+ℓ-1},  L = q_P ≍ Z^ℓ,  ‖βX‖² ≍ 1.

Lemma `hybrid-energy` gives E ≪ (K + N) Z^ε. The wall excess is w := M + ℓ − 1 = log_Z(N/K).
Note that ã(P) contains the cubic Gauss-sum phases τ⁺_{p,2} of the marked primes (P78 l.1823,
ω_{p,0}). So X_i is **not** a Hecke-character prime sum: τ⁺_{p,2}χ_p(n)^{-2} is the shifted Gauss
sum g(χ_p², n)/√p. Equivalently, β(n)X(n) = Σ_P a(P)·(theta coefficient at Pn), up to normalisation.

## 1. What W(w) is worth (model)

| hypothesis | θ (r4 model) | geometry (l_x, l_y, ℓ) |
|---|---|---|
| wE = 0 (paper) | 0.874957 | |
| wE = 0.06 | 0.871206 | |
| wE = 0.064 | **0.870955 < 27/31** | (0.3662, 0.4936, 0.2042) |
| wE = 0.07 | 0.870580 | (0.3673, 0.4950, 0.2077), M = 0.862 |
| wE = 0.10 | 0.868705 | |

Empirically θ(wE) = 0.874957 − wE/16 on [0, 0.1]. So 27/31 needs **wE ≥ 0.0638**.

## 2. Sharp quantitative obstruction: only an *offset* helps, never a slope

Near the wall the shadow price of the wall constraint is 1/16, since θ(wE) = θ₀ − wE/16. Cauchy–Schwarz
charges half the energy excess on the low side. Consider an energy bound whose excess grows at rate f
beyond the wall:

    E ≪ Z^{M + f·(M+ℓ−1)+ε}.

Locally such a bound pays only if f/2 < 1/16, that is f < 1/8. Globally the LP can jump to a far
vertex (ℓ ≈ 0.29, M + ℓ ≈ 1.19). Checked with `model_wall.py` (`wall_factor` = f, wE = 0):

| f | θ | optimum (l_x, l_y, ℓ) |
|---|---|---|
| 0.50, 0.25, 0.22, 0.20, 0.18, 0.16 | 0.874957 (no gain) | paper vertex |
| 0.14 | 0.873822 | (0.384, 0.522, 0.295) |
| 0.12 | 0.872232 | (0.384, 0.522, 0.291) |
| 0.105 | 0.871070 | (0.385, 0.522, 0.289) |
| 0.10 | 0.870689 < 27/31 | (0.385, 0.522, 0.288) |
| 0 (wall piece deleted) | 0.863636 = 19/22 | (0.400, 0.509, 0.273) |

So a slope-type bound K^{1−f}N^f is worthless unless f ≲ 0.15, and 27/31 needs f ≲ 0.103.
Every transfer or interpolation mechanism in §3 gives f ≥ 1/2, including the hypothetical
decoupled one. What is needed is an **offset**: full diagonality E ≪ (K + N^{1−η})Z^ε with
η ≥ w/(2M+ℓ−1) ≈ 0.068 at the target geometry. That is Lindelöf on average for a family of
twisted cubic-theta sums of length N = K·Z^{0.064}. This size of saving is a GRH-type
(Lindelöf-on-average) statement, not a technical loss.

## 3. Mechanisms tried, and exactly where each fails

**(O1) Phase-blind barrier (PROVED).** Fix the true marks X and keep only |β(n)| = N^{−1/2}.
Then β′(n) = N^{−1/2}·conj(X(n)ψ_{k0}(n))/|X(n)| has the same moduli, and

    |S(k0)|² = (Σ|β′X|)² ≥ N (E|X|)² ≥ N (E|X|²)³/E|X|⁴ ≫ N Z^{−ε},

using bounded 4th moments of X, which follow from the cubic large sieve since L² ≤ N. So any
argument that majorises |β| in the column variable n, or uses the marks only through moments or
large values, cannot beat the column term. This kills Hölder with high moments of the slot
polynomials, and "moments of different lengths" applied to X. Cancellation must come from the
phases of the theta coefficients β against Xψ_k.

**(O2) Exchange-rate law (`exchange_rate.py`; proved as exponent bookkeeping, even when the
optimal hybrid large sieve is assumed for every enlarged family).** The only access to β's phases is:
- reflection R, which is an involution: dual excess w, primal excess ℓ − w at best;
- twisted multiplicativity (TM) β(n₁n₂) = β(n₁)β(n₂)(n₁/n₂)₃, which splits off a factor;
- Cauchy–Schwarz of a coherent variable of log-length v into the rows.

Excess after transfer, at M = 0.862, ℓ = 0.208, w = 0.070:
- a mark slot moved to the rows: max(0, w − v) + v ≥ w. This is break-even, because the columns do not shrink.
- a TM factor n₂ of the dual index: excess ≥ w/2 + ℓ = 0.243 (0.139 with an optimistic sparse sieve). Every prime of n is coupled to the marks through χ_P(n₂)^{−2} = (P²/n₂)₃, so the column index becomes n₁P² and absorbs the marks.
- a factor decoupled from both the marks and ψ_k: excess w/2, which is still a slope of 1/2 and useless by §2. The only such variable, the cube part b of nb³ (χ_P(b³)^{−2} = 1), exists only on dyads that are already cheaper (ℓ_b > 0).
- the primal (unreflected) side: excess max(1 + ℓ − M, (2 + 2ℓ − M)/3) = 0.518 with the BGL sextic sieve, and 1 − M = 0.138 even if its sparsity were usable.

All transfers satisfy excess ≥ w as soon as ℓ ≥ w, which is equivalent to M ≤ 1 and always holds here.

**(O3) Rough part.** TM needs n to have a divisor in a fixed middle range. The complement, for
example n prime or n = (large prime)·(small), has density ≫ 1, not a power of Z. Bounding it needs
Type-I information: twisted theta sums at conductor ≫ length, where even convexity is worse than
trivial (N < K√L, the effective AFE point). A partial positive fact, also proved: for each single
row k₀, the prime-weighted part Σ_n β(n)X(n)ψ_{k₀}(n)Λ_I(n), with Λ_I(n) = Σ_{p|n, p∈I} log p and
I = [N^{0.3}, N^{0.7}], is ≪ (L N^{0.7})^{1/2} Z^ε. The proof writes n = pm, applies TM,
Cauchy–Schwarz over P (cost L) and the cubic large sieve in (p, m). At the optimum L N^{0.7} ≤ K
(0.860 ≤ 0.862). This excludes "one row carries the column term" on that part. It does nothing against
the generic bad configuration, in which K rows each have |S|² ≈ N/K.

**(O4) The a-priori zero-free half-plane β* ≤ 0.8749 does not help.**
- On the low side it is useless: the dual marks carry Gauss sums τ⁺_{p,2}, so they are not Hecke prime sums, and (O1) applies anyway.
- On the high side the mark amplitude P^{1/2+xδ} is capped only when (1+δ)/2 > β*, that is δ > 3/4. The binding δ is about 0.39, so there is no gain.

**(O5) Mark redesigns.**
- A prime mark (j = 0 zero mask) costs conductor p² for primal length p, a net factor ×p on the dual.
- A cube-part mark 1_{p|b} is the κ = 3 rescaling. It gives B_m ↦ (m/p)₂·B_m at scale Z, which is the m-twist dead end (NOTES §3, §5).
- Marks on s (the A_m side) multiply the Gram exceptional term P_a^{1/6}, which is the signal's shadow, by the same P^{1/6} as the signal.

## 4. Numerics (`dual_energy.py`, N₀ = 40000, 14127 columns, 14110 rows)

Per-row energy E(K)/(#rows·‖c w‖²), dual side, quadratic rows, at K/N₀ = 1/256 … 1:

| coefficient | range | wall (provable up to) |
|---|---|---|
| θ (unmarked) | 0.69–0.94 | everywhere |
| θ·X, slot [31, 62) | 0.90–1.17 | K ≤ N₀/31 |
| θ·X, slots 7·61 | 0.70–0.99 | K ≤ N₀/427 |
| θ·X, slots 7·19·43 | 0.75–0.99 | K ≤ N₀/5719 (whole range beyond) |
| θ·χ_P^{−2}, P = 7, 31, 97, 211 | 0.74–1.05 | K ≤ N₀/P |
| random phases·\|X\| | 0.76–0.95 | |
| adversarial \|θX\|ψ_{k0} | 1.81 at K/N₀ = 1/2; max\|S\|²/diag = 6270 | |

The maximum of |S(k)|²/diag is 8–11 over 14110 rows, which is the size of a maximum of
exponentials. So W(w) holds numerically up to N/K ≈ 5700, with Gaussian row statistics. An
excess would need the coefficients βX to correlate at level (N/K)^{1/2} with a positive proportion
of the quadratic characters; nothing of the kind is visible.

## 5. Status

- PROVED: (O1); the bookkeeping in (O2); the partial single-row bound in (O3); the model thresholds (wE ≥ 0.0638; slope f ≲ 0.15 for any gain, ≲ 0.103 for 27/31).
- HEURISTIC/NUMERICAL: W(w) is true.
- OPEN (and needed): a Lindelöf-on-average offset for quadratic twists of the marked cubic theta beyond the diagonal, E ≪ (K + N^{1−η})Z^ε with η ≈ 0.07. Within {R, P, TM, cubic/quadratic/hybrid large sieves, Cauchy–Schwarz, Hölder} this is unreachable by (O1)–(O3).
- A route would need a genuinely spectral input, such as amplification plus Kuznetsov on the cubic cover (Proskurin/Louvel), giving individual-twist savings for theta, uniformly in the level K²L. Even that must be upgraded to an average statement, because of the generic K-row configuration.
