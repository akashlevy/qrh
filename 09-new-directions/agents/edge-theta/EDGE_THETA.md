# edge-theta: truth, provability, and the affine-edge moments behind W(w) and (H2)

Agent `edge-theta`, 10 Oct 2026.

**Verdict.**
1. **The power saving is TRUE, and it is large.** Measured with exact numerics on both sides of the Poisson duality, the relative off-diagonal has random-matrix size (family)^{-1/2}. The exponent is 0.50 ± 0.08 on the primal side and 0.53 ± 0.02 on the Möbius dual side. There is no dependence on the family/length ratio F and no secondary term for P78's theta β = γ₂ᾱ. The needs are much smaller: W needs only O(Z^ε) on the primal side, i.e. 1/F on the dual side. At the 27/31 geometry the margin is about Z^{0.43}. (H2) holds with saving K^{-1/2}, against a need of L^{-0.93}, a margin of about Z^{0.18}. **No truth obstruction exists:** W(w) is not "log-only true".
2. **No route to a proof.** Our extra inputs cannot replace the step that is missing:
   - β* and GRH are pointwise inputs, and pointwise inputs reproduce exactly the hybrid large sieve.
   - The theta structure is the reflection R, and R stops at the wall.
   - The marks can be chosen freely, but every choice behaves randomly.
3. **The sharpest obstruction is not "affine edge ⇒ logs".** The unmarked Möbius-sextic family *does* have a power-saving asymptotic, thanks to R. W sits deep inside the diagonal regime of its Eulerian dual (F_m ≈ L²Z^{-w}). For Eulerian GL2-type coefficients that regime is unconditional, because Perron reduces it to at-edge L-value moments. Here the only substitute for that Perron step is theta reflection, and reflection ends at M + ℓ = 1.

Labels: **[P]** proved (given the P78 identities and NOTES §6), **[B]** proved as exponent bookkeeping, **[C]** conditional, **[N]** numerical, **[H]** heuristic, **[L]** literature (abstracts and search snippets only; arXiv is blocked).

## Files

| file | content |
|---|---|
| `edge_lib.py` | lookup-table sextic symbols, β = γ₂ᾱ by twisted multiplicativity, G(n) and the Möbius partner a = μḠ = conj(γ₁β), checked on primes |
| `test_lib.py` | checks for `edge_lib.py` |
| `primal_energy.py` → `primal_N{10,20,40,80}000.{out,json}` | reflected (quadratic-row) energy with exact n-diagonal, 12 angular row sectors, cross-mark Gram for (H2) |
| `dual_mobius.py` → `dual_unmarked_N*.{out,json}`, `dual_marked_N*_L*.{out,json}` | Möbius-sextic family D(H) in the diagonal regime, unmarked and marked (exact coincidence diagonal) |
| `pairs_mobius.py` → `pairs_N*.{out,json}` | per-pair twisted moments T(P1,P2): a direct test of Lemma M (§5) |
| `bias_check.py` → `.out`, `.json` | Patterson bias of Σβ, Σγ₂, Σγ₂α, Σμ |
| `analyze.py` → `analysis.txt`, `fits.json` | all exponent fits |
| `model_et.py` (copy of wall-spectral `model_ws.py`) → `model_checks.out` | model checks |

## 0. Objects and dictionary

**Primal.** E(K) = Σ_{k~K} |S_k|², with S_k = Σ_{n~N} c(n) w ψ_k(n).
- k and n are squarefree, odd, primary; ψ_k = χ_k³ is the quadratic row character.
- The coefficients are c = β, β·X (marks X = Σ_P c_P χ_P^{-2}), or random phases.
- ρ = E/Diag − 1, with the exact diagonal d_k = Σ_n |cw|² 1_{(n,k)=1}.
- W(w) means ρ ≪ Z^ε for F = N/K ≤ Z^w.

**Dual (NOTES §6, §20).** D(H) = Σ_{h~H} |Π_h|² |M_h|², where
- M_h = Σ_{n~Z} a(n) w (h/n)_6;
- Π_h = Σ_P c_P (h/P)_6 is the mark polynomial;
- the Möbius partner is a = μ conj G. It satisfies a·γ₁ = conj β (this follows from γ₁γ₂ = μαG, verified at primes in `test_lib.py`), so Poisson in h returns exactly the theta energy.

Write F_m = H/(ZL), and δ = D/Diag − 1, where Diag includes the coincidences n₁P₁ = n₂P₂. In this language:
- W(w) ⟺ δ ≪ Z^ε/F_m on F_m ∈ [L²Z^{-w}, L²], i.e. H ∈ [Z^{1−w}L³, ZL³]. **[P]**, NOTES §6.
- The wall M + ℓ = 1 corresponds to H = ZL³, i.e. F_m = L². So W lives deep in the diagonal regime of an *Eulerian* family, not at its edge.

## 1. Truth check [N]

### 1.1 Primal (`analysis.txt`, 30 windows, N₀ = 10⁴…8·10⁴, F = 1.2…312)

Fit: log ρ_rms = a + b·log(rows per sector) + c·log F.

| coefficient | b | c |
|---|---|---|
| β | −0.497 ± 0.078 | +0.013 ± 0.077 |
| β·X, coherent slots | −0.567 ± 0.072 | −0.071 |
| β·X, random slots | −0.505 ± 0.071 | +0.02 |
| random control | −0.493 ± 0.027 | −0.001 |

The ratio of theta-type to random sector rms is 0.94 on average (60 windows).

**Cross-mark correlation (H2).** For the normalised Gram entries c_{P₁P₂} = G_{P₁P₂}/√(G_{P₁P₁}G_{P₂P₂}):
- theta: rms = e^{0.357}·rows^{−0.500±0.012}·F^{−0.03};
- random control: e^{0.301}·rows^{−0.494}.

### 1.2 Dual Möbius-sextic family

Unmarked runs: N₀ = 1000, 2000, 4000 with F ≤ 256. Marked runs: L₀ = 13 and 31 with F_m ≤ 128, entirely beyond the wall.

| family | δ_rms ∝ (#h)^b (theta partner) | random control | theta/random | sector δ / allowed (median, max) |
|---|---|---|---|---|
| unmarked | −0.529 ± 0.019 | −0.502 ± 0.007 | 0.93 | 0.02, 0.13 |
| marked | −0.533 ± 0.018 | −0.533 ± 0.006 | 0.94 | 0.03, 0.18 |

After the Möbius family's own variable #h is accounted for, the residual F_m-exponent is −0.05 ± 0.07 (unmarked) and −0.13 ± 0.06 (marked). Both are in the favourable direction; the random control gives −0.06 ± 0.02.

### 1.3 Per-pair twisted moments

For P₁ ≠ P₂ (`pairs_*.out`):
- |T(P₁,P₂) − C(P₁,P₂)| = (0.27–0.45)·#cols·√#h, which equals the random control.
- C is the explicit coincidence main term.

### 1.4 Patterson bias (`bias_check.out`)

- Σ_{n~N} β(n) = O(√count): **β = γ₂ᾱ is bias-free.**
- Σ γ₂ ≈ 0.37·count^{5/6}: the Patterson pole is present.

So the only secondary main term is the one wall-spectral found: the Gauss-sum marks produce the (NL)^{2/3} term. It is harmless at 27/31 and caps W at 19/22.

### 1.5 Extrapolation to the binding geometries

| lever (geometry) | needed | truth (random-matrix law) | margin |
|---|---|---|---|
| W(0.064), (l_x, l_y, ℓ) = (0.366, 0.494, 0.204), M = 0.860 | ρ ≤ Z^ε | ρ ≈ K^{-1/2} = Z^{-0.43} | Z^{0.43} |
| same lever, dual form: H = Z^{1.549}, F_m = Z^{0.345} | δ ≤ Z^{-0.345} | H^{-1/2} = Z^{-0.774} | Z^{0.43} |
| (H2), (0.385, 0.522, 0.288), M = 0.907 | L^{-0.93} = Z^{-0.268} | 1.4·K^{-1/2} = Z^{-0.453} | Z^{0.185} |

If the truth held, the model would give 19/22 = 0.863636 (`model_checks.out`, Hoff = 0, bias kept).

**Caveats.** The numerics only reach rows ≤ 9·10³, so they cannot exclude log factors or secondary terms with tiny constants. W tolerates Z^ε, so neither would matter.

## 2. Function-field comparison [L]

- **Predicted shape.** For the affine D4^(1) case (4th moment of quadratic Dirichlet L-functions), Diaconu–Goldfeld–Hoffstein predict infinitely many lower-order terms of size X^{(j+1)/(2j)+ε}, j ≥ 2, accumulating at X^{1/2} (as quoted in Florea's paper). **So power savings are true in the affine class,** but with a tower of secondary main terms. Our analogue of those terms is the Patterson-bias (NL)^{2/3} term.
- **Fixed q, analytic method.** Florea gets only 3 of the 11 terms of the degree-10 polynomial: a saving of a power of g = log_q X. This happens even though RH (Weil) is a theorem there.
- **What does give power savings.**
  - Diaconu–Paşol–Popa (arXiv 2111.11062 and a sequel): an exact formula for a weighted 4th moment over F_q(t), with infinitely many terms. They continue the five-variable affine series to its optimal domain via "a new type of functional equation".
  - Bergström–Diaconu–Petersen–Westerland (2302.07664): all moments for q large, using stable homology of braid groups with symplectic coefficients plus the stability range of Miller–Patzt–Petersen–Randal-Williams.
  - Whitehead and Sawin–Whitehead: affine Weyl group series characterised axiomatically, continued to the maximal (Tits-cone) domain.
- **Lesson.** Power savings at an affine edge come from *family-level* inputs: rationality plus affine functional equations, or homological stability. Pointwise RH for the members does not supply them.
- **Our needs in this language.** In wall-spectral's coordinates, W needs I = ½ for diagonal mark pairs and I ≤ 0.263 for cross pairs. Both lie strictly inside the Tits cone I > 0, which the function-field continuation reaches. In number fields, Bochner hulls are stuck at inf I of the starting region (wall-spectral, [P]).

## 3. Why SY needs GRH, why Li gets logs, and why neither step is ours

### 3.1 Anatomy [L]

Soundararajan–Young (JEMS 2010), on GRH:
1. Remove μ²(d) by Σ_{a²|d} μ(a).
2. Poisson in d. The non-square dual frequencies k ≠ □ produce shifted second moments of the *same* family at the edge, Σ_k L(½+u, f⊗χ_k) L(½+v, f⊗χ_k).
3. Bound these shifted moments, and the large-a tail, by Soundararajan's GRH method (log|L| ≤ a short prime sum), which is sharp up to (log)^ε.
4. The main term is ≍ X log X. The error from this method is X(log X)^{1/2+ε}: Li's paper attributes this rate to SY's method, and states that its own unconditional error is of the same quality, improvable to X(log X)^ε.

Li (Invent. 2024) makes step 3 unconditional with an "optimal" large sieve for the special coefficients λ_f(n)/√n, using prime-square inflation (Forti–Viola, Heath-Brown). Shen–Stucky do the same for the 4th moment, with a large sieve adapted to fourth moments and cross-cancellation between diagonal and off-diagonal terms.

**Why only logs [H, structural].** At the edge, Poisson maps the off-diagonal to a moment of the same family, again at the edge. This is the fixed point of the involution, |ε−δ| = 1−κ = 0. It can only be *bounded*, and sharp-up-to-log bounds give an asymptotic with log savings. A power saving would require evaluating that dual moment asymptotically, which is the same problem again.

### 3.2 Our problem is not at that fixed point [P, given P78 and §6]

W is an *upper bound* beyond the edge by F = Z^w. Its Poisson partner is a *diagonal-regime asymptotic* at F_m ≈ L²Z^{-w}.

- **For GL2 quadratic twists,** the beyond-edge Lindelöf-on-average bound is unconditional: Perron reduces partial sums of any length to L(½+it, f⊗χ_k), and Heath-Brown's at-edge second moment bounds those. Neither GRH nor SY is needed.
- **The theta analogue of that Perron step is reflection R,** and it works exactly up to the wall.
- **Proposition U [C on P78 Lemma reflected/hybrid-energy and §6].** For the *unmarked* Möbius-sextic family, Z ≤ H ≤ Z^{3/2} (M ∈ [½, 1]),
  Σ_{h~H} |M_h(Z)|² = Diag + O(Z^{2+ε}),
  a power-saving asymptotic with relative error Z^{1+ε}/H. For generic coefficients the sextic large sieve gives *no* saving in this range, since its (HZ)^{2/3}Z term exceeds Diag·Z/H. **So theta structure does beat "affine ⇒ logs" for the unmarked family.** Beyond Z^{3/2}, P78's frozen powerful-row term (1−2M)/4 enters.
- **What breaks with marks.** Each mark is a twist of conductor P, which costs level P² per length P. I re-checked the half-integral-weight intuition [H]: U_P should map level N to level NP with a cubic nebentypus. Even so, the reflected length is Q²L/Z either way (conductor Q²L at length Z, or (QL)² at length ZL), so the U_P picture is a reparametrisation (NOTES §3).

### 3.3 Our three extra inputs, with exponents

**(i) β* ≤ 0.875 and GRH [B].**
- On the Eulerian (Möbius) side, GRH gives |M_h| ≪ Z^{1/2+ε}. That yields δ ≪ Z^ε, which under Poisson is exactly the hybrid large sieve, so wE = 0.
- β* gives δ ≪ Z^{2β*−1+ε} = Z^{0.75}, which is worse than the unconditional bound.
- On the theta side there is no Euler product, so β* has nothing to act on.
- W ⟺ cancellation in Σ_h (|Π_h|²|M_h|² − E), a joint statement over h that no pointwise bound sees.
- The bootstrap map β* ↦ wE is therefore identically 0. Even the GL2 analogue with marks, Σ_k L(f⊗χ_kψ₁) conj L(f⊗χ_kψ₂) with saving L^{-1} for distinct cubic ψ_i mod P_i, is open on GRH: Soundararajan's method controls |L|-moments, not signed cross-correlations.

**(ii) Theta structure [P/N].**
- Twisted multiplicativity and Kazhdan–Patterson periodicity generate R and the identity (∗). Cube marks amplify nothing (NOTES §9).
- The Kubota residue only produces secondary *main* terms: the (NL)^{2/3} term, Σγ₂ ~ count^{5/6}.
- For β = γ₂ᾱ the twisted series have no pole (§1.4).
- R yields Proposition U and nothing past M + ℓ = 1, by the move calculus (OBSTRUCTIONS).

**(iii) Freedom to choose the marks [B/N].**
- In the dual, marks become family weights: |Π_h|² = Σ c_{P₁} conj(c_{P₂}) ψ_{P₁P₂}(h). So W ⟸ per-pair bounds T(P₁,P₂) − C ≪ Z^{2+ε} (§5).
- Poisson plus the large sieve gives per pair only a relative error max(1/F_m, 1/L), against a need of 1/(F_m L). R is needed even to reach the wall.
- Numerically every T(P₁,P₂) is random-sized for coherent, random and Kubota-aligned (conj γ₁) slots. No choice of c_P creates cancellation beyond √(#pairs), and coherence is forced by the signal.

## 4. Sharpest obstruction (what is proved)

- **O1 [P] Exact location.** W(w) ⟺ the Möbius-sextic family with prime-slot family weights has relative error ≤ Z^ε/F_m for F_m ∈ [L²Z^{-w}, L²]. The wall is F_m = L², deep in the diagonal regime, not at the affine edge F_m = 1.
- **O2 [B] Pointwise inputs give exactly the large sieve.** On the Eulerian side, GRH gives δ ≪ Z^ε and β* gives Z^{0.75}. Both have zero payoff for every β* > ½.
- **O3 [N] No truth obstruction.** The off-diagonal is (#family)^{-1/2} on both sides of the duality and per pair, with W's need exceeded by Z^{0.43} and (H2)'s by Z^{0.18}.
- **O4 [L/H] What a proof would need.** A family-level input placing Z(s₁,s₂,w) at I = ½ and I ≤ 0.263 inside the Tits cone. In function fields that input is rationality or homological stability; it has no known number-field substitute. Equivalently, a beyond-the-wall analogue of "Perron + at-edge L-values" for the marked theta.

## 5. Weakest sufficient statement (open) and payoff

**Lemma M(w).** Let a be the Poisson partner of P78's completed theta row, and let ψ_P(h) be the mark character, (h/P)₃^{±1} in P78's normalisation (sextic in my numerics). Then for primes P₁ ≠ P₂ ~ L and Z^{1−w}L³ ≤ H ≤ ZL³:

  Σ_{h~H} ψ_{P₁}(h) conj ψ_{P₂}(h) |Σ_{n~Z} a(n)(h/n)₆|² = C(P₁,P₂) + O(Z^{2+ε}),

where C is the coincidence term over n₁P₁ = n₂P₂.

- **What it implies.** Together with the P₁ = P₂ asymptotic (Proposition U with a coprimality condition), it implies W(w), since L'²Z² ≤ L'·Z²L.
- **Status.** It is an *Eulerian* statement: Möbius sums of Hecke characters, twisted by a character of conductor L² in the family variable. In the random model its error is ≈ Z·H^{1/2} ≤ Z^{1.78}, a margin Z^{0.225} below the Z² target. Numerically the error equals the random control (§1.3). OPEN; it is equivalent in difficulty to (H2).

**Payoff** (`model_checks.out`, bias kept):

| hypothesis | θ |
|---|---|
| W(0.0638) | 0.8709680 (27/31 threshold ≈ 0.06381) |
| W(0.064) | 0.8709554 |
| DIAG + W(0.0196) | 27/31 (diag agent) |
| W for all w, i.e. (H1) + (H2) full | 19/22 = 0.8636364 |
