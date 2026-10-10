# 11 — The energy wall program (started 10 Oct 2026, evening)

## Why this lever
In the weighted-numerator architecture (our 10499/12000, now the best framework-verified tracker
entry), a **lossless** relaxation of the energy wall by w ≈ 0.041 (≈ 0.036 if the same saving also
feeds lem:marked) closes θ < 41/47 by itself (`10-below-41-47/NOTES.md` §15–16). Among all open
inputs it is the least hopeless: it collides with no famous barrier (unlike zero density, which
would beat ζ, or the Gram pairing, which is signal-circular), it is numerically true with margin
≈ Z^{0.43}, and the *unmarked* family already has a power-saving asymptotic past the analogous
point (Proposition U). It is still an open beyond-the-edge problem. This workstream attacks it in
three steps:
1. **Numerics at the weighted geometry** (confirm the margin where we need it).
2. **Simplest case first: one marked prime** with any lossless w > 0 — an extension of Prop. U
   from the unmarked family to one level-P² twist. If that fails, the full statement fails too.
3. **Function-field analogue**, where affine moments of this type have power-saving secondary
   terms; look for a mechanism.

Labels: **[P]** proved (given P78's lemmas), **[B]** exponent bookkeeping, **[C]** conditional,
**[N]** numerical, **[H]** heuristic, **[L]** literature.

## Sources (read the relevant parts first)
- `10-below-41-47/agents/wall4/NOTES.md` — **the reduction of W(w) to one core sum** (§1–2), why the
  Lemma-D split fails (§3: single-mark energies are self-dual and each has the full wall),
  Proposition T (rank-one wall W_θ), model thresholds (§8).
- `09-new-directions/agents/edge-theta/EDGE_THETA.md` — truth numerics (random-matrix law),
  **Proposition U** (§3.2), Lemma M(w) (§5), function-field literature pointers (§2).
- `09-new-directions/agents/energy-wall/ENERGY_WALL.md` — phase-blind barrier (O1), exchange-rate
  law (O2), rough part (O3), numerics (§4).
- `09-new-directions/agents/wall-spectral/WALL_SPECTRAL.md` — identity (∗)
  β(n)ψ_k(n)χ_P(n)^{−2} = μ(n)G(n)𝔯(n,k)·conj γ₁(n; k³P⁴); (H1)/(H2); conditional 19/22; FE groups.
- `09-new-directions/agents/{h2,hls,mark-design}/` — (H2), the hybrid large sieve (HLS), mark numerics.
- P78: `/home/user/openai/math/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex`:
  lem:reflected-energy (l.7657ff), lem:hybrid-energy (l.7510), marked completion/reflection
  (l.7249–7485), prop:probe-low (l.8061–8650).
- Numerics library: `09-new-directions/agents/edge-theta/edge_lib.py` (exact sextic symbols in Z[ω],
  β = γ₂ᾱ, G, Möbius partner), `primal_energy.py`, `dual_mobius.py`, `pairs_mobius.py`;
  `09-new-directions/agents/cB-numerics/cubic_data.c` (fast exact cubic Gauss sums to norm 10⁷).
- Exponent models: `10-below-41-47/agents/wall4/wm4.py`, `10-below-41-47/agents/combo/combo_model.py`.

## The objects (from wall4 §1–2; P78 normalisations)
Reflected marked energy at the binding dyad:

  𝒩 = Σ_{k sf, q_k ≍ K} | Σ_P c_P 1_{(P,k)=1} Σ_{n sf ~ N} β(n) ψ_k(n) χ_P(n)^{−2} |²,

with K = Z^M (rows), L = q_P = Z^ℓ (marks), N = K²L/Z = K·Z^w (reflected column length),
w = M + ℓ − 1 (wall excess), c_P = a(P)q_P^{−1/2} (|a| ≤ 1, coherent: amplitude matters), β = γ₂ᾱ the
cubic-theta coefficient on squarefree n, ψ_k = χ_k³ quadratic. lem:hybrid-energy: 𝒩 ≪ (K + N)N Z^ε.
**W(w):** 𝒩 ≪ K N Z^ε for N ≤ K Z^w. Only the "core" dyad matters (wall4 §2): all slots active,
no dual cube part, rows in the top Z^{w/2}, powerful row part O < 2w, dual length within Z^{2w/3} of T_d.
Weighted geometry: ℓ ≈ 0.18–0.20, b ≈ 0.06–0.15, M ≈ 0.82–0.86, w ≈ 0.04; note L = Z^{0.19} ≫ N/K = Z^{0.04}.

**Single-mark problem SM(w):** for one prime P (q_P = L),

  E_P = Σ_{k sf ~ K} | Σ_{n sf ~ N} β(n) ψ_k(n) χ_P(n)^{−2} |² ≪ K N Z^ε  for  N ≤ K Z^w,

i.e. quadratic twists of the single twisted theta θ ⊗ χ_P^{−2} (level ~P²). Each twist has conductor
(q_k q_P)² and AFE length KL > N, so S_{k,P} is shorter than its AFE length; the k-family is beyond
the edge by L. wall4 §3.1(c): E_P is self-dual under Poisson + μ-absorption with the χ_P⁴-twisted
Möbius family Σ_{h ≍ Z²/Q}|Σ_n μ(n)χ_P(n)⁴(h/n)₆|², so every Poisson/reflection/large-sieve method stops
at the wall. A proof needs a genuinely new input (spectral, amplification, Hecke structure, or a
family-level functional equation).

**Proposition U** (edge-theta §3.2, [C on P78 reflected/hybrid energy]): for the unmarked
Möbius-sextic family and Z ≤ H ≤ Z^{3/2}, Σ_{h~H}|M_h(Z)|² = Diag + O(Z^{2+ε}) (relative error Z^{1+ε}/H).

## Working rules
- Write only under `/home/user/qrh/11-energy-wall/<your-name>/`. Do not commit or push.
- Be rigorous and skeptical; label every claim. Do not rediscover the known obstructions above —
  build on them.
- Final message ≤ 600 words: what you did, what you proved (precise statements), numbers, and the
  precise obstruction if it failed; point to your files.
