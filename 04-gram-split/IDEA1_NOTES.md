# Idea 1: removing the Gram P_a^{1/6} loss (working notes, 9 Oct 2026)

Sources:
- 7/8 paper: `formalize_11_12/papers/paper78.tex`. Relevant parts:
  - Prop `probe-gram`, lines ~8340–8560;
  - Prop `probe-low`, lines ~8564–8640;
  - the low separation `eq:low-separated`, line 3470;
  - the probe `eq:general-probe`;
  - the Poisson identity `eq:probe-poisson-identity`, line 3831;
  - the reflection `eq:reflection`, line 1763.
- LP model: `gram_split/theta_split_model.py`, built on `six_sevenths/routes/r4/model.py`.

## 1. What the P_a^{1/6} term is

Low separation (paper eq. low-separated):

> I = Q^{-1/2} Σ_σ ∫ Ŵ_0(iv) Σ_m Ω ξ(m) (q_m/Q)^{-iv} A_{m,σ,v}(Y) B_{m,σ}(Z) dv.

Here A_m = Y^{-3/2} Σ_s c(s) W_1(r_s) r_s^{-1+iv} g(s,-m), a sum of sextic Gauss sums of modulus s ~ Y, and B_m is the completed theta row.

The paper applies Cauchy–Schwarz in m with Σ_m |A_m|² ≪ (Q/Y)(1 + P^{1/6} + P²/Y), where P = Y²/Q. The P^{1/6} comes from the dual frequencies k = a^6 e, where the column correlation does not oscillate.

**Claim (theta component).** A_m(Y) = Y^{-1/3} Ŵ_1(1/6+iv) ρ(m) + A^rest_m(Y). Here ρ(m) is the residue at w = 2/3 of the sextic Kubota series D(w;m) = Σ_s c(s) g̃(s,-m) N(s)^{-w}, that is, a Fourier coefficient of the sextic theta function θ_6 on the 6-fold cover over Q(ω). The excess (Q/Y)P^{1/6} = Q^{5/6} Y^{-2/3} is ‖A^sig‖².

**Numerical confirmation** (`exp_rank1.py`, `summarize_rank1.py`): exact sextic Gauss sums in Z[ω], with s on disjoint dyadic ranges [Y_i, 2Y_i], Y_i = 256·2^i.
- The cross-Gram satisfies G_ij (Y_i Y_j)^{1/3} ≈ κ for i ≠ j. There are no common s, so this term is pure excess.
- The diagonal fits G_ii·Y = c·Q + κ' Y^{1/3}, with κ' ≈ κ:
  - Q = 3000: κ' = 0.258, off-diagonal κ = 0.27 (large Y);
  - Q = 400: κ' = 0.0055, off-diagonal κ = 0.050 (noisy at small Q).
- κ scales like Q^{5/6}: measured ratio 5.50, predicted 5.36.
- **So the excess is rank one in the scale variable: a single theta vector.**

## 2. Negative result: the excess cannot be cancelled by changing weights or scales

The paper's Mellin weight is 𝒲 = X^{1/2-z} Z^{x+z-1} Y^{w-1} Φ M(z) Ŵ_1(w). The principal-row signal sits at w = 1, z = 1/6, and is proportional to X^{1/3} Ŵ_1(1) M(1/6), with M(1/6) ∝ Ŵ_0(5/6).

At fixed Q = qXY, both the signal and the theta term scale as Y^{-1/3}. More precisely:
- for each m, the theta term's weight is K(μ) = ∫ W_1(t) t^{1/6} W_0(μ/t) dt/t, with Mellin transform K̂(s) = Ŵ_1(s+1/6) Ŵ_0(s);
- the signal coefficient is K̂(5/6).

So any signed combination of probes with different weights or scales (X, Y) that kills the theta term (K ≡ 0) also kills the signal. The "two-scale" trick I(X,Y,Z) − λ^{1/3} I(X/λ, λY, Z) cancels both.

**The signal is the s = 5/6 pole of Σ_m ρ(m) B_m q_m^{-s}.**

## 3. Positive route: bound the theta pairing directly, without Cauchy–Schwarz

|I| ≤ Q^{-1/2}( ‖A^rest‖·‖B‖ + Y^{-1/3}|Ŵ_1| |𝒯| ), where 𝒯 = Σ_m Ω ξ (q_m/Q)^{-iv} ρ(m) B_m.

Use the reflected form of B: B_m = Σ_{q_μ ≲ L} β_μ ψ_m(μ) V♯, with L = q_{rad m}²/Z ≈ Q²/Z < Q and β_μ = d(μ) α(μ) ϑ / √q_μ (cubic theta cusp coefficients). For squarefree m, ψ_m(μ) is the quadratic symbol (μ/m)_2 up to a fixed ray phase.

Then 𝒯 = Σ_μ β_μ F(μ), with F(μ) = Σ_m ρ(m) (μ/m)_2 Ω.

**Heuristic (V).** By Voronoi (automorphy of θ_6 at the cusps d/μ), F(μ) has dual length q_μ²/Q. The d-sum combines the quadratic twist with the sextic Kubota multiplier into a cubic Gauss sum, so

> F(μ) ≈ ε(μ) (Q/q_μ) Q^{-1/12} Σ_{n ≤ q_μ²/Q} a'(n) (n/μ)_3^{±1}, with |ε| = 1 and a' the unitary-normalised theta coefficients at another cusp,

plus main terms only for μ ∈ cubes, which total O(Q^{5/6}).

Then Heath-Brown's cubic large sieve gives

> |𝒯| / (‖ρ‖·‖B‖) ≲ [L/Q + L²/Q² + L²/Q^{5/3}]^{1/2} = Z^{-s(M')},  where s(M') = min((1-M')/2, 1 - 7M'/6).

**LP model** (`theta_split_model.py`): Gram exponent g = max(0, b/6 − 2s(M'), 2b − l_y + d).

| model | θ | geometry (l_x, l_y, ℓ) |
|---|---|---|
| paper (no split) | 0.874957 | — |
| split, cubic large sieve | **0.864017** | (0.339, 0.506, 0.155) |
| split, ideal ("random") large sieve | 0.859796 | — |
| split, ideal large sieve, P²/Y term also dropped | 0.851918 | — |

**Numerics for (V)** (`exp_twist.py`, `exp_voronoi.py`), measuring coherence of twisted theta sums against random size:
- quadratic twist (j = 3): ratio ≈ 0.3–1.4 at all μ sizes, so it is incoherent;
- sextic twist (j = 1): mildly coherent, ratio 3–9 at moderate N(μ)²/Q, decaying to 1;
- cubic twist (j = 2): ratio ≈ 1.

These results are consistent with the heuristic but not conclusive at small μ (noise).

## 4. What remains to prove (for an unconditional θ ≈ 0.864)

- **(V)** A Voronoi summation formula for Σ_m ρ(m) χ_μ(m) W(q_m/Q), with χ_μ the quadratic character mod μ, uniform in μ, for the finitely many twisted sextic theta functions. The dual must have the cubic-character structure, with dual coefficients satisfying Σ_{q_n≤N} |a'(n)|² ≪ N^{1+ε} and main terms only for cube μ.
- **(D)** Σ_m Ω |A_m − A^sig_m|² ≪ (Q/Y)(1 + P²/Y) Z^ε. In words: the exceptional-frequency main term of the Gram equals ‖A^sig‖² to relative accuracy P^{-1/6}. This is a Rankin–Selberg unfolding of a truncated Eisenstein series against θ_6.
- **(P)** Compatibility with Part II, covering:
  - the marked/compensated rows B^J (slots);
  - non-squarefree rows m, where the twist is ψ_m(μ) = conj χ_{m·rad(m)²}(μ);
  - the v-integral;
  - the restriction to the ray class s ∈ σ.

## 5. Updates (later the same day)

**Marked rows (P).** In Part II's low side, the B-row is the marked row B^J. Its marked slot primes are active in the reflection, so the reflected length is L_J = Q²·Z^{ℓ'−1}, where ℓ' = ℓ − d is the marked length; this follows from paper eq. old-eq:4.4, T_d − H = M'+ℓ'−1. The saving is then

> s = −½·max(M'+ℓ'−1, 7M'/3 + 2ℓ' − 2, −M'/6).

| route | θ |
|---|---|
| model `LSmarked` | **0.870719** |
| model `randommarked` (ideal large sieve) | 0.863101 |

**Exact certificate** (`exact_cert_theta.py`, the r4 checker with a theta-split Gram):
- geometry (l_x, l_y, ℓ) = (2671/10000, 2671/5000, 107/2500), σ0 = 871/1000;
- result: ALL CHECKS PASS (40623 high-side boxes); low side = 52243/60000;
- **so (V)+(D)+(P) ⇒ ζ(s) ≠ 0 for Re s > 0.871.**

**(D) via Lindelöf-on-average.** Shift the w-contour of A_m(Y) = Y^{-1}(2πi)^{-1} ∫ Ŵ_v(w) Y^w D(w;m) dw past the theta pole at w = 2/3 to Re w = 1/2+ε. Then

> A^rest_m = Y^{-1} ∫_{(1/2+ε)} Ŵ_v Y^w D dw,  so ‖A^rest‖² ≪ Y^{-1+2ε} ∫ |Ŵ_v| Σ_m |D(1/2+ε+it; m)|² dt.

Hence (D) ⟸ (LoA6): Σ_{N(m)≤Q} |D(1/2+ε+it; m)|² ≪ Q^{1+ε}(1+|t|)^R.

This is the n=6 analogue of [DdFDS24 = David–de Faveri–Dunn–Stucky, arXiv 2410.03048, Lemma 6.6] (cubic), which is proved by an approximate functional equation plus a large sieve. For n=6 the sextic large sieve is sufficient: BGL gives (AB)^{2/3} = Q at A = Q, B = Q^{1/2}, and de Faveri's Thm 1.1 gives Ξ_6(Q, Q^{1/2}) ≪ Q^{1+ε}. The functional equation of Kubota series for general n is in Brubaker–Bump, "On Kubota's Dirichlet series", Crelle 2006.

No GRH is needed: the zeros of ζ(6s−5), with s = w + 1/2, do not meet 1/2+ε ≤ Re w ≤ 1.

Bonus: this bounds all of A^rest by (Q/Y)Z^ε, so the paper's P²/Y Gram term also disappears. With it dropped, model `LSmarked` gives **0.870307**.

**(V) reformulated.** (V) is equivalent to a functional equation for quadratic twists of the sextic theta L-series L(z, θ_6 ⊗ χ_μ), with:
- conductor ≍ N(μ)²;
- the dual twist being the cubic character (·/μ)_3^{±1}, by the Kubota multiplier computation; numerics agree (j=3 incoherent, j=1 coherent);
- polynomial growth.

It can be derived from the residue at w = 2/3 of the twisted A2 double Dirichlet series Σ_s g̃(s) L(z, χ̄_s χ_μ) N(s)^{-w} (Brubaker–Bump–Friedberg–Hoffstein type), but uniformity in μ is not in the literature.

Remark: centre-point Lindelöf-on-average plus CS over μ is WORSE. The dual side of the approximate functional equation carries cubic characters, and the cubic large sieve costs L^{4/3}. The Voronoi form, with dual length q_μ²/Q < q_μ, is what wins.

## 6. Check: row counts at small or zero slot supply

**Concern.** Prop detector-counts assumes prime supply > 7/37. Its zero-slot remark gives only t = 1 ⇒ U^{1−2δ/3}. With that weaker R = 1 − 2δ/3 the slot-free split would give 0.8751–0.8779, no gain.

**Resolution.** The witness proposition (Two saturated witnesses, line 4510) holds for every t ∈ [1, 3/2] with no slots. In the count proof, every case has a slot-free bound:
- r ≥ 1: no-slot inverse amplification, eq. no-slot-inverse-count, gives ≤ L(t);
- m ≥ 1/2: zero-slot plain gives ≤ 1−δ;
- r < 1 and m < 1/2: the unweighted (zero-capacity) estimates 1−δr (no-slot-inverse-count for r ≤ 1) and 1−2δm (zero-slot plain, m ≤ 1/2).

Hence R_0(t) = max(L(t), 1 − 2δt/3), minimised over t. This is exactly the R4 model at S = 0; for example δ = 0.3 gives 0.7727, against 0.8 at t = 1. With positive supply below 7/37 the capacities are min(requested, supply), which is the R4 model's rule. The 7/37 in the proposition only guarantees full capacity.

**Validity of the uncapped formula.** The split's g-piece should be min(b/6, b/6+φ). θ_true = min(LP with b/6, LP with b/6+φ) by distributivity of max over min. At both certified geometries φ < 0, so the certificates are valid as stated:
- slot-free: φ = −0.073;
- slotted: φ = −0.045 at d = 0 and −0.119 at d = ℓ.

## 7. Reducing gap G3 (theta mean squares at other cusps) to paper-type inputs

Claim: Σ_{N(n)≍N} |ρ_𝔞(n)|² ≪ N^{5/6+ε} at every cusp 𝔞. Proof sketch:
- Let A^𝔞(Y) be the additive polynomial built from the Kubota Gauss sums at 𝔞 (the same Gauss sums, with shifted classes).
- By (LoA6 at 𝔞), Y^{-1/3}|Ŵ| ‖ρ_𝔞‖ = ‖A^{sig,𝔞}‖ ≤ ‖A^𝔞‖ + ‖A^{rest,𝔞}‖.
- The Prop. probe-gram argument (Poisson in m) applies verbatim at 𝔞 and gives ‖A^𝔞‖² ≪ (N/Y)(1 + P^{1/6} + P²/Y), with P = Y²/N.
- Take Y = N^{0.6}. Then:
  - Y^{2/3}(N/Y)P^{1/6} = N^{5/6};
  - Y^{2/3}(N/Y) = N·Y^{-1/3} ≤ N^{5/6};
  - Y^{2/3}(N/Y)P²/Y = Y^{8/3}/N ≤ N^{5/6}, since Y ≤ N^{11/16}.

So G3 ⇐ (Gram at 𝔞) + (LoA6 at 𝔞), with no Rankin–Selberg regularisation needed. G4, the non-squarefree dual indices, remains. de Faveri's Θ_3 already allows cube-free columns, so only the cube-full parts n_1³ need a growth bound on ρ_𝔞(n_0 n_1³). A similar Gram argument at level n_1³ is the natural route.
