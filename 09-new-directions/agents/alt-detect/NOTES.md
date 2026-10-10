# alt-detect — breaking links 1–3 of the structural chain (10 Oct 2026)

Agent `alt-detect`. Labels: **PROVED** (finite algebra / standard), **R-** (rigorous modulo standard
theory not written out), **HEUR**, **MODEL** (exponent bookkeeping; LP exact for the stated pieces),
**NUM** (numerics). Files are listed in sec 8.

## 0. Verdict

No candidate survives with rigorous inputs. One structure is genuinely new, the **cubic family with a
double-Gauss-sum B-row**. It is the only known way to make the low-side dual *full-index* (no residual
twist). Its value depends on one unresolved number: the conductor exponent c_B in the row h of the cubic
dual B_h. All the evidence I have points to c_B ≥ 3, which kills it.

| candidate | residual twist after reflection | conductor exponent c_B | theta (P2 shape) | theta (P78 shape, MODEL) | status |
|---|---|---|---|---|---|
| sextic → cubic theta (P78/P2) | quadratic (**forced**, sec 1) | 2 | 11/12 | 7/8 (GRH-rows floor 13/15) | known |
| (a) cubic → theta(x)theta / sextic Kubota index h^4 | **none (full-index)** | 4 rigorous; heuristics 3; 2 = optimistic | 1 / 23/24 / 37/42 | 10/11 / 10/11 / **0.797** | open only if c_B = 2 (sec 3) |
| (a) quartic, quintic, sextic double rows | — | ≥ 3 or no known object | ≥ 1 | — | dead (sec 2) |
| (b) Lambda / log L / L'/L | no automorphic dual | — | ≥ 1 (cubic), 9/8 (sextic) | — | dead (sec 4) |
| (c) |theta|^2, theta^3, Shintani zeta | absorption impossible | — | — | — | dead (sec 5) |
| (d) other orders and fields | quadratic **forced** for every single-Gauss-sum row | — | — | — | PROVED (sec 1) |

## 1. Jacobi triples: the quadratic residual is forced (PROVED)

**The triple identity.** Let chi = chi_n have order N on squarefree primary n, and let gamma_j = gamma(chi^j). By Stickelberger/Weil, J(chi^a, chi^b) is −(Hecke character) when a, b, a+b ≢ 0. Hence

  **mu(n) · Hecke(n) = gamma_a gamma_b gamma_c   whenever a + b + c ≡ 0 (mod N), with a, b, c ≢ 0.**

**How a probe uses it.** In any probe the Poisson step produces one factor gamma_a, where chi^a is the twist in the probe variable. The B-row must supply the other two factors, gamma_b gamma_c.

**Single-Gauss-sum B-rows.** A degree-2 theta has a single Gauss sum, and by KP this exists only on GL_2^{(3)}.
- The product gamma_b gamma_a is mu·Hecke only if gamma_c is a Hecke/ray quantity.
- If c ≡ 0, the product is gamma_b gamma_{-b} = ±1, which detects L, not 1/L.
- So c ≡ N/2, and gamma_c is the quadratic Gauss sum, which is a ray quantity.
- The twist on the row is then chi^a = chi^{−b}·chi^{N/2}. The part chi^{−b} collapses into an index shift. **The residual chi^{N/2} has order exactly 2.**

This holds for every N and every field containing mu_N, and also for KP thetas on GL_r^{(k)}, whose reflected degree is r ≥ 3. So link 3 cannot be broken by changing the order or the field while keeping a single-Gauss-sum B-row. This answers (d).

**Collapse-only twists cannot absorb mu.** A pure collapse twist gives c_B = 1. That is the only way to get a finite-type energy MDS: the Cartan matrix [[2,0,−1],[0,2,−1],[−c_B,−c_B,2]] has det = 8 − 4c_B. But pure collapse is exactly the case c ≡ 0 above. **Finite type and mu-absorption are therefore incompatible** for single-Gauss-sum rows.

## 2. Classification of double-Gauss-sum B-rows (PROVED algebra; HEUR on objects)

These are the B-rows with gamma_b gamma_c and b, c ≢ N/2.

| N | triples | B-row object |
|---|---|---|
| 6 | (1,1,4) | DH turns it into the P78 row; nothing new |
| 6 | (2,5,5) | mu × theta, or the RS square of the sextic theta, which has no Gauss-sum Whittaker model; dead |
| 4 | (2,3,3) | B-row = mu·Hecke·ray, i.e. Mobius itself; dead |
| 5 | any | no automorphic object with Gauss-sum coefficients (KP thetas only at degree ≥ 4); dead |
| **3** | (1,1,1) | **gamma_1^2 = Hecke·mu·conj(gamma_1) = coefficient of theta_3⊗theta_3** |

For N = 3 the row coefficient gamma_1^2 has three realizations:
- the Rankin–Selberg square theta_3⊗theta_3, via ∫theta^2 E_kappa;
- the cubic A_2 WMDS residue (the cross-check in `04-gram-split/wf2` §3.1);
- by DH, a sextic Gauss sum: mu·g_3(psi) = Hecke·g_6(phi) with phi^2 = psi-bar.

The twist is fully absorbed (index shifts in both theta factors), so the low-side dual is again the cubic family, *full-index* (Kubota index h, not h^3). Copies are u = a^3, so the copy exponent is c = 1/3.

## 3. The N = 3 architecture in detail

**3.1 Which dual actually occurs (PROVED).** Absorption pairs the Gauss sum of the twist character with that same twist. The Poisson dual is therefore always B_h = Σ mu nu(n) g~(psi_n) psi-bar_n(h), with numerator h.
- In theta language, **both** factors carry the collapse twist: B_h ∝ Σ α(n) conj(G~(n h))^2.
- In sextic language, this is the coprime part of the twisted sextic Kubota series at **index h^4**. Kubota's rank-one theory over Q(omega) is rigorous here since mu_6 ⊂ F, and it gives an FE with conductor N(h)^4.
- The conductor-N(h)^2 object (sextic index h^2, twist on only one factor) has the opposite orientation. It is *not* a Poisson dual, in either P2 or P78 shape.

**3.2 Conductor exponent c_B (open).** The Kubota bound is 4. Heuristics:
- "minimal sextic index with the same coprime coefficients" gives 4;
- a local-type model calibrated on P2's table (j=4 collapse → 1, j=2 → 2) gives 3;
- the cubic A_2 residue gives 1 after the KP collapse. But its non-squarefree completion contains Σ_p g(1,p)·(·), a prime sum of Gauss sums, so B is not that residue divided by an L-function.

NUM: `conductor_test2.py` uses controls of known conductor. At the reachable sizes, u = log X/log N(m) ≤ 5, it is inconclusive: the predicted decay beyond the conductor is only (C/X)^{1/6}.

**3.3 Normaliser trap and bypass (R-).** KP periodicity gives T(Ab^3)^2 = N(b)T(A)^2.
- Every small-conductor completion (RS, A_2) therefore equals B·L(3s−1,(alpha nu)^3). Inverting this costs D^{(1+β')/3} per row, coherently in h.
- With this cost: P2 bootstrap fixed point 13/14; P78 with r4 counts gives θ(7/8) = 7/8 exactly, i.e. no gain (`runs/cubic_beta_*.json`).
- Bypass: B is holomorphic in Re s > 1/2 (sextic Kubota representation), so contour shifts cross no zeros of L(3s−1). After reflection the dual is again squarefree with harmless N(d)^{-2}-type cube removal. This makes the trap irrelevant, but only if c_B ≤ 2.

**3.4 Bookkeeping (MODEL; `p2_orbit.py`, `p78_cubic_lp.py`, `model_cubic.py`).**
- **P2 shape.** I computed the orbit of Poisson and reflection moves, with the cubic large sieve as base. The moves preserve |family − length|.

  | c_B | 1 | 2 | 3 | 4 |
  |---|---|---|---|---|
  | θ (P2 shape) | 5/6 | 37/42 | 23/24 | 1 |

- **P78 shape, GRH for rows.** This includes the slot Patterson-bias piece (M'+4l'−2)/3. Without it the cubic template runs away to huge geometry; the sextic does not.

  | c_B | 2 | 3 | 4 |
  |---|---|---|---|
  | GRH-rows floor | 13/17 = 0.765 | 0.8596 | 10/11 |

- **P78 shape with r4 counts, transplanted.**
  - Copy term 1/6 → 1/3. Amplification slope (2/3)·h_min, where h_min is from the P2 orbit.
  - Results: c_B = 2 → **0.797**; c_B = 3 or 4 → **10/11**.
  - Sensitivity at c_B = 2: removing the prime-gain lemmas gives 0.798; h_min = 1.25 gives 0.819; h_min = 1.5 gives 0.872.

**Precise obstruction.** The chain breaks at link 3 (full-index dual), but link 2 comes back in a new form. The absorbed twist lands on both theta factors, so the dual's conductor exponent is ≥ 3 (heuristic) and ≤ 4 (proved). The FE group of the energy MDS is then hyperbolic (det < 0), not merely affine.

**The single sufficient new input:** a conductor-N(h)^2 functional equation (with polynomial growth) for h ↦ Σ_{A sqfree,(A,h)=1} mu(A) g~(A) (h/A)_3^{-1} N(A)^{-s}. Even then, the P78 value needs cubic analogues of the count lemmas (HEUR).

## 4. (b) Lambda, log L, L'/L

- Lambda is not twisted-multiplicative, so Lambda·gamma is not an automorphic coefficient.
- Vaughan or Heath-Brown decompositions give two kinds of term:
  - Type I terms contain mu, which brings back absorption;
  - Type II bilinear terms in the order-N family have only the coefficient-blind BGL bound F + L + (FL)^{2/3}.
- In P2 shape with no reflection this gives h ≥ 3/2, hence θ = 1 (cubic) and 9/8 (sextic) (`p2_orbit.py`).
- The first-moment paper's dual.tex gains only Y^{-η}, with tiny ineffective η, at two boundary configurations. NOTES §14 shows that savings vanishing at the edge need η ≳ 0.85 to matter.
- A power-saving Patterson asymptotic would not detect zeros: Σ g(p) N(p)^{-s} has no Euler product and no explicit formula.

## 5. (c) |theta|^2, theta^3, prehomogeneous zeta functions

- Absorption needs Gauss sums whose **modulus is the theta index**. So the probe must sum characters mod A over their argument, which gives one character per modulus.
- Full families sum over characters to a fixed modulus. Orthogonality then gives congruences A ≡ A' (mod q), i.e. shifted convolutions (Fourier coefficients of |theta|^2, an ordinary function), but never gamma(A).
- So **fullness at the probe level and mu-absorption are mutually exclusive**. Fullness can only appear after reflection, and sec 1 and 2 cover that case.
- theta^3: its coefficients are additive ternary convolutions, so there is no Gauss sum at the index.
- Shintani/prehomogeneous zeta functions have nonnegative orbit-count coefficients. Their cubic content enters only through h_3(d), which is again one character per quadratic field. They are L-type, never 1/L-type.

## 6. Other observations

- **Coefficient-blind regime.** In the cubic template, the trivial "Poisson node" (energy diagonal once ZP ≤ Q^{2/3}) drives θ_GRH to 27/40 at huge geometry. This happens only if the slot-bias piece is dropped. The sextic does not move.
- **A curiosity.** In the cubic P78 model with the normaliser trap, the bootstrap fixed point is 7/8 to 2·10^-7.
- **Checks (NUM, `check_gauss.py`).**
  - G~(pi)^3 = −pi/|pi| exactly. This is mu-absorption for cubic twists.
  - Twisted multiplicativity matches a direct sum.

## 7. What would change the picture

- **Settle c_B.** One route is to compute the p-part polynomials of the coprime twisted sextic Kubota series at index h^4: does the excess N(h)^2 drop out, as in the n = 2 example of the A_2 referee? Another is to run the conductor test at X ≫ 10^6, which needs a fast cubic Gauss-sum routine.
- If c_B = 2 can be proved, the cubic architecture is the strongest lead found this session.

## 8. Files

| file | content |
|---|---|
| `p2_orbit.py` | P2-shape orbit bookkeeping (Poisson/reflection moves + order-k LS); A_2 bootstrap |
| `p78_cubic_lp.py` | P78 template with pluggable energy nodes (validated: sextic → 13/15) |
| `p78_compare.py`, `p78_cubic_bootstrap.py` | Poisson-node comparisons; normaliser bootstrap (GRH rows) |
| `model_cubic.py` (+ `r4_model_copy.py`) | heuristic cubic transplant of r4 counts (count function verified identical to r4 in sextic mode) |
| `r4_poisson_node.py`, `cubic_bigbox.py` | geometry-box / Poisson-node validity checks |
| `conductor_test.py`, `conductor_test2.py`, `check_gauss.py` | Z[omega] Gauss sums, identities, conductor scaling test |
| `runs/*.json` | all model outputs quoted above |
