# R3: the zero-density reformulation of M_kappa (verdict: OBSTRUCTION)

*Working notes for route R3, 8 Oct 2026.*

**Scripts** (run with python3 / venv):
- `r3_verify.py`: exact `fractions`; every assertion passes.
- `r3_coeffblind.py`: numerics with exact sextic symbols in Z[omega].

**Paper references.**
- **P** = 7/8 paper `paper.tex` (Sep 30).
- **P2** = companion `paper2.tex` (Oct 5).
- Line numbers refer to those files.

---------------------------------------------------------------------------------------------------

## 0. Verdict in one paragraph

The reformulation is exact in one direction and has **no slack** in the other.

**What a multiplicity-beating density estimate is.** A density estimate for the sextic family {L(s, nu psi_u) : N(u) <= H} that beats the H^{1/6} multiplicity of the target (A(1-sigma) < 1/6) is not a density statement at all. Up to epsilon, it is *equivalent* to the zero-free half-plane Re s > sigma for every member of the family (Theorem A).

**Why the standard tools cannot reach it.** No large-values inequality that depends only on the size of the coefficients can give it (Theorem B). This covers the large sieve (Heath-Brown cubic, BGL sextic, and even an optimal sextic sieve), Halasz–Montgomery, Huxley subdivision and Guth–Maynard-type estimates, with or without the "duplicates are exact copies" structure.

**Where the target's zero is detected.** At the target copies u = a^6, saturated plain witnesses can only have length U^{o(1)}: they are smooth sums of the fixed target character, which are O(1) beyond its analytic conductor. A zero of the target is therefore detected only by an inverse (Mobius) witness of length >= U^{1-o(1)}, where U is the row norm. So the density problem at the copies *is* the sub-diagonal Mobius large-values problem, that is, M_kappa itself (Proposition C).
- The only Mobius-specific mean values available are diagonal-limited (P2 Prop thm:ms; P eq. raw-moment). They give exactly A = 2, i.e. theta = 11/12.
- Adding P's prime slots (Lemma marked) to this simple extraction makes things worse: theta(z) = 11/12 + z/3.

**Inside the Part II probe.** Here the copies are the signal and the density input concerns the non-principal rows. Two facts:
- The **floor bin** (rows whose zeros all lie near Re s = 1/2, essentially all rows) cannot be removed by any density estimate. Together with the low-side lemmas it forces theta >= 13/15 for every geometry (exact LP certificate).
- The density hypothesis for rows (A = 2), with detector floor -> 1/2, already *attains* 13/15 at l_x = l_y = 2/5, ell = 1/5 (model level).

So zero-density inputs stronger than DH (A < 2, up to GRH) are worthless in Part II. In particular, **no zero-density estimate of any strength yields 6/7 in either existing architecture.**

**Best provable exponents near sigma in [6/7, 7/8].**
- *Distinct sixth-power-free rows*, from Part II's Prop. detector-counts with no prime information:
  - R(delta) = 1 - delta + delta(5-6delta)/(3(5-2delta)), with delta = 2sigma-1;
  - **A(6/7) = 7/3, A(7/8) = 16/7**;
  - DH (A = 2) is not proven in this range.
  - For comparison, BGL Cor. 1.6 gives A ≈ 4.16 and 4.31.
- *Multiplicity-counted (the sense in which A = 7/6 <=> 6/7):*
  - DH strength (A = 2) at the copies, from the diagonal Mobius mean square, giving 11/12;
  - A = 0 for sigma > 7/8, from the 7/8 theorem itself (Part II's "effective" A = 4/3 is just a restatement of its theta).
- **Resulting theta: no improvement on 34999/40000 = 7/8 - 1/40000.**

---------------------------------------------------------------------------------------------------

## 1. Setting and dictionary

**Notation.** Fix a finite-order Hecke character nu over F = Q(sqrt(-3)). Then
- psi_u(n) = nu(n) chi_n(u), where chi_n(u) = (u/n)_6 with zero extension;
- A_u(D) = sum_n mu(n) nu(n) chi_n(u) W(N n / D).

**The target copies.** For primary a with (a, S) = 1,
- psi_{a^6}(n) = nu(n) 1_{(n,a)=1}.

So L(s, nu psi_{a^6}) = L(s, nu) times prod_{p | a} (1 - nu(p) N p^{-s}), and the finite product is nonzero on Re s > 0. There are ≍ H^{1/6} such copies with N(u) <= H. Distinct primary a give distinct a^6, because the primary generator is unique.

**Dictionary.** A family density estimate N_mult(sigma; H) << H^{A(1-sigma)+eps}, counting rows with multiplicity, contradicts a target zero at beta as soon as A(1-beta) < 1/6. So

>  theta = 1 - 1/(6 A_eff)  (exact, `r3_verify.py` §1):

| theta | A_eff |
|---|---|
| 11/12 | 2 (DH) |
| 7/8 | 4/3 |
| 34999/40000 | 20000/15003 ≈ 1.33307 |
| 13/15 | 5/4 |
| **6/7** | **7/6** |

**M_kappa as large values (rigorous, trivial).** M_kappa says sum_{N(u)<=H}|A_u(D)|^2 << D^{1+kappa+eps} with H = D^kappa. It gives

  #{u : N(u) <= H, |A_u(D)| >= D^sigma} << D^{1+kappa-2sigma+eps} = H^{A_kappa(sigma)(1-sigma)},

with A_kappa(sigma) = (1+kappa-2sigma)/(kappa(1-sigma)).
- A_{6/7}(6/7) = 7/6. This is the number quoted for the route.
- A_1 ≡ 2 (DH).
- Restricted to the copies, it gives beta <= (6+5kappa)/12.
- A single row gives |A_u(D)|^2 <= D^{1+kappa}, so every row L-function is zero-free in Re s > (1+kappa)/2 (= 13/14 at kappa = 6/7).

---------------------------------------------------------------------------------------------------

## 2. Theorem A: the multiplicity-beating density statement has no slack

**Hypothesis LV(sigma_1, kappa, eta).** For all large D, with H = D^kappa,

  #{u in O : 0 < N(u) <= H, max_{D^{1-kappa/6} <= D' <= D} |A_u(D')| >= D^{sigma_1}} <= H^{1/6-eta}.

**Theorem A.**
- (i) LV(sigma_1, kappa, eta) implies that for every fixed u_0, A_{u_0}(D) << D^{sigma_1}. Hence L(s, nu psi_{u_0}) has no zero in Re s > sigma_1.
- (ii) Conversely, if every L(s, nu psi_u) is zero-free in Re s > sigma_1 - eps', then LV holds with count 0 for large D.
  - This uses the standard bound 1/L << (q(|t|+2))^{eps} inside a zero-free region, e.g. P Lemma logarithmic-control.

**Proof of (i).** The exact column identity P eq. (sixth-power-column-identity) (line 12426) with a = p prime, p ∤ u_0 S, gives directly from the definition

  A_{u_0}(D) = A_{u_0 p^6}(D) - psi_{u_0}(p) A_{u_0 p^6}(D / N p).

This is exact: split n by p | n, use mu(p) = -1, and use psi_{u_0 p^6}(m) = psi_{u_0}(m) 1_{p∤m}.

Suppose |A_{u_0}(D)| >= 2 D^{sigma_1}. Then for every such prime p with N(p) <= (H/N(u_0))^{1/6}, at least one of A_{u_0p^6}(D) and A_{u_0p^6}(D/Np) has modulus >= D^{sigma_1}. Here D/Np >= D^{1-kappa/6}.
- The rows u_0 p^6 are distinct and satisfy N(u_0 p^6) <= H.
- There are >> (H/N u_0)^{1/6} / log H > H^{1/6-eta} of them for large D.

This contradicts LV. So |A_{u_0}(D)| < 2D^{sigma_1} for all large D. Mellin inversion over a smooth dyadic partition then makes 1/L(s, nu psi_{u_0}) holomorphic on Re s > sigma_1. ∎

LV is to be read uniformly over the profiles W(x) x^{-s} with s in compact sets. This is the same rowwise-parameter Sobolev device as in P Lemma inverse-amplification, so it costs nothing.

**Consequence.** The statement "N_mult(sigma;H) << H^{1/6-eta}" needed by R3 is, up to epsilons, *the theorem itself* for the family. It is not a weaker density statement that could be proved first.
- M_kappa is strictly stronger: it also controls moderate values.
- Its interest is only that it is a *family mean value* that might be provable by duality (Poisson) without singling out the target.

---------------------------------------------------------------------------------------------------

## 3. Theorem B: coefficient-blind large-values estimates cannot beat the multiplicity

**Definition.** Call B(H, D, V) a coefficient-blind large-values bound for the family if, for every coefficient vector with |c_n| <= 1 supported on primary n with N n in (D, 2D] and (n, S) = 1, and every fixed weight w,

  #{u : 0 < N(u) <= H, |sum_n c_n chi_n(u) w(N n/D)| >= V} <= B(H, D, V).

The weight w may be complex, for example the detector profile x^{-sigma-i gamma}.

This class covers every bound that depends on the coefficients only through their sizes:
- Montgomery/Huxley large values;
- the Heath-Brown and BGL large sieves (P Lemma sextic-large-sieve, line 4707);
- an optimal sextic sieve;
- Halasz–Montgomery with any off-diagonal input, including exact vanishing of smooth sums beyond the conductor;
- Huxley subdivision;
- Guth–Maynard-type bounds (their Thm 1.1 is for |b_n| <= 1).

**Theorem B.** For H, D large, V <= c_1 D / log log D, and H^{1/6} <= D^{O(1)},

  B(H, D, V) >= c_0 H^{1/6}.

**Proof.** Take c_n = conj(w(N n/D)) / |w(N n/D)| on the support, and 0 elsewhere. For u = a^6 with a primary, (a, S) = 1 and N(a)^6 <= H, we have chi_n(a^6) = 1_{(n,a)=1}. Hence

  sum_n c_n chi_n(a^6) w = sum_{(n,a)=1} |w(N n/D)| = c_w D phi(a)/N(a) + O(2^{omega(a)} D^{1/2+eps}),

by Mobius inversion over d | a and the lattice-point count in each residue class. Since phi(a)/N(a) >> 1/log log N(a), all ≍ H^{1/6} copies are >= V. ∎

**Corollaries.**
1. If the target L(s, nu) has a zero at beta > 1/2, the detector (P Prop. detector-witness, line 4510) makes every copy row carry a witness of size D^{beta - eps}, which is <= V. A coefficient-blind bound can contradict this only if B < #copies, and Theorem B forbids that. So **coefficient-blind tools give no information at all at the target**, not even the trivial-looking 11/12.
2. Restricting to sixth-power-free rows (as Part II and BGL Cor. 1.6 do) removes the copies, so such bounds say nothing about the target either. The target can then enter only as a *signal* (the Part II probe; see §5).
3. For distinct rows, the coefficient-blind optimum is still worse than DH.
   - de Faveri–Dunn–Hoffstein (arXiv 2607.07911, Thm 1.1) prove that the (QN)^{2/3} term in the cubic large sieve is unavoidable.
   - With their conjectural optimal sextic sieve, the induced density exponent is A(6/7) ≈ 3.17 and A(7/8) ≈ 3.33 (`r3_verify.py` §3).

**Numerical sanity check** (`r3_coeffblind.py`, exact sextic symbols, w = sin^2). This is the per-row energy divided by the diagonal, split into copies (u = a^6) and the rest.

| D | kappa | rows | #copies | mu: total / copies | c = 1: total / copies | random ±1: total / copies |
|---|---|---|---|---|---|---|
| 3000 | 0.6 | 438 | 1 | 1.072 / 0.000 | 2.045 / 1.284 | 1.067 / 0.000 |
| 3000 | 6/7 | 3462 | 2 | 0.914 / 0.000 | 1.236 / 0.325 | 0.974 / 0.000 |
| 10^4 | 0.6 | 912 | 1 | 0.963 / 0.000 | 2.735 / 2.059 | 0.934 / 0.000 |
| 10^4 | 6/7 | 9720 | 2 | 0.936 / 0.000 | 1.270 / 0.386 | 0.940 / 0.000 |

- For c = 1, each copy carries ≈ 0.020·D, compared with ≈ 0.10 per row on the diagonal. So the copy part grows like D^{1 - 5kappa/6}: 1.28 -> 2.06 at kappa = 0.6, exactly as predicted, and 0.33 -> 0.39 at kappa = 6/7. Once there are ≍ D^{1/7} copies, M_{6/7} must fail for c = 1.
- For mu and random signs the copies are negligible.
- M_kappa is therefore a property of the *Mobius* coefficients and cannot follow from any coefficient-blind inequality.
- Caveat: this is a sanity check of a trivial construction. It is not evidence that M_{6/7} is provable.

---------------------------------------------------------------------------------------------------

## 4. Proposition C: the copies are detected only by long inverse witnesses

The detector (P Lemma buffered-bins, Lemma detector-dyads, Prop. detector-witness, lines 4281–4687) produces a product M_r S_m with:
- r + m >= t, t ∈ [1, 3/2];
- |M_r|^2 >= U^{delta r - eps} and |S_m|^2 >= U^{delta m - eps} (both saturated).

**Plain witnesses at a copy u = a^6.** Write U = N(u) = N(a)^6, let nu be nonprincipal of fixed conductor q_nu, and let N = U^m. The profile is W_S(y) = W_2(y) y^{-sigma - i(gamma-nu')}, with height |gamma - nu'| << T_1 = Z^tau <= U^{1/100}. Mobius inversion over d | a gives

  sum_{(n,a)=1} nu(n) W_S(N n/N) = sum_{d | a} mu(d) nu(d) sum_k nu(k) W_S(N(dk)/N).

- If N/N(d) >= (q_nu T_1)^{1+eps}, the inner sum is a smooth sum of a nonprincipal character beyond its analytic conductor, so it is O_A(N^{-A}).
- Otherwise it has O((q_nu T_1)^{1+eps}) terms of size O(1).

Hence the unnormalised sum is O(tau(a) (q_nu T_1)^{1+eps}). After the normalisation N^{-1/2},

  |S_m|^2 << U^{-m + 2tau/d + eps}.

Saturation requires |S_m|^2 >= U^{delta m - eps}. Together these force m(1+delta) <= 2tau/d + 2eps, so m = O(tau/d + eps), which is arbitrarily small. Then r >= t - O(tau/d + eps) >= 1 - o(1).

If nu is principal (target zeta_F), the plain sums at the copies are dominated by the pole's main term and carry no zero information. P's detector removes principal presentations explicitly (lines 4258–4267) and uses that L_orig is entire (line 4567). Either way the copies are witnessed only by inverse polynomials.

So at a copy the saturated plain witness has m = o(1), and the inverse witness has r = t - o(1) >= 1 - o(1). In particular, **the zero of the target is detected at the copies only by Mobius sums whose length is at least the row norm**: rows U, columns U^t with t >= 1.
- At r <= 1/(1+c) < 1 this is the diagonal regime of P2 Prop. thm:ms (line 679) and P eq. raw-moment (line 12353). The count there is U^{1 - delta r}, which beats U^{1/6} iff delta r > 5/6. Taking t -> 1 gives **exactly 11/12**.
- r > 1 is the sub-diagonal regime. The count needed there is precisely the large-values form of M_{1/r}, with kappa = 1/r < 1.

So R3 at the target reduces *exactly* to M_kappa (large-values form). The zero-density language adds no leverage.

**Lemma marked as an amplified simple extraction** (`r3_verify.py` §4). P Lemma marked (line 9250) sums over *all* u in O. Its sum therefore contains the copies, and slots with a_i(p) = conj(nu(p)) are large at the copies.
- Restricting it to the copies gives |M_1(Z^r)|^2 <= Z^{5m/6 - z}.
- The width conditions are r + 2z <= m - c_1 and 2r + 8z <= 3m - c_2.
- With r = 1 and m = max(1+2z, (2+8z)/3), this gives theta(z) = 1/2 + (5m/6 - z)/2 = 11/12 + z/3 for z <= 1/2. The minimum is 11/12 at z = 0.
- Slots cost 2z of row length but amplify by only Z^z. Part II gets below 11/12 only through the probe's signal normalisation (the zeta(6z) residue and C(s)), not through counting.

---------------------------------------------------------------------------------------------------

## 5. Theorem D: zero-density inputs saturate at DH and stop at 13/15 inside Part II

**High side.** The general-geometry high exponent is P eq. (common-high-exponent) (line 15726, Lemma high-bin). At d = h:

  sigma_0 >= (1+delta)(1-l_y)/2 - ell/2 + q ell + h(R + delta/2 - 1/6),   with h = 1 - l_x + ell.

This was checked against the paper's first line: a - sigma_0 + h(17/50-1/6) - a l_y - (1-a)ell - (delta/2-q)ell + d(R+delta/2-17/50).

**Low side.** From P Prop. probe-low (line 8564; exponent l_x/2 + Gram/2 + energy/2) against C(s) = l_x/2 + s - 1 + h/6 (eq. 5.11, line 8644):

  sigma_0 >= 5/6 + l_x/6 - ell/6 + g/2 + e/2,

with g >= max(0, b/6) (Gram P_a^{1/6}) and e >= max(0, M+ell-1) (energy wall).

**Floor bin.** The floor bin a = 51/100 (or 1/2 + eta after the audited floor lowering) needs no zero witness. Its count is R = 1, and this count cannot be reduced by any density estimate:
- every row has zeros with real part >= 1/2 within the height window;
- under GRH every row is in the floor bin.

Its per-row bound is the pointwise one. Only cancellation *across rows* could improve it, and that is not a density input.

**Theorem D (exact LP certificate, `r3_verify.py` §5a; also `barriers.py` part B).** Take the combination with weights:
- 4/5 · low side;
- 1/5 · (g >= 0);
- 1/5 · (g >= b/6);
- 1/3 · (e >= 0);
- 1/15 · (e >= M+ell-1);
- 1/5 · floor-bin high side at delta_0 = q = 0.

Every geometric variable cancels and the combination yields

  **sigma_0 >= 13/15 for every geometry (l_x, l_y, ell) and every row-count input.**

For delta_0 > 0 and q >= 0, the floor row equals the delta_0 = 0 row plus delta_0[(1-l_y)/2 + h/2] + q ell >= 0. So the bound persists.

**Theorem D' (DH saturates; model level, `r3_verify.py` §5b).** Assume DH for rows, R(delta) = 1 - delta for every delta in (0, 5/6] and every amplitude x, with floor -> 1/2. Take l_x = l_y = 2/5, ell = 1/5, h = 4/5. Then:
- the low side is exactly 13/15 (b = 0, M + ell = 1, case-A energy (1+3ell-2M)/4 = 0);
- the high side is s_high(delta, x) = 13/15 + delta(x/5 - 1/10) <= 13/15 on [0,5/6] x [0,1/2];
- the side constraints hold: principal margin 23/375, l_y - ell = 1/5, M - 2ell = 2/5, 5ell - h = 1/5.

So DH alone reaches the GRH-for-rows barrier. Every zero-density estimate beyond DH, up to GRH, adds nothing in Part II.
- Caveat: D' is at the exponent-bookkeeping level. Turning it into a theorem would require proving DH for rows *and* re-verifying every Part II side lemma at the new geometry, which I have not done.
- The *obstruction* direction (Theorem D) only uses the paper's own exponents and is exact.

---------------------------------------------------------------------------------------------------

## 6. Best provable exponents (the task's question)

**(a) Distinct sixth-power-free rows** (rows of norm ≍ U with a zero in the bin a = (1+delta)/2).

P Prop. detector-counts (line 15185) bins rows by prime amplitude x. The count is R_*(delta, x), largest at x = 0, so the *pure* density estimate is R_*(delta, 0) with t optimised. The two counts cross at t* = 5/(5-2delta), giving

  N(sigma; U) << U^{R(delta)+eps},  R(delta) = 1 - delta + delta(5-6delta)/(3(5-2delta)),  delta = 2sigma-1.

Exponents A(sigma) = R/(1-sigma) (`r3_verify.py` §3):

| sigma | R | Part II (proved) | DH | BGL Cor 1.6 | H-M, exact vanishing (N >= H^{2+eps}) | H-M, convexity off-diag. | optimal sextic LS (conj.) |
|---|---|---|---|---|---|---|---|
| 6/7 | 1/3 | **7/3 = 2.333** | 2 | 4.162 | 4 | 14/3 | ≈3.17 |
| 13/15 | 49/159 | 2.311 | 2 | 4.237 | 4 | 30/7 | — |
| 7/8 | 2/7 | **16/7 = 2.286** | 2 | 4.308 | 4 | 4 | ≈3.33 |
| 11/12 | 1/6 | 2 (= DH) | 2 | 4.727 | 4 | 3 | — |

Halasz–Montgomery is worse than the large sieve for this thin family:
- Every non-duplicate off-diagonal character psi_u conj(psi_{u'}) = psi_{u u'^5} has conductor ≍ N(rad(u u' / gcd^2)), typically ≍ H^2. So exact vanishing of smooth sums needs length N >= H^{2+eps}, and then the count is R <= N^{2-2sigma}, i.e. A >= 4.
- With convexity on L(1/2, psi_{uu'^5}) one needs N^{2sigma-3/2} > H^{1/2}, i.e. A = 2/(4sigma-3).
- Burgess gives A = 3/(8sigma-6) (7/2 at 6/7).
- Only Lindelof for the off-diagonal L-values (GRH strength) makes H–M competitive.
- At the copies H–M is vacuous: R copies of one character give R^2 V^2 <= ||a||^2 N R^2, i.e. Cauchy–Schwarz.

There is also no "distance" in the conductor aspect for Huxley subdivision to exploit:
- In the t-aspect, subdivision gains because sum_n n^{i(t-t')} decays with |t-t'| and zeta has subconvexity.
- Here all non-duplicate pairs are equally far apart, with conductor ≍ H^2 generically. Only the near-duplicates u' = u v with N(v) small have small ratio conductors, and those are different L-functions.

**(b) Multiplicity-counted** (the sense in which A = 7/6 <=> 6/7):
- At the copies, the provable strength is DH (A = 2 + O(c)). This comes from the Mobius-specific diagonal mean square (Proposition C) and gives theta = 11/12.
- For sigma > 7/8 the count is 0 by the 7/8 theorem.
- On [6/7, 7/8] nothing better than DH (11/12-strength) is provable as a *density* statement. Part II's "effective" A = 4/3 (or 20000/15003) only restates theta = 7/8 (- 1/40000), which is obtained through the probe and not through counting.

**Resulting theta from R3: 34999/40000 (no improvement).**

---------------------------------------------------------------------------------------------------

## 7. Notes on the remaining tools in the brief

- **Guth–Maynard.** Their large-values theorem, R << T^{o(1)}(N^{2-2sigma} + N^{18/5-4sigma} + T N^{12/5-4sigma}), never improves the mean-value term N^{2-2sigma}, which is DH level. They improve only the other terms, and essentially only near sigma ≈ 3/4.
  - The mean-value term is sharp for |b_n| <= 1. For example, b_n = R^{-1/2} sum_{j<=R} n^{-i t_j} with R = N^{2-2sigma} has R large values.
  - Their key inputs are the Fourier-analytic structure of the phases n^{it} (stationary phase, intervals, additive energy of the set of t's). The conductor-aspect analogue is Gauss-sum phases on K*/K*^6, which have no ordering or stationary phase.
  - Their bounds are coefficient-blind, so Theorem B applies.
  - Nothing transfers that could beat DH, let alone the multiplicity.
- **Reflected-energy lemma (P Lemma reflected-energy).** This is a low-side input. The 13/15 barrier of Theorem D is a low-side barrier (Gram P_a^{1/6} plus energy wall M + ell <= 1). Improving it is a different route (R1/R2-type), not a density input.
- **Inverse plus plain witnesses (Lemmas marked and plain).**
  - The plain fourth moment (P Lemma plain, case z = 0) is Lindelof-on-average for products of two plain polynomials of *any* lengths. It is very strong, but plain witnesses vanish at the copies (§4) and the count is 1 - 2 delta m with m <= 1/2.
  - DH for distinct rows would need a *mixed* moment, sum_u |M_r S_m|^2 << U for r + m <= 1. Lemma marked cannot supply it: its marks must be prime sums costing 2z of row length, giving r + 2m <= 1.
  - The inverse fourth moment is not diagonal for r > 1/2: mu*mu admits no Gauss–Jacobi absorption.
  - So DH for rows (which would give 13/15) is itself open with these tools.
- **"Duplicates are exact copies."**
  - The exact copy identity is used *positively* only in the simple extraction (P Lemma inverse-amplification, P2 eq. A_u-vs-A_1), and there it yields exactly theta = (6+5kappa)/12 from M_kappa.
  - In every counting argument the copies count as one character with weight H^{1/6} (large sieve) or make H–M vacuous.
  - A coefficient-blind argument cannot "use" them to beat themselves (Theorem B).

---------------------------------------------------------------------------------------------------

## 8. What a proof of 6/7 would need, and why it lies outside R3

The needed input is a *Mobius-specific sub-diagonal* mean value: M_{6/7}, or at least its large-values form for the copies. In R3 language this is a multiplicity-beating density estimate, and by Theorem A it must be proved as a mean value, not counted.

The only known mechanism that is coefficient-specific is P2's chain:
- Poisson in u;
- Gauss–Jacobi absorption mu gamma_{-1} = cubic Gauss sum;
- theta reflection;
- quadratic large sieve (Goldmakher–Louvel).

For H = D^kappa < D, Poisson is an exact identity whose dual has D^2/H > D rows.
- Its positive-mean-square bound (P2 eq. intro-poisson-comparison, line 339) has heuristic size D^2 > D H.
- The true dual is the fluctuation, because its main term carries the factor ∫ hat(weight) = weight(0) = 0. This is prior analysis point (3), and the reformulation of M_kappa as a dual asymptotic is restated from there, not re-derived here.
- So M_kappa for kappa < 1 is equivalent to a power-saving *asymptotic* for the dual mean square of twisted cubic-Gauss sums beyond its diagonal, with error D^{2+eps} against a main term D^{3-kappa}.

Literature status (cf. R5 notes):
- The quadratic analogue is known only *above* the diagonal and only under GRH (Gao–Zhao, arXiv 2212.09241).
- The asymptotic large sieve of Conrey–Iwaniec–Soundararajan needs an extra family average that this family lacks; this is my understanding, not re-checked.
- No such result exists for the sextic/cubic family.

This is the content of other routes, not a density estimate.

---------------------------------------------------------------------------------------------------

## 9. Sharpest conditional statements in R3 language

1. **LV(6/7 + eps, 6/7, eta)** of §2 implies zero-free Re s > 6/7 + eps for every L(s, nu psi_u). Combined with the continuation framework over all nu, this gives beta_* <= 6/7.
   - By Theorem A, LV is *equivalent* to the conclusion for the family.
   - It is weaker than M_{6/7}: LV follows from M_kappa for all kappa in [6/7, 1], with eta = 7eps/3, after a Sobolev step in log D'.
   - It is not weaker than the theorem, so it is not a useful hypothesis beyond M_{6/7}.
2. **M_{6/7}** (prior analysis): implies 6/7 through either extraction. The only nontrivial family-level hypothesis, and not a density statement.
3. **DH for the row family in Part II's binned form, plus floor 1/2 + eta**: implies 13/15 + O(eta) at l_x = l_y = 2/5, ell = 1/5. This is model level (Theorem D'). It is the best any zero-density input can do in Part II (Theorem D), and 13/15 > 6/7.

---------------------------------------------------------------------------------------------------

## 10. Files

| File | Contents |
|---|---|
| `r3_verify.py` | Exact checks: dictionary; A_kappa; Part II pure-density R(delta) (closed form vs brute-force t-grid); BGL / H–M / optimal-LS exponents; Lemma-marked amplified extraction theta(z) = 11/12 + z/3; floor-bin LP certificate (13/15, any delta_0, q); DH attains 13/15 at (2/5, 2/5, 1/5) with side constraints. |
| `r3_coeffblind.py` | Mobius vs plain vs random coefficients over rows u <= D^kappa, with the sixth-power copies separated (exact sextic symbols; reuses `strengthening/wall_experiment.py` read-only). |
