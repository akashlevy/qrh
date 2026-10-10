# Route R2 — self-improving (bootstrap) argument toward Re s > 6/7

Working notes, 8 Oct 2026. Papers: `paper.tex` (7/8, line numbers below refer to it unless marked P2) and
`paper2.tex` (11/12, "P2"). Scripts in this directory: `exact_certificates.py` (stdlib, exact),
`wbeta_consistency.py` (stdlib, exact), `grh_barrier_lp.py`, `grh_combos.py`, `bootstrap_fixed_point.py` (venv).

## 0. Verdict: OBSTRUCTION (plus one small rigorous lemma and a sharp conditional statement)

A bootstrap hypothesis H(β) ("every primitive finite-order Hecke L-function over F = Q(√−3) is zero-free in
Re s > β") cannot by itself carry either architecture to 6/7:

* **Part II architecture.** H(β) only acts on non-principal rows. Even the strongest conceivable row information,
  GRH for every non-principal row, leaves the floor-bin constraint (rows with no detected zero, count R = 1).
  Together with the zero-blind low side this forces σ₀ ≥ **13/15** (exact LP certificate). With the paper's detector
  floor 51/100 the floor is **167/192 ≈ 0.86979**. The self-consistent dynamic bootstrap inside the paper's own lemma set
  (bin ceiling δ ≤ 2β−1 and Lemma plain with κ = 2β−1, even extended below 3/4) does not move θ: float model
  θ(β) ≈ 0.87501 for every β ∈ [0.86, 0.875]; the exact barrier for that lemma set is 0.874957.
* **Extraction architecture (P2 §3).** The only input that reaches 6/7 is a below-diagonal Möbius mean square M_κ with κ ≤ 6/7.
  H(β) cannot supply it. The counterfactual configuration W_β is consistent with H(β), with every standard consequence of H(β),
  and with the conclusion of every mean-value lemma of both papers. In W_β one target has a zero at β ∈ (6/7, 11/12],
  and every other primitive finite-order Hecke L-function satisfies GRH. M_κ fails in W_β for every κ < 6(2β−1)/5,
  solely because of the principal sixth-power class H^{1/6}|A₁|². That is the target itself, so the circularity is intrinsic.
* **Family zero-density.** Density-hypothesis-strength counting for the sextic row family gives exactly 11/12.
  The sixth-power multiplicity H^{1/6} must be beaten. At σ = 6/7 the required exponent is < 1/6, while DH gives 2/7.
  The "DH-world" is consistent with all bootstrap consequences.

The rigorous by-product: H(β) removes the D²/Y² term of P2's prime extraction. Hence M_κ alone gives β* ≤ (6+5κ)/12
for every κ ∈ (0,1] (§2–3). Without the bootstrap, 6/7 is exactly the *floor* of P2's prime extraction for any
mean-square input. With it, M_{6/7} gives 6/7 and any κ < 6/7 would give less.

## 1. What H(β) provides (all already in paper.tex)

* Continuation criterion, Prop. `lem:continuation-criterion` (l.400–502). The high side is measured against C(β*),
  and the low side against C(σ₀)+ω with ω < β*−σ₀. Equivalent bookkeeping: contradiction iff max(σ_low, σ_high) < β*.
  So the C(β*) feature does not by itself relax anything (§4).
* Lemma `lem:logarithmic-control` (l.1531–1600), global form. For b ≥ β* and v > 0:
  |L(s,ψ)| + |L(s,ψ)⁻¹| ≪_{v,ε} (Q(3+|t|)²)^ε on Re s ≥ b+v, and |L'/L| ≪ log. This is uniform in the conductor Q
  (Borel–Carathéodory plus Hadamard three circles, from Lemma `lem:hecke-strip-growth`).
  Deleted Euler factors cost q_R^ε (Lemma `lem:deleted-euler-factors`).
* Consequences:
  * Möbius sums: Σ_n μ(n)ψ(n)W(N n/X) ≪ X^{β+v}(QN(rad))^ε.
  * Plain sums: X^{β+v}Q^ε, a Lindelöf-type saving for short sums.
  * Zero-free rectangles for the detector (Lemma `lem:buffered-bins`, l.4281).
  * The bin ceiling δ ≤ 2β*−1 (l.4379).
* Part II already uses the bootstrap. Lemma `lem:plain` assumes β* ≤ (1+κ)/2, κ ∈ [3/4,1], and is applied with
  κ = 2β*−1 (l.12564, 12581, 15097).

## 2. Lemma B (bootstrapped prime extraction) — rigorous

**Lemma B.** Assume H(β) with β > 1/2. Fix ν, W ∈ C_c^∞((0,∞)) and v > 0. For every primary prime p ∉ S,
uniformly in p,

  |A₁(D) − A_{p⁶}(D)| ≪_{ν,W,v} (D/N p)^{β+v}.

*Proof.* χ_n(p⁶) = 1_{p∤n}. So A₁ − A_{p⁶} = Σ_{p|n} μ(n)ν(n)W(Nn/D) = −ν(p)·A^{(p)}(D/Np), where
A^{(p)}(X) = Σ_{(m,pS)=1} μ(m)ν(m)W(Nm/X) = (2πi)⁻¹∫_{(2)} Ŵ(s)X^s /L^{S∪{p}}(s,ν) ds.
Here 1/L^{S∪{p}} = (1/L^S)(1 − ν(p)Np^{−s})⁻¹. The second factor is ≤ (1−2^{−1/2})⁻¹ on Re s ≥ 1/2. The first is
≪ (3+|t|)^ε on Re s ≥ β+v by Lemma `logarithmic-control` (fixed conductor), and Ŵ decays rapidly. Shift to
Re s = β+v. ∎

**Corollary (extraction).** Let Y = (H/2)^{1/6} and J ≍ Y/log Y (Landau). Then

  |A₁(D)|² ≤ (2/J)Σ_{Y/2<Np≤Y}|A_{p⁶}(D)|² + O((D/Y)^{2β+2v}),

which replaces the term D²/Y² of P2 eq. `prime-extract` (P2 l.711–716).

## 3. Conditional results

**Theorem C.** Let κ ∈ (0,1]. Suppose (M_κ): for all primitive ν, W and ε > 0,
Σ_{0<N u≤D^κ}|A_u^ν(D)|² ≪ D^{1+κ+ε}. Then β* ≤ (6+5κ)/12. In particular M_{6/7} gives Re s > 6/7.

*Proof.* Suppose β* > θ_κ := (6+5κ)/12 and put Δ = β*−θ_κ. H(β*) holds by definition of β*. By the corollary with
H = D^κ:
|A₁|² ≪ D^{1+5κ/6+2ε} + D^{(1−κ/6)(2β*+2v)}.
The first exponent is 2β* − 2Δ + 2ε. The second is ≤ 2β* − κβ*/3 + 2v.
Take ε = Δ/2 and v = κβ*/12. This gives |A₁^ν(D)| ≪_{ν,W} D^{β*−c} with c = min(Δ/2, κβ*/12), independent of ν.
The Mellin argument of P2 §3 (l.726–757) then makes 1/L^S(s,ν) holomorphic on Re s > β*−c for every ν
(for a zero ρ, take W = y^{−ρ}φ). This contradicts the definition of β*. ∎

(The same exponent follows unconditionally from the identity eq. `sixth-power-column-identity` in Lemma
`inverse-amplification`, l.12426, with raw-moment replaced by M_κ. Lemma B is an alternative five-line proof.)

**Weakest form usable by R2 (BMS).** It suffices that for every β ∈ (6/7, 7/8] the following holds:
under H(β) there are κ < 6(2β−1)/5 and c > 0 with

  Σ_{N u≤D^κ}|A_u(D)|² ≪ D^{1+κ+ε} + D^{2β−c+κ/6}.

Then Lemma B forces |A₁| ≪ D^{β−c'}, and so β* ≤ 6/7. BMS is strictly weaker than M_{6/7}, in three ways:
* it is conditional;
* it has an error allowance just below the principal-class contribution;
* κ need only approach 6/7 as β ↓ 6/7.

But BMS is false in W_β (§5), so a proof of BMS must contain zero-blind information about the principal class.

**Remark (prime-extraction floor).** Without Lemma B, P2's extraction gives θ(κ) = max((6+5κ)/12, 1−κ/6). Its minimum
over all κ is exactly 6/7, attained at κ = 6/7. So any argument feeding a below-diagonal mean square into
P2 eq. `prime-extract` stops at exactly 6/7, however strong the input.

## 4. Obstruction I — Part II bookkeeping (exact)

The high exponent, eq. `common-high-exponent` (l.15726), at d = h is
E = (1+δ)(1−l_y)/2 − σ₀ − ℓ/2 + qℓ + h(R + δ/2 − 1/6), with h = 1−l_x+ℓ.

The floor bin (l.15843) has R = 1, q ≤ δ₀/2 and δ₀ = 1/50. It contains every row without a detected zero, which is
most rows in any world. It is therefore present whatever row information H(β) supplies. GRH for rows corresponds to
δ₀ → 0 with no intermediate bins.

| configuration (low side as in paper: Prop `probe-low`, `probe-gram`, Lemma `reflected-energy`) | exact optimum | certificate |
|---|---|---|
| GRH for all non-principal rows (floor at 1/2) | **13/15** | ⅘·low + ⅕·floor + ⅕·[g ≥ b/6] + ⅕·[g≥0] + 1/15·[e ≥ M+ℓ−1] + ⅓·[e≥0] |
| detector floor 51/100, rows otherwise perfect | **167/192** = 0.869792 | 103/128·low + 25/128·floor + 53/256·[g≥b/6] + … |
| GRH rows + energy wall relaxed by w (W(w)) | 13/15 − w/60 | vertex search |
| GRH rows + Gram bound **without** the P_a^{1/6} term | **6/7 exactly**, at (l_x,l_y,ℓ) = (2/7,4/7,1/7) | 17/21·low + 4/21·floor + 1/42·[g ≥ 2b−l_y] + 1/14·[e ≥ M+ℓ−1] + slack |

* The 13/15 certificate uses only the Gram exceptional term b/6 (sixth-power frequencies, l.8534–8556) and the energy wall
  M+ℓ−1. Both are zero-blind statements about theta coefficients and additive Gauss sums. H(β) has nothing to act on there.
* The last row shows what Part II would need for 6/7, beyond the bootstrap. It would need GRH-quality rows. At the 6/7
  vertex these must satisfy R(δ,x) ≤ 1 − δ(3/4 + x/6) in every intermediate bin; DH counts suffice and the paper's
  R* ≈ 1−2δ/3 does not. It would also need removal of the structural P_a^{1/6} term.
* Dynamic bootstrap inside the paper's lemma set (`bootstrap_fixed_point.py`, corrected float model) gives
  θ(β) = 0.87501 for β ∈ {0.875, 0.87, 0.8667, 0.86}. DH counts or GRH rows with floor 51/100 give 0.869808.
  GRH with floor → 1/2 gives 0.866683.
* Why the bootstrap does not even give DH counts. The detector (Prop. `detector-witness`, l.4510) gives a split
  r + m ≥ t, m ≤ 1/2, which is not controlled.
  * For t > 1, the inverse witness can have r > 1, where e(r) = (1+5r)/6 loses α = 5/6. That loss is a below-diagonal
    inverse moment for sixth-power-free rows, false in the DH-world (§6).
  * The mixed split needs Σ_u|M_r S_m|² ≪ U at r+m = 1, i.e. the sextic large sieve without the (QN)^{2/3} term
    (open; zero-blind).

## 5. Obstruction II — extraction architecture

(a) **Circularity.** Any bound for S(H,D) = Σ_{N u≤H}|A_u|² contains the principal class Σ_a|A_{a⁶}|² ≍ H^{1/6}|A₁|².
H(β) bounds it only by H^{1/6}D^{2β}, which is consistent with a zero at β. The large sieve, and Halász–Montgomery with
H(β)-bounds (|Σ_n χ_n(uu'^5)W| ≪ D^{β}), give nothing on this class. Its Gram block is the all-D block of size H^{1/6}.

(b) **W_β is consistent with every input** (`wbeta_consistency.py`, exact). Every row outside the principal class is at
diagonal size in W_β; only that class matters.
* thm:ms / raw-moment: needs β ≤ 11/12 + 5θ/12.
* inverse-amplification at U ≍ 1: needs β ≤ 11/12 + 5c/12.
* Lemma marked, all admissible (m,r,z), principal class with aligned slot coefficients so that |Q_u|² = Z^z:
  max r(2β−1) + z − 5m/6 over the vertices (m,0), (m/2,m/4), (0,3m/8) is ≤ 0 iff β ≤ 11/12.
* Lemma plain: the principal class is excluded from R_z for z > 0, and plain sums do not see zeros.
* Detector and row counts: there are no intermediate rows.
* Part II probe inequalities at every geometry: consistent iff β ≤ 13/15 (§4).
* Low-side lemmas: zero-blind theorems.

Hence no exponent-level combination of these inputs with H(β) proves β* < β for any β ≤ 13/15. Since 6/7 < 13/15,
R2 fails.

(c) **Where the obstruction lives, in all three pictures.** Take H = D^κ with κ < 1.
* Möbius picture: rows H < columns D, so the problem is below the diagonal.
* Poisson dual (P2 Step 3): cubic-Gauss-sum sums B_h(D), rows calH = D^{2−κ} > D. The positive-weight dual contains the dual
  diagonal D·calH = D^{3−κ}, giving S ≲ D², which is P2's lossy step. With the exact signed weight Φ̂ (Φ vanishes near 0),
  every polynomial-in-N(h) main term cancels, the dual diagonal included. The remainder is the dual off-diagonal, and Poisson
  back turns it into the original off-diagonal (circular).
* Theta reflection (P2 Step 4): rows calH, reflected columns calH²/D = D^{3−2κ} > calH iff κ < 1, so again below the diagonal.

So every picture is beyond the large sieve exactly when κ < 1. What is needed is an *asymptotic* evaluation, with power
saving past the diagonal, of a quadratic-twisted mean square of cubic theta coefficients. That is metaplectic and zero-blind.

The only Hecke L-functions on the theta side enter through cube removal (P2 l.426–611), as 1/L(3s−1/2, w_h) with
w_h = ᾱ³ξ³(h/·)₂. Since α(n) = n/|n|, these are Grössencharacters of infinity type 3, and they are not covered by β*,
which is defined over finite-order characters. Cube removal is in any case handled unconditionally by P2's recursion.

## 6. Obstruction III — family density

The explicit formula under H(β) gives S(H,D) ≲ DH + max_σ N_fam(σ;H)D^{2σ}. M_κ is therefore the density statement
N_fam(σ;H) ≪ H^{1−(2σ−1)/κ}. The principal class forces N_fam(β;H) ≥ H^{1/6}.

* DH for the family (H^{2−2σ}; the Part II count U^{1−δ} at r = 1) gives only σ ≤ 11/12.
* 6/7 needs N_fam(σ;H) < H^{1/6} as σ ↓ 6/7, against DH's H^{2/7}: a gap of H^{5/42}.
* BMS even needs the number of bad sixth-power-free rows u ~ U to be ≤ U^{1/6}D^{−c}.
* The DH-world has U^{2−2β} bad rows per scale. It satisfies H(β), all large sieves, thm:ms (β ≤ 11/12) and
  inverse-amplification (with equality at r = 1). In it M_κ fails for every κ < 1, with excess exponent
  κ(2−2β) + 2β − 1 − κ > 0.

H(β) gives no edge-improved density. Three-lines between Re s = 1/2 and β+ε does not detect zeros within o(1) of the
edge (checked: the condition reduces to δ ≥ ε'(2β−2), which is always true).

## 7. Summary of numbers

* 0.874957: barrier for the paper's lemma set, the bootstrap already included (prior analysis). 7/8 − 1/40000 is certified.
* 167/192: any row improvement with detector floor 51/100.
* 13/15: any row improvement including GRH for rows.
* 11/12: extraction with known mean squares (thm:ms, κ > 1), with or without the bootstrap.
* (6+5κ)/12: extraction with M_κ (Theorem C). Equals 6/7 iff κ = 6/7.
* 6/7: GRH rows plus the Gram bound without P_a^{1/6}, in Part II bookkeeping (exact). Equivalently, M_{6/7} in the
  extraction. Both inputs are beyond the bootstrap.

## 8. Status of each claim (what is proved vs exponent-level)

* **Rigorous, from statements proved in the papers:**
  * Lemma B and Theorem C. These use Lemma `logarithmic-control`, Landau's prime ideal theorem, and the Mellin argument of P2 §3.
  * The prime-extraction floor of 6/7 (elementary algebra on P2 eq. `prime-extract`).
* **Exact but model-level:** the LP certificates in §4 (`exact_certificates.py`). They are exact statements about the
  Part II exponent bookkeeping as transcribed in QRH-analysis.md and the audit:
  * low side: Prop. `probe-low`/`probe-gram`, Lemma `reflected-energy` with the paper's case-A term;
  * floor bin: eq. `floor-bound`, `common-high-exponent`.
  The 13/15 certificate uses only the robust rows (energy wall M+ℓ−1 and Gram b/6), so it does not depend on the
  disputed energy modelling.
* **Exponent-level relativization (heuristic in the usual sense):**
  * W_β consistency (§5b, `wbeta_consistency.py`);
  * the DH-world (§6);
  * the "three pictures below the diagonal" description (§5c).
  These show that no argument using only the listed inputs plus H(β) can reach below 13/15 (Part II) or 11/12
  (extraction). They do not exclude a genuinely new zero-blind input. That input is exactly what is missing:
  M_{6/7}, or a Part II low side without the P_a^{1/6} Gram term together with GRH-quality rows.
* **Numerics** (float, coarse grid): `bootstrap_fixed_point.py` (dynamic bootstrap θ(β) ≈ 0.87501, constant in β).
