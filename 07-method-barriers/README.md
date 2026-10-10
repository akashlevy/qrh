# 07 — Method barriers and low-side no-go results (Oct 9, 2026)

These exact barrier certificates and audits come from the GPT-6 Astra research session. Each one bounds what a stated **package of sufficient estimates** can certify. None is a statement about where zeros actually are.

| File | Barrier | Status |
|---|---|---|
| `feedback_exponent_barrier.txt`, `feedback_barrier_certificate.py` | **κ-feedback barrier.** With the unchanged exponent package (geometry plus feeding proved bounds back into κ), no target below **θ\* ≈ 0.874957019420099** can be certified. θ\* = 11/12 − λ\*/4 with 657λ³ − 954λ² + 21λ + 20 = 0; equivalently θ\* is the root in [13/15, 7/8] of 7884θ³ − 18819θ² + 14643θ − 3686 = 0. The bound is sharp for the three-constraint relaxation. This is the same cubic as the public "cubic-root" endpoint; it was derived independently here. | Exact; passes |
| `method_barriers.py` | **13/15 low-majorant barrier.** The unchanged direct-Gram and reflected-row low majorants force B ≥ 13/15 via an exact non-negative dual identity. This blocks 6/7 and 1/2 without a low-side change; the Lean lemma is `Parameters.physical_low_barrier` in `../08-…/lean`. Also: **prime-sixth-power extraction has exact minimum 6/7**; the critical-class substitution; dual weights giving θ ≥ (7+18δ)/(9+18δ). | Exact; passes. It imports `bootstrap_certificate` from `../06-kappa-bootstrap/certificates`. |
| `independent_low.txt` | **Low-side no-gos:** an auxiliary residue character with ξ⁶ ≠ 1 removes the z = 1/6 pole; varying-Q kernel cancellation kills the signal; sums of separately bounded probes cannot beat min_j(L_j − c_j); a direct reciprocal-L per-column strategy needs λ < 2(1 − B) ≈ 0.25 but has λ ≈ 4/5. **Inconclusive ideas:** projection onto ρ, angular modulation, Gauss-phase cancellation (τ₆(m³) unknown). | Proved, scoped |
| `independent_signed_dual.txt` | **Signed dual of paper2's mean square.** Proves exact common-divisor cancellation, removes the fake diagonal (O(D^{1+ε})), and gives a conductor-uniform large-j tail ≪ D³/(B²J^{17/8−δ}). **Gap:** the unit j-block needs a coefficient-sensitive O(D^{2+ε}) bound for the coprime-pair form; absolute or positive completion loses D/H. So M_κ for κ < 1 does not follow. | Partial |

The weighted-numerator result in `../08-weighted-numerator-10499-12000` is not covered by the first barrier: it changes the high-side numerator estimate, which is outside the "unchanged exponent package". The 13/15 low-side barrier still applies to it.

```bash
python3 feedback_barrier_certificate.py
python3 method_barriers.py
```
