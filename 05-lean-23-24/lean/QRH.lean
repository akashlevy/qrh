import QRH.Spine.ZeroFree
import QRH.Spine.Extraction
import QRH.Spine.MeanSquare
import QRH.Spine.Milestones

/-!
# A quasi-Riemann hypothesis: `ζ(s) ≠ 0` for `Re s > 23/24`

The proof follows the paper `paper/qrh.tex`. The *spine* (`QRH/Spine/*`) states each step
of the paper as a Lean theorem; this file only chains them:

* `QRH.meanSquare` (Proposition 3.1): the sextic family has square-root cancellation on
  average over rows `u` of norm `≤ D^{11/10}`;
* `QRH.saving_of_meanSquare` (Proposition 3.2): averaging the rows `u = b⁶` extracts a power
  saving `A^φ(D) ≪ D^{23/24+ε}` for the totient-weighted Möbius sum;
* `QRH.riemannZeta_ne_zero_of_saving` (Proposition 3.3): a power saving for these sums
  continues `1/ζ_K` holomorphically to `Re s > 23/24`.

The remaining `QRH.*` modules are the machinery behind these steps: a declaration-level
tree-shaking (scripts/treeshake.py) of the closure of `OAI.ShortDraft.riemannZeta_ne_zero_of_re_gt`
in OpenAI's lean-base (openai/math @ adc7f124).
-/

namespace QRH

open scoped ContDiff

/-- **Theorem 1.1.** The Riemann zeta function has no zeros in the half-plane `Re s > 23/24`. -/
theorem riemannZeta_ne_zero_of_twentyThree_twentyFourths_lt_re
    {s : ℂ} (hs : (23 / 24 : ℝ) < s.re) : riemannZeta s ≠ 0 := by
  refine riemannZeta_ne_zero_of_saving (fun W hW hWs ε hε => ?_) hs
  have hsupp : Function.support W ⊆ Set.Icc 1 2 := fun x hx =>
    ⟨le_of_not_gt fun h => hx (hWs x (Or.inl h.le)),
      le_of_not_gt fun h => hx (hWs x (Or.inr h.le))⟩
  exact saving_of_meanSquare W 1 2 hsupp hW ε hε
    (meanSquare W 1 2 one_pos hsupp hW (2 * ε) (by positivity))

end QRH
