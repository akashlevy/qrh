import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

set_option autoImplicit false

namespace QRHResearch.KappaThirteenEighteenths

-- Algebraic parameter lemmas only; this file does not assert an L-function theorem.
theorem terminal_slot (M A z κ : ℝ)
    (hκ : (13/18:ℝ) ≤ κ) (hz : 0 ≤ z)
    (hAz : z ≤ A) (hcap : A + (6*κ-1)*z ≤ M) :
    κ*z ≤ M/6 ∧ z ≤ 3*M/13 := by
  have hprod := mul_le_mul_of_nonneg_right hκ hz
  constructor <;> nlinarith

theorem slot_budgets (M total ell κ : ℝ)
    (hell : 0 ≤ ell) (hκ : (13/18:ℝ) ≤ κ)
    (hlarge : 5*M/6 ≤ total+ell)
    (hcap : total+ell+(6*κ-1)*ell ≤ M) :
    ell ≤ M/20 ∧ (6*κ-1)*ell ≤ M/6 := by
  have hkell : (10/3:ℝ)*ell ≤ (6*κ-1)*ell :=
    mul_le_mul_of_nonneg_right (by linarith) hell
  constructor <;> nlinarith

theorem robust_slot_budgets (M total ell κ defect : ℝ)
    (hell : 0 ≤ ell) (hκ : (13/18:ℝ) ≤ κ)
    (hlarge : 5*M/6-defect ≤ total+ell)
    (hcap : total+ell+(6*κ-1)*ell ≤ M) :
    ell ≤ M/20+3*defect/10 ∧ (6*κ-1)*ell ≤ M/6+defect := by
  have hkell : (10/3:ℝ)*ell ≤ (6*κ-1)*ell :=
    mul_le_mul_of_nonneg_right (by linarith) hell
  constructor <;> nlinarith

theorem robust_balanced_long_nonneg (M total ell κ xi defect : ℝ)
    (hM : 0 ≤ M) (hell : 0 ≤ ell) (hκ : (13/18:ℝ) ≤ κ) (hxi : 0 ≤ xi)
    (hlarge : 5*M/6-defect ≤ total+ell)
    (hcap : total+ell+(6*κ-1)*ell ≤ M)
    (hreserve : xi+8*defect/5 ≤ M/15) :
    8*M/15-13*defect/10 ≤ total-M/4 ∧ 0 ≤ total-M/4 := by
  obtain ⟨he, _⟩ := robust_slot_budgets M total ell κ defect hell hκ hlarge hcap
  constructor <;> linarith

theorem robust_reflected_low (M short along ell z κ xi defect reflected : ℝ)
    (hM : 0 ≤ M) (hsM : short ≤ M/4)
    (hze : z ≤ ell) (hκ : (13/18:ℝ) ≤ κ)
    (he : ell ≤ M/20+3*defect/10)
    (hbudget : (6*κ-1)*ell ≤ M/6+defect)
    (hxi : 0 ≤ xi) (hdefect : 0 ≤ defect) (hdefectM : defect ≤ M)
    (hhigh : 5*M/6 < short+max 0 along+z)
    (href : reflected ≤ max 0 (M-along+xi)) :
    0 < along ∧ reflected+short+z ≤ 23*M/30+xi+3*defect/5 ∧
      reflected+short+6*κ*z ≤ 14*M/15+xi+8*defect/5 := by
  have haz : 0 < along := by
    by_contra hn
    have ha : along ≤ 0 := le_of_not_gt hn
    rw [max_eq_left ha] at hhigh
    linarith
  rw [max_eq_right haz.le] at hhigh
  have hzbudget : (6*κ-1)*z ≤ M/6+defect :=
    (mul_le_mul_of_nonneg_left hze (by linarith)).trans hbudget
  have hzcap : z ≤ M/20+3*defect/10 := hze.trans he
  refine ⟨haz, ?_⟩
  by_cases hr : 0 ≤ M-along+xi
  · rw [max_eq_right hr] at href
    constructor <;> nlinarith
  · rw [max_eq_left (le_of_not_ge hr)] at href
    constructor <;> nlinarith

theorem robust_reflected_low_admissible (M xi defect short reflected z κ : ℝ)
    (hd : 0 ≤ defect) (hreserve : xi+8*defect/5 ≤ M/15)
    (ht : reflected+short+z ≤ 23*M/30+xi+3*defect/5)
    (hc : reflected+short+6*κ*z ≤ 14*M/15+xi+8*defect/5) :
    reflected+short+z ≤ 5*M/6 ∧ reflected+short+6*κ*z ≤ M := by
  constructor <;> linarith

#print axioms terminal_slot
#print axioms slot_budgets
#print axioms robust_slot_budgets
#print axioms robust_balanced_long_nonneg
#print axioms robust_reflected_low
#print axioms robust_reflected_low_admissible
end QRHResearch.KappaThirteenEighteenths
