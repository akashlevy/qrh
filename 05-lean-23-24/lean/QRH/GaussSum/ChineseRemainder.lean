import QRH.Arithmetic.EisensteinCoordinates
namespace OAI

open scoped BigOperators
open MulChar AddChar
open Filter Asymptotics MeasureTheory
open scoped Topology
open Real
open scoped FourierTransform SchwartzMap
noncomputable section
open Finset Complex

open scoped Classical

namespace GaussCRT

end GaussCRT

namespace IdealGaussCRT

variable {T R S : Type*} [CommRing T] [CommRing R] [CommRing S]

noncomputable def quotientAddChar {A : Type*} [CommRing A]
    (I : Ideal A) (ψ : AddChar A ℂ)
    (hper : ∀ a : A, a ∈ I → ψ a = 1) : AddChar (A ⧸ I) ℂ where
  toFun := fun x => Quotient.liftOn' x ψ (by
    intro a b hab
    have hsub : a - b ∈ I := I.quotientRel_def.mp hab
    have h := ψ.map_add_eq_mul (a - b) b
    rw [sub_add_cancel] at h
    rw [hper _ hsub, one_mul] at h
    exact h)
  map_zero_eq_one' := by
    change ψ 0 = 1
    exact AddChar.map_zero_eq_one ψ
  map_add_eq_mul' := by
    intro x y
    refine Quotient.inductionOn₂' x y ?_
    intro a b
    exact AddChar.map_add_eq_mul ψ a b

@[simp] theorem quotientAddChar_mk {A : Type*} [CommRing A]
    (I : Ideal A) (ψ : AddChar A ℂ)
    (hper : ∀ a : A, a ∈ I → ψ a = 1) (a : A) :
    quotientAddChar I ψ hper (Ideal.Quotient.mk I a) = ψ a := rfl

theorem bezout_of_principal_coprime {A : Type*} [CommRing A]
    {a b : A} (h : IsCoprime (Ideal.span {a}) (Ideal.span {b})) :
    ∃ u v : A, u * a + v * b = 1 := by
  obtain ⟨x, hx, y, hy, hxy⟩ := Ideal.isCoprime_iff_exists.mp h
  obtain ⟨u, hu⟩ := Ideal.mem_span_singleton.mp hx
  obtain ⟨v, hv⟩ := Ideal.mem_span_singleton.mp hy
  refine ⟨u, v, ?_⟩
  rw [hu, hv] at hxy
  convert hxy using 1; ring

noncomputable def rawTraceChar {A : Type*} [CommRing A]
    (ι : A →+* ℂ) (ψ : AddChar ℂ ℂ) (lam : ℂ) (c : A) : AddChar A ℂ where
  toFun x := ψ (ι x / (ι c * lam))
  map_zero_eq_one' := by simp
  map_add_eq_mul' x y := by
    rw [map_add, add_div, AddChar.map_add_eq_mul]

theorem rawTraceChar_period {A : Type*} [CommRing A]
    (ι : A →+* ℂ) (ψ : AddChar ℂ ℂ) (lam : ℂ) (c : A)
    (hc : ι c ≠ 0) (hlam : lam ≠ 0)
    (hbase : ∀ x : A, ψ (ι x / lam) = 1)
    (x : A) (hx : x ∈ Ideal.span {c}) :
    rawTraceChar ι ψ lam c x = 1 := by
  obtain ⟨y, hy⟩ := Ideal.mem_span_singleton.mp hx
  rw [hy]
  change ψ (ι (c * y) / (ι c * lam)) = 1
  rw [map_mul]
  have hfrac : ι c * ι y / (ι c * lam) = ι y / lam := by
    field_simp
  rw [hfrac]
  exact hbase y

noncomputable def traceModChar {A : Type*} [CommRing A]
    (ι : A →+* ℂ) (ψ : AddChar ℂ ℂ) (lam : ℂ) (c : A)
    (hc : ι c ≠ 0) (hlam : lam ≠ 0)
    (hbase : ∀ x : A, ψ (ι x / lam) = 1) :
    AddChar (A ⧸ Ideal.span {c}) ℂ :=
  quotientAddChar (Ideal.span {c}) (rawTraceChar ι ψ lam c)
    (rawTraceChar_period ι ψ lam c hc hlam hbase)

@[simp] theorem traceModChar_mk {A : Type*} [CommRing A]
    (ι : A →+* ℂ) (ψ : AddChar ℂ ℂ) (lam : ℂ) (c : A)
    (hc : ι c ≠ 0) (hlam : lam ≠ 0)
    (hbase : ∀ x : A, ψ (ι x / lam) = 1) (x : A) :
    traceModChar ι ψ lam c hc hlam hbase
      (Ideal.Quotient.mk (Ideal.span {c}) x) =
      ψ (ι x / (ι c * lam)) := rfl

theorem bezout_inverse_left {A : Type*} [CommRing A]
    {a b u v : A} (hbez : u * a + v * b = 1) :
    (Ideal.Quotient.mk (Ideal.span {a}) b) *
      (Ideal.Quotient.mk (Ideal.span {a}) v) = 1 := by
  rw [← map_mul]
  apply Ideal.Quotient.eq.mpr
  change b * v - 1 ∈ Ideal.span {a}
  rw [Ideal.mem_span_singleton]
  refine ⟨-u, ?_⟩
  calc
    b * v - 1 = b * v - (u * a + v * b) := by rw [hbez]
    _ = a * -u := by ring

theorem gauss_shift_bezout_left {A : Type*} [CommRing A]
    {a b u v : A} [Fintype (A ⧸ Ideal.span {a})]
    (hbez : u * a + v * b = 1)
    (χ : MulChar (A ⧸ Ideal.span {a}) ℂ)
    (ψ : AddChar (A ⧸ Ideal.span {a}) ℂ) :
    gaussSum χ (ψ.mulShift (Ideal.Quotient.mk (Ideal.span {a}) v)) =
      χ (Ideal.Quotient.mk (Ideal.span {a}) b) * gaussSum χ ψ := by
  let mb := Ideal.Quotient.mk (Ideal.span {a}) b
  let mv := Ideal.Quotient.mk (Ideal.span {a}) v
  have hprod : mb * mv = 1 := bezout_inverse_left hbez
  have hu : IsUnit mv := IsUnit.of_mul_eq_one_right mb hprod
  have hχ : χ mb * χ mv = 1 := by rw [← map_mul, hprod, map_one]
  have hshift := gaussSum_mulShift χ ψ hu.unit
  have hshift' : χ mv * gaussSum χ (ψ.mulShift mv) = gaussSum χ ψ := by
    simpa only [hu.unit_spec] using hshift
  calc
    gaussSum χ (ψ.mulShift mv) =
        (χ mb * χ mv) * gaussSum χ (ψ.mulShift mv) := by rw [hχ, one_mul]
    _ = χ mb * (χ mv * gaussSum χ (ψ.mulShift mv)) := by ring
    _ = χ mb * gaussSum χ ψ := by rw [hshift']

end IdealGaussCRT

end

end OAI
