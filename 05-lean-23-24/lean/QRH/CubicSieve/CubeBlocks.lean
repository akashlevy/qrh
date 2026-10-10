import QRH.QuadraticSieve.UniformTailBounds
namespace OAI

noncomputable section

open scoped BigOperators
open MulChar AddChar
open Filter Asymptotics MeasureTheory
open scoped Topology
open Real
open scoped FourierTransform SchwartzMap
open Finset Complex
open scoped Classical
open ActualEisensteinCubic
open RationalPrimeExtraction ShortDraftLatticeCount
open EisensteinEmbedding ConcreteTraceCRT
open scoped LSeries.notation ArithmeticFunction.Moebius
open Ideal NumberField RingOfIntegers UniqueFactorizationMonoid
open scoped ContDiff
open Set

namespace SecondPassArithmetic

section

open scoped BigOperators Classical SchwartzMap ContDiff
open MeasureTheory
open ActualEisensteinCubic
open FirstPassCubeLabels (primeProduct firstLogDensity)
open ConcreteTraceCRT (eisEmbedding)
open RayFourExpansion (RayCharacter)

theorem globalFirstCoreRows_integrated_transfer
    (ε : ℝ) (hε : 0<ε) (g V : 𝓢(ℝ,ℂ)) (M : ℝ) (hM : 0≤M)
    (hgM : ∀ t,g t ≠ 0 → |t|≤M) (side : Bool) (A J N : ℕ) :
    ∃ (windows : Fin 7 → ℝ → ℂ) (C Ct : ℝ),0≤C ∧ 0<Ct ∧
      (∀ i,HasCompactSupport (windows i)) ∧ (∀ i,ContDiff ℝ ∞ (windows i)) ∧
      (∀ i t,windows i t ≠ 0 → |t|≤M+6+1) ∧
      ∀ {ι : Type*} [DecidableEq ι]
      (p : ι → O) (hp : ∀ i,p i ≠ 0) [∀ i,(Ideal.span {p i}).IsMaximal]
      (hcop : Pairwise (Function.onFun IsCoprime (fun i => Ideal.span {p i})))
      (hg : ∀ i,lambda ∉ Ideal.span {p i})
      (_hinj : Function.Injective (fun i => Ideal.span {p i}))
      (_hc : ∀ i,ringChar (O ⧸ Ideal.span {p i}) ≠ 2) (_hpr : ∀ i,lambda^2 ∣ p i-1)
      (pool : Finset ι) (s : Finset (GlobalFirstData ι))
      (w : GlobalFirstData ι → ℝ) (Γ K ell B F H : ℝ)
      (Ψ : O →* ℂ) (m : O) (χ : RayCharacter) (r : FirstCoreIndex)
      (T : GlobalFirstData ι → Finset O),
      0≤Γ → (∀ b∈s,0≤w b ∧ w b≤Γ) → 0<K → 0<ell → 0<B → 0<F → 0≤H →
      (∀ b∈s,GlobalFirstAdmissible b) → (∀ a,‖Ψ a‖≤1) →
      (∀ b∈s,‖eisEmbedding (primeProduct p b.cube.support b.cube.leftExponent)‖^2≤B) →
      (∀ b∈s,‖eisEmbedding (primeProduct p b.cube.support b.cube.rightExponent)‖^2≤B) →
      (∀ b∈s,∀ z∈T b,(Ideal.absNorm (Ideal.span {z}) : ℝ)≤globalFirstPooledRow p K ell B F side b) →
      (∀ b∈s,∀ z∈T b,z≠0) →
      let source := globalFirstSupportedPool p pool s
        (fun b G E => (globalFirstConcreteCutoff p K ell B F M H side b G E).erase 0) ell M side
      let weight := fun b => w b*globalFirstCoefficient p K ell B F side b
      let Ψmode := firstCoreTwist side χ Ψ r
      (∫ t : ℝ,firstLogDensity (2*(J+2)) t*
        globalFirstFiniteEnergy p hg pool s weight Ψ m g V ell side χ r T t) ≤
      (∫ t : ℝ,firstLogDensity 0 t)*
        (globalFirstDiagonalCost p hg pool s (fun b => weight b) Ψmode m
            (firstCoreBaseProfile g V side) rowMajorant ell side (globalFirstPooledRow p K ell B F side) +
          globalChildBudget p hp hcop hg pool source Ψmode m windows (Γ*C) ε 0 0 K ell B F M H A J side +
          Ct*globalFirstTailCost p s (fun b => weight b) K ell B F M H N side) := by
  obtain ⟨windows,C,hC,hwc,hws,hwb,hrows⟩ := globalFirstCoreRows_actual_transfer ε hε g V M hM hgM side A J
  obtain ⟨Ct,hCt,htail⟩ := globalFirstSourceTail_bound N g V rowMajorant side
  refine ⟨windows,C,Ct,hC,hCt,hwc,hws,hwb,?_⟩
  intro ι _ p hp _ hcop hg hinj hc hpr pool s w Γ K ell B F H Ψ m χ r T hΓ hw hK hell hB hF hH hs hΨ hb₁ hb₂ hT hT0
  dsimp only
  let source := globalFirstSupportedPool p pool s
    (fun b G E => (globalFirstConcreteCutoff p K ell B F M H side b G E).erase 0) ell M side
  let weight := fun b => w b*globalFirstCoefficient p K ell B F side b
  let Ψmode := firstCoreTwist side χ Ψ r
  let d := globalFirstDiagonalCost p hg pool s (fun b => weight b) Ψmode m
    (firstCoreBaseProfile g V side) rowMajorant ell side (globalFirstPooledRow p K ell B F side)
  let c := globalChildBudget p hp hcop hg pool source Ψmode m windows (Γ*C) ε 0 0 K ell B F M H A J side
  let q := Ct*globalFirstTailCost p s (fun b => weight b) K ell B F M H N side
  have hΨmode : ∀ a,‖Ψmode a‖≤1 := fun a => (firstCoreTwist_norm_le side χ Ψ r a).trans (hΨ a)
  have hweight (b : GlobalFirstData ι) (hb : b∈s) : 0≤weight b :=
    mul_nonneg (hw b hb).1 (globalBinFirstCoefficient_nonneg K ell B F hK.le hell hB hF _)
  have hd : 0≤d := globalFirstDiagonalCost_nonneg p hp hg pool s _ Ψmode m _ rowMajorant ell side _
  have hq : 0≤q := mul_nonneg hCt.le (globalFirstTailCost_nonneg p hp s _ K ell B F M H N side hK hell hB hF hH)
  have hpoint (t : ℝ) : globalFirstFiniteEnergy p hg pool s weight Ψ m g V ell side χ r T t ≤
      1*(d+((1+‖t‖)^(J+2)*(1+‖t‖)^(J+2))*c+q) := by
    have hrow := hrows p hp hcop hg hinj hc hpr pool s w Γ K ell B F H Ψ m χ r t T
      hΓ hw hK hell hB hF hH hs hΨ hb₁ hb₂ hT hT0
    dsimp only at hrow
    have hsmall : globalFirstFiniteEnergy p hg pool s weight Ψ m g V ell side χ r T t ≤
        ∑ b∈s,weight b*((∑ z∈T b,‖globalFirstCoreRow p hg pool Ψ m g V ell side χ r t b z‖^2)+
          globalFirstCoreCorrection p g V ell side t b) := by
      apply Finset.sum_le_sum
      intro b hb
      exact mul_le_mul_of_nonneg_left (le_add_of_nonneg_right (sq_nonneg _)) (hweight b hb)
    have hz := globalFirstSourceZero_bound p hp hg hinj pool s (fun b => weight b) Ψmode m
      g V rowMajorant K ell B F M H side t
    have ht := htail p hp hg hinj hc pool s (fun b => weight b) Ψmode m K ell B F M H t
      hΨmode hK hell hB hF hH hgM
    have hupper := add_le_add (add_le_add (le_refl
      (globalChildBudget p hp hcop hg pool source Ψmode m windows (Γ*C) ε t t K ell B F M H A J side)) hz) ht
    have hfinal := hsmall.trans (hrow.trans hupper)
    rw [globalChildBudget_height] at hfinal
    simpa only [d,c,q,one_mul,add_comm,add_left_comm,add_assoc] using hfinal
  have hint := FirstPassCubeLabels.integral_pair_moment_bound J
    (globalFirstFiniteEnergy p hg pool s weight Ψ m g V ell side χ r T)
    (globalFirstFiniteEnergy_integrable p hg pool s weight Ψ m g V ell side χ r T (2*(J+2)))
    1 d c q (by norm_num) hd hq hpoint
  simpa only [mul_one] using hint

end
section

open scoped BigOperators Classical SchwartzMap
open MeasureTheory
open ActualEisensteinCubic
open FirstPassCubeLabels (primeProductNorm normalizedColumn columnLog firstLogDensity)
open RayFourExpansion (RayCharacter crossCoeff)

@[ext] structure GlobalCubeBlock (ι : Type*) where
  cube : CubeCoordinates ι
  common : Finset ι
  firstDivisor : Finset ι

namespace GlobalCubeBlock
variable {ι : Type*}
def withCommon (b : GlobalCubeBlock ι) (A : Finset ι) : GlobalFirstData ι :=
  ⟨b.cube,A,b.common,b.firstDivisor⟩
end GlobalCubeBlock

namespace GlobalFirstData
variable {ι : Type*}
def block (b : GlobalFirstData ι) : GlobalCubeBlock ι := ⟨b.cube,b.common,b.firstDivisor⟩
@[simp] theorem block_withCommon (b : GlobalCubeBlock ι) (A : Finset ι) : (b.withCommon A).block=b := rfl
@[simp] theorem withCommon_block (b : GlobalFirstData ι) : b.block.withCommon b.firstCommon=b := by cases b; rfl
end GlobalFirstData

variable {ι : Type*} [DecidableEq ι]

def globalCubePool (pool : Finset ι) (b : GlobalCubeBlock ι) : Finset ι := pool\(b.cube.support∪b.common)

def globalFirstFamilySet (pool : Finset ι) (blocks : Finset (GlobalCubeBlock ι)) : Finset (GlobalFirstData ι) :=
  blocks.biUnion (fun b => (globalCubePool pool b).powerset.image b.withCommon)

theorem sum_globalFirstFamilySet {M : Type*} [AddCommMonoid M]
    (pool : Finset ι) (blocks : Finset (GlobalCubeBlock ι)) (f : GlobalFirstData ι → M) :
    ∑ d∈globalFirstFamilySet pool blocks,f d =
      ∑ b∈blocks,∑ A∈(globalCubePool pool b).powerset,f (b.withCommon A) := by
  unfold globalFirstFamilySet
  rw [Finset.sum_biUnion]
  · apply Finset.sum_congr rfl
    intro b hb
    rw [Finset.sum_image]
    intro A hA B hB he
    exact congrArg GlobalFirstData.firstCommon he
  · intro b hb c hc hne
    apply Finset.disjoint_left.mpr
    intro x hx hy
    obtain ⟨A,hA,rfl⟩ := Finset.mem_image.mp hx
    obtain ⟨B,hB,he⟩ := Finset.mem_image.mp hy
    exact hne (congrArg GlobalFirstData.block he).symm

theorem mem_globalFirstFamilySet (pool : Finset ι) (blocks : Finset (GlobalCubeBlock ι)) (d : GlobalFirstData ι) :
    d∈globalFirstFamilySet pool blocks ↔ d.block∈blocks ∧ d.firstCommon⊆globalCubePool pool d.block := by
  simp only [globalFirstFamilySet,Finset.mem_biUnion,Finset.mem_image,Finset.mem_powerset]
  constructor
  · rintro ⟨b,hb,A,hA,rfl⟩;exact ⟨hb,hA⟩
  · rintro ⟨hb,hA⟩;exact ⟨d.block,hb,d.firstCommon,hA,d.withCommon_block⟩

variable (p : ι → O) [∀ i,(Ideal.span {p i}).IsMaximal]
  (hg : ∀ i,lambda ∉ Ideal.span {p i})

def globalCubeColumnScale (ell : ℝ) (side : Bool) (b : GlobalCubeBlock ι) : ℝ :=
  globalFirstBlockColumnScale p ell side (b.withCommon ∅)

def globalCubeCoefficient (K ell B F : ℝ) (side : Bool) (b : GlobalCubeBlock ι) : ℝ :=
  globalFirstCoefficient p K ell B F side (b.withCommon ∅)

omit [∀ (i : ι), (span {p i}).IsMaximal] in
theorem globalFirstCoefficient_block (K ell B F : ℝ) (side : Bool) (d : GlobalFirstData ι) :
    globalFirstCoefficient p K ell B F side d=globalCubeCoefficient p K ell B F side d.block := by
  simp [globalFirstCoefficient,globalCubeCoefficient,globalBinFirstCoefficient,globalFirstIndex,
    globalScaleIndex,globalScaleVector,globalLogRep,GlobalFirstData.append,
    GlobalCubeBlock.withCommon,GlobalFirstData.block,globalCubeActiveRoot]

def globalCubeInputFamilyEnergy (pool : Finset ι) (b : GlobalCubeBlock ι) (Ψ : O →* ℂ) (m : O)
    (g V : 𝓢(ℝ,ℂ)) (ell : ℝ) (side : Bool) (J : ℕ) (T : Finset O) : ℝ :=
  firstInputFamilyEnergy p hg (globalCubePool pool b) b.cube.support
    (fun i => b.cube.leftExponent i+b.cube.rightExponent i) b.cube.leftBit b.cube.rightBit side Ψ m
    (fun A => g (columnLog p (globalCubeColumnScale p ell side b) A)) V
    (globalCubeColumnScale p ell side b) (primeSubsetGenerator (fun i => Ideal.span {p i}) b.common)
    (primeSubsetGenerator (fun i => Ideal.span {p i}) b.firstDivisor) J T

def globalFirstFamilyWeight (w : GlobalCubeBlock ι → ℝ) (Ψ : O →* ℂ) (m : O)
    (side : Bool) (r : FirstCoreIndex) (d : GlobalFirstData ι) : ℝ :=
  w d.block*‖firstCoreOuter p hg d.cube.support
    (fun i => d.cube.leftExponent i+d.cube.rightExponent i) d.cube.leftBit d.cube.rightBit
    side Ψ m d.firstCommon r‖

theorem globalFirstFiniteEnergy_integral_sum
    (pool : Finset ι) (s : Finset (GlobalFirstData ι)) (w : GlobalFirstData ι → ℝ)
    (Ψ : O →* ℂ) (m : O) (g V : 𝓢(ℝ,ℂ)) (ell : ℝ) (side : Bool)
    (χ : RayCharacter) (r : FirstCoreIndex) (T : GlobalFirstData ι → Finset O) (J : ℕ) :
    (∫ t,firstLogDensity J t*globalFirstFiniteEnergy p hg pool s w Ψ m g V ell side χ r T t) =
      ∑ d∈s,w d*∫ t,firstLogDensity J t*∑ z∈T d,‖globalFirstCoreRow p hg pool Ψ m g V ell side χ r t d z‖^2 := by
  have hi (d : GlobalFirstData ι) : Integrable (fun t => firstLogDensity J t*
      ∑ z∈T d,‖globalFirstCoreRow p hg pool Ψ m g V ell side χ r t d z‖^2) := by
    simpa only [globalFirstFiniteEnergy,Finset.sum_singleton,one_mul] using
      globalFirstFiniteEnergy_integrable p hg pool {d} (fun _ => 1) Ψ m g V ell side χ r T J
  have he (t : ℝ) : firstLogDensity J t*globalFirstFiniteEnergy p hg pool s w Ψ m g V ell side χ r T t =
      ∑ d∈s,w d*(firstLogDensity J t*∑ z∈T d,‖globalFirstCoreRow p hg pool Ψ m g V ell side χ r t d z‖^2) := by
    simp only [globalFirstFiniteEnergy,Finset.mul_sum,mul_left_comm]
  simp_rw [he]
  rw [integral_finsetSum _ (fun d hd => (hi d).const_mul (w d))]
  simp only [integral_const_mul]

theorem globalCubeInputFamilyEnergy_eq_integrated_rows
    (pool : Finset ι) (blocks : Finset (GlobalCubeBlock ι)) (w : GlobalCubeBlock ι → ℝ)
    (Ψ : O →* ℂ) (m : O) (g V : 𝓢(ℝ,ℂ)) (K ell B F : ℝ) (side : Bool) (J : ℕ)
    (T : GlobalCubeBlock ι → Finset O) :
    (∑ b∈blocks,w b*globalCubeCoefficient p K ell B F side b*
      globalCubeInputFamilyEnergy p hg pool b Ψ m g V ell side J (T b)) =
    ∑ a : RayCharacter×RayCharacter,‖crossCoeff a.1 a.2‖*∑ r : FirstCoreIndex,
      ∫ t,firstLogDensity J t*globalFirstFiniteEnergy p hg pool (globalFirstFamilySet pool blocks)
        (fun d => globalFirstFamilyWeight p hg w Ψ m side r d*globalFirstCoefficient p K ell B F side d)
        Ψ m g V ell side (if side then a.1 else a.2) r (fun d => T d.block) t := by
  simp_rw [globalFirstFiniteEnergy_integral_sum, sum_globalFirstFamilySet]
  simp only [globalCubeInputFamilyEnergy,firstInputFamilyEnergy,Finset.mul_sum]
  conv_lhs => rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro a ha
  conv_lhs => arg 2; ext b; rw [Finset.sum_comm]
  conv_lhs => rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro r hr
  apply Finset.sum_congr rfl
  intro b hb
  apply Finset.sum_congr rfl
  intro D hD
  rw [globalFirstCoefficient_block]
  simp only [globalFirstFamilyWeight,GlobalCubeBlock.withCommon,
    GlobalFirstData.block,globalFirstCoreRow,globalCubePool,globalCubeColumnScale,globalFirstBlockColumnScale]
  ring

end

open ActualEisensteinCubic
open FirstPassCubeLabels (primeProduct)
open ConcreteTraceCRT (eisEmbedding)
open SecondPassIntegration (densityChildEnergy)

theorem globalChild_large_radial_bound (deltaLoss : ℝ) (hδ : 0<deltaLoss) (η : ℝ) (hη : 0<η) :
    ∃ Cbox : ℝ,0<Cbox ∧ ∀ {ι : Type*} [DecidableEq ι]
      (p : ι → O) (hp : ∀ i,p i ≠ 0) [∀ i,(Ideal.span {p i}).IsMaximal]
      (hcop : Pairwise (Function.onFun IsCoprime (fun i => Ideal.span {p i})))
      (hg : ∀ i,lambda ∉ Ideal.span {p i})
      (_hc : ∀ i,ringChar (O ⧸ Ideal.span {p i}) ≠ 2)
      (_hinj : Function.Injective (fun i => Ideal.span {p i}))
      (pool : Finset ι) (Ψ : O →* ℂ) (m : O) (V₁ V₂ : ℝ → ℂ)
      (C K ell B F M H U Z θ A₀ Vmax : ℝ) (N J : ℕ)
      (source : Finset (GlobalSecondData ι)) (side : Bool),
      0≤C → 1≤U → 0<K → 0<ell → 1≤B → 1≤F → 0≤H →
      K≤U → ell≤U → B≤U → H≤U → Real.exp M≤U → 0<Z →
      0≤A₀ → 0≤Vmax → (∀ t,‖V₁ t‖≤Vmax) → (∀ t,‖V₂ t‖≤Vmax) →
      (∀ t,V₁ t≠0 → |t|≤A₀) → (∀ t,V₂ t≠0 → |t|≤A₀) →
      (∀ a,‖Ψ a‖≤1) → (∀ x∈source,x.source.frequency≠0) →
      (∀ x∈source,‖eisEmbedding (primeProduct p x.cube.support x.cube.leftExponent)‖^2≤B) →
      (∀ x∈source,‖eisEmbedding (primeProduct p x.cube.support x.cube.rightExponent)‖^2≤B) →
      (∑ j∈(globalLogBox (globalNormCaps ell B (globalRowScaleFloor K ell B F M) M H)).filter
          (fun j => Z^θ≤globalBinRadial ell (globalPooledRowScale K ell B F j) j),
        ∑ ray : SecondRayIndex,∑ q∈globalBinTriples p source side j,
          globalDescentWeight C η 0 0 K ell B F (N+1) (2*J) j ray q *
          (densityChildEnergy p hp hcop hg pool (secondRayMinus Ψ ray) (secondRayPlus Ψ ray)
            (fixedTripleMask m q) (globalArithmeticTargets p source side q j) V₁ V₂
            (globalPooledColumnScale ell j) (globalPooledColumnScale ell j) (2*J) /
            (globalPooledColumnScale ell j*globalPooledLabelScale j))) ≤
      C*Cbox*(globalTerminalConstant A₀ Vmax*globalRayTripleConstant)*K*((ell*B^3)*F)*
        U^(deltaLoss+8*η)/Z^(θ*N) := by
  obtain ⟨Cbox,hCbox,hdecay⟩ := globalLogBox_radial_decay deltaLoss hδ η hη.le
  refine ⟨Cbox,hCbox,?_⟩
  intro ι _ p hp _ hcop hg hc hinj pool Ψ m V₁ V₂ C K ell B F M H U Z θ A₀ Vmax N J source side
    hC hU hK hell hB hF hH hKU hellU hBU hHU heU hZ hA₀ hVmax hV₁ hV₂ hVs₁ hVs₂ hΨ hk hb₁ hb₂
  let Q := C*(globalTerminalConstant A₀ Vmax*globalRayTripleConstant)*K*((ell*B^3)*F)
  have hQ : 0≤Q := by
    have hc0 := globalTerminalConstant_nonneg A₀ Vmax
    have hB0 := zero_lt_one.trans_le hB
    have hF0 := zero_lt_one.trans_le hF
    unfold Q globalRayTripleConstant
    positivity
  calc
    _ ≤ ∑ j∈(globalLogBox (globalNormCaps ell B (globalRowScaleFloor K ell B F M) M H)).filter
        (fun j => Z^θ≤globalBinRadial ell (globalPooledRowScale K ell B F j) j),
        Q*((globalBinLossScale B j)^η/(1+globalBinRadial ell (globalPooledRowScale K ell B F j) j)^N) := by
      apply Finset.sum_le_sum
      intro j hj
      have ht := globalRadialBin_terminal_bound_all_bins p hp hcop hg hc hinj pool Ψ hΨ m V₁ V₂
        C η 0 0 K ell B F A₀ Vmax N J hC hη hK hell hB hF hA₀ hVmax hV₁ hV₂ hVs₁ hVs₂
        source side j hk hb₁ hb₂
      apply ht.trans_eq
      simp only [norm_zero,add_zero,one_pow,mul_one,Q]
      ring
    _ = Q*(∑ j∈(globalLogBox (globalNormCaps ell B (globalRowScaleFloor K ell B F M) M H)).filter
        (fun j => Z^θ≤globalBinRadial ell (globalPooledRowScale K ell B F j) j),
        (globalBinLossScale B j)^η/(1+globalBinRadial ell (globalPooledRowScale K ell B F j) j)^N) := by
      rw [Finset.mul_sum]
    _ ≤ Q*(Cbox*U^(deltaLoss+8*η)/Z^(θ*N)) := mul_le_mul_of_nonneg_left
      (hdecay K ell B F M H U Z θ N hU hK hell hB hF hH hKU hellU hBU hHU heU hZ) hQ
    _ = _ := by dsimp [Q];ring

end SecondPassArithmetic

namespace InitialMeanSquare
open ActualEisensteinCubic SecondPassArithmetic
open FirstCauchyArithmetic (supportMobius supportConjugateSum)
open ConcreteTraceCRT (eisEmbedding)

def conjugateMonoid (Ψ : O→*ℂ) : O→*ℂ where
  toFun x := star (Ψ x)
  map_one' := by simp
  map_mul' x y := by simp [map_mul]

variable {ι : Type*} [DecidableEq ι]
  (p : ι→O) [∀i,(Ideal.span {p i}).IsMaximal]
  (hg : ∀i,lambda∉Ideal.span {p i})

def mobiusRow (F : Finset ι) (Ψ : O→*ℂ) (m : O) (H : Finset ι→ℂ) (z : O) : ℂ :=
  ∑S∈F.powerset,supportMobius (fun i => Ideal.span {p i}) S*
    inputColumn p Ψ m 1 H S*finiteSquarefreeRow (fun i => Ideal.span {p i}) hg S z

lemma star_mobiusRow (F : Finset ι) (Ψ : O→*ℂ) (m : O) (H : Finset ι→ℂ) (z : O) :
    star (mobiusRow p hg F Ψ m H z)=
      inputConjugateRow p hg F (conjugateMonoid Ψ) m 1 1 (fun S => star (H S)) z := by
  unfold mobiusRow inputConjugateRow supportConjugateSum
  rw [star_sum]
  apply Finset.sum_congr rfl
  intro S hS
  have hm : star (rowCoprimeMask (fun i => Ideal.span {p i}) S m)=rowCoprimeMask (fun i => Ideal.span {p i}) S m := by
    unfold rowCoprimeMask
    split_ifs <;> simp
  simp only [star_mul,FirstCauchyArithmetic.star_supportMobius,
    inputColumn,secondInputCoefficient,conjugateMonoid,MonoidHom.coe_mk,OneHom.coe_mk,
    mul_one,hm]
  simp only [finiteSquarefreeRow,Finset.prod_const_one,one_pow,map_one,mul_one]
  ring

lemma mobiusRow_norm (F : Finset ι) (Ψ : O→*ℂ) (m : O) (H : Finset ι→ℂ) (z : O) :
    ‖mobiusRow p hg F Ψ m H z‖=
      ‖inputConjugateRow p hg F (conjugateMonoid Ψ) m 1 1 (fun S => star (H S)) z‖ := by
  rw [←star_mobiusRow p hg F Ψ m H z,norm_star]

variable (hp : ∀i,p i≠0) (hinj : Function.Injective (fun i => Ideal.span {p i}))
  (hc : ∀i,ringChar (O⧸Ideal.span {p i})≠2)
include hc

theorem mobiusRow_smoothed_source (F : Finset ι) (Ψ : O→*ℂ) (m : O)
    (H : Finset ι→ℂ) (W : 𝓢(ℝ,ℂ)) (Y : ℝ) (hY : 0<Y)
    (K : Finset ι→Finset ι→Finset O) :
    (∑'z:O,W (‖eisEmbedding z‖^2/Y)*(‖mobiusRow p hg F Ψ m H z‖^2:ℝ))=
      truncatedSecondSource p hp hg hinj F (conjugateMonoid Ψ) m 1 1
        (fun S => star (H S)) W Y (fun G E => (K G E).erase 0)+
      truncatedSecondZero p hg F (conjugateMonoid Ψ) m 1 1
        (fun S => star (H S)) W Y K+
      secondSourceTail p hp hg hinj F (conjugateMonoid Ψ) m 1 1
        (fun S => star (H S)) W Y K := by
  simp_rw [mobiusRow_norm p hg]
  rw [inputConjugateRow_eq_truncated_add_tail p hp hg hinj hc F
    (conjugateMonoid Ψ) m 1 1 (fun S => star (H S)) W Y hY K]
  rw [truncatedSecondSource_eq_nonzero_add_zero p hp hg hinj]

theorem mobiusRow_finite_source (F : Finset ι) (Ψ : O→*ℂ) (m : O)
    (H : Finset ι→ℂ) (Y : ℝ) (hY : 0<Y) (T : Finset O)
    (hT : ∀z∈T,(Ideal.absNorm (Ideal.span {z}):ℝ)≤Y) (hT0 : ∀z∈T,z≠0)
    (K : Finset ι→Finset ι→Finset O) :
    (∑z∈T,‖mobiusRow p hg F Ψ m H z‖^2)≤
      ‖truncatedSecondSource p hp hg hinj F (conjugateMonoid Ψ) m 1 1
        (fun S => star (H S)) rowMajorant Y (fun G E => (K G E).erase 0)‖+
      (truncatedSecondZero p hg F (conjugateMonoid Ψ) m 1 1
        (fun S => star (H S)) rowMajorant Y K).re+
      ‖secondSourceTail p hp hg hinj F (conjugateMonoid Ψ) m 1 1
        (fun S => star (H S)) rowMajorant Y K‖-‖H ∅‖^2 := by
  have hb := finite_nonzero_row_energy_le_majorant
    (inputConjugateRow p hg F (conjugateMonoid Ψ) m 1 1 (fun S => star (H S))) T Y hY hT hT0
  rw [inputConjugateRow_nonzero_correction p hg F (conjugateMonoid Ψ) m 1 1
    (fun S => star (H S)) Y hY] at hb
  have he := mobiusRow_smoothed_source p hg hp hinj hc F Ψ m H rowMajorant Y hY K
  simp_rw [mobiusRow_norm p hg] at he ⊢
  change (∑'z:O,rowMajorant (‖eisEmbedding z‖^2/Y)*
    (‖inputConjugateRow p hg F (conjugateMonoid Ψ) m 1 1 (fun S => star (H S)) z‖^2:ℝ))=_ at he
  rw [he,Complex.sub_re,Complex.add_re,Complex.add_re,Complex.ofReal_re,norm_star] at hb
  exact hb.trans (by gcongr <;> exact Complex.re_le_norm _)

end InitialMeanSquare

namespace CompletedGauss

section
open ActualEisensteinCubic CanonicalQuadraticSieve QuadraticSquarefreeKernel UniqueFactorizationMonoid

def PowerfulIdeal (I : Ideal O) : Prop :=
  I≠0 ∧ ∀P∈normalizedFactors I, 2≤(normalizedFactors I).count P

theorem powerful_squarefree_dvd_squarePart (I : Ideal O) (hI : PowerfulIdeal I) :
    squarefreePart I ∣ squarePart I := by
  have ha := squarePart_ne_zero hI.1
  have hb := (squarefree_squarefreePart I).ne_zero
  apply (dvd_iff_normalizedFactors_le_normalizedFactors hb ha).mpr
  apply Multiset.le_iff_count.mpr
  intro P
  have he := congrArg normalizedFactors (squarePart_sq_mul_squarefreePart I)
  rw [normalizedFactors_mul (pow_ne_zero _ ha) hb,normalizedFactors_pow] at he
  have hc := congrArg (Multiset.count P) he
  simp only [Multiset.count_add,Multiset.count_nsmul] at hc
  have hsf := Multiset.nodup_iff_count_le_one.mp
    ((squarefree_iff_nodup_normalizedFactors hb).mp (squarefree_squarefreePart I)) P
  by_cases hzero : (normalizedFactors (squarefreePart I)).count P=0
  · omega
  · have hm : P∈normalizedFactors I := Multiset.count_pos.mp (by omega)
    have hp := hI.2 P hm
    omega

theorem powerful_exists_square_cube (I : Ideal O) (hI : PowerfulIdeal I) :
    ∃A : Ideal O, A≠0 ∧ A^2*squarefreePart I^3=I := by
  obtain ⟨A,hA⟩ := powerful_squarefree_dvd_squarePart I hI
  have ha : A≠0 := by
    intro hz
    have h := squarePart_ne_zero hI.1
    apply h
    rw [hA,hz,mul_zero]
  refine ⟨A,ha,?_⟩
  calc
    _ = (squarefreePart I*A)^2*squarefreePart I := by ring
    _ = squarePart I^2*squarefreePart I := by rw [←hA]
    _ = I := squarePart_sq_mul_squarefreePart I

def powerfulBase (I : Ideal O) : Ideal O :=
  if hI : PowerfulIdeal I then (powerful_exists_square_cube I hI).choose else 0

theorem powerfulBase_spec (I : Ideal O) (hI : PowerfulIdeal I) :
    powerfulBase I≠0 ∧ powerfulBase I^2*squarefreePart I^3=I := by
  rw [powerfulBase,dite_eq_left hI]
  exact (powerful_exists_square_cube I hI).choose_spec

theorem powerfulBase_injective_on_fiber (I J : Ideal O)
    (hI : PowerfulIdeal I) (hJ : PowerfulIdeal J)
    (hb : squarefreePart I=squarefreePart J) (ha : powerfulBase I=powerfulBase J) : I=J := by
  rw [←(powerfulBase_spec I hI).2,←(powerfulBase_spec J hJ).2,ha,hb]

end

open ActualEisensteinCubic CanonicalQuadraticSieve QuadraticSquarefreeKernel

end CompletedGauss

end

end OAI
