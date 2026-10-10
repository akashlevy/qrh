import WeightedQRH.LowDetector.LowNormalized
import WeightedQRH.Parameters

/-! Scope and axiom inspection of the completed low-side analytic port. -/
namespace WeightedQRH

theorem low_signal_identity : (2327 : ℝ) / 12000 = signal theta := by
  norm_num [signal, theta, b]

end WeightedQRH

#check OAI.SevenEighths.WeightedProbeLow.original_normalized_compensatedPhysicalProbe_low
#print axioms OAI.SevenEighths.WeightedProbeLow.original_normalized_compensatedPhysicalProbe_low
#print axioms OAI.SevenEighths.WeightedLowReflection.low_original_completed_energy
#print axioms OAI.SevenEighths.WeightedProbeLow.lowUnselectedMass_decay
