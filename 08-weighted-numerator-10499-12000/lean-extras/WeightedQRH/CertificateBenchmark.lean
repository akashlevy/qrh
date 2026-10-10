import WeightedQRH.Parameters

namespace WeightedQRH

/-- Exact comparison with the lower endpoint of the updated public cubic-root
bracket. This arithmetic statement does not certify the external QRH claim. -/
theorem theta_lt_public_cubic_lower :
    theta < (8749570194 : ℝ) / 10000000000 := by norm_num [theta]

/-- Exact comparison with the newer rational public benchmark supplied in the
research audit. This is only a comparison of real numbers. -/
theorem theta_lt_updated_public_rational :
    theta < (874957019420098946128603850561452983 : ℝ) /
      1000000000000000000000000000000000000 := by norm_num [theta]

end WeightedQRH

#print axioms WeightedQRH.theta_lt_public_cubic_lower
#print axioms WeightedQRH.theta_lt_updated_public_rational
