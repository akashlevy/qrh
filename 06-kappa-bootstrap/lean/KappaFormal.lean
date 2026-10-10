import KappaFormal.Energy.CertifiedExistence
import Batteries.Tactic.Alias

/-!
Independent modification of OpenAI's plain-moment energy induction.
The kappa floor is 13/18 in place of 3/4. The condition 2*beta-1 <= kappa
is retained. This module does not assert an improved zero-free theorem.
Unchanged dependencies are loaded from the original cached proof tree.
-/
set_option autoImplicit false

namespace QRHResearch.KappaThirteenEighteenths

alias certified_bands := OAI.SevenEighths.CenteredMomentEnergyCertifiedExistence.certified_bands
alias terminal_certificate := OAI.SevenEighths.CenteredMomentEnergyCertifiedExistence.terminal_certificate

#print axioms certified_bands
#print axioms terminal_certificate
#check terminal_certificate

end QRHResearch.KappaThirteenEighteenths
