#!/usr/bin/env python3
"""lit-edge calibration: translate the open inputs into literature units.

1. (H1) => pointwise bound D_P(1/2,k) << K^{1/2}.  Each twist has analytic conductor
   C ~ (KL)^2 (wall-spectral §0), so the implied exponent is x = M/(4(M+l)) in C^x
   (convexity 1/4, Burgess-type 3/16, Weyl 1/6).
2. W(w) as a first moment over the family of a product polynomial of length N^2,
   N = K Z^w, i.e. N^2 = K^eta with eta = 2(1 + w/M).  Conrey-Rodgers / Cech:
   0/1-swap terms only for eta < 2; eta > 2 needs 2-swap terms.
3. Conrey-Rodgers error N1^{1/4} D^{1/2} (N1 = polynomial length, under LH) vs main D:
   relative error exponent (eta/4 - 1/2) in units of log D.
"""

geoms = {
    "27/31 geometry (wall-spectral mds_reach)": dict(M=0.860, l=0.204, w=0.064),
    "19/22 optimum (wall-spectral, (H1)+(H2))": dict(M=0.385 + 0.522, l=0.288, w=None),
}
for name, g in geoms.items():
    M, l, w = g["M"], g["l"], g["w"]
    x = M / (4 * (M + l))
    print(f"{name}: M={M:.3f} l={l:.3f}  (H1) => D_P << C^{x:.4f}"
          f"  [convexity 0.25, Burgess 0.1875, Weyl 0.1667]; L = K^{l/M:.3f}")
    if w is not None:
        eta = 2 * (1 + w / M)
        print(f"   W({w}) as a first moment: eta = {eta:.4f} (>2: 2-swap regime);"
              f" CR error/main = D^{eta/4 - 0.5:+.4f}")
