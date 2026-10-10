# Lemma LoA6: Lindelöf on average for sextic Kubota series, and the clean Gram bound (D)

*Working draft, 9 Oct 2026. The status labels are those of `V_derivation.md`:*
- *[R] rigorous given the cited theorems;*
- *[R-] routine but not fully written;*
- *[G] gap.*

## Statement

Notation:
- K = Q(ω);
- c is a fixed finite-order ray-class function (the paper's c(s) = χ_s(b_*)/ξ(s), restricted to a class σ);
- D(w;m) = Σ_{s primary, (s,6)=1} c(s) g̃(s,−m) N(s)^{-w}, with g̃ = g/√N(s), for 0 ≠ m ∈ O;
- D converges absolutely for Re w > 1.

**Lemma LoA6.**
1. D(w;m) continues meromorphically to Re w > 1/2. Its only pole there is a simple pole at w = 2/3, with residue ρ_c(m).
2. For every ε > 0 there are A and R with
   Σ_{0 < N(m) ≤ Q} |D(1/2+ε+it; m)|² ≪_ε Q^{1+ε}(1+|t|)^R.

**Corollary D.** Let A_m(Y) = Y^{-1} Σ_s c(s) W_1(N s/Y)(N s/Y)^{-1/2+iv} g̃(s,−m), the paper's additive polynomial up to normalisation. Put A^sig_m = Y^{-1/3} Ŵ_1(1/6+iv) ρ_c(m). Then for every ε > 0,

  Σ_{N(m) ≍ Q} |A_m(Y) − A^sig_m|² ≪_ε (1+|v|)^{R'} (Q/Y) Q^ε Y^{2ε},

uniformly in Y ≥ 1. In particular the P_a^{1/6} and P_a²/Y terms of Prop. probe-gram both disappear for A^rest := A − A^sig.

*Proof of the corollary from LoA6 [R].* Mellin inversion gives A_m(Y) = Y^{-1}(2πi)^{-1} ∫_{(2)} Ŵ_{1,v}(w) Y^w D(w;m) dw, where Ŵ_{1,v}(w) = Ŵ_1(w − 1/2 + iv) has rapid decay in vertical strips.
1. Shift the contour to Re w = 1/2+ε. The residue at 2/3 is exactly A^sig. The horizontal segments vanish, by polynomial growth of D in vertical strips, which comes from the functional equation and Phragmén–Lindelöf.
2. Cauchy–Schwarz in t gives |A^rest_m|² ≤ Y^{-1+2ε} ‖Ŵ_{1,v}‖_{L¹} ∫ |Ŵ_{1,v}(1/2+ε+it)| |D(1/2+ε+it;m)|² dt.
3. Sum over m and apply LoA6. ∎

## Proof outline of LoA6

**1. Meromorphic continuation, the pole and the functional equation [R-, by citation].**
- D(w;m) is (up to Γ-factors and a power of |m|) the m-th Fourier coefficient at ∞ of the Kubota Eisenstein series on the 6-fold cover of SL_2 over K, twisted by the ray-class data. References: Kubota 1966/69 and the Kyoto lectures; Patterson 1977, §§1–6 for the general n-th order scheme; Brubaker–Bump, *On Kubota's Dirichlet series*, Crelle 598 (2006), whose explicit functional equations hold for all n.
- Continuation and the functional equation w ↔ 1−w, which mixes the finitely many cusps and the finitely many m mod K^{×6}·(units), are standard.
- Poles in Re w > 1/2:
  - The constant-term normaliser is ζ_K(6w−2)/ζ_K(6w−1), up to S-factors. Its only pole in Re w > 1/2 is at 6w − 2 = 2, i.e. w = 2/3 (the theta pole).
  - The point 6w − 2 = 1, i.e. w = 1/2, lies on the boundary.
  - The zeros of ζ_K(6w−1) lie in Re w < 1/2.
- So no GRH is needed (contrast Dunn–Radziwiłł, who must go left of Re = 5/6 in the s-normalisation).

**2. Approximate functional equation [R-].**
- In the m-aspect the analytic conductor is ≍ N(m)(1+|t|)^{O(1)}.
- Evidence: the cubic convexity bound ζ(3s−2)ψ^{(3)}(a,s) ≪ N(a)^{3/4−Re s/2}, Dunn–Radziwiłł / [DdFDS (6.13)], has slope ½·log N(a) per unit of Re s. The sextic case has the same shape.
- Hence, for Re w = 1/2+ε,

  D(w;m) = Σ_s c(s) g̃(s,−m) N(s)^{-w} V(N(s)/X) + Σ_𝔞 ε_𝔞(w,m) Σ_s c_𝔞(s) g̃_𝔞(s,−m) N(s)^{w−1} V_𝔞(N(s)/X'),

  with X X' ≍ N(m)(1+|t|)^{O(1)}. Take X = X' = (N(m)(1+|t|)^{O(1)})^{1/2}.
- Here g̃_𝔞 are the Kubota Gauss sums at the other cusps: the same Gauss sums with a fixed shift of m and fixed ray twists. The factors |ε_𝔞| are ≪ (1+|t|)^{O(1)} on Re w = 1/2+ε.

**3. Mean square of each piece via the sextic large sieve [R-].**
- For squarefree s coprime to m, g̃(s,−m) = χ̄_s(−m) g̃(s). Write m = m_0 m_1^6 with m_0 sixth-power-free; then χ_s(m) = χ_s(m_0)·1_{(s,m_1)=1}.
- The pieces with (s,m) ≠ 1 or s non-squarefree are reduced to this case by the triangular relation of `V_derivation.md` (G1). It costs only divisor factors.
- Then

  Σ_{N(m)≤Q} |Σ_{N(s)≤X} b_s χ̄_s(m)|² ≤ Σ_{N(m_1)≤Q^{1/6}} Θ_6(Q/N(m_1)^6, X) Σ |b_s|² ≪ Q^{1+ε} Σ_{N(s)≤X} N(s)^{-1-2ε},

  with |b_s| ≤ N(s)^{-1/2-ε}. This uses de Faveri, *Optimal large sieve for fixed order characters*, arXiv:2610.04045, Thm 1.1, with n = 6:

  Θ_6(A,B) ≪ (AB)^ε (A + B + A^{5/6}B^{1/3} + A^{1/3}B^{5/6}).

- At A = Q, B = Q^{1/2} every term is ≤ Q: A^{5/6}B^{1/3} = Q and A^{1/3}B^{5/6} = Q^{3/4}.
- Blomer–Goldmakher–Louvel's older bound A + B + (AB)^{2/3} = A + B + Q also suffices.
- The dual pieces are handled identically.

**Conclusion.** LoA6 holds with R = O(1). This is the same argument as David–de Faveri–Dunn–Stucky (arXiv:2410.03048, Lemma 6.6) for n = 3. De Faveri–Dunn–Hoffstein (arXiv:2607.07911, §2.2) say they "expect similar arguments to work for all n ≥ 5".

## Gaps
- **[G-a]** Writing out the sextic functional equation with explicit cusp data. Brubaker–Bump give it in the WMDS language, and it has to be translated into Kubota-series form. Exponents are unaffected.
- **[G-b]** The (s,m) ≠ 1 and non-squarefree reductions; divisor-bounded losses only.
- **[G-c]** The ray-class restriction s ∈ σ is a finite linear combination of ray characters c. Each gives an Eisenstein series with a character. This is fine.

## Numerical support
`exp_rank1.py`, at Q = 3000 and Y = 256 … 65536:
- the diagonal Gram fits G_ii·Y = 53.9 + 0.258·Y^{1/3};
- the cross-scale estimate gives κ ≈ 0.27, so the residual ‖A − A^sig‖²·Y ≈ 54 is flat;
- this holds even where the paper's P²/Y term would be up to 10^7 times the diagonal.
