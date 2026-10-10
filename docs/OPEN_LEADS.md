# Open leads

Leads are ranked roughly by expected value per unit of effort within each section. None of them is established. "Cheap" means days, not weeks.

## Toward a smaller θ (mathematics)

1. **Change the low side.**
   - Every route here, including the weighted numerator, is capped at **13/15** by the unchanged Gram and reflected-row low majorants ([OBSTRUCTIONS](OBSTRUCTIONS.md) §1).
   - Going below needs one of: an arithmetic cancellation that survives the reflection phase π = ρ·ζ_m; a non-radial or angular probe with a newly derived high-side identity; or a different probe altogether.
   - Probes outside the template classified in `04-gram-split/vls/D_probe_redesign.md` are unexplored. Examples: several Poisson variables, Eisenstein or Kubota columns, coupled high sides, signals not from a Davenport–Hasse collapse.
2. **Push the weighted numerator further** ([08](../08-weighted-numerator-10499-12000)).
   - On its own it appears saturated near 0.87489, from a coarse exploratory search.
   - **Untried lever** (`audits/independent_count.txt` item 4): factor the row-independent Euler term D·W = η(p)Np^{−s−w} to extract L_F^S(s + w, η), extending the Euler region across s + w = 1.
   - Combining the weighted numerator with κ-feedback (c > 2/9) is undocumented. **Caution:** κ-feedback is publicly Li Xiang's (and Guo's and Shunia's) mechanism; per the project rule, do not copy Li Xiang's work.
3. **Prove M_κ for some κ < 1** ([02](../02-strengthening), [03](../03-six-sevenths)).
   - M_κ is square-root cancellation on average for Möbius sums over a family of sextic twists *smaller than the sum length*.
   - κ = 6/7 gives 6/7 by Theorem C1; κ = 3/5 gives 3/4.
   - Any of the following would do:
     - a column functional equation for μ;
     - a power-saving second moment of quadratic twists of the cubic theta function at unbalanced length 𝓗^{7/8};
     - a non-positive (sieve or spectral) treatment of the coprime-pair form Q_cop, whose cancellation is visible numerically;
     - a coefficient-sensitive O(D^{2+ε}) bound for the unit j-block of the signed dual ([07](../07-method-barriers)).
4. **A lossless marked energy beyond the wall:** W(w) with w ≳ 0.285, plus the proved refined count, gives 6/7 ([03](../03-six-sevenths) R4). Numerically the energies stay flat to 8–128× past the wall.
5. **A lossless one-slot marked moment** at (r, z, m) = (7/12, 5/12, 1) gives 6/7. It is true under GRH ([03](../03-six-sevenths) R1, C4).
6. **A power saving s₀ > 0 for (V′)** in the Gram split ([04](../04-gram-split)).
   - The saving must hold in the d = 0 block at M ∈ [0.80, 0.86].
   - (V′_LS) with η > 0.066, or LoA_π with θ_π ≲ 1.31, would beat 0.874957.
   - Smaller items there: write out Lemma C at all cusp pairs, so that Corollary D is unconditional; this gains structure, not θ. Try the finite-scale (V_X) double-Poisson route, which avoids metaplectic Voronoi.
7. **A non-Cauchy–Schwarz treatment of the Gram sixth-power-frequency component**, plus M_{8/9}-strength amplification (vertex B in [03](../03-six-sevenths)).

## Formalization

1. **15/16 short proof (cheap).**
   - The formal mean square already holds for every σ > 1/20.
   - Make the literals 11/60, 11/10, 21/10, 137/120 and 23/24 in `QRH/Spine/{Setup,Extraction,ZeroFree}.lean` parametric in ϑ, then take ϑ ↓ 1/20.
   - Estimated 100–150 Lean lines; not yet compiled. [05](../05-lean-23-24)
2. **Parametric descent toward 11/12.** Rebuild on the dead `ShortDraftDescent.integral_geometric_descent`, seeded with `initial_invariant_power_general`. Thread κ through `canonicalRankMargin`, `fixedDepthRank_step` and `canonicalPadding_budget`. MAP.md estimates 4.5–7.5k lines. [05](../05-lean-23-24)
3. **Remove the exclusion mask from the descent state** with a Lean analogue of paper2's Lemma remove-exclusions. This is a cheap test.
4. **Close paper2's SL₂(ℤ)-invariance gap.** paper2 asserts SL₂(ℤ)-invariance of Kubota's θ without proof. The gap reduces to the finite Gauss-sum residue identity (9√3/2)·principalArithmeticResidue h = unramifiedGaussResidue(9h). If that fails, the failure is itself a finding about paper2.
5. **Shrink the theta construction** (30% of the 23/24 code) with a generic Fredholm/Laurent "simple-pole residue" lemma, then re-run reachability.
6. **Formalize what is honest to formalize from the 6/7 work:**
   - Theorem C1 as a conditional statement on the existing Hecke/Mellin infrastructure;
   - the exact rational certificates: 34999/40000, the LP floors, W(4/7), W(3/10), the GRH 6/7 duals, the energy duals.
7. **Shorten the 7/8 proof:** bundle the about 40 repeated hypotheses in the ladder modules; use the generic hybrid sieve; try the decentered majorant on A ≤ 5M/6. Agent-driven proof rewriting: pilot on 2–3 modules before spending a large token budget.
8. **Upstream generic Mellin criteria** (`dirichlet_zero_free_of_mellin_bound`, `inverse_mellin_analytic`, `pair_zero_free_of_identity_gt_two`) to Mathlib or PNT+.

## Verification and housekeeping

1. **Independent sandboxed Linux rerun** of the 10499/12000 three-kernel check (Comparator + NanoDa + con-ron), to the tracker's standard. Human expert review of the weighted manuscript.
2. Delete the 253 stale pruned `.olean` objects in the original fork workspace before any reuse.
3. Re-certify W(2/7) with a real margin, or headline W(3/10). Run `certify.py 0 17/48 23/48 1/6 7/8` once without a timeout; the w = 0 sanity check is unverified.
4. Fix the stale defaults in `02-strengthening/conditional_certificate.py`. Apply the errata listed in [RESULTS](RESULTS.md) §F to the working files, if those files are to be cited.
5. Monitor [tracker PR #6](https://github.com/GTimG/QuasiRiemannTracker/pull/6) for maintainer review.
