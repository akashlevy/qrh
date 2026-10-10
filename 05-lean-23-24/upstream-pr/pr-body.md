> Opened on the fork because openai/math does not accept pull requests from outside collaborators. The base branch is the fork's `main`, which is identical to openai/math@fd4aeeb, so this diff is the diff against upstream. Compare view on upstream: https://github.com/openai/math/compare/main...akashlevy:math:qrh-23-24-pruned

This PR is a set of diffs for discussion, not for merging. It shows how the quasi-Riemann hypothesis (QRH) material could be made smaller and easier to read, using things the repository already contains:

1. **Lean:**
   - a minimal formal proof of a QRH for every Dirichlet L-function, L(s, χ) ≠ 0 for Re s > 23/24 (hence ζ(s) ≠ 0 there), stated step by step in the paper's numbering;
   - an in-place pruning of the QRH Lean development: every declaration that no QRH result uses is deleted, about 114k code lines (26%).
2. **Paper:** `paper2.tex` (the 11/12 companion) is edited in place into the 23/24 argument for Dirichlet L-functions that the Lean proof formalizes. The diff of `paper2.tex` is the complete list of changes.

## Why

**The Lean proof of a QRH can be much smaller than the 7/8 development.** The library already contains unconditional 23/24 theorems, for every Dirichlet L-function (`OAI.ShortDraft.dirichlet_LFunction_ne_zero`) and for ζ (`OAI.ShortDraft.riemannZeta_ne_zero_of_re_gt`). Nothing used them, so they were easy to miss.
- It is paper2's argument at the fixed parameter ϑ = 1/10, extracting from integer rows b⁶ instead of prime rows.
- Its import closure is 286 modules, against 2,924 for the 7/8 theorem.
- It shares no module with the 7/8 argument chain (probe, detector, row counts, moment inductions), only the common infrastructure.

Sizes of the local (non-Mathlib) import closures, in code lines (non-blank, non-comment), before and after this PR:

| | modules | code lines before | code lines after | change |
|---|---:|---:|---:|---:|
| 7/8, ζ (`DirichletL.Nonvanishing`) | 2,924 → 2,671 | 441,096 | 327,627 | −25.7% |
| 23/24, Dirichlet and ζ (`DirichletL.QRH`, new) | 292 | 173,227 | 122,571 | −29.2% |
| the whole `DirichletL` directory | 2,932 → 2,679 | 441,401 | 327,932 | −25.7% |

After pruning, the 23/24 proof is 37% of the size of the 7/8 proof, and 28% of the size of the 7/8 proof as it is on `main`. It is also a different, shorter argument. Most of what pruning deletes are superseded variants of live lemmas (alternative transfer and log-separation lemmas, unused extraction routes) and whole files of the 7/8 chain that its final proof no longer uses: 95 modules in `Moments`, 67 in `Descent` and 42 in `Reflection` become empty and are removed.

**The paper and the formal proof should line up.** paper2 proves 11/12 for all Hecke L-functions, with prime rows and ϑ → 0. The formal 23/24 proof is a special case organized differently at the top level, so a reader cannot match the two. In this PR:
- The paper is specialized to what the Lean proof proves: the Hecke characters ν = χ∘N that come from Dirichlet characters χ, ϑ = 1/10, and extraction from all rows u = b⁶ with a weight φ(𝔫)/N𝔫.
- Section 3 is rewritten as three propositions that are exactly the three Lean theorems the main theorem chains.
- `Milestones.lean` names the results of the paper's §§4–7 in the paper's numbering. The numbers cited in the docstrings were checked against the compiled paper.
- Two remarks, and Appendix C, record where the Lean proof uses different bookkeeping:
  - height-twist expansion of test functions;
  - exclusion masks;
  - a compactly supported row majorant;
  - a fixed cube cutoff D^{1/1000} with at most 3000 descent levels.

## Lean, part 1: the 23/24 proof (commits 2 and 4)

- `OAI.DirichletCharacter.LFunction_ne_zero_of_twenty_three_twenty_fourths_lt_re`, for every Dirichlet character χ of every modulus. It has exactly the shape of the 7/8 Dirichlet theorem, with 23/24 in place of 7/8: `(23/24 : ℝ) < s.re → ¬(χ = 1 ∧ s = 1) → L(s, χ) ≠ 0`. Comparator challenge: `ComparatorChallenges/DirichletTwentyThreeTwentyFourths.{lean,json}`.
- `OAI.riemannZeta_ne_zero_of_twenty_three_twenty_fourths_lt_re : (23/24 : ℝ) < s.re → riemannZeta s ≠ 0`, the corollary for ζ, with no `s ≠ 1` hypothesis. Comparator challenge: `ComparatorChallenges/QuasiRiemannTwentyThreeTwentyFourths.{lean,json}`.
- Commit 2 proved the ζ case. Commit 4 generalizes the same steps to an arbitrary Dirichlet character. Every library lemma involved already took a general character, so this costs 23 code lines and needs nothing that the pruning deleted.
- `OAI/NumberTheory/DirichletL/QRH/`:
  - `Setup`: A_u(D), A^φ(D), S, and the two bounds passed between steps.
  - `MeanSquare` (Prop 3.1 = Prop 4.5 ∘ Prop 5.1).
  - `Extraction` (Prop 3.2).
  - `ZeroFree` (Prop 3.3).
  - `Milestones`: aliases for Lemma 4.1, Lemma 4.2, Props 5.2–5.4, Prop 6.2, Lemmas 6.4–6.5 and Lemma 7.2.
- The main theorem is the same three-line composition as the paper's proof of Theorem 1.1.
- `formalization.yaml` and `docs/003.md` list the new result.

## Lean, part 2: pruning (commit 3)

The pruning is computed inside Lean:
1. **Roots:** the results that must keep working — 7/8 for ζ, Dirichlet and Hecke, and the new 23/24 theorems — plus the whole module `DirichletL/CubicDyadicDecay`, the only one in this directory that code outside it imports (`NumberTheory/CubicMoment`).
2. **Reachability:** compute every constant reachable from the roots through types and proof terms.
3. **Deletion:** delete every top-level declaration in `OAI/NumberTheory/DirichletL` that is not reachable: 8,648 of 29,407. The 253 modules left empty are removed.
4. **Kept anyway:** closed `rfl`/`simp` lemmas, which tactics use without leaving a trace in proof terms, and notation.
5. **Cleanup:** remove what the deletions leave behind:
   - `variable`, `include`, `omit` and `local notation` commands that nothing can see anymore;
   - `open` atoms naming namespaces that no longer exist;
   - pruned names in `open NS (..)`.

Five modules are kept byte-identical even though they contain unreachable declarations: `ThetaProduct`, `EisensteinTheta`, `Hecke.Theta`, `Hecke.Family` and `Hecke.IdealBridge`. The reason:
- `ComparatorChallenges/HeckeSevenEighths.lean` restates parts of them verbatim, and Comparator requires the restated constants to equal the library's exactly, proof terms included.
- Deleting an earlier unreachable declaration renames the auxiliary lemmas (`foo._proof_1`) that Lean shares between proofs. The proofs stay valid but no longer match: the first Comparator run of `HeckeSevenEighths` failed on exactly this, for `ThetaProduct.rectangularPair_k`.

The old `OAI.ShortDraft` 23/24 assembly is unreachable once the new theorems replace it, so it is among the deletions. The new theorems cover the same statements: every Dirichlet L-function, and ζ. `lean/scripts/treeshake/` contains the tools and the exact commands, so the deletion can be regenerated and audited.

## Paper changes (commits 1 and 5)

The paper is still paper2: 89.1% of the new text is paper2 verbatim, 2.9% of paper2 is removed, and `paper2.tex` changes by +284 / −143 lines. Kept:
- the introduction;
- prior work;
- the consequences;
- Appendices A (arithmetic identities, theta transformation) and B (smooth weights), unchanged.

Changed:
- **Dirichlet characters, ϑ fixed:** Theorem 1.1 is stated for every Dirichlet L-function, and ν = χ∘N throughout, so paper2's ν text is unchanged. ϑ = 1/10 is fixed, and the exponent is 11/12 + 5ϑ/12 = 23/24.
- **Extraction (Step 2 of the outline, §3):** it uses the rows b⁶, so no prime ideal theorem is needed. The Mellin step uses the base change L_K(s, ν) = L(s, χ)·L(s, χχ₋₃), with its Euler factors checked prime by prime, and an Euler product E_φ(s, ν) that is holomorphic for Re s > 0. The identity theorem is applied away from s = 1.
- **Consequences:** Corollary 1.2 (primes in progressions) and the further consequences are paper2's, with 23/24 in place of 11/12.
- **New:**
  - §1.2 (the formalization);
  - Remarks 4.6 and 5.5 (the formal bookkeeping);
  - Appendix C (paper ↔ Lean table).

## Verification

- **Build:** a standalone Lake workspace containing only the QRH-family modules, with the repository's Lean (4.34.1) and Mathlib pins (the rest of the library's dependencies were not fetched). The whole pruned directory builds: all 2,679 modules, and the Comparator challenges.
- **Outside importers:** no module outside `DirichletL` imports it except `NumberTheory/CubicMoment`, through `CubicDyadicDecay`, which is untouched.
- **Comparator** (same pinned commits as the repository; on macOS, run with an unsandboxed landrun shim) accepts all five QRH-family challenges after pruning, with only the standard axioms:
  - `QuasiRiemannHypothesis` (7/8, ζ);
  - `DirichletSevenEighths`;
  - `HeckeSevenEighths`;
  - `DirichletTwentyThreeTwentyFourths` (new);
  - `QuasiRiemannTwentyThreeTwentyFourths` (new).
- **No escape hatches:** there is no `sorry`, `axiom`, `native_decide` or `implemented_by` in `DirichletL`.
- **Import closure:** the 23/24 closure contains no module of the 7/8 argument chain, checked by script.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
