# Handoff: a Lean proof of the 11/12 zero-free half-plane, following OpenAI's companion paper

*Prepared 8 Oct 2026 at the end of a long research session. This file stands alone; the rest of the parent folder is not needed.*

## 1. Goal

Produce a Lean 4 proof of

```lean
theorem riemannZeta_ne_zero_of_eleven_twelfths_lt_re
    {s : ℂ} (hs : (11 / 12 : ℝ) < s.re) : riemannZeta s ≠ 0
```

The proof must follow the argument of OpenAI's companion paper *The Quasi-Riemann Hypothesis* (5 Oct 2026, `papers/paper2.tex`), and it should keep the complexity reasonably small.

Requirements:
- Use only the standard axioms (`propext`, `Quot.sound`, `Classical.choice`).
- No `sorry`.
- Check it Comparator-style, the same way as the existing 7/8 challenge.

**Import boundary (essential).** OpenAI has already formalized Re s > 7/8, and 11/12 follows from it in one line. So the new proof must not import:
- the 7/8 end theorems: `OAI.NumberTheory.DirichletL.Nonvanishing`, or anything whose proof depends on them;
- the 7/8 paper's own argument chain: its probe, detector, row counts, and marked/plain moment inductions.

Infrastructure modules may be reused: Eisenstein integers, residue symbols, reciprocity, Gauss sums, Poisson summation, theta reflection, the quadratic large sieve, Hecke L-functions, Mellin transforms, prime counting. Write a script that computes the transitive import closure of the new development and fails on any denylisted module, and record each allow/deny decision.

**Why bother, given 7/8 is already formal?** The aim is an independent, much shorter formal proof of a different argument. That argument is easier to audit and to read. The companion paper has *not* been refereed or formalized, so the formalization may also surface gaps. If it does, record them honestly and don't work around them silently.

## 2. Materials in this folder

| Path | What it is |
|---|---|
| `papers/paper2.tex`, `paper2.pdf` | The companion 11/12 paper: 3988 lines, 49 pages. **This is the argument to formalize.** |
| `papers/paper78.tex`, `paper78.pdf` | The 7/8 paper, 199 pages. Reference only; the existing Lean code is organised around it. |
| `lean-base/` | Extraction of OpenAI's Lean formalization of the 7/8 theorem. Repo `CoolRmal/quasi-riemann-seven-eighths` @ `6ff3bd5`; its sources are byte-identical to `openai/math` @ `adc7f1241b42e322a6451854ab7e4b4c146bf78a` (see `provenance.json`). |

Details of `lean-base/`:
- Lean `v4.34.1`, with Mathlib pinned in `lake-manifest.json`.
- Dependencies `PrimeNumberTheoremAnd` and `rellich-kondrachov`, carrying the upstream patches in `patches/`.
- `scripts/` holds the setup, dependency-preparation and Comparator scripts. They were written for Linux runners.
- 2924 modules and about 486k lines, all under `OAI/NumberTheory/DirichletL/`.
- The existing challenge is `ComparatorChallenges/QuasiRiemannHypothesis.lean`. Mirror it for 11/12.

The upstream repository is `https://github.com/openai/math` (`lean/` directory).

## 3. The argument of paper2

Each step below gives line ranges in `paper2.tex`. Its outline (lines 206–660, about 10 pages) explains the whole idea; read it first.

| Step | Lines | Content |
|---|---|---|
| Setup | 184–206, 660–680 | K = ℚ(√−3) and O = ℤ[ω]. Primary means ≡ 1 mod 3. χ_n(u) = (u/n)₆ is extended by zero. A_u(D) = Σ_{(n,S)=1} μ(n)ν(n)χ_n(u)W(Nn/D). |
| **Reduction** (sec:reduction) | 660–773 | Prop. `thm:ms`, the mean square Σ_{N(u)≤D^{1+ϑ}} \|A_u(D)\|² ≪ D^{2+ϑ+ε}, implies the extraction bound via rows u = p⁶ with Y/2 < N(p) ≤ Y and Y = H^{1/6}. This needs about Y/log Y primes (Landau's prime ideal theorem). The result is A_1(D) ≪ D^{11/12+5ϑ/12+ε}. The Mellin argument then shows Ŵ(s)/L_K^S(s,ν) is holomorphic on Re s > 11/12, so there is no zero. For ζ, use L_K(s,1) = ζ(s)L(s,χ₋₃) together with L(1,χ₋₃) > 0. |
| **Poisson** (sec:initialization) | 773–1254 | Lemma `arithmetic`: Gauss–Jacobi identities `gj`, `convert1`/`convert2` (μ-absorption), `recip`, `crt-a`, `quotient`. Then Lemma `poisson` (Poisson summation with excluded primes), the dual energy (eq. `energy`), Lemma `remove-exclusions`, and Prop. `poisson-reduction`. |
| **Descent** (sec:descent) | 1254–1632 | Prop. `canonical` (finite descent), the completed sums T, Prop. `R` (completed mean square), Lemma `cube-reduction` (cube inversion), Prop. `transfer`. Ends with the exponent 11/12 at sec:completion. |
| **Reflection** (sec:reflection) | 1632–2217 | Lemma `theta-realization` (cubic theta coefficients); Prop. `lem:reflection` (the theta transformation); Lemmas `reflection-uniformity` and `theta-bounds`; Lemma `quadratic` (Goldmakher–Louvel quadratic large sieve); squarefree rows (Lemma `squarefree-completed`). |
| **Transfer** (sec:transfer-proof) | 2217–2687 | First Poisson transform (Lemma `first-transfer`), the second Poisson identity (Lemma `arithmetic-poisson`), Lemma `second-transfer`. |
| **Appendix** | 2687–3988 | Gauss-sum and reciprocity identities; the cubic theta transformation with fixed ray-class twists: cusp expansions, Patterson's coefficients, following Dunn–Radziwiłł §5 and Appendix A. |

**External inputs that need Lean proofs.** Check `lean-base` first; much of this probably exists there in some form.
- Patterson's cubic theta coefficients and its transformation law. Likely locations: `Reflection/`, `Cusp/`, `Reciprocity/KubotaCharacter.lean`, `Eisenstein/`.
- The quadratic large sieve over K (Goldmakher–Louvel): `QuadraticSieve/`.
- Continuation of Hecke L-functions over K: `Hecke/`.
- Prime ideal counting for K, lower bounds suffice: `PrimeCounting/`.
- Mellin and contour tools: `Mellin/`, `Continuation*.lean`.
- Poisson summation over the Eisenstein lattice: `Poisson/`.
- Gauss sums: `GaussSum/`.
- Dirichlet's L(1,χ₋₃) ≠ 0: in Mathlib.

**Caution.** The 7/8 paper's reflection, Prop. `completed-reflection`, is not stated the same way as paper2's Prop. `lem:reflection`; normalizations and twists differ. Reusing it may need an adapter lemma. Verify exact statements before relying on them.

## 4. Ways to reduce complexity (decide these early)

1. **Target ζ only, so ν is the trivial character over K.**
   - Paper2's argument works one character at a time. The Mellin step needs only a fixed ν, with no supremum over a family of characters (unlike the 7/8 paper's β\*).
   - The auxiliary twists ξ = νη still range over a fixed finite set of ray-class characters. That part cannot be dropped.
2. **Fixed-exponent milestone first.**
   - Prove the mean square with a fixed ϑ = 1/10 and a fixed small ε. That already gives a fixed half-plane: ϑ = 1/10 yields 11/12 + 1/24 = 23/24, plus ε.
   - Only afterwards add the "for every ϑ, ε" quantifiers needed to reach 11/12. This postpones most of the ε-management.
3. **Extraction rows.** Two options; pick whichever has more existing infrastructure.
   - *Prime rows p⁶.* Needs a lower bound ≫ Y/log Y (or even Y^{1−ε}) on the number of primes of K in (Y/2, Y].
   - *Integer rows a⁶.* Uses the exact identity in the 7/8 paper's Lemma `inverse-amplification`, which needs only lattice counting, but adds a supremum-over-scales (Sobolev) step.
4. **Fix the test function.** The contradiction uses one W = y^{−ϱ}φ(y) per putative zero ϱ. So most seminorm and uniformity bookkeeping over W can be specialized.
5. **Constants.** Constants may depend on every fixed datum. Uniformity in the conductor is not needed.
6. **Reuse infrastructure aggressively**, but check each statement against paper2's exact need.

## 5. Suggested milestones

- **M0 — Environment.**
  - No Lean toolchain is installed on this Mac, and about 30 GB of disk is free.
  - Ask the user before downloading: elan, Lean 4.34.1, and the Mathlib cache (several GB).
  - Build the infrastructure subset of `lean-base`.
  - Optionally run the existing 7/8 Comparator challenge to confirm the base is sound.
- **M1 — Dependency map** (`MAP.md`).
  - For every lemma or proposition of paper2, record the matching `lean-base` declarations or mark it NEW, with an estimated size.
  - Mark the import boundary.
  - This is a broad search; fan it out one paper2 section at a time.
- **M2 — Blueprint skeleton.**
  - Write Lean statements of all paper2 lemmas with `sorry`, wired into the final theorem, so the top level compiles.
  - Keep a dependency graph, blueprint-style.
- **M3 — Fill the proofs bottom-up.**
  - Start with the reduction and extraction plus the Mellin step.
  - Then do the Poisson identities, then descent and transfer.
  - Do the theta reflection last, or reuse it.
- **M4 — Verification.**
  - Run the Comparator check of the 11/12 challenge.
  - Run the import-closure audit.
  - Write a short report: total line count, what was reused, and any gaps found in paper2.

## 6. Risks

- Paper2 is unrefereed and unformalized, so expect to find small gaps.
- The ray-class-twisted cubic theta reflection (the appendix) is the largest potential *new* piece if `lean-base` does not already provide exactly what is needed.
- Build time and disk space on macOS; the Comparator scripts assume Linux.
- Keep the import boundary honest throughout. It is easy to pull in a 7/8 module transitively without noticing.
