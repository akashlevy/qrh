

# KERNELS

## summary
Every accepted "framework-verified" kernel dossier in the tracker came from the same pipeline. The source was built with Lean 4.34.1 and exported with lean4export 076e8e5. Those exports were then checked on an x86_64 Linux "verification worker" by PalomarSubmission's unmodified scripts/verify_submission.py at d4e41c1 (MIT). The checker and all three kernels come from the official Lean v4.35.0-rc2 Linux release: `lake comparator`, `leanchecker` (Lean's fresh kernel, logged as "Lean default"), `nanoda_bin` and `con-ron`. There are three such dossiers: public/proofs/palomar-20261009/{qrh,openai}, public/proofs/nielstron-20261009-kernels (merged PR #1), and public/proofs/nielstron-algebraic-20261009-kernels (open PR #2, head 0284b10d, branch nielstron-algebraic-endpoint-20261009; the fork nielstron/QuasiRiemannTracker has only that branch, main and the PR #1 branch).

The tracker's own Dockerfile/worker pipeline is a different thing. It uses comparator d03acab, nanoda_lib 3a24072 and landrun, is described as "uncommissioned", and produced none of the accepted evidence.

What CONTRIBUTING.md requires for "framework-verified":
- The exact claimed statements and their definitions pass Comparator.
- The exported proofs are accepted by Lean, NanoDa and con-ron.
- Only propext, Classical.choice and Quot.sound are permitted.
- Sources, checker pins, exact statements and acceptance reports are kept for maintainer review.
- core/catalogue.mjs also demands a parseable first_verified_at for this status. The per-entry validators require first_verified_at == result.json verified_at_utc and a matching source_commit.

PR #2 sets the precedent for a new entry: it adds its own validator (core/algebraic-kernel-evidence.mjs, plus scripts/check-algebraic-kernels.mjs and tests/algebraic-kernel-evidence.test.mjs). That validator hard-codes the Linux x86_64 tool hashes and the '/tools/lean-4.35.0-rc2-linux/bin/...' kernel paths.

Faithful reproduction on this Mac is blocked:
- `lake comparator` refuses to run on non-Linux. Its source says it "sandboxes the code it checks with `bwrap`, which needs Linux namespaces".
- The Palomar scripts and the prior drivers assert sys.platform == 'linux'. They need bubblewrap 0.12.0 and cgroup-v2 delegation (systemd-run, or `sudo -n` with cgroup_delegate.py).
- `docker info` currently returns a 500 Internal Server Error, so the Docker daemon is not usable.
- Only 26 GiB of disk is free (94% used). The host has 36 GiB RAM and an M3 Pro (12 cores).

Two big risks for the akashlevy proof:
1. In the fork c388e88e6, Hecke/Family.lean and Hecke/IdealBridge.lean match the pinned upstream hashes (9871d175..., 6c1fc33b...). But 300 of the 306 OAI modules imported (directly or indirectly) by IdealBridge are modified compared with openai/math fd4aeeb2, by commit 817ddb5e6 "Prune the quasi-Riemann development". Earlier dossiers tied the statement definitions to pristine pinned upstream files, so the Challenge should be compiled against pristine openai/math fd4aeeb2 plus Mathlib d13f23b7. That way Comparator checks the fork's definitions against upstream; if the pruning changed any definition, Comparator will reject.
2. The earlier local Comparator runs in formalize_11_12/logs used comparator d03acab with the landrun shim (no sandbox) and ran only the Lean default kernel. NanoDa and con-ron were not run.

Copies of the files I fetched for reference: /private/tmp/claude-501/-Users-akashlevy-Documents-riemann/5258226a-75a9-431d-830a-875c5394bd2d/scratchpad/pr2/ (PR #2 dossier and validator) and .../scratchpad/palomar/ (PalomarSubmission@d4e41c1: README.md, scripts/verify_submission.py with sha256 c575759d... matching the pin, install_bwrap.sh, profiles).

## nanoda
What the accepted dossiers used: the nanoda_bin bundled in the official Lean v4.35.0-rc2 release, not a separate build.
- Release archive: lean-4.35.0-rc2-linux.tar.zst, sha256 b67d0e0ba959fc74e186f1b44aebca32d6ecb6a28f8e9a029fc4e9f73398a049, 595 MB. It was fetched through the GitHub release API and checked against the asset digest (palomar_prepare_remote.py).
- Binary sha256: 8241c5e6baa29490aa118d382f97961b208ae730fdb96545c2abed74c9ece5a8. Pinned in tool-pins.json, checker-binary-pins.json and the TOOLS table of core/native-kernel-evidence.mjs.

Where that binary comes from: Lean v4.35.0-rc2 src/CMakeLists.txt has ExternalProject_Add(nanoda) with GIT_REPOSITORY https://github.com/ammkrn/nanoda_lib and GIT_TAG 4c544ed4099c8227f07d5de77ad1e69fb0740a27 (2026-09-09, merge of PR #31 "flt_perf"; Cargo.toml says version 0.4.17).
- Build: `cargo build --release --locked`, then target/release/nanoda_bin is copied into bin/.
- On standalone Linux builds it is linked through leanc via CARGO_TARGET_<arch>_LINKER.
- Release CI step "Build nanoda and con-ron" runs `make -C build/$TARGET_STAGE nanoda con-ron` on every non-wasm release platform, so the darwin_aarch64 and linux_aarch64 archives also contain it.

How it is invoked: Palomar's protected comparator.json registers it as "nanoda": ["<judge_prefix>/bin/nanoda_bin"] (PROTECTED_KERNELS in verify_submission.py). Any `external_kernels` or `enable_nanoda` in a submitted config is ignored or rejected.
- `lake comparator` (Lake/CLI/Check.lean runExternalKernel) writes a temporary config: {"use_stdin":false, "export_file_path":<solution export>, "permitted_axioms":["Classical.choice","Quot.sound","propext"], "unpermitted_axiom_hard_error":false, "num_threads":4, "nat_extension":true, "string_extension":true}.
- It appends that config path as the last argument and runs nanoda_bin inside a nested bwrap. The log shows "Running nanoda kernel on solution" followed by "nanoda kernel accepts the solution".
- Comparator itself does the axiom check.

Not used for any acceptance:
- The tracker's uncommissioned worker pins nanoda_lib 3a2407216ee84a75f9e1aead6803d0578be06ae7 (v0.4.19, 7 commits ahead of 4c544ed). It is built in rust:1.90.0-bookworm@sha256:3914072ca0c3b8aad871db9169a651ccfce30cf58303e5d6f2db16d1d8a7e58f with `git checkout --detach ... && cargo build --release --locked` and run via comparator d03acab with COMPARATOR_NANODA=/opt/bin/nanoda_bin (worker/Dockerfile, verifier/pins.json).
- The standalone comparator d03acab passes nanoda the same config but with "unpermitted_axiom_hard_error": true.

## con_ron
What the accepted dossiers used: the con-ron bundled in Lean v4.35.0-rc2. Linux x86_64 binary sha256 4e5616d94374cae37324dad2594f6230c2bb2297240249fc78ff508f3bc5acef.
- Upstream: https://github.com/leanprover/con-ron (Apache-2.0, "A rust port of con-leche"), pinned in Lean v4.35.0-rc2 src/CMakeLists.txt at GIT_TAG cb54c9a4e5f604a678b0230f96b8358b9fd73c82 (2026-09-13).
- Build: GIT_SUBMODULES "" (the Aeneas submodule is not needed), `cargo build --release --locked` (a workspace of crates/con-ron-core, con-ron-dump and con-ron; release profile has overflow-checks=true; uses mimalloc), then target/release/con-ron is copied into bin/.

How it is invoked: Palomar registers "con-ron": ["<judge_prefix>/bin/con-ron", "--jobs=2"] (CON_RON_JOBS=2, chosen for memory: the comparator phase used about 10.7 GiB on a 734 MB export). `lake comparator` appends the solution export path and runs it in a nested bwrap.
- Acceptance log lines: "con-ron: accepted 122217 declarations (--verified)" (N24 run) and 122304 (algebraic run), then "con-ron kernel accepts the solution".
- The positive control logged "accepted 466 declarations".
- On the ill-typed control it printed "rejected: type mismatch in declaration [at thm palomar_preflight ...]" and exited with 1.
- Every run also prints the informational line "con-ron: 1 inductive blocks modelled in-process: Lean.Syntax ...".

The tracker's own uncommissioned worker does not run con-ron at all.

## lean4export_and_comparator
Source compiler and exporter:
- Lean 4.34.1, Linux x86_64, commit 5045d0056413266e57c625dcd7c365b10e377c52, lean binary sha256 e8baaa71855a616dc351028f3ad2200051b0671f423a1696a100e809302d5550.
- Exporter: https://github.com/leanprover/lean4export @ 076e8e57707e813375e8f9da8bf989799ace9680. Built by `git init; git fetch --depth 1 https://github.com/leanprover/lean4export.git 076e8e5...; git checkout --detach FETCH_HEAD; <lean-4.34.1>/bin/lake build lean4export`. Binary .lake/build/bin/lean4export, sha256 8c5d64ba68f4d3a68b3bcb7e1c5f188e35ae573b29470b8d4530a62324c9ddd1.
- The bundled 4.35.0-rc2 `leanexport` is used only for the preflight controls, because the export format must match the compiler that produced the .olean files.

Export procedure (the drivers verify_candidate.py and verify_algebraic.py, through Palomar's v.sandboxed_run and v.export_module):
1. Write the wrappers: src/PalomarCanonical<24hex>/Challenge.lean and src/Solution.lean.
2. Compile each one: `<lean-4.34.1>/bin/lean -j16 -R <src> -o <lib>/<Module>.olean <src>/<Module>.lean`.
   - The Challenge's LEAN_PATH is <lib> plus a symlinked dependency view with every QRH* entry removed.
   - The Solution's LEAN_PATH is <lib> plus the full candidate build.
3. Export each module: `lean4export <Module> -- <targets>`, with stdout captured to a file outside the sandbox.
   - The targets (export-targets.json) are Quot, Quot.mk, Quot.lift, Quot.ind, the theorem names, Classical.choice, Quot.sound, propext, then the primitives Nat.add, Nat.sub, Nat.mul, Nat.pow, Nat.gcd, Nat.div, Nat.mod, Nat.beq, Nat.ble, Nat.land, Nat.lor, Nat.xor, Nat.shiftLeft, Nat.shiftRight, String.ofList, Char.ofNat, List, eagerReduce, Nat, String, String.mk, Char, optParam, autoParam, semiOutParam, outParam.
   - The primitive list is read from 4.35.0-rc2 src/lean/lake/Lake/CLI/Check.lean primitiveTargets.
4. v.verify_export requires the first line to start with {"meta":.
- Earlier export sizes: challenge 382 MB, solution 1.70 GB. The exports are not published; their sha256 and size go into result.json.

Judge:
- Lean v4.35.0-rc2, linux x86_64 release sha256 b67d0e0b... Binary sha256s: lake 8ba83c98f76562d03b8dbc33f1057151ba6cbb121679dde76b2260ccc90495e1, lean bf8d54e4714cc4b03d3f6bb34c83b7202b87e49c8bfcbff6895c085bb90ceb38, leanexport c5bc1a10a21e22cbb3671d65987cdef2076833b14182d823ccb2153adf2fa954, leanchecker 9e36e955a251b1313c8eb5b7fa101282d51df109d34d8bb7c1cfdf51985aba29.
- Command: `<4.35.0-rc2>/bin/lake comparator --config comparator.json --challenge-from-export challenge.export --solution-from-export solution.export`.
  - Environment: PATH=<4.35.0-rc2>/bin:/usr/bin:/bin, LEAN_ABORT_ON_PANIC=1, COMPARATOR_BWRAP=<bwrap>, fresh HOME and TMPDIR.
  - It runs under Palomar sandboxed_run(nested_sandbox=True). The kernels run in order con-ron, nanoda, then Lean default (`leanchecker --silent --from-export`).
- Success requires exit code 0, the "<kernel> kernel accepts the solution" line for all three kernels, and "Your solution is okay!". Verdicts are parsed by v.comparator_verdict.

Protected config, written by v.protected_comparator_config from requested-comparator.json:
- Keys: challenge_module "PalomarCanonical<secrets.token_hex(12)>.Challenge", solution_module "Solution", theorem_names [...], definition_names [], permitted_axioms ["Classical.choice","Quot.sound","propext"], external_kernels {con-ron: [..., "--jobs=2"], nanoda: [...]}.

Palomar pins:
- https://github.com/PalomarRegistry/PalomarSubmission @ d4e41c1d5b0d114c4859e6e5831dc6d3ad1d0d44 (2026-10-05). File sha256s: verify_submission.py c575759d82c506f181b3736a110950fd32931c08a1be454262027d6205f178a2, supervise_cgroup.py 3cf9af4fa405e0e2e1bf1e06abe4b53887c1e4a7d5ebe9a6f20b5cf6a66f0e58, cgroup_delegate.py c8a44ea58ed9a94db2caf21c66cee3f664c86370fabd50ba56e22fa9d9345824.
- It is imported unmodified (`from scripts import verify_submission as v`). The drivers only call v.configure_bwrap and override VERIFICATION_LIMITS (memory_high_percent 40, memory_max_percent 48, tasks_max 4096) and the execution deadline (19800 s).

bubblewrap:
- v0.12.0, built from https://github.com/containers/bubblewrap/releases/download/v0.12.0/bubblewrap-0.12.0.tar.xz (sha256 9760d007363e3abba7c747489910f9f82d9fca53ba3bd3282e396fa3c97a3314) with meson and ninja (Palomar scripts/install_bwrap.sh). Binary sha256 bb807d18eaee5dad15afd5cf48c8de1c3206ff97e852af693d679d62394eb5a2.
- Root copied it to /opt/qrh-verifier-20261009/bwrap and loaded a temporary AppArmor profile 'qrh-pinned-bwrap' (abi 4.0, flags=(unconfined), `userns,`), because Ubuntu 24.04 restricts unprivileged user namespaces.
- The earlier ProofCouncil/OpenAI runs used the system /usr/bin/bwrap (sha256 6da06f15...).

Standalone comparator:
- https://github.com/leanprover/comparator @ d03acab154d269c06e60e4de7e4cc85deebff94b, with lean-toolchain changed from 4.34.0 to 4.34.1, and lean4export 076e8e5 as a Lake package.
- This is only the tracker's own signed-admission pin (verifier/pins.json, worker/Dockerfile: `lake build comparator lean4export`, run as `lake env /opt/bin/comparator /input/comparator.json` with COMPARATOR_LANDRUN, COMPARATOR_LEAN4EXPORT and COMPARATOR_NANODA set). Landrun is Zouuup/landrun @ 811cfff51ceaf3d9843708aa6d22e9b84ccac8b4.
- It matches the local build at /Users/akashlevy/Documents/riemann/formalize_11_12/.tools/comparator (HEAD d03acab; lean4export package at 076e8e5; that working tree shows ` M lean-toolchain`).
- The tracker says this path is uncommissioned, and it produced none of the accepted evidence.

## palomar_canonical_challenge
"PalomarCanonical<hex>" is not a fixed published file name. It is a random per-run module alias: protected_challenge_module returns f"PalomarCanonical{secrets.token_hex(12)}.Challenge". Palomar compiles the trusted Challenge under that namespace so that a Solution cannot capture it.

The statement text itself is the fixed three-target template in public/proofs/palomar-20261009/qrh/src/Challenge.lean (sha256 eee59775c5f4f892199f7456f85f7664a14ef3010fc1b3df97ed2c0f8a1db4e6). The OpenAI 7/8 copy is identical except that the literal is (7 / 8 : ℝ). Exact text:

```lean
import Mathlib.NumberTheory.LSeries.DirichletContinuation
import OAI.NumberTheory.DirichletL.Hecke.IdealBridge

namespace QRHPalomar
open scoped _root_.DirichletCharacter

theorem allDirichlet {q : ℕ} [NeZero q] (χ : _root_.DirichletCharacter ℂ q) {s : ℂ}
    (hs : (874957019421 / 1000000000000 : ℝ) < s.re)
    (hpole : ¬ (χ = 1 ∧ s = 1)) : _root_.DirichletCharacter.LFunction χ s ≠ 0 := by
  sorry

theorem zeta {s : ℂ} (hs : (874957019421 / 1000000000000 : ℝ) < s.re) : riemannZeta s ≠ 0 := by
  sorry

theorem allHecke (χ : OAI.SevenEighths.HeckeFamily.Character) (s : ℂ)
    (hs : (874957019421 / 1000000000000 : ℝ) < s.re)
    (hpole : s ≠ 1 ∨ χ.residue ≠ 1) : OAI.SevenEighths.HeckeFamily.LFunction χ s ≠ 0 := by
  sorry
end QRHPalomar
```

How prior entries used it:
- N24 (public/proofs/nielstron-20261009-kernels/src/PalomarCanonical852af299eac8f1607f3d35b3/Challenge.lean) is the template with the literal replaced in all three places. core/native-kernel-evidence.mjs enforces this byte for byte ("Independent Challenge differs from the inherited fixed three-target template").
- PR #2 (src/PalomarCanonical3d11d053040a6a513a35e5ef/Challenge.lean) uses the same three statements with a new literal, adds four more theorems (allDirichletExact, zetaExact, allHeckeExact, existsUniqueRoot) and a doc comment, and is checked by its own validator.

Solution wrapper: `import <entry module>`, then the identical statement text with each `sorry` replaced by `exact <decl> ...`. The OpenAI Solution turns its Hecke pole hypothesis into `s ≠ 1 ∨ χ.residue ≠ 1` with `apply ...; intro h; rcases hpole with hne | hne`.

For theta = 10499/12000 the literal would be `(10499 / 12000 : ℝ)`. The fraction is already reduced; the template's three statements cover zeta, all Dirichlet L-functions and the OAI finite-order Eisenstein Hecke family, matching the claimed result. Replacing the literal in all three places in the template gives text with sha256 4624fc67bbd5792408525558a48481a5cd3acd258c78caf1ec4d29528703c287; I computed this in memory and wrote no file.

Must a submission use it? No rule applies it automatically to a new ID. But CONTRIBUTING requires that "the exact claimed statements and their definitions pass Comparator", and both precedents reuse this template with only the literal changed. Doing the same is the safe and reviewable choice.

Statement-provenance checks in earlier dossiers (challenge-source-audit.json / verify_algebraic DEFINITION_PINS):
- mathlib d13f23b723b8a846827a245b89c10fc7d3f11612: DirichletContinuation.lean 7ef38fb53f462e19dfd6746fc1826423e55bb5e2672c2b5667e2fd4889f401b6 and RiemannZeta.lean 8d24b729a30f1fef0820f0f7c9cdbb5910df30dfac2c597eb666794bead08b1c.
- openai/math fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb: Hecke/Family.lean 9871d175749b34ec40cf66ccc14b4faa0a9ef08831cb3fa901c944e57c7b9ec6 and Hecke/IdealBridge.lean 6c1fc33b7372108323bf3af3d4ab6f6e454ea9896b89e78964463a8e35df735f, each matched against its pinned Git blob.
- No QRH module may appear anywhere in the Challenge's imports (4386 modules in those runs).

Not the same thing: the tracker's signed-track canonical challenge (verifier/upstream-challenge.lean, sha256 8fb13ad9...; verifier/challenge.mjs rewrites it to namespace QRH, theorem canonical_bound, `import Mathlib`, all-Dirichlet only) belongs to the uncommissioned admission path.

## controls
Palomar's v.comparator_preflight runs before any candidate code, on the same runner, with the judge toolchain (lean, leanexport and lake comparator from 4.35.0-rc2).

Setup:
- Three modules compiled against Init only:
  - PalomarPreflightChallenge and PalomarPreflightSolution: `theorem palomar_preflight : True := trivial`.
  - PalomarPreflightWrong: `theorem palomar_preflight : True ∧ True := ⟨trivial, trivial⟩`.
- A fourth export, PalomarPreflightIllTyped, is made by v.ill_typed_export. It rewrites the NDJSON of the Solution export so the theorem's `value` is replaced by its `type`: well-formed, but ill-typed.
- All of them are judged with the same protected kernel config: permitted_axioms are the standard three, external kernels nanoda and con-ron.

Expected results:
1. Positive (PalomarPreflightSolution): exit 0, "Your solution is okay!", and all three "<kernel> kernel accepts the solution" lines. con-ron accepted 466 declarations.
2. Mismatch (PalomarPreflightWrong): exit 1 with "error: Challenge and solution theorem statement do not match: 'palomar_preflight'". It is rejected before any kernel runs, and must be classified as code comparator.rejected.
3. Ill-typed (PalomarPreflightIllTyped): exit 1. con-ron, nanoda and "Lean default" all print "kernel rejected the solution". The run fails unless Lean's own kernel is among the rejecters. nanoda panicked at src/tc.rs:955 ("assertion failed: self.def_eq(u, v)").

How the drivers recorded them: they temporarily wrap v.judge_exports (observe_preflight) to save each control's stdout+stderr and exit code without changing arguments or outcomes. The saved files are control-results.json ({"cases":[{name, exit_code, log}]} in the order PalomarPreflightSolution, PalomarPreflightWrong, PalomarPreflightIllTyped) and controls/matching.log, controls/mismatched.log, controls/ill-typed.log. result.json then carries "palomar_preflight": "passed" and preflight_cases.

The validators (core/native-kernel-evidence.mjs, PR #2 algebraic-kernel-evidence.mjs) re-check those exact names, exit codes and log substrings.

Other integrity checks in the drivers:
- Before compiling: every source and .olean hash in the frozen build is rechecked, the Challenge imports must contain no QRH modules, the pinned Palomar and lean4export Git HEADs must be clean, and every tool hash is checked.
- After the judge: v.verify_tool_snapshot runs and the source hashes are re-read.
- PR #2 additionally checks that the published source Git blobs equal the checked bytes, that LEAN_PATH resolves to the checked .olean files, and that the Solution wrapper contains no sorry and no added axiom.

## evidence_files_expected
["public/proofs/<id>-kernels/collection.json (schema_version 1, kind 'local-mechanical-kernel-evidence', source_commit, plus a files map of sha256 for every other file; the file inventory must match exactly)", "public/proofs/<id>-kernels/result.json (status PASS, comparator_exit_code 0, judge_exit_code, theta 'N/D', declarations, kernels ['Lean default','nanoda','con-ron'], allowed_axioms, source_commit, source_manifest_sha256, source_archive_sha256, palomar_commit d4e41c1..., exporter_commit, palomar_preflight 'passed', preflight_results, preflight_cases, statement_definitions_compared true, source_compiler '4.34.1', judge_toolchain '4.35.0-rc2', verified_at_utc, phases, resource_limits, native_input_validation, acceptance_markers, artifacts{path: sha256}, exports{challenge/solution sha256+bytes}, export_sizes, comparator_config, judge_log, challenge_source, solution_source, tool_pins, runner, host_policy_changed, global_userns_restriction_unchanged, registration_limits, publication_note)", "public/proofs/<id>-kernels/comparator.json (protected: PalomarCanonical<24hex>.Challenge, solution_module 'Solution', definition_names [], 3 axioms, external_kernels nanoda [nanoda_bin] and con-ron [con-ron, '--jobs=2'])", 'public/proofs/<id>-kernels/requested-comparator.json', "public/proofs/<id>-kernels/judge.log (all three 'kernel accepts the solution' lines plus 'Your solution is okay!')", 'public/proofs/<id>-kernels/control-results.json', 'public/proofs/<id>-kernels/controls/matching.log', 'public/proofs/<id>-kernels/controls/mismatched.log', 'public/proofs/<id>-kernels/controls/ill-typed.log', 'public/proofs/<id>-kernels/src/PalomarCanonical<24hex>/Challenge.lean', 'public/proofs/<id>-kernels/src/Solution.lean', 'public/proofs/<id>-kernels/export-targets.json', 'public/proofs/<id>-kernels/tool-pins.json (nested record of the run: palomar_commit, exporter_commit, source_compiler_version, judge_tools sha256s incl. bwrap, source_lean_sha256, source_exporter_sha256, driver_sha256, verifier_script_sha256, apparmor profile sha256 if used)', 'public/proofs/<id>-kernels/checker-binary-pins.json (flat input pins: {lake, lean, leanexport, leanchecker, nanoda_bin, con-ron: {sha256, matches_reference}}, exporter {sha256, source_commit})', 'public/proofs/<id>-kernels/checker-source-audit.json (verify_submission.py, supervise_cgroup.py, cgroup_delegate.py sha256 each matching its pinned Git blob; checker_code_modified false)', 'public/proofs/<id>-kernels/challenge-source-audit.json (imports of the Challenge, QRH_modules_in_challenge_closure [], source and .olean hashes, 4 protected definition files matching pinned Git blobs at mathlib d13f23b7 and openai/math fd4aeeb2)', 'public/proofs/<id>-kernels/build-inputs.json (native build rehash: modules_checked, qrh_modules_checked, challenge_imports_no_QRH true, sources_and_oleans_rehashed true)', 'public/proofs/<id>-kernels/source-inputs.json (manifest_sha256 and per-file sha256 of the checked source)', 'public/proofs/<id>-kernels/Challenge-compile.log, Challenge-export.log, Solution-compile.log, Solution-export.log', 'public/proofs/<id>-kernels/namespace-provisioning.json and namespace-cleanup.json (plus temporary-apparmor-profile.txt if an AppArmor profile was used; the earlier validators demand these)', 'public/proofs/<id>-kernels/original-report.json and publication-transformations.json (original vs published sha256 for each file whose workspace paths were neutralized)', 'public/proofs/<id>-kernels/PalomarSubmission-LICENSE.txt (MIT notice) and README.txt (method, pins, scope limits, reproduction command)', 'public/proofs/<id>-kernels/verify_<...>.py (the exact driver; its sha256 = tool-pins driver_sha256)', 'Proof exports themselves are omitted (too large); their sha256 and byte sizes are recorded in result.json', 'proofs/<id>/ (immutable source snapshot: LICENSE, NOTICE, README.md, PUBLICATION.md, formalization with lakefile, lake-manifest and lean-toolchain, audit/*.json build graph and records, scripts, inputs incl. manuscript, verification/src wrappers) and proofs/<id>.sha256.json', 'public/proofs/<id>/ (README.txt, result.json, source-public.tar.gz + SHA256SUMS.txt, collection.json, token-counts.json, etc.)', 'catalogue/results.json and public/catalogue.json record: status framework-verified, first_verified_at == result.verified_at_utc, source_commit, timeline_at, date_label, verification_note, attribution_note, references to the proofs/<id>-kernels/* files', 'core/<name>-kernel-evidence.mjs validator + scripts/check-<name>-kernels.mjs + tests/<name>-kernel-evidence.test.mjs + scripts/test-deployment.mjs entries (PR #2 precedent); must pass npm run check and npm run test:browser']

## feasible_on_macos
Not as a faithful reproduction of the precedent on native macOS. A partial version is possible natively; the faithful one needs a Linux VM or container.

Hard blockers on native macOS:
- Lean 4.35.0-rc2 `lake comparator` (Lake/CLI/Check.lean mkContext) returns "cannot run" on non-Linux. The only way around it is `--inadvisably-no-sandbox`, which prints "WARNING: Sandbox disabled, this run is not trustworthy." Lean's own comparator-from-export test is skipped on non-Linux.
- PalomarSubmission's verify_submission.py and both prior drivers assert Linux. They need bubblewrap 0.12.0 with user namespaces (nested bwrap in the judge phase), a seccomp filter fd, and cgroup-v2 delegation. Without systemd as PID 1, that delegation goes through passwordless `sudo -n python cgroup_delegate.py`.

Current machine state:
- `docker info` returns "500 Internal Server Error", so the daemon is not usable. Docker Desktop would have to be (re)started by the user.
- 26 GiB disk free (94% used).
- 36 GiB RAM, M3 Pro with 12 cores.
- cargo/rustc 1.88.0 (Homebrew).
- Only Lean v4.34.1 is installed. It has no leanexport, nanoda_bin or con-ron.

Option A, closest to precedent (needs downloads and the user's approval; I did none of this):
- Run in a privileged linux/arm64 container (ubuntu:24.04), or a Linux VM. Precedent: Palomar's own namespace profile sets container.privileged=true.
- Downloads:
  - lean-4.34.1-linux_aarch64.
  - lean-4.35.0-rc2-linux_aarch64.tar.zst (sha256 b702aca1e0bceb4d902ae1d269f514dc571195d9f6912590f7a81c983bf45324, 596 MB).
  - PalomarSubmission @ d4e41c1, lean4export @ 076e8e5, and the bubblewrap 0.12.0 tarball (sha256 9760d007...).
  - Packages: git, python3 ≥ 3.11 (for hashlib.file_digest), meson, ninja, gcc, sudo.
- Rebuild the Mathlib d13f23b7 + OAI + 228-module QRH closure inside Linux with Lean 4.34.1.
- Write a driver modeled on verify_candidate.py. It may import v unmodified and only override limits; preflight and kernels must not be disabled.
- Give the VM roughly 28–30 GiB RAM. The earlier run's limits had minimum_host_memory_bytes 30064771072, and con-ron --jobs=2 needs about 10–18 GiB on large exports. Disk space must be freed first.
- The binary hashes will be aarch64 ones, not the x86_64 pins (8241c5e6..., 4e5616d9..., bb807d18..., 8c5d64ba..., e8baaa71...). The existing validators hard-code the x86_64 hashes and '/tools/lean-4.35.0-rc2-linux/bin/...' paths, so the new dossier's validator must pin the arm64 values. The difference must be stated openly; it must not be relabeled.
- x86_64 emulation through Rosetta or QEMU to match the old hashes is impractical for a Mathlib-sized build and gigabyte-scale kernel replay.

Option B, macOS-native but weaker (no sandbox, not Palomar's judge):
- Use the already-built comparator d03acab and lean4export 076e8e5 with the landrun shim. Add "external_kernels": {"nanoda": ["<abs>/nanoda_bin"], "con-ron": ["<abs>/con-ron", "--jobs=2"]} to the config (d03acab supports external_kernels; enable_nanoda cannot be combined with it).
- Get nanoda_bin and con-ron one of two ways:
  - Extract bin/nanoda_bin and bin/con-ron from lean-4.35.0-rc2-darwin_aarch64.tar.zst (sha256 e721cf30671d3830fa6120c7b1a4e6b22624650ce1d618c0ca28a4889f5ef324, 572 MB).
  - Build `cargo build --release --locked` from ammkrn/nanoda_lib @ 4c544ed4099c8227f07d5de77ad1e69fb0740a27 and leanprover/con-ron @ cb54c9a4e5f604a678b0230f96b8358b9fd73c82 without submodules. Compatibility with rustc 1.88 is unverified.
- Prior runs showed that 4.35.0-rc2's nanoda/con-ron accept exports from lean4export 076e8e5 of 4.34.1-built files.
- The three controls would have to be recreated by hand to mirror Palomar's preflight.
- This meets the letter of CONTRIBUTING (Comparator + Lean + NanoDa + con-ron) but not the precedent's sandboxed Palomar judge. TRUST.md also says the macOS attempt "only compiled Comparator/exporter", so maintainers may treat this as weaker or "verification-pending".

Under either option, compile the Challenge against pristine openai/math fd4aeeb2 + Mathlib d13f23b7 rather than the pruned fork c388e88e6. In the fork, 300 of the 306 OAI modules imported (directly or indirectly) by IdealBridge differ from upstream, even though Family.lean and IdealBridge.lean themselves match. The source audit can then truthfully say the statement definitions match the pinned upstream Git blobs.

# REQUIREMENTS
 I read the tracker clone (main @ 40844a7, unchanged) and merged PR #1, and ran the checks in scratch copies. Main is still 40844a7. PR #2 from nielstron is open and touches the same files.

**Baseline run** (copy at /private/tmp/claude-501/-Users-akashlevy-Documents-riemann/5258226a-75a9-431d-830a-875c5394bd2d/scratchpad/tracker/check-copy, Node v24.17, Python 3.13):
- `npm ci --ignore-scripts` passed in 1.6 s. It reports 2 audit warnings; leave them alone.
- `npm run check` passed with exit 0 in 10.2 s:
  - the publication policy authenticated 8740 files and the archive contents;
  - 25 of 25 unit tests passed;
  - the build authenticated 8332 snapshot files and an empty signed registry, and tsc and vite succeeded.
- The tracked files `public/catalogue.json` and `public/site.json` were rebuilt byte-identical.
- `npm run test:browser` failed all 14 tests at browser launch. Playwright 1.56.1 needs `chromium_headless_shell-1194`, and only build 1200 is installed. Fixing this needs `npx playwright install chromium`, which is a download I did not run. CI runs this suite anyway.

**Draft record test** (copy check-copy2): I appended the template record below and added placeholder public files.
- Only `tests/catalogue.test.mjs` failed, because it hard-codes 5 records.
- After changing 5 to 6, `npm run check` passed again.
- Validator probes showed:
  - a framework-verified record with no `first_verified_at`, a pending record with a `first_verified_at`, the status `verified`, the rational 20998/24000, an `http://` link, a link containing a space, a predecessor dated later, an uppercase id, an empty `date_label` and any `exact_bound` are all rejected;
  - a link to a missing `public/` file fails the build.

**Gaps to cover in our PR:**
- No existing check scans a new `proofs/<id>/` folder. I put a `/Users/akashlevy/` path in a new `proofs/<id>/leak.txt` and the publication check still passed. The validator only covers the two existing packages and `public/proofs`.
- In `public/proofs`, a personal path is rejected. Binary files there are skipped only if they end in `.gz`, `.png`, `.jpg` or `.pdf`; a `.olean` there crashes the check with a UTF-8 decode error.
- So PR #1 and PR #2 each added their own validator for their package, wired into `check:publication` and the build. We need to do the same.

**What the normal result path requires** (CONTRIBUTING.md "Add a result or improve the website"):
- Fork the repo and branch. No `akashlevy/QuasiRiemannTracker` fork exists yet; `gh` is logged in as akashlevy.
- Append one record, last in the file, to `catalogue/results.json`. Status must be `framework-verified` or `verification-pending`.
- Add an immutable snapshot `proofs/<id>/`, its hash manifest `proofs/<id>.sha256.json`, and public downloads in `public/proofs/<id>/` (archive, checksum file, a hash list of the downloads, README/NOTICE copies, the manuscript PDF and TeX, and the Comparator evidence).
- Add `-text` lines to `.gitattributes` for the new folders.
- Add the new package validator plus tests, and update the hard-coded result counts in the unit, browser and deployment tests.
- Commit the rebuilt `public/catalogue.json`.
- Open a PR using the template. It needs GTimG's review, and a maintainer must approve CI for a first-time contributor. The site redeploys automatically after merge.

**Verification label:**
- `framework-verified` needs Comparator plus acceptance by Lean, NanoDa and con-ron (CONTRIBUTING.md:6).
- Our evidence so far is Comparator with the Lean kernel only. As it stands the record must be `verification-pending` with an empty `first_verified_at`.
- If we later run NanoDa and con-ron for real, we can propose `framework-verified`, using the actual completion timestamp from that run, with a separate immutable evidence folder as in PR #1.

**Signed admission path (`submissions/`): not usable now.**
- CONTRIBUTING.md:13-15 says it is not running and not commissioned.
- `verifier/trust-policy.json` has a null repository and worker image, no keys and no acceptance record, so the CI check and the signing step both fail closed.
- `verifier/acceptance.json` does not exist, and local checking needs Linux (exit code 78 on this Mac).
- It accepts only the all-Dirichlet statement, from a self-contained `src/Candidate` tree under strict limits (at most 256 files, 2 MiB each, 10 MiB total, no lakefiles). Dependencies must come from the pinned openai/math fd4aeeb image, which rules out a proof built on the modified akashlevy/math fork.
- The registry is empty, so `builds_on` would have to be `[]`.

**The record template below passes the validator, but some values are placeholders:**
- `timeline_at` 2026-10-09T23:00:00Z is not a real event time. Replace it with the actual Comparator completion time in UTC from the run log, or with the PR's `created_at` and the label "Public submission". It must be later than PR #2's 2026-10-09T19:51:12.662856+00:00.
- `source_commit`, `pr`, `submitted_at`, both `entrypoint` fields and the bracketed part of `method` must be filled with real values.

- REQ: Use the normal result path: fork GTimG/QuasiRiemannTracker, create a branch, add the result to catalogue/results.json keeping exact constants, authorship, publication dates and evidence links, then open a PR. Site checks plus maintainer review and merge trigger automatic deployment to GitHub Pages.
  SRC: CONTRIBUTING.md:3-9; docs/GITHUB.md:3; docs/DEPLOYMENT.md:15-17; .github/workflows/ci.yml:56-72
  HOW: Fork to akashlevy/QuasiRiemannTracker (it does not exist yet; gh is logged in as akashlevy). Branch from main@40844a7614b67840801e65b836b3da2e1f7ff029, e.g. akashlevy-weighted-numerator-20261009, and open the PR against GTimG:main. Get the user's go-ahead before pushing or opening the PR, since both are public.

- REQ: Append exactly one new record, LAST in records[]. The unit test requires the array to be in timeline order with every record strictly lowering the bound.
  SRC: tests/catalogue.test.mjs:11-24; core/rational.mjs:126-149; core/catalogue.mjs:4,14-16
  HOW: Add the record after the nielstron-20261009-tightening record, or after nielstron-algebraic-20261009 if PR #2 merges first. 10499/12000 is about 0.8749167, below every listed bound, so it advances the frontier. Use a UTC timeline_at (Z or +00:00) later than every earlier record, because the code compares these timestamps as strings.

- REQ: Record id: lowercase /^[a-z0-9][a-z0-9-]+$/, unique, equal to the proof package id.
  SRC: core/catalogue.mjs:23-25; precedent catalogue/results.json:242 (id == proofs/nielstron-20261009-tightening)
  HOW: id: "akashlevy-20261009-weighted-numerator". Use the same string for proofs/<id>/, proofs/<id>.sha256.json and public/proofs/<id>/.

- REQ: theta: numerator and denominator as canonical decimal strings, positive denominator, at most 200 digits, gcd 1.
  SRC: core/rational.mjs:2-25; core/catalogue.mjs:26
  HOW: {"numerator":"10499","denominator":"12000"} (gcd = 1; I checked it). Do not set exact_bound or bound_interval: the validator rejects any exact_bound other than Liu's expression.

- REQ: status is exactly 'framework-verified' or 'verification-pending'. 'verified' is reserved for signed registry receipts and is rejected. framework-verified means Comparator checked the exact statements and their definitions and Lean, NanoDa AND con-ron accepted the exports, with only propext, Classical.choice and Quot.sound allowed. pending covers results not reproduced in that full framework and claims no defect.
  SRC: CONTRIBUTING.md:6,13; core/catalogue.mjs:5-13,27-28; README.md:7; docs/TRUST.md:3
  HOW: Comparator with the Lean kernel alone gives status 'verification-pending'. To propose 'framework-verified', run the full Palomar-style judge (Comparator plus Lean, NanoDa and con-ron, with positive, mismatched-statement and ill-typed-proof controls) on Linux and retain the source, checker pins, exact statements and acceptance reports. Maintainers decide either way.

- REQ: Date fields: timeline_at must parse, with a non-empty date_label explaining what the date is. verification_note is required. Pending records must have an empty first_verified_at; framework-verified records need a parseable genuine check timestamp. Never invent a local verification date.
  SRC: CONTRIBUTING.md:6; core/catalogue.mjs:29-41; docs/TRUST.md:5; .github/PULL_REQUEST_TEMPLATE.md:19-20; core/native-kernel-evidence.mjs:393-398 (precedent: first_verified_at === report.verified_at_utc)
  HOW: Pending: first_verified_at "". timeline_at = the real Comparator completion UTC from the run log with date_label "Author-side Comparator check", or the real PR created_at with date_label "Public submission". Framework-verified: first_verified_at = timeline_at = the dossier's verified_at_utc, date_label "Independent kernel verification". Leave submitted_at, merged_at and published_at as "" and pr as null, or fill submitted_at and pr from the real PR after opening it.

- REQ: builds_on lists only actual mathematical reuse. Each id must exist in the catalogue with a strictly earlier timeline_at.
  SRC: core/catalogue.mjs:79-84; CONTRIBUTING.md:22; .github/PULL_REQUEST_TEMPLATE.md:13
  HOW: ["openai-baseline"] (the proof sits on the akashlevy/math fork of openai/math). Do not list proofcouncil, nielstron or liu unless their mathematics or code is actually reused; cite them in references or the method text instead.

- REQ: references: a non-empty array of {label, url}. url is either https:// or matches ^proofs/[a-zA-Z0-9_./-]+$ with no '..', and every non-https url must exist at public/<url>. URLs must be unique within the record because the page uses them as React keys.
  SRC: core/catalogue.mjs:42-51; scripts/catalogue.mjs:42-45; src/CatalogueDetail.tsx:81-92; src/urls.ts:1-4
  HOW: Point local links at public/proofs/<id>/ files (README.txt, manuscript.pdf, manuscript.tex, Challenge.lean, Solution.lean, comparator.json, judge.log, result.json, source-public.tar.gz) and add an https link to akashlevy/math@c388e88e6519d139ad2ca5571ef0c7fbf0e245b5. A url into the proofs/<id>/ snapshot only resolves if the file is also copied into public/proofs/<id>/. No spaces, unicode or '+' in file names.

- REQ: Attribution: authors (a non-empty list of {name, url?}) and title are required. The optional attribution_note is shown as 'Credits'. Precedent lists only humans or organisations as authors and describes AI assistance in attribution_note.
  SRC: core/catalogue.mjs:42-43; src/CatalogueDetail.tsx:29-54; catalogue/results.json:248-253,313; .github/PULL_REQUEST_TEMPLATE.md:11
  HOW: authors: [{"name":"Akash Levy","url":"https://github.com/akashlevy"}]. title: "Akash Levy · weighted numerator". attribution_note: state that GPT-6 Astra (Ultra) was used to discover the proof and perform partial formalization, and that Opus 5.5 Ultracode completed the formalization and submitted the result. Put the same text as a footnote (\thanks or \footnote) in the manuscript and repeat it in NOTICE, README and the PR body.

- REQ: The other fields used by existing records must be present: submitted_at, published_at, merged_at, pr, repository, source_commit, first_verified_at, entrypoint{module, declaration}, license, method, verification_note, date_label.
  SRC: src/types.ts:2-29; catalogue/results.json:225-314; src/CatalogueDetail.tsx:55-79
  HOW: Copy the template's field set. source_commit is the real 40-character commit holding the exact checked sources (tracker-fork snapshot commit or akashlevy/math commit), with repository to match. entrypoint holds the all-Dirichlet theorem's module and full declaration name from the Solution wrapper, or "" if not yet fixed. license: "Apache-2.0 (retained upstream notices)".

- REQ: Proof sources go under proofs/<stable-id>/ with license, attribution and reproduction instructions. Existing snapshots are immutable; revisions go in a new directory. Heavy Lean checks stay out of website CI.
  SRC: CONTRIBUTING.md:7; proofs/README.txt:14-19; docs/DEPLOYMENT.md:19
  HOW: Create proofs/akashlevy-20261009-weighted-numerator/ holding the 228-module Lean development (lean-toolchain 4.34.1, lake-manifest pinning openai/math or akashlevy/math@c388e88e6 and Mathlib d13f23b7), the manuscript TeX source, audit/target files and reproduction scripts. Do not add a Lean build to CI. Never edit proofs/qrh-20261009 or proofs/nielstron-*.

- REQ: License files: LICENSE (Apache-2.0 text), NOTICE (new attribution and derivative statement), UPSTREAM-NOTICE (precedent: a verbatim copy of the tracker's root NOTICE), and third-party notices that do not collide by case.
  SRC: proofs/nielstron-20261009-tightening/{LICENSE,NOTICE,UPSTREAM-NOTICE}; scripts/package-native-tightening.py:171-172,194-204; NOTICE:1-6; proofs/qrh-20261009/audit/LICENSES.md; maintainer fix 19b061b (license/ -> third-party-notices/)
  HOW: Copy the root LICENSE (identical to the existing package LICENSE files) into proofs/<id>/LICENSE. NOTICE should name Akash Levy, carry the AI footnote text, state the derivation from OpenAI/math fd4aeeb via akashlevy/math c388e88e6 and Mathlib d13f23b7, summarise modifications and say no endorsement is implied. Put the dependency LICENSE copies (openai-math, mathlib, PrimeNumberTheoremAnd, rellich-kondrachov if imported) in third-party-notices/, never in a 'license/' directory next to 'LICENSE'.

- REQ: README.md and PUBLICATION.md in the snapshot. README covers scope, exact bound, exact theorem names, pins, status, attribution and reproduction commands. PUBLICATION.md covers what is included or excluded, sanitisation, the allowlist and manifest, and how to regenerate.
  SRC: proofs/nielstron-20261009-tightening/README.md:1-99; proofs/nielstron-20261009-tightening/PUBLICATION.md:1-38; proofs/qrh-20261009/PUBLICATION.md:1-11
  HOW: README states 'Catalogue status: verification pending', what Comparator checked and that NanoDa and con-ron were not run, plus exact pins and runnable reproduction steps. PUBLICATION.md lists the exclusions (oleans, binaries, toolchains, caches, upstream checkouts, credentials, operational records) and any path neutralisation.

- REQ: publication-files.json: a non-empty, sorted, unique JSON array of safe relative POSIX paths that lists itself and matches the files on disk exactly. No symlinks, no case or NFC collisions, no path that is both a file and a directory.
  SRC: scripts/publication_policy.py:31-68
  HOW: Generate it with sorted() over the reviewed file list and include 'publication-files.json'. Remove .DS_Store and __pycache__ before checking (set PYTHONDONTWRITEBYTECODE=1). Check it with publication_policy.approved_files(Path('proofs/<id>')).

- REQ: proofs/<id>.sha256.json maps every allowlisted path to its sha256, with the same key set as the allowlist (including publication-files.json), written as json.dumps(indent=2) plus a trailing newline.
  SRC: scripts/package-proof.py:17,30; scripts/publication_policy.py:162-166; scripts/catalogue.mjs:10-22
  HOW: Generate it with a packager modelled on scripts/package-proof.py. Do not hand-edit it after packaging, and keep public/proofs/<id>/source-sha256.json identical to it (precedent).

- REQ: If any bytes are changed for publication (for example to neutralise a workspace path), record each change in publication-transformations.json. Lean sources, manuscripts and license notices keep their original bytes.
  SRC: proofs/nielstron-20261009-tightening/PUBLICATION.md:10-17; proofs/nielstron-20261009-tightening/publication-transformations.json; docs/DEPLOYMENT.md:25
  HOW: Replace /Users/akashlevy/..., /private/tmp/claude-501/... and the scratchpad paths in logs and JSON with /work/qrh-proof. Record {path, original_sha256, published_sha256, change} for each, plus source_policy and omitted notes.

- REQ: public/proofs/<id>/source-public.tar.gz: gzip; every member a regular file named <id>/<allowlisted path> with a sha256 equal to the manifest; no directory, symlink or duplicate entries; complete; uid/gid 0 and empty uname/gname; each member at most 10 MiB and the total at most 100 MiB uncompressed; every member valid UTF-8 that passes the privacy scan.
  SRC: scripts/publication_policy.py:71-91; scripts/package-proof.py:19-29; scripts/package-native-tightening.py:80-88
  HOW: Build it with Python tarfile adding files only (mode 0644, mtime 0, gzip mtime=0 and filename=''), not 'tar czf', which adds directory entries and owner names. Check it with publication_policy.check_archive(archive, proof, manifest).

- REQ: public/proofs/<id>/SHA256SUMS.txt is a single line '<sha256>  source-public.tar.gz' (two spaces) plus a newline.
  SRC: scripts/publication_policy.py:122-124,172-174; scripts/package-proof.py:31-32
  HOW: Write it from the final archive bytes. scripts/test-deployment.mjs (PR #2 pattern) re-downloads the archive and compares its hash.

- REQ: Public downloads: a collection.json mapping every other file in public/proofs/<id>/ to its sha256. Precedent also includes status and snapshot_exports, which map public copies back to snapshot paths with identical bytes.
  SRC: scripts/publication_policy.py:107-119; scripts/package-native-tightening.py:36-49,89-94; public/proofs/nielstron-20261009-tightening/collection.json; public/proofs/nielstron-20261009-kernels/collection.json
  HOW: Publish README.txt, NOTICE.txt, PUBLICATION.txt, manuscript.pdf, manuscript.tex, Challenge.lean, Solution.lean, comparator.json, judge.log, result.json, source-sha256.json, SHA256SUMS.txt, source-public.tar.gz and collection.json with status 'verification-pending'. If NanoDa and con-ron are run later, add a separate immutable public/proofs/<id>-kernels/ evidence folder.

- REQ: .gitattributes gets '-text' entries for the new snapshot and public folders so the hashed bytes are not changed by end-of-line conversion.
  SRC: .gitattributes:1-4
  HOW: Append 'proofs/akashlevy-20261009-weighted-numerator/** -text' and 'public/proofs/akashlevy-20261009-weighted-numerator/** -text' (plus a -kernels line if added). Confirm with: git check-attr text -- <file> (expect 'unset').

- REQ: Forbidden content: compiled oleans/ileans, compiler binaries, toolchains, .lake, build/upstream/downloads dirs, caches, credentials, raw agent session or event streams, operational files. Every file in the snapshot must be UTF-8 text.
  SRC: proofs/README.txt:9-12; docs/DEPLOYMENT.md:23; .gitignore:7-15; proofs/nielstron-20261009-tightening/PUBLICATION.md:10-11; scripts/publication_policy.py:18,21-28
  HOW: Ship sources, logs and JSON only; reproduction rebuilds compiled objects from pins. Avoid file names matching PROGRESS.txt, FORMALIZATION_PROMPT.txt, auth.json, .env*, *-events.jsonl, 8-digit+T timestamp names, preflight-*, launcher.log, continuation.log, final-resources.json or space-cleanup-*. Rename Palomar preflight logs, e.g. to controls/matching.log.

- REQ: Privacy scan rejects: /Users/<x>/ or /scratch/userdata/<x>/ paths, SSH key paths and private keys, CODEX_HOME assignments, GitHub or OpenAI/Anthropic tokens, and the words session/thread/resume followed within 160 characters by a UUID.
  SRC: scripts/publication_policy.py:11-17,176-181
  HOW: Before packaging, grep for /Users/, /private/tmp, claude-501, scratchpad, gh[pousr]_, github_pat_, sk-ant-, 'BEGIN .*PRIVATE KEY' and session UUIDs; neutralise and record each change. PDFs are not scanned, so also check manuscript.pdf with strings for /Users/ or /private/ paths.

- REQ: No symlinks anywhere in proofs/<id>, public/proofs or evidence.
  SRC: scripts/publication_policy.py:54-62,177-179; core/native-kernel-evidence.mjs:59-72
  HOW: Find symlinks with 'find proofs/<id> public/proofs/<id> -type l' and replace each with a real file copy or remove it.

- REQ: Manuscript placement: TeX source belongs in the snapshot and tarball; a PDF may only go in public/proofs, where .pdf/.png/.jpg/.gz skip the UTF-8 scan. Any other binary there crashes the check (I verified this with a .olean).
  SRC: scripts/publication_policy.py:180-181; catalogue/results.json:214-217; public/proofs/qrh-20261009/accepted-answer.tex (identical to proofs/qrh-20261009/inputs/.../accepted-answer.tex)
  HOW: Put the clean 10499/12000 paper as proofs/<id>/manuscript/<name>.tex (+ .bib, text figures), copy it to public/proofs/<id>/manuscript.tex, and place the compiled PDF only at public/proofs/<id>/manuscript.pdf. Include the GPT-6 Astra (Ultra) / Opus 5.5 Ultracode footnote and Akash Levy as author.

- REQ: Add an integrity and privacy validator for the new package, wired into the existing checks with tests. No current check covers a new proofs/<id>/: I placed a /Users/ path there and the check still passed. Precedent: PR #1 and PR #2 each added one; README says the build checks the complete snapshot hashes.
  SRC: scripts/publication_policy.py:158-182 (only qrh + nielstron); scripts/catalogue.mjs:6-22; README.md:22; PR #1 files (core/native-kernel-evidence.mjs, scripts/publication_policy.py:94-155,175); PR #2 files (scripts/check-algebraic-publication.py, scripts/publication_policy.py +2 lines)
  HOW: Add e.g. scripts/check-akashlevy-publication.py, which runs approved_files, the manifest hash check, collection.json, check_archive and SHA256SUMS. It should bind result.json's theta to the catalogue record and require PASS, the three declarations and allowed axioms only, and require that pending status keeps first_verified_at empty. Call it from check_publication via subprocess (PR #2 pattern), optionally also from scripts/catalogue.mjs, and add a node:test test. Do not modify the existing package validators.

- REQ: Comparator evidence must use the fixed three-target template (allDirichlet, zeta, allHecke) with only the rational literal replaced, no definition holes, permitted axioms limited to the standard three, and Hecke definitions from unmodified pinned OpenAI/math.
  SRC: public/proofs/palomar-20261009/qrh/src/Challenge.lean:1-19; public/proofs/palomar-20261009/qrh/comparator.json; core/native-kernel-evidence.mjs:149-177,294-308,351-358
  HOW: Publish Challenge.lean as the ProofCouncil template with '874957019421 / 1000000000000' replaced by '10499 / 12000', plus Solution.lean, comparator.json (definition_names []) and the judge log. Show that the challenge's definition files (e.g. OAI/NumberTheory/DirichletL/Hecke/IdealBridge closure) match openai/math fd4aeeb git blobs, since the fork pruned DirichletL in place. Record the exact Comparator, Lean and exporter pins.

- REQ: Update the tests that hard-code result counts and the newest record.
  SRC: tests/catalogue.test.mjs:11-20; tests/browser/dashboard.spec.ts:40-48,122-126; tests/browser/chart-labels.spec.ts:49-58,86; scripts/test-deployment.mjs:69-87,127
  HOW: Counts go 5 to 6 (7 if PR #2 merges first). Verified-id list unchanged if pending. Teal tags 3, amber 3 (pending). The bound-sorted first row becomes 'Akash Levy'. Verified filter 3 and pending filter 3. Chart dots 6 and verified dots 3. Make the chart-labels 'newest' check use the new id. Add the new public evidence paths and archive hash check to test-deployment.mjs.

- REQ: Commit the regenerated public/catalogue.json (the build writes it from catalogue/results.json).
  SRC: scripts/catalogue.mjs:105-109; PR #1 file list (public/catalogue.json +90)
  HOW: Run npm run build or npm run check after editing the catalogue and commit public/catalogue.json. public/site.json should stay unchanged.

- REQ: Existing proof snapshots, historical evidence and the registry stay unchanged.
  SRC: CONTRIBUTING.md:7; .github/PULL_REQUEST_TEMPLATE.md:16; proofs/README.txt:18-19; verifier/audit-history.mjs:7,21; .github/workflows/ci.yml:28-36
  HOW: Confirm 'git diff origin/main --stat -- proofs/qrh-20261009 proofs/nielstron-20261009-tightening* public/proofs/qrh-20261009 public/proofs/nielstron-* public/proofs/palomar-20261009 registry evidence' is empty. Revert evidence/screenshots/*.png rewritten by browser tests (PR #1 did not commit them).

- REQ: Required local commands: npm ci --ignore-scripts, npm run check (publication policy, unit tests, build), npm run test:browser (install Chromium if needed). CI also runs the production-path deployment test.
  SRC: CONTRIBUTING.md:8; package.json:12,14,18,21; .github/workflows/ci.yml:23,37-45; README.md:15
  HOW: Baseline confirmed: npm ci passes in 1.6 s and npm run check passes in 10.2 s (25/25 tests). Browser tests need 'npx playwright install chromium' first (only build 1200 is installed; 1194 is needed). Then run QRH_BASE_PATH=/QuasiRiemannTracker/ npm run build && QRH_BASE_PATH=/QuasiRiemannTracker/ node scripts/test-deployment.mjs.

- REQ: The PR description follows the template: type, summary, stable id and exact bound, source/paper/proof links, verification status and evidence, authors and AI collaboration attribution, license and upstream notices, mathematical predecessors, and the three checkboxes.
  SRC: .github/PULL_REQUEST_TEMPLATE.md:1-20; CONTRIBUTING.md:9; PR #1 body
  HOW: Type 'result metadata / native Lean evidence'. Id akashlevy-20261009-weighted-numerator, θ = 10499/12000. Say honestly that Comparator with the Lean kernel passed and NanoDa/con-ron were not run, and that the label is verification pending. Include the AI footnote text and predecessor openai-baseline. Tick only what is true. End with the Claude Code attribution line; commits end with 'Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>'.

- REQ: Status and evidence changes need maintainer review; no author can grant their own verified badge. All paths are owned by @GTimG. A first-time contributor's CI run needs maintainer approval.
  SRC: CONTRIBUTING.md:13; .github/CODEOWNERS:2; docs/DEPLOYMENT.md:17; PR #1 body (CI approval note)
  HOW: Present framework-verified only as a proposal, and only with full three-kernel evidence. Expect GTimG to approve the CI run and review; the site deploys automatically from main after merge.

- REQ: Signed Comparator + NanoDa admission (submissions/) is not usable now.
  SRC: CONTRIBUTING.md:13,15-26; verifier/trust-policy.json:1-10; verifier/ci-check.mjs:11-16; verifier/sign-result.mjs:12-15; verifier/runner.mjs:49; verifier/manifest.mjs:64-84; registry/README.md; docs/TRUST.md:53; docs/REPRODUCTION.md:16,43
  HOW: Do not use it. The trust policy has a null repository and worker image, no keys and no acceptance record, verifier/acceptance.json is absent, and local runs exit 78 on non-Linux. It also admits only the all-Dirichlet statement from a self-contained src/Candidate tree (at most 256 files, 2 MiB each, 10 MiB total, no lakefiles, image-pinned openai/math fd4aeeb dependencies only), with builds_on drawn from an empty registry.

- REQ: Coordinate with open PR #2 (nielstron-algebraic-20261009, timeline_at 2026-10-09T19:51:12.662856+00:00), which edits the same catalogue, test, deployment-test and publication-policy files.
  SRC: PR #2 metadata and files (gh api repos/GTimG/QuasiRiemannTracker/pulls/2); tests/catalogue.test.mjs:19-24
  HOW: Keep our timeline_at later than 19:51:12Z UTC. Whichever PR merges second must rebase: put our record after theirs, bump counts to 7, add their id to any verified lists, and re-run npm run check.

# RECORD TEMPLATE
 {
  "submitted_at": "",
  "published_at": "",
  "merged_at": "",
  "pr": null,
  "repository": "akashlevy/QuasiRiemannTracker",
  "source_commit": "",
  "first_verified_at": "",
  "entrypoint": {
    "module": "",
    "declaration": ""
  },
  "license": "Apache-2.0 (retained upstream notices)",
  "builds_on": ["openai-baseline"],
  "id": "akashlevy-20261009-weighted-numerator",
  "title": "Akash Levy · weighted numerator",
  "theta": {
    "numerator": "10499",
    "denominator": "12000"
  },
  "authors": [
    {
      "name": "Akash Levy",
      "url": "https://github.com/akashlevy"
    }
  ],
  "method": "Unconditional nonvanishing of the Riemann zeta function, every Dirichlet L-function and the finite-order Hecke L-functions of Q(ζ₃) for Re s > 10499/12000. Lean 4.34.1 development of 228 modules on the akashlevy/math fork (c388e88e6) of OpenAI/math fd4aeeb, with Mathlib d13f23b7. [Replace with a one-paragraph description of the weighted-numerator argument.]",
  "references": [
    {
      "label": "Manuscript (PDF)",
      "url": "proofs/akashlevy-20261009-weighted-numerator/manuscript.pdf"
    },
    {
      "label": "Manuscript (TeX)",
      "url": "proofs/akashlevy-20261009-weighted-numerator/manuscript.tex"
    },
    {
      "label": "Exact checked statements",
      "url": "proofs/akashlevy-20261009-weighted-numerator/Challenge.lean"
    },
    {
      "label": "Solution wrapper",
      "url": "proofs/akashlevy-20261009-weighted-numerator/Solution.lean"
    },
    {
      "label": "Comparator configuration",
      "url": "proofs/akashlevy-20261009-weighted-numerator/comparator.json"
    },
    {
      "label": "Comparator acceptance log",
      "url": "proofs/akashlevy-20261009-weighted-numerator/judge.log"
    },
    {
      "label": "Check report and artifact hashes",
      "url": "proofs/akashlevy-20261009-weighted-numerator/result.json"
    },
    {
      "label": "Scope, attribution and reproduction",
      "url": "proofs/akashlevy-20261009-weighted-numerator/README.txt"
    },
    {
      "label": "Complete Lean source, manuscript source and evidence",
      "url": "proofs/akashlevy-20261009-weighted-numerator/source-public.tar.gz"
    },
    {
      "label": "Lean base: akashlevy/math fork of OpenAI/math",
      "url": "https://github.com/akashlevy/math/tree/c388e88e6519d139ad2ca5571ef0c7fbf0e245b5"
    }
  ],
  "status": "verification-pending",
  "timeline_at": "2026-10-09T23:00:00Z",
  "date_label": "Author-side Comparator check",
  "verification_note": "Author-side check: Comparator compared the three exact all-Dirichlet, zeta and finite-order Eisenstein Hecke statements at 10499/12000 and their definitions, and Lean's kernel replayed the proof, with only propext, Classical.choice and Quot.sound permitted. NanoDa and con-ron have not yet checked the exported proof, so this result is listed as verification pending. This is not Palomar registration, editorial review or a signed tracker receipt.",
  "attribution_note": "Akash Levy (akashlevy) is the human contributor. GPT-6 Astra (Ultra) was used to discover the proof and perform partial formalization; Opus 5.5 Ultracode completed the formalization and submitted the result. Builds on the Apache-2.0 OpenAI/math and Mathlib developments; upstream notices are retained."
}

# LAYOUT
 Changes in the PR branch, relative to the tracker repo root. Nothing else changes; existing snapshots and registry/ stay untouched.

.gitattributes                                   (+2 or 3 lines)
    proofs/akashlevy-20261009-weighted-numerator/** -text
    public/proofs/akashlevy-20261009-weighted-numerator/** -text
    [public/proofs/akashlevy-20261009-weighted-numerator-kernels/** -text   only if a 3-kernel dossier is added]
catalogue/results.json                           (append the record LAST)
public/catalogue.json                            (regenerated by npm run build; commit it)
proofs/akashlevy-20261009-weighted-numerator.sha256.json   ({path: sha256} for every allowlisted file)
proofs/akashlevy-20261009-weighted-numerator/    (immutable snapshot; UTF-8 text only; no symlinks/oleans/binaries/.DS_Store/__pycache__)
    LICENSE                    Apache-2.0 (byte-identical to repo LICENSE)
    NOTICE                     Akash Levy attribution + GPT-6 Astra (Ultra) / Opus 5.5 Ultracode statement, derivation from openai/math fd4aeeb via akashlevy/math c388e88e6, Mathlib d13f23b7, modifications
    UPSTREAM-NOTICE            (precedent) verbatim upstream notice
    README.md                  scope, θ = 10499/12000, exact theorem names, pins, status (verification pending), what Comparator checked, reproduction commands
    PUBLICATION.md             inclusions/exclusions, sanitisation, allowlist/manifest, regeneration command
    publication-files.json     sorted unique allowlist (includes itself)
    publication-transformations.json   (only if any bytes were neutralised)
    formalization/             lean-toolchain, lakefile, lake-manifest.json, QRH.lean (or equivalent root) + the 228 modules
    manuscript/<paper>.tex     clean 10499/12000 manuscript source with the AI footnote (+ .bib, text figures)
    audit/                     exact Challenge/target statements, axiom report, build status, dependency-licenses.json, patch of fork changes vs openai/math fd4aeeb
    third-party-notices/       openai-math/mathlib/PNT/rellich LICENSE copies (never a 'license/' dir beside LICENSE)
    scripts/                   reproduction scripts + packager (no bytecode)
public/proofs/akashlevy-20261009-weighted-numerator/
    source-public.tar.gz       deterministic gzip tar of exactly the allowlist under akashlevy-20261009-weighted-numerator/
    SHA256SUMS.txt             '<sha256>  source-public.tar.gz'
    source-sha256.json         == proofs/<id>.sha256.json
    collection.json            {status: "verification-pending", files: {name: sha256}, snapshot_exports: {...}}
    README.txt, NOTICE.txt, PUBLICATION.txt   (byte copies of the snapshot files)
    manuscript.pdf             (PDF allowed only here; scrub absolute paths from its metadata)
    manuscript.tex             (copy of the snapshot TeX)
    Challenge.lean, Solution.lean, comparator.json, judge.log, result.json   (Comparator evidence, neutralised paths)
[public/proofs/akashlevy-20261009-weighted-numerator-kernels/]   optional later immutable dossier with Comparator + Lean + NanoDa + con-ron results and controls (needed for framework-verified)
scripts/check-akashlevy-publication.py           new validator (allowlist, manifest, archive, SHA256SUMS, collection, theta/status binding, privacy)
scripts/publication_policy.py                    +1-2 lines calling the new validator (PR #2 pattern)
[scripts/catalogue.mjs]                          optional +1-2 lines to validate the new evidence during build
tests/catalogue.test.mjs                         counts 5 -> 6 (or 7 after PR #2)
tests/<new-validator>.test.mjs                   regression test for the new validator
tests/browser/dashboard.spec.ts                  row/dot/tag/filter counts; bound-sorted first row 'Akash Levy'
tests/browser/chart-labels.spec.ts               dot counts; newest-label id
scripts/test-deployment.mjs                      dot/row counts; new evidence paths + archive hash check

# COMMANDS
 gh repo fork GTimG/QuasiRiemannTracker --clone=false   # public action: only after the user confirms
git clone https://github.com/akashlevy/QuasiRiemannTracker && cd QuasiRiemannTracker && git checkout -b akashlevy-weighted-numerator-20261009 40844a7614b67840801e65b836b3da2e1f7ff029
npm ci --ignore-scripts
PYTHONDONTWRITEBYTECODE=1 python3 -c "import sys,json,hashlib; from pathlib import Path; sys.path.insert(0,'scripts'); from publication_policy import approved_files, check_archive; p=Path('proofs/akashlevy-20261009-weighted-numerator'); n=approved_files(p); m=json.load(open('proofs/akashlevy-20261009-weighted-numerator.sha256.json')); assert set(n)==set(m); assert all(hashlib.sha256((p/k).read_bytes()).hexdigest()==v for k,v in m.items()); check_archive(Path('public/proofs/akashlevy-20261009-weighted-numerator/source-public.tar.gz'), p, m); print('snapshot OK', len(n))"
python3 -c "import hashlib;d,n=open('public/proofs/akashlevy-20261009-weighted-numerator/SHA256SUMS.txt').read().strip().split('  ');assert n=='source-public.tar.gz' and hashlib.sha256(open('public/proofs/akashlevy-20261009-weighted-numerator/'+n,'rb').read()).hexdigest()==d;print('SHA256SUMS OK')"
grep -rInE '/Users/|/private/tmp|/scratch/userdata|claude-501|scratchpad|\.ssh/|BEGIN [A-Z ]*PRIVATE KEY|gh[pousr]_[A-Za-z0-9]{30,}|github_pat_|sk-(proj-|ant-)' proofs/akashlevy-20261009-weighted-numerator public/proofs/akashlevy-20261009-weighted-numerator || echo 'no private strings'
find proofs/akashlevy-20261009-weighted-numerator public/proofs/akashlevy-20261009-weighted-numerator \( -type l -o -name '*.olean' -o -name '*.ilean' -o -name '*.c' -o -name '.DS_Store' -o -name '__pycache__' -o -size +10M \) -print
strings public/proofs/akashlevy-20261009-weighted-numerator/manuscript.pdf | grep -E '/Users/|/private/' || echo 'pdf metadata clean'
git check-attr text -- proofs/akashlevy-20261009-weighted-numerator/README.md public/proofs/akashlevy-20261009-weighted-numerator/source-public.tar.gz   # expect 'unset'
python3 scripts/check-akashlevy-publication.py   # the new validator
npm run check:publication
npm test
npm run build
npm run check
npx playwright install chromium   # download; get user OK first (local chromium_headless_shell-1194 is missing)
npm run test:browser
QRH_BASE_PATH=/QuasiRiemannTracker/ npm run build && QRH_BASE_PATH=/QuasiRiemannTracker/ node scripts/test-deployment.mjs
git checkout -- evidence/   # drop screenshots rewritten by browser tests
git diff --stat origin/main -- proofs/qrh-20261009 proofs/nielstron-20261009-tightening proofs/nielstron-20261009-tightening.sha256.json proofs/qrh-20261009.sha256.json public/proofs/qrh-20261009 public/proofs/nielstron-20261009-tightening public/proofs/nielstron-20261009-kernels public/proofs/palomar-20261009 registry   # must be empty
git status --short && git diff --stat origin/main
gh pr create --repo GTimG/QuasiRiemannTracker --base main --head akashlevy:akashlevy-weighted-numerator-20261009 --title 'Add 10499/12000 weighted-numerator QRH result (Akash Levy)' --body-file pr-body.md   # public action: only after the user confirms; body follows .github/PULL_REQUEST_TEMPLATE.md and ends with the Claude Code attribution line

# PITFALLS
 - Label honesty: Comparator with the Lean kernel only does NOT meet framework-verified, which also needs NanoDa and con-ron (CONTRIBUTING.md:6). Submit as verification-pending with first_verified_at "" unless a real three-kernel run exists. validateCatalogue accepts any framework-verified record that has a date (I checked), so only maintainer review prevents a self-granted badge; do not rely on that.
- PLACEHOLDER: timeline_at 2026-10-09T23:00:00Z in the template is not a real event time. Replace it with the actual Comparator completion UTC from the run log (date_label 'Author-side Comparator check') or the real PR created_at (date_label 'Public submission'). Never invent a date. Use UTC only, because timeline_at values are compared as strings; it must be later than 2026-10-09T19:51:12.662856+00:00 (PR #2).
- Also fill these in the template: source_commit (real 40-character commit of the checked sources; set repository to akashlevy/math if that is where it lives), entrypoint module/declaration (from the Solution wrapper's all-Dirichlet theorem), the bracketed part of method, and pr/submitted_at once the PR exists.
- No existing check scans a new proofs/<id>/ folder: I put a /Users/akashlevy/ path there and npm run check:publication still passed. Write and wire a validator (as PR #1 and PR #2 did) and run the privacy grep yourself; otherwise path leaks and hash drift go unnoticed.
- Any non-UTF-8 file in proofs/<id>/ or inside the tarball makes check_content raise UnicodeDecodeError. Keep PDFs, images and oleans out of the snapshot. In public/proofs only .pdf/.png/.jpg/.gz are skipped; I confirmed a .olean there crashes the check. PDFs are not privacy-scanned, so scrub absolute paths from PDF metadata manually.
- Building the tarball with 'tar czf' fails check_archive: it adds directory entries and uid/uname. Use Python tarfile with only file members, mode 0644, mtime 0, uid/gid 0, empty uname/gname and gzip mtime=0. Keep each member at most 10 MiB and the total at most 100 MiB; never ship proof exports (the precedent export was 1.7 GB), record their size and sha256 instead.
- Path rules: no case-insensitive or NFC collisions. The maintainer had to rename compressed/license/ to third-party-notices/ in commit 19b061b. Also no path that is both a file and a directory, and no .DS_Store or __pycache__ on disk (approved_files compares the real files to the allowlist).
- Neutralise /Users/akashlevy/... and /private/tmp/claude-501/.../scratchpad paths, plus 'session|thread|resume' followed by a UUID, in logs and JSON, and record the original/published hashes in publication-transformations.json. Lean sources, the TeX manuscript and license texts must keep their exact bytes. Rename Palomar 'preflight-*' logs and avoid 'YYYYMMDDT...' file names, which the operational-file regex rejects.
- The Challenge must be the fixed three-target ProofCouncil template with only the literal changed to 10499 / 12000. The Hecke target imports OpenAI definitions, and the akashlevy/math fork pruned DirichletL in place. Show the challenge's definition files are byte-identical to openai/math fd4aeeb git blobs, or the Comparator result proves a different statement than the tracker's.
- Open PR #2 adds a sixth record and edits the same tests and files. Whichever PR merges second must rebase, keep chronological record order, bump counts to 7 and re-run checks. Expect conflicts in catalogue/results.json, public/catalogue.json, tests/catalogue.test.mjs, the browser specs, scripts/test-deployment.mjs and scripts/publication_policy.py.
- Unit tests and browser tests hard-code result counts. Only tests/catalogue.test.mjs fails under npm run check, but CI also runs Playwright and test-deployment, so update dashboard.spec.ts (counts and the bound-sorted first row becoming 'Akash Levy'), chart-labels.spec.ts (dot counts and newest-label id) and test-deployment.mjs (counts and new paths).
- I could not run the browser suite locally: all 14 tests failed because chromium_headless_shell-1194 is missing (build 1200 is installed). Run 'npx playwright install chromium' (a download, so get user approval) before claiming 'npm run test:browser passes' in the PR. Revert evidence/screenshots/*.png rewritten by those tests.
- builds_on is for actual mathematical reuse only (openai-baseline). Cite ProofCouncil, Nielstron or Liu as references, not predecessors, unless their mathematics or code is reused; method similarity is not dependence.
- Do not modify any existing snapshot or evidence folder (proofs/qrh-20261009, proofs/nielstron-*, public/proofs/{qrh-20261009,nielstron-*,palomar-20261009}), registry/, verifier/trust-policy.json or package-lock.json. Do not run 'npm audit fix' despite the 2 reported vulnerabilities.
- Do not try the submissions/ signed path. It is uncommissioned (null policy, no keys, no acceptance record, Linux only, exit 78 here). It accepts only the all-Dirichlet statement from a self-contained src/Candidate tree within strict size limits, so the fork-based 228-module proof does not fit, and its registry is empty.
- List only humans as authors (Akash Levy, https://github.com/akashlevy). Put the required AI statement (GPT-6 Astra (Ultra) discovered the proof and did partial formalization; Opus 5.5 Ultracode completed the formalization and submitted) as a footnote in the manuscript, and in attribution_note, NOTICE, README and the PR body.
- Forking, pushing and opening the PR are public actions. The user's request authorises submission once everything is confirmed, so confirm the final package with the user before running them. A maintainer must approve the first-time contributor's CI run, and GTimG must review; deployment then happens automatically on merge.
