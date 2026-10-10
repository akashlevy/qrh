#!/usr/bin/env python3
"""Import-boundary audit for the 11/12 development.

Computes the transitive local import closure of the 11/12 roots (every module under
ElevenTwelfths/ plus the 11/12 Comparator challenge) and fails if the closure touches
any module denied by audit/policy.json.

Policy rules (see audit/policy.json, every rule carries a recorded reason):
  deny_modules        exact module names that may not be imported
  deny_prefixes       module-name prefixes that may not be imported
  deny_dependents_of  modules M such that any module whose closure contains M is denied
                      (this catches "anything whose proof depends on" an end theorem)
  allow_prefixes      documented allow decisions; purely informational except that the
                      report lists every used lean-base module not covered by any allow rule

Exit status: 0 = clean, 1 = boundary violation, 2 = usage/config error.
Usage: python3 scripts/audit_imports.py [--roots MOD ...] [--report audit/closure-report.md]
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importgraph import ROOT, closure, module_path  # noqa: E402

POLICY = ROOT / "audit" / "policy.json"


def default_roots() -> list[str]:
    roots = []
    for base in ("QRH",):
        for p in sorted((ROOT / base).rglob("*.lean")):
            roots.append(".".join(p.relative_to(ROOT).with_suffix("").parts))
    challenge = ROOT / "ComparatorChallenges" / "QRH.lean"
    if challenge.is_file():
        roots.append("ComparatorChallenges.QRH")
    top = ROOT / "QRH.lean"
    if top.is_file():
        roots.append("QRH")
    return roots


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--roots", nargs="*", help="root modules (default: all 11/12 modules)")
    ap.add_argument("--report", type=Path, help="write a markdown closure report here")
    args = ap.parse_args()

    policy = json.loads(POLICY.read_text())
    roots = args.roots or default_roots()
    if not roots:
        print("no root modules found", file=sys.stderr)
        return 2
    missing = [r for r in roots if module_path(r) is None]
    if missing:
        print(f"unknown root modules: {missing}", file=sys.stderr)
        return 2

    mods, external = closure(roots)

    violations: list[tuple[str, str]] = []
    deny_exact = {r["module"]: r["reason"] for r in policy.get("deny_modules", [])}
    for m in sorted(mods):
        if m in deny_exact:
            violations.append((m, deny_exact[m]))
        for r in policy.get("deny_prefixes", []):
            if m.startswith(r["prefix"]):
                violations.append((m, r["reason"]))
    for r in policy.get("deny_dependents_of", []):
        target = r["module"]
        if target in mods:
            violations.append((target, "in closure; " + r["reason"]))
    for prefix in policy.get("deny_external_prefixes", []):
        for e in sorted(external):
            if e.startswith(prefix["prefix"]):
                violations.append((e, prefix["reason"]))

    base = [m for m in mods if m.startswith("OAI.")]
    allow = policy.get("allow_prefixes", [])
    unclassified = sorted(m for m in base if not any(m.startswith(a["prefix"]) for a in allow))
    lines = sum(len(module_path(m).read_text().splitlines()) for m in base)

    print(f"roots: {len(roots)}; local closure: {len(mods)} modules "
          f"({len(base)} from lean-base, {lines} lines)")
    print("external roots:", dict(Counter(e.split(".")[0] for e in external)))
    if unclassified:
        print(f"note: {len(unclassified)} lean-base modules not covered by an allow rule")
    if args.report:
        by_dir = Counter(".".join(m.split(".")[3:5]) if m.count(".") > 4 else m for m in base)
        out = ["# Import-closure report", "", f"Roots: {', '.join(roots)}", "",
               f"Local closure: {len(mods)} modules; lean-base part: {len(base)} modules, {lines} lines.", "",
               "| lean-base directory | modules |", "|---|---|"]
        out += [f"| {d} | {n} |" for d, n in sorted(by_dir.items())]
        out += ["", "## Violations", ""] + ([f"- `{m}`: {why}" for m, why in violations] or ["None."])
        out += ["", "## lean-base modules without an allow rule", ""] + ([f"- `{m}`" for m in unclassified] or ["None."])
        args.report.write_text("\n".join(out) + "\n")
    if violations:
        print(f"BOUNDARY VIOLATION: {len(violations)} denied modules in closure", file=sys.stderr)
        for m, why in violations[:50]:
            print(f"  {m}: {why}", file=sys.stderr)
        return 1
    print("import boundary: OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
