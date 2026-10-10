#!/usr/bin/env python3
"""Probe which theorems of a built module can be re-proved by a one-line tactic.

For every theorem defined in --module (with a source range of at least --min-lines lines,
optionally restricted to the reachable names in --reach), elaborates its statement as a
fresh goal and tries each tactic in --tactics (as `intros; TAC`) under a heartbeat budget.
Nothing is edited; output is a TSV: name, start, end, lines, winning tactic (or -), ms.

Usage: python3 scripts/golf_probe.py --module MOD --out FILE [--reach TSV] [--min-lines 3]
           [--heartbeats 200] [--tactics "simp" "grind" ...]
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_TACTICS = ["simp", "fun_prop", "positivity", "omega", "norm_num", "linarith", "grind"]

TEMPLATE = r'''import {mod}
open Lean Elab Command Meta

def tryTac (type : Expr) (tac : Syntax) (hb : Nat) : TermElabM Bool := do
  let s ← saveState
  -- heartbeat timeouts are runtime exceptions: plain `try/catch` does not see them
  let ok ← tryCatchRuntimeEx
    (withTheReader Core.Context (fun c => {{ c with maxHeartbeats := hb * 1000 }}) do
      Core.withCurrHeartbeats do
        let mvar ← mkFreshExprMVar type
        let goals ← Tactic.run mvar.mvarId! (Tactic.evalTactic tac)
        Term.synthesizeSyntheticMVarsNoPostponing
        let v ← instantiateMVars mvar
        pure (goals.isEmpty && !v.hasMVar && !v.hasSorry))
    (fun _ => pure false)
  s.restore
  return ok

elab "#golf_probe" : command => liftTermElabM do
  let env ← getEnv
  let target := `{mod}
  let allow : Option NameSet ← do
    if "{reach}" == "" then pure none else
      let mut set : NameSet := {{}}
      for r in (← IO.FS.lines "{reach}") do
        let cols := r.splitOn "\t"
        if cols.length ≥ 5 && cols[4]! == "1" then set := set.insert cols[1]!.toName
      pure (some set)
  let tacs : List (String × Syntax) := ({tactics} : List String).filterMap fun t =>
    match Parser.runParserCategory env `tactic s!"(intros; {{t}})" with
    | .ok stx => some (t, stx)
    | .error _ => none
  let mut rows : Array String := #[]
  for (n, ci) in env.constants.toList do
    let some idx := env.getModuleIdxFor? n | continue
    unless env.header.moduleNames[idx.toNat]! == target do continue
    let .thmInfo ti := ci | continue
    if let some a := allow then unless a.contains n do continue
    let some rg ← findDeclarationRanges? n | continue
    let lines := rg.range.endPos.line - rg.range.pos.line + 1
    if lines < {minlines} then continue
    let t0 ← IO.monoMsNow
    let mut win := "-"
    for (name, stx) in tacs do
      if ← tryTac ti.type stx {hb} then
        win := name
        break
    let t1 ← IO.monoMsNow
    rows := rows.push s!"{{n}}\t{{rg.range.pos.line}}\t{{rg.range.endPos.line}}\t{{lines}}\t{{win}}\t{{t1 - t0}}"
  IO.FS.writeFile "{out}" (String.intercalate "\n" rows.toList ++ "\n")
  logInfo m!"probed {{rows.size}} theorems with {{tacs.length}} tactics"

#golf_probe
'''


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reach", default="")
    ap.add_argument("--min-lines", type=int, default=3)
    ap.add_argument("--heartbeats", type=int, default=200)
    ap.add_argument("--tactics", nargs="*", default=DEFAULT_TACTICS)
    a = ap.parse_args()
    tactics = "[" + ", ".join('"' + t.replace('"', '\\"') + '"' for t in a.tactics) + "]"
    src = TEMPLATE.format(mod=a.module, reach=str(Path(a.reach).resolve()) if a.reach else "",
                          tactics=tactics, minlines=a.min_lines, hb=a.heartbeats,
                          out=str(Path(a.out).resolve()))
    with tempfile.NamedTemporaryFile("w", suffix=".lean", delete=False, dir=os.environ.get("TMPDIR")) as f:
        f.write(src)
        path = f.name
    env = dict(os.environ, PATH=f"{Path.home()}/.elan/bin:" + os.environ.get("PATH", ""))
    r = subprocess.run(["lake", "env", "lean", path], cwd=ROOT, env=env, capture_output=True, text=True)
    os.unlink(path)
    print("\n".join(l for l in (r.stdout + r.stderr).splitlines() if not l.startswith("warning:")))
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
