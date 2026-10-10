#!/usr/bin/env bash
# Comparator check of the QRH challenge on macOS.
# Differences from lean-base/scripts/run-comparator-service.sh (Linux):
#   * Landrun (Linux Landlock) is replaced by .tools/landrun-shim, which runs commands
#     unsandboxed. This is fine for a solution we wrote; for an adversarial-grade check run
#     the sandboxed Linux workflow.
#   * Comparator/lean4export are built from the same pinned commits by hand
#     (comparator d03acab1, lean4export 076e8e57) with the project's Lean toolchain.
set -euo pipefail
task_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$task_root"
export PATH="$HOME/.elan/bin:$PATH"
export COMPARATOR_LANDRUN="$task_root/.tools/landrun-shim"
export COMPARATOR_LEAN4EXPORT="$task_root/.tools/comparator/.lake/packages/lean4export/.lake/build/bin/lean4export"
config=${1:-ComparatorChallenges/QRH.json}
mkdir -p verification
date -u +'%Y-%m-%dT%H:%M:%SZ' > verification/comparator-started.txt
printf 'Comparator starting: %s\n' "$config"
exec lake env "$task_root/.tools/comparator/.lake/build/bin/comparator" "$config"
