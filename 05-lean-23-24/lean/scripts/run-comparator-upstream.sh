#!/usr/bin/env bash
# Comparator check of a challenge of the openai/math fork (math/lean), built in the upstream/
# workspace. Same macOS caveats as run-comparator-macos.sh (unsandboxed landrun shim).
# Usage: bash scripts/run-comparator-upstream.sh QuasiRiemannHypothesis
set -euo pipefail
task_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$task_root/upstream"
export PATH="$HOME/.elan/bin:$PATH"
export COMPARATOR_LANDRUN="$task_root/.tools/landrun-shim"
export COMPARATOR_LEAN4EXPORT="$task_root/.tools/comparator/.lake/packages/lean4export/.lake/build/bin/lean4export"
config="../math/lean/ComparatorChallenges/${1:?challenge name}.json"
printf 'Comparator starting: %s\n' "$config"
exec lake env "$task_root/.tools/comparator/.lake/build/bin/comparator" "$config"
