#!/usr/bin/env bash
# M0: install elan (without touching shell profiles), the pinned Lean toolchain,
# the patched git dependencies, and Mathlib's prebuilt cache.
set -euo pipefail
task_root=$(cd "$(dirname "$0")/.." && pwd)
cd "$task_root"
export PATH="$HOME/.elan/bin:$PATH"

if ! command -v elan >/dev/null; then
  curl -sSfL https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -o /tmp/elan-init.$$.sh
  sh /tmp/elan-init.$$.sh -y --default-toolchain none --no-modify-path
  rm -f /tmp/elan-init.$$.sh
fi
elan toolchain install "$(cat lean-toolchain)"
lean --version
bash scripts/prepare-dependencies.sh
lake exe cache get
df -h "$task_root" | tail -1
