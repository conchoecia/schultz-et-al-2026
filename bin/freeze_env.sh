#!/bin/bash
# Capture the exact resolved dep versions in the currently-active
# conda env. Writes `environment.lock.yml` (full conda list) and
# `requirements.lock.txt` (pip freeze) to the repo root. Commit these
# alongside a paper submission to ensure bit-exact reproducibility.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: no conda env activated (\$CONDA_PREFIX is empty)." >&2
    echo "       Activate egt-repro first, then re-run." >&2
    exit 1
fi

if command -v micromamba >/dev/null 2>&1; then
    MAMBA=micromamba
elif command -v mamba >/dev/null 2>&1; then
    MAMBA=mamba
else
    MAMBA=conda
fi

echo "==> Freezing env at $CONDA_PREFIX with $MAMBA"

$MAMBA env export -p "$CONDA_PREFIX" > "$REPO_ROOT/environment.lock.yml"
python -m pip freeze                 > "$REPO_ROOT/requirements.lock.txt"

echo "==> Wrote:"
echo "    $REPO_ROOT/environment.lock.yml"
echo "    $REPO_ROOT/requirements.lock.txt"
