#!/bin/bash
# Create (or update) the `egt-repro` conda env from environment.yml.
#
# Works with conda, mamba, or micromamba; prefers mamba/micromamba if
# available for speed. On HPCs that provide conda via lmod (e.g. LiSC's
# Conda/Miniforge3 module), load that module first.
#
# Usage:
#   bash bin/setup_env.sh [env-path-or-name]
#
# With no argument, creates a named env `egt-repro`. With an argument, uses
# it as the env prefix (conda -p) or name (conda -n, if no slashes).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_SPEC="$REPO_ROOT/environment.yml"

if [ ! -f "$ENV_SPEC" ]; then
    echo "ERROR: $ENV_SPEC not found." >&2
    exit 1
fi

TARGET="${1:-egt-repro}"

# Pick the fastest available resolver.
if command -v micromamba >/dev/null 2>&1; then
    MAMBA=micromamba
elif command -v mamba >/dev/null 2>&1; then
    MAMBA=mamba
elif command -v conda >/dev/null 2>&1; then
    MAMBA=conda
else
    echo "ERROR: neither conda, mamba, nor micromamba is on PATH." >&2
    echo "       On LiSC: module load Conda/Miniforge3" >&2
    exit 1
fi

# Name vs prefix: if target looks like a path, use -p; else -n.
if [[ "$TARGET" == */* ]]; then
    FLAG="-p"
else
    FLAG="-n"
fi

echo "==> Using $MAMBA"
echo "==> Target env: $FLAG $TARGET"
echo "==> Spec:       $ENV_SPEC"
echo

# Check if env already exists; update vs create.
if $MAMBA env list | awk '{print $1, $2}' | grep -qE "(^| )${TARGET}( |$)"; then
    echo "==> Env exists — updating"
    $MAMBA env update "$FLAG" "$TARGET" -f "$ENV_SPEC" --prune
else
    echo "==> Env does not exist — creating"
    $MAMBA env create "$FLAG" "$TARGET" -f "$ENV_SPEC"
fi

echo
echo "==> Done."
echo
echo "Activate with:"
if [ "$FLAG" = "-p" ]; then
    echo "  conda activate $TARGET"
else
    echo "  conda activate $TARGET"
fi
echo
echo "Verify egt installed:"
echo "  egt --help"
