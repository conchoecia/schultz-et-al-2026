#!/bin/bash
#SBATCH --job-name=go_sweep_new
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Run the GO enrichment sweep against the post-rebuild unique_pairs.tsv.gz.
# Just delegates to the existing go_enrichment_sweep/run.sh which
# auto-prefers the rebuilt table when present.

STEP_NAME=go_enrichment_sweep
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"
bash run.sh
