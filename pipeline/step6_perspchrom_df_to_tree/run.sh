#!/bin/bash
#SBATCH --job-name=tree_step6
#SBATCH --cpus-per-task=25
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# Monte-Carlo-supported mapping of perspective-chromosome changes onto the tree.
# Produces: statsdf.tsv (+ per-clade simulation PDFs unless --skip-traces).
egt perspchrom-df-to-tree \
    "../step4_persp_chr/per_species_ALG_presence_fusions.tsv" \
    "$ALG_RBH" \
    --num-simulations "$NUM_SIMULATIONS" \
    --num-processes "$CORES_STEP6" \
    --skip-traces
