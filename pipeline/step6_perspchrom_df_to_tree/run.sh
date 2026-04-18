#!/bin/bash
#SBATCH --job-name=tree_step6
#SBATCH --cpus-per-task=25
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
STEP_NAME=step6_perspchrom_df_to_tree
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/pipeline/$STEP_NAME"
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
