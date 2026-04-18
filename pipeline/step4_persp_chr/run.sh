#!/bin/bash
#SBATCH --job-name=tree_step4
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
STEP_NAME=step4_persp_chr
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/pipeline/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# Perspective-chromosome / ALG-fusion inference across species.
# Produces: perspchrom.tsv, per_species_ALG_presence_fusions.tsv, tree1.tsv.gz,
#           locdf.tsv, changestring_checkpoints/, unique_changes_summary.tsv
egt alg-fusions \
    -d "$RBH_DIR" \
    -a "$ALG_NAME" \
    -r "$ALG_RBH" \
    -t "../step2_divergence_times/${TREE_PREFIX}.node_information.tsv" \
    --ncores "$CORES_STEP4" \
    --parallel
