#!/bin/bash
#SBATCH --job-name=verify_zero_distances
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=0-00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Cross-check: find stored-zero entries in the COO and look up both
# orthologs' positions in the per-species RBH file. If the distances
# computed from positions are 0 for the reported cells, the COO's
# "observed zero" cells are what they claim to be (exact-neighbor loci).

STEP_NAME=defining_features_pairs
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

COO_ROOT="${COO_ROOT:-/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP}"

cd "$SCRIPT_DIR"
python verify_zero_distances.py \
    --coo      "$COO_ROOT/allsamples.coo.npz" \
    --sampledf "$COO_ROOT/sampledf.tsv" \
    --combo    "$COO_ROOT/combo_to_index.txt" \
    --rbh-dir  "$RBH_DIR" \
    --max-report 20
