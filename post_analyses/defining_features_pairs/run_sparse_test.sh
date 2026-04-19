#!/bin/bash
#SBATCH --job-name=defining_features_sparse_test
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Correctness + resource-use test for the sparse-native rewrite of
# `egt defining-features`. Runs a single clade (Annelida, taxid 6340)
# into a scratch directory, so its output can be diff'd against the
# reference Annelida_6340 TSV produced by the original (slow) code path.
# Memory budget intentionally much smaller than the 400 GB the dense
# path needs.

STEP_NAME=defining_features_pairs
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

OUT="$SCRIPT_DIR/sparse_test"
mkdir -p "$OUT"
cd "$OUT"

COO_ROOT="${COO_ROOT:-/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP}"
COO_PATH="$COO_ROOT/allsamples.coo.npz"
SAMPLE_DF="$COO_ROOT/sampledf.tsv"
COO_COMBO="$COO_ROOT/combo_to_index.txt"

echo "=== egt defining-features (sparse single-clade test) ==="
echo "  output dir:  $(pwd)"
echo "  COO:         $COO_PATH"
echo

/usr/bin/time -v egt defining-features \
    --coo_path             "$COO_PATH" \
    --sample_df_path       "$SAMPLE_DF" \
    --coo_combination_path "$COO_COMBO" \
    --taxid_list           "6340"

ls -la
