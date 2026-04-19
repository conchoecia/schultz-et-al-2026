#!/bin/bash
#SBATCH --job-name=rebuild_coo_202509
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

STEP_NAME=defining_features_pairs
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

PROD=/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP
OUT="$PROD/allsamples.coo.REBUILD_20260419.npz"

echo "[rebuild] input sampledf    : $PROD/sampledf.tsv"
echo "[rebuild] input algcomboix  : $PROD/combo_to_index.txt"
echo "[rebuild] input gb.gz dir   : $PROD/distance_matrices/"
echo "[rebuild] output            : $OUT"
echo "[rebuild] egt version       : $(egt --version 2>/dev/null || echo '?')"
echo "[rebuild] egt git HEAD      : $(git -C /lisc/data/scratch/molevo/dts/git/egt log -1 --format='%h %s' 2>/dev/null)"
echo ""

cd "$SCRIPT_DIR"
time egt phylotreeumap combine-distances \
    --sampledf "$PROD/sampledf.tsv" \
    --algcomboix "$PROD/combo_to_index.txt" \
    --output "$OUT"

echo ""
echo "[rebuild] output npz size:"
ls -la "$OUT"
