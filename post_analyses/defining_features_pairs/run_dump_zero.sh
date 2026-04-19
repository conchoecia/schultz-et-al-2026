#!/bin/bash
#SBATCH --job-name=dump_zero_pairs
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=0-02:00:00
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

COO_ROOT="${COO_ROOT:-/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP}"

OUT_TSV="$SCRIPT_DIR/zero_pairs_chrom_entries_annelida.tsv"
cd "$SCRIPT_DIR"

python dump_zero_pairs.py \
    --coo          "$COO_ROOT/allsamples.coo.npz" \
    --sampledf     "$COO_ROOT/sampledf.tsv" \
    --combo        "$COO_ROOT/combo_to_index.txt" \
    --rbh-dir      "$RBH_DIR" \
    --out-tsv      "$OUT_TSV" \
    --clade-taxid  6340

ls -la "$OUT_TSV"
head -5 "$OUT_TSV"
