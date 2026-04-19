#!/bin/bash
#SBATCH --job-name=defining_features_new
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Rerun defining-features against the newly-verified canonical COO.
# Runs after grid_verify_new_coo has promoted allsamples.coo.REBUILD to
# allsamples.coo.npz. Uses the sparse-native rewrite so mem is ~64G
# (was 400G OOMing on the densified path).

STEP_NAME=defining_features_pairs
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT"
cd "$OUT"

COO_ROOT="${COO_ROOT:-/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP}"
COO_PATH="$COO_ROOT/allsamples.coo.npz"
SAMPLE_DF="$COO_ROOT/sampledf.tsv"
COO_COMBO="$COO_ROOT/combo_to_index.txt"

for f in "$COO_PATH" "$SAMPLE_DF" "$COO_COMBO"; do
    [ -s "$f" ] || { echo "ERROR: missing input $f" >&2; exit 1; }
done

TAXIDS=(
    # Original 20 (from SupplementaryTable_16).
    6340 33213 6544 6073 6606 10197 215450 33511 7147 7586
    6448 50557 6447 6231 33340 6040 33317 32584 2697495 7742
    # New 8.
    6656 6605 7711 40674 32443 42113 6960 6157
)
TAXID_CSV=$(IFS=,; echo "${TAXIDS[*]}")

echo "=== defining-features (post-rebuild) ==="
echo "  COO:         $COO_PATH"
echo "  sample_df:   $SAMPLE_DF"
echo "  combo_index: $COO_COMBO"
echo "  taxids:      ${#TAXIDS[@]} entries"
echo

egt defining-features \
    --coo_path             "$COO_PATH" \
    --sample_df_path       "$SAMPLE_DF" \
    --coo_combination_path "$COO_COMBO" \
    --taxid_list           "$TAXID_CSV"

if [ -f unique_pairs.tsv.gz ]; then
    echo "produced: $(readlink -f unique_pairs.tsv.gz)"
    zcat unique_pairs.tsv.gz | head -1
    ncla=$(zcat unique_pairs.tsv.gz | awk -F'\t' 'NR>1 {print $1}' | sort -u | wc -l)
    nrow=$(zcat unique_pairs.tsv.gz | wc -l)
    echo "clades:$ncla  rows:$nrow"
fi
