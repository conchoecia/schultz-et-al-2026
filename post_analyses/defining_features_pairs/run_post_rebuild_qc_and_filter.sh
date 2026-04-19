#!/bin/bash
#SBATCH --job-name=qc_and_unique_pairs
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Post-rebuild cleanup: for any *_unique_pair_df.tsv.gz missing a
# matching *_qc.pdf (i.e. the 9 we quarantined), regenerate the
# detailed QC PDF. Then rerun summarize_qc on all 28 (regenerates
# all *_qc_summary.pdf + the cross-clade summary). Then build
# the filtered unique_pairs.tsv.gz.

STEP_NAME=defining_features_pairs
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR/out"

echo "=== [1/3] detailed QC PDFs for any TSV missing one ==="
MISSING=()
for t in *_unique_pair_df.tsv.gz; do
    pdf="${t%.tsv.gz}_qc.pdf"
    [ ! -s "$pdf" ] && MISSING+=("$t")
done
echo "  ${#MISSING[@]} missing detailed _qc.pdf"
if [ ${#MISSING[@]} -gt 0 ]; then
    printf '%s\n' "${MISSING[@]}" | xargs -n1 -P8 -I{} python -c "
from egt.defining_features_qc_plots import write_qc_plots_from_tsv
out = write_qc_plots_from_tsv('{}')
print(f'[ok] wrote {out}')
"
fi

echo ""
echo "=== [2/3] summarize_qc on all 28 ==="
cd "$SCRIPT_DIR"
python summarize_qc.py out

echo ""
echo "=== [3/3] build unique_pairs.tsv.gz ==="
python build_unique_pairs_tsv.py \
    --clade-stats-dir out \
    --rbh-file /lisc/data/scratch/molevo/dts/manifold/post_decay_202401/BCnSSimakov2022.rbh \
    --pair-combination-path /lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP/combo_to_index.txt \
    --output out/unique_pairs.tsv.gz \
    --sigma 2

echo ""
echo "=== done ==="
ls -la out/unique_pairs.tsv.gz out/all_clades_selection_summary.pdf
