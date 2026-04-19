#!/bin/bash
#SBATCH --job-name=defining_features_qc_plots
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-01:00:00
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

cd "$SCRIPT_DIR/out"

# Parallel QC PDF generation: 8 concurrent plots x ~2 min each ~= ~7 min total.
ls *_unique_pair_df.tsv.gz | xargs -n1 -P8 -I{} python -c "
from egt.defining_features_qc_plots import write_qc_plots_from_tsv
import sys
out = write_qc_plots_from_tsv('{}')
print(f'[ok] wrote {out}')
"

echo ""
echo "=== PDFs produced ==="
ls -la *_unique_pair_df_qc.pdf | wc -l
ls -la *_unique_pair_df_qc.pdf
