#!/bin/bash
#SBATCH --job-name=tree_step7_vmt
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
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

# Vertebrata (taxid 7742) minus Teleostei (taxid 32443) — clade subset.
egt branch-stats-vs-time \
    -e "../step2_divergence_times/${TREE_PREFIX}.edge_information.tsv" \
    -n "../step2_divergence_times/${TREE_PREFIX}.node_information.tsv" \
    -s "../step6_perspchrom_df_to_tree/statsdf.tsv" \
    -i "${REPO_ROOT}/${EXTINCTION_INTENSITY}" \
    -t "$CORES_STEP7" \
    --analyze_single_clade 7742 \
    --exclude_subclades 32443
