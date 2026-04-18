#!/bin/bash
#SBATCH --job-name=tree_step8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# Plot the calibrated tree annotated with per-branch rate statistics.
egt branch-stats-tree \
    -n "../step7_tree_analysis_branchstatsvtime/branch_stats_output/modified_node_list.tsv" \
    -e "../step7_tree_analysis_branchstatsvtime/branch_stats_output/modified_edge_list.tsv" \
    -p "../step4_persp_chr/perspchrom.tsv"
