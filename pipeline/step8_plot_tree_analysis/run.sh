#!/bin/bash
#SBATCH --job-name=tree_step8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
STEP_NAME=step8_plot_tree_analysis
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/pipeline/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# Plot the calibrated tree annotated with per-branch rate statistics.
egt branch-stats-tree \
    -n "../step7_tree_analysis_branchstatsvtime/branch_stats_output/modified_node_list.tsv" \
    -e "../step7_tree_analysis_branchstatsvtime/branch_stats_output/modified_edge_list.tsv" \
    -p "../step4_persp_chr/perspchrom.tsv"
