#!/bin/bash
#SBATCH --job-name=tree_step9_pmc
#SBATCH --cpus-per-task=1
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

# Fourier analysis on the Protostomia minus Clitellata custom clade subset
# produced by step7/run_protost_minus_clitellata.sh.
CUSTOM_DIR="../step7_tree_analysis_branchstatsvtime/branch_stats_output/custom_clade_analyses/Protostomia_33317_minus_Clitellata_42113"
OUTDIR="${CUSTOM_DIR}/fourier_analysis"
mkdir -p "$OUTDIR"

WEIGHTED_FILE="${CUSTOM_DIR}/Protostomia_minus_Clitellata_changes_vs_age_weighted.tsv"
if [ ! -f "$WEIGHTED_FILE" ]; then
    echo "ERROR: weighted file not found: $WEIGHTED_FILE" >&2
    exit 1
fi

egt fourier-of-rates \
    --rates     "$WEIGHTED_FILE" \
    --agecol    age \
    --ratecol   fusion_rate_at_this_age_mean \
    --polynomial 3 \
    --outprefix "${OUTDIR}/Protostomia_minus_Clitellata_fusions"

egt fourier-of-rates \
    --rates     "$WEIGHTED_FILE" \
    --agecol    age \
    --ratecol   dispersal_rate_at_this_age_mean \
    --polynomial 3 \
    --outprefix "${OUTDIR}/Protostomia_minus_Clitellata_dispersals"
