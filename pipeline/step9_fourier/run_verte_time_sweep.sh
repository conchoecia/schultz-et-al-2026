#!/bin/bash
#SBATCH --job-name=tree_step9_vsweep
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

# Vertebrata-minus-Teleostei Fourier analysis with a max_time sweep
# (250..550 Mya in 10-Mya steps). Produces the time-window-support plots.
CUSTOM_DIR="../step7_tree_analysis_branchstatsvtime/branch_stats_output/custom_clade_analyses/Vertebrata_7742_minus_Teleostei_32443"
OUTDIR="${CUSTOM_DIR}/fourier_analysis"
mkdir -p "$OUTDIR"

WEIGHTED_FILE="${CUSTOM_DIR}/Vertebrata_minus_Teleostei_changes_vs_age_weighted.tsv"
if [ ! -f "$WEIGHTED_FILE" ]; then
    echo "ERROR: weighted file not found: $WEIGHTED_FILE" >&2
    exit 1
fi

for MAX_TIME in $(seq 250 10 550); do
    echo "=== max_time = ${MAX_TIME} Mya ==="
    egt fourier-of-rates \
        --rates     "$WEIGHTED_FILE" \
        --agecol    age \
        --ratecol   fusion_rate_at_this_age_mean \
        --min_time  1 \
        --max_time  "$MAX_TIME" \
        --polynomial 3 \
        --outprefix "${OUTDIR}/Vertebrata_minus_Teleostei_fusions${MAX_TIME}"

    egt fourier-of-rates \
        --rates     "$WEIGHTED_FILE" \
        --agecol    age \
        --ratecol   dispersal_rate_at_this_age_mean \
        --min_time  1 \
        --max_time  "$MAX_TIME" \
        --polynomial 3 \
        --outprefix "${OUTDIR}/Vertebrata_minus_Teleostei_dispersals${MAX_TIME}"
done
