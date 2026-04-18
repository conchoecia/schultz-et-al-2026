#!/bin/bash
#SBATCH --job-name=tree_step9
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# (a) Collapsed-tree visualization.
egt collapsed-tree \
    -n "../step7_tree_analysis_branchstatsvtime/branch_stats_output/modified_node_list.tsv" \
    -e "../step7_tree_analysis_branchstatsvtime/branch_stats_output/modified_edge_list.tsv" \
    -o collapsed_tree.pdf

# (b) Fourier analysis of rates, per clade.
INDIR="../step7_tree_analysis_branchstatsvtime/branch_stats_output/per_clade_analyses"
OUTDIR="${INDIR}/fourier_analysis"
mkdir -p "$OUTDIR"

CLADES=(Mollusca Spiralia Deuterostomia Vertebrata Metazoa Myriazoa Porifera Cnidaria Protostomia Arthropoda)

for CLADE in "${CLADES[@]}"; do
    echo "=== Fourier: ${CLADE} ==="
    WEIGHTED_FILE=$(find "$INDIR" -name "${CLADE}_*_changes_vs_age_weighted.tsv" -type f | head -1)
    if [ -z "$WEIGHTED_FILE" ]; then
        echo "  (no weighted file for ${CLADE}, skipping)"
        continue
    fi

    egt fourier-of-rates \
        --rates     "$WEIGHTED_FILE" \
        --agecol    age \
        --ratecol   fusion_rate_at_this_age_mean \
        --polynomial 3 \
        --outprefix "${OUTDIR}/${CLADE}_fusions"

    egt fourier-of-rates \
        --rates     "$WEIGHTED_FILE" \
        --agecol    age \
        --ratecol   dispersal_rate_at_this_age_mean \
        --polynomial 3 \
        --outprefix "${OUTDIR}/${CLADE}_dispersals"
done
