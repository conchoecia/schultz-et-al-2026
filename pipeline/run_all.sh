#!/bin/bash
# Run the full pipeline, stages 1→9, in order. Respects set -e: on any
# stage failure the driver halts. Each stage logs to its own runlog.txt.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$PIPELINE_DIR/.." && pwd)"

if [ ! -f "$REPO_ROOT/config.yaml" ]; then
    echo "ERROR: $REPO_ROOT/config.yaml not found." >&2
    echo "       Copy config.template.yaml to config.yaml and edit it first." >&2
    exit 1
fi

STAGES=(
    step1_generate_newick
    step2_divergence_times
    step3_dispersal_characterization
    step4_persp_chr
    step5_perspchangeplot
    step6_perspchrom_df_to_tree
    step7_tree_analysis_branchstatsvtime
    step8_plot_tree_analysis
    step9_fourier
)

for stage in "${STAGES[@]}"; do
    echo "================================================================"
    echo "  ▶ ${stage}"
    echo "================================================================"
    bash "${PIPELINE_DIR}/${stage}/run.sh" 2>&1 | tee "${PIPELINE_DIR}/${stage}/runlog.txt"
done

echo
echo "================================================================"
echo "  ✓ pipeline complete"
echo "================================================================"
