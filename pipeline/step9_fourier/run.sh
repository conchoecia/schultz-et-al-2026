#!/bin/bash
#SBATCH --job-name=tree_step9
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
STEP_NAME=step9_fourier
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/pipeline/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# (a) Collapsed-tree visualization.
echo "=== step 9 (a): collapsed-tree ==="
egt collapsed-tree \
    -n "../step7_tree_analysis_branchstatsvtime/branch_stats_output/modified_node_list.tsv" \
    -e "../step7_tree_analysis_branchstatsvtime/branch_stats_output/modified_edge_list.tsv" \
    -o collapsed_tree.pdf

# (b) Fourier analysis of rates across 10 standard clades (step 7 main output).
INDIR="../step7_tree_analysis_branchstatsvtime/branch_stats_output/per_clade_analyses"
OUTDIR="${INDIR}/fourier_analysis"
mkdir -p "$OUTDIR"

CLADES=(Mollusca Spiralia Deuterostomia Vertebrata Metazoa Myriazoa Porifera Cnidaria Protostomia Arthropoda)

for CLADE in "${CLADES[@]}"; do
    echo "=== step 9 (b): Fourier ${CLADE} ==="
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

# (c) Fourier on the Protostomia-minus-Clitellata custom clade subset
#     produced by step 7 (b).
PMC_DIR="../step7_tree_analysis_branchstatsvtime/branch_stats_output/custom_clade_analyses/Protostomia_33317_minus_Clitellata_42113"
PMC_OUT="${PMC_DIR}/fourier_analysis"
PMC_WEIGHTED="${PMC_DIR}/Protostomia_minus_Clitellata_changes_vs_age_weighted.tsv"
if [ -f "$PMC_WEIGHTED" ]; then
    echo "=== step 9 (c): Fourier Protostomia-minus-Clitellata ==="
    mkdir -p "$PMC_OUT"
    egt fourier-of-rates \
        --rates     "$PMC_WEIGHTED" \
        --agecol    age \
        --ratecol   fusion_rate_at_this_age_mean \
        --polynomial 3 \
        --outprefix "${PMC_OUT}/Protostomia_minus_Clitellata_fusions"
    egt fourier-of-rates \
        --rates     "$PMC_WEIGHTED" \
        --agecol    age \
        --ratecol   dispersal_rate_at_this_age_mean \
        --polynomial 3 \
        --outprefix "${PMC_OUT}/Protostomia_minus_Clitellata_dispersals"
else
    echo "  (skipping step 9 (c): $PMC_WEIGHTED not present)"
fi

# (d) Fourier on the Vertebrata-minus-Teleostei custom clade subset
#     produced by step 7 (c), plus a max-time window sweep (250..550 Mya
#     in 10-Mya steps) used for the time-window-support figures.
VMT_DIR="../step7_tree_analysis_branchstatsvtime/branch_stats_output/custom_clade_analyses/Vertebrata_7742_minus_Teleostei_32443"
VMT_OUT="${VMT_DIR}/fourier_analysis"
VMT_WEIGHTED="${VMT_DIR}/Vertebrata_minus_Teleostei_changes_vs_age_weighted.tsv"
if [ -f "$VMT_WEIGHTED" ]; then
    echo "=== step 9 (d): Fourier Vertebrata-minus-Teleostei ==="
    mkdir -p "$VMT_OUT"
    egt fourier-of-rates \
        --rates     "$VMT_WEIGHTED" \
        --agecol    age \
        --ratecol   fusion_rate_at_this_age_mean \
        --polynomial 3 \
        --outprefix "${VMT_OUT}/Vertebrata_minus_Teleostei_fusions"
    egt fourier-of-rates \
        --rates     "$VMT_WEIGHTED" \
        --agecol    age \
        --ratecol   dispersal_rate_at_this_age_mean \
        --polynomial 3 \
        --outprefix "${VMT_OUT}/Vertebrata_minus_Teleostei_dispersals"

    echo "=== step 9 (e): Vertebrata-minus-Teleostei max-time sweep 250..550 Mya ==="
    for MAX_TIME in $(seq 250 10 550); do
        echo "--- max_time = ${MAX_TIME} Mya ---"
        egt fourier-of-rates \
            --rates     "$VMT_WEIGHTED" \
            --agecol    age \
            --ratecol   fusion_rate_at_this_age_mean \
            --min_time  1 \
            --max_time  "$MAX_TIME" \
            --polynomial 3 \
            --outprefix "${VMT_OUT}/Vertebrata_minus_Teleostei_fusions${MAX_TIME}"
        egt fourier-of-rates \
            --rates     "$VMT_WEIGHTED" \
            --agecol    age \
            --ratecol   dispersal_rate_at_this_age_mean \
            --min_time  1 \
            --max_time  "$MAX_TIME" \
            --polynomial 3 \
            --outprefix "${VMT_OUT}/Vertebrata_minus_Teleostei_dispersals${MAX_TIME}"
    done
else
    echo "  (skipping step 9 (d,e): $VMT_WEIGHTED not present)"
fi
