#!/bin/bash
#SBATCH --job-name=tree_step7
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
STEP_NAME=step7_tree_analysis_branchstatsvtime
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/pipeline/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# Branch-wise ALG-change rates against geologic time + extinction intensity.
# Produces: branch_stats_output/{modified_node_list.tsv, modified_edge_list.tsv,
# per_clade_analyses/…, custom_clade_analyses/…} consumed by steps 8 and 9.

# (a) Main run — no clade filter.
echo "=== step 7 (a): main branch-stats-vs-time ==="
egt branch-stats-vs-time \
    -e "../step2_divergence_times/${TREE_PREFIX}.edge_information.tsv" \
    -n "../step2_divergence_times/${TREE_PREFIX}.node_information.tsv" \
    -s "../step6_perspchrom_df_to_tree/statsdf.tsv" \
    -i "$EXTINCTION_INTENSITY" \
    -t "$CORES_STEP7"

# (b) Clade variant: Protostomia (33317) minus Clitellata (42113).
#     Annelida ancestral ALG dispersals are dominated by Clitellata; this
#     subset probes whether the protostome signal holds without them.
echo "=== step 7 (b): Protostomia minus Clitellata ==="
egt branch-stats-vs-time \
    -e "../step2_divergence_times/${TREE_PREFIX}.edge_information.tsv" \
    -n "../step2_divergence_times/${TREE_PREFIX}.node_information.tsv" \
    -s "../step6_perspchrom_df_to_tree/statsdf.tsv" \
    -i "$EXTINCTION_INTENSITY" \
    -t "$CORES_STEP7" \
    --analyze_single_clade 33317 \
    --exclude_subclades 42113

# (c) Clade variant: Vertebrata (7742) minus Teleostei (32443).
#     Teleosts dominate the vertebrate fusion rate signal; this subset
#     probes the non-teleost contribution.
echo "=== step 7 (c): Vertebrata minus Teleostei ==="
egt branch-stats-vs-time \
    -e "../step2_divergence_times/${TREE_PREFIX}.edge_information.tsv" \
    -n "../step2_divergence_times/${TREE_PREFIX}.node_information.tsv" \
    -s "../step6_perspchrom_df_to_tree/statsdf.tsv" \
    -i "$EXTINCTION_INTENSITY" \
    -t "$CORES_STEP7" \
    --analyze_single_clade 7742 \
    --exclude_subclades 32443
