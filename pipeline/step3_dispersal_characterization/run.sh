#!/bin/bash
#SBATCH --job-name=tree_step3
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
STEP_NAME=step3_dispersal_characterization
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/pipeline/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# Pairwise evolutionary-decay curves against divergence times.
egt decay-pairwise \
    --divergence_file "../step2_divergence_times/${TREE_PREFIX}.divergence_times.txt" \
    --ALG_rbh "$ALG_RBH" \
    --ALG_rbh_dir "$RBH_DIR" \
    --ALGname "$ALG_NAME" \
    -c "$GENOME_CONFIG_YAML" \
    -d "$SYNTENY_DIR" \
    -T "$TARGET_SPECIES" \
    --min_scaf_size "$MIN_SCAF_SIZE" \
    --bin_size "$BIN_SIZE" \
    --cache_dir ./cache

# Optional: ALG conservation plot (commented out in the original).
# Uncomment after step 4 has produced per_species_ALG_presence_fusions.tsv.
#
# egt alg-dispersion \
#     -d "$RBH_DIR" \
#     -a "$ALG_RBH" \
#     -n "$ALG_NAME" \
#     --metadata "../step4_persp_chr/per_species_ALG_presence_fusions.tsv" \
#     -o ./alg_dispersion_plots \
#     --lineage_col taxidstring \
#     -m 0.05
