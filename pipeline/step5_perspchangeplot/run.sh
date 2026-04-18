#!/bin/bash
#SBATCH --job-name=tree_step5
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# Chromosome number vs ALG-change count scatter.
egt chrom-number-vs-changes \
    "../step4_persp_chr/per_species_ALG_presence_fusions.tsv" \
    "$CHROM_COUNTS"
