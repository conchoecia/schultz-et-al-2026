#!/bin/bash
#SBATCH --job-name=pigeonhole_check
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT"

# Clade set is paper-specific; edit as needed.
CLADES="Metazoa:33208,Porifera:6040,Cnidaria:6073,Bilateria:33213,Protostomia:33317,Deuterostomia:33511,Nematoda:6231,Arthropoda:6656,Mollusca:6447,Annelida:6340,Vertebrata:7742"

egt pigeonhole-check \
    --presence-fusions "../../pipeline/step4_persp_chr/per_species_ALG_presence_fusions.tsv" \
    --chrom-counts     "$CHROM_COUNTS" \
    --clade-groupings  "$CLADES" \
    --n-simulations    10000 \
    --out-dir          "$OUT"
