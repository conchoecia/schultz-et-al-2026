#!/bin/bash
#SBATCH --job-name=divergence_vs_dispersal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
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

DIV_TSV="$SCRIPT_DIR/out/divergence.tsv"
if [ ! -s "$DIV_TSV" ]; then
    echo "ERROR: $DIV_TSV is missing. Produce it first via compute_divergence.sh." >&2
    exit 1
fi

CLADES="Metazoa:33208,Porifera:6040,Cnidaria:6073,Bilateria:33213,Protostomia:33317,Deuterostomia:33511,Nematoda:6231,Arthropoda:6656,Mollusca:6447,Annelida:6340,Vertebrata:7742"

egt divergence-vs-dispersal \
    --divergence-tsv    "$DIV_TSV" \
    --divergence-column median_bitscore \
    --presence-fusions  "../../pipeline/step4_persp_chr/per_species_ALG_presence_fusions.tsv" \
    --clade-groupings   "$CLADES" \
    --out-dir           "$OUT"
