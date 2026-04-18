#!/bin/bash
#SBATCH --job-name=family_naming_map
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

# Find the human per-species RBH file.
HUMAN_RBH=$(ls "${RBH_DIR}"/BCnSSimakov2022_Homosapiens-9606-*.rbh 2>/dev/null | head -1)
if [ -z "$HUMAN_RBH" ]; then
    echo "ERROR: no Homo sapiens RBH under $RBH_DIR" >&2
    exit 1
fi

# Pass 1: direct RBH join. Fast.
egt build-family-naming-map \
    --alg-rbh     "$ALG_RBH" \
    --human-rbh   "$HUMAN_RBH" \
    --output      "$OUT/bcns_family_to_human_gene.tsv"

# To run the optional pass 2 (HMM-consensus DIAMOND for families with no
# direct human hit), see pass2_hmm_consensus.sh in this directory.
