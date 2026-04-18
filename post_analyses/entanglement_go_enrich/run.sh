#!/bin/bash
#SBATCH --job-name=entanglement_go_enrich
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

GAF="$OUT/goa_human.gaf.gz"
if [ ! -s "$GAF" ]; then
    echo "ERROR: $GAF missing. Run download_go.sh first." >&2
    exit 1
fi

FAMILY_MAP="../build_family_naming_map/out/bcns_family_to_human_gene.tsv"
ENTANGLED="../entanglement_browse/out/entangled_pairs_per_clade.tsv"

for f in "$FAMILY_MAP" "$ENTANGLED"; do
    if [ ! -s "$f" ]; then
        echo "ERROR: $f missing. Run the upstream post-analysis first." >&2
        exit 1
    fi
done

egt entanglement-go-enrich \
    --alg-rbh         "$ALG_RBH" \
    --family-gene-map "$FAMILY_MAP" \
    --entangled-pairs "$ENTANGLED" \
    --human-go        "$GAF" \
    --fdr             0.05 \
    --out-dir         "$OUT"
