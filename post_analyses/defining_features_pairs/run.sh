#!/bin/bash
#SBATCH --job-name=defining_features
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Regenerate the per-clade unique-pairs table (the basis of
# SupplementaryTable_16) against the full 5,821-species 202509 dataset.
# Runs `egt defining-features` over the pre-built 5.9 GB COO + its
# companions and emits unique_pairs.tsv.gz into this step's out/ dir.
#
# Clade list: the 20 original SupplementaryTable_16 clades plus 8 new
# ones (Arthropoda, Cephalopoda, Chordata, Mammalia, Teleostei,
# Clitellata, Hexapoda, Platyhelminthes). Downstream plots (go_enrichment
# _sweep) pick this up automatically via --supp-table override.

STEP_NAME=defining_features_pairs
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT"
cd "$OUT"

# Pre-built COO + companions for the 202509 full dataset.
COO_ROOT="${COO_ROOT:-/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP}"
COO_PATH="$COO_ROOT/allsamples.coo.npz"
SAMPLE_DF="$COO_ROOT/sampledf.tsv"
COO_COMBO="$COO_ROOT/combo_to_index.txt"

for f in "$COO_PATH" "$SAMPLE_DF" "$COO_COMBO"; do
    [ -s "$f" ] || { echo "ERROR: missing input $f" >&2; exit 1; }
done

# 28-taxid list: 20 originals + 8 additions.
TAXIDS=(
    # Original 20.
    6340        # Annelida
    33213       # Bilateria
    6544        # Bivalvia
    6073        # Cnidaria
    6606        # Coleoidea
    10197       # Ctenophora
    215450      # Decapodiformes
    33511       # Deuterostomia
    7147        # Diptera
    7586        # Echinodermata
    6448        # Gastropoda
    50557       # Insecta
    6447        # Mollusca
    6231        # Nematoda
    33340       # Neoptera
    6040        # Porifera
    33317       # Protostomia
    32584       # Scaphopoda
    2697495     # Spiralia
    7742        # Vertebrata
    # New 8.
    6656        # Arthropoda
    6605        # Cephalopoda
    7711        # Chordata
    40674       # Mammalia
    32443       # Teleostei
    42113       # Clitellata
    6960        # Hexapoda
    6157        # Platyhelminthes
)
TAXID_CSV=$(IFS=,; echo "${TAXIDS[*]}")

echo "=== egt defining-features ==="
echo "  COO:         $COO_PATH"
echo "  sample_df:   $SAMPLE_DF"
echo "  combo_index: $COO_COMBO"
echo "  taxids:      ${#TAXIDS[@]} entries"
echo

egt defining-features \
    --coo_path             "$COO_PATH" \
    --sample_df_path       "$SAMPLE_DF" \
    --coo_combination_path "$COO_COMBO" \
    --taxid_list           "$TAXID_CSV"

# Rename in place — egt writes unique_pairs.tsv.gz to cwd.
if [ -f unique_pairs.tsv.gz ]; then
    echo "produced: $(readlink -f unique_pairs.tsv.gz)"
    zcat unique_pairs.tsv.gz | head -1 >&2
    ncla=$(zcat unique_pairs.tsv.gz | awk -F'\t' 'NR>1 {print $1}' | sort -u | wc -l)
    nrow=$(zcat unique_pairs.tsv.gz | wc -l)
    echo "clades:$ncla  rows:$nrow"
fi
