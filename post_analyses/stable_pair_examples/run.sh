#!/bin/bash
#SBATCH --job-name=stable_pair_examples
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Build the stable-pair-examples figure: distance distribution across all
# qualifying stable pairs, plus per-clade locus tracks for the top-1 pair
# per (clade, distance bin) for the two focal clades (Vertebrata, Diptera).
# Stage 0 caches ReMap + cCRE + REDfly regulatory BED tracks.

STEP_NAME=stable_pair_examples
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
cd "$SCRIPT_DIR"

GENOME_YAML="$REPO_ROOT/genome_database/genome_list.yaml"
RBH="$REPO_ROOT/reference_data/BCnSSimakov2022.rbh"

# Where the annotated_genomes tree lives in place. Overridable from the env
# or a sibling symlink (see genome_database/annotated_genomes_link.placeholder).
if [ -z "${GENOMES_DIR:-}" ]; then
    if [ -L "$REPO_ROOT/genome_database/annotated_genomes_link" ]; then
        export GENOMES_DIR="$REPO_ROOT/genome_database/annotated_genomes_link"
    else
        export GENOMES_DIR="/lisc/data/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes"
    fi
fi
echo "GENOMES_DIR=$GENOMES_DIR"

echo "=== Stage 0: probe + epigenomics fetch ==="
python loci_io.py --probe --genome-yaml "$GENOME_YAML" || true
python epigenomics_fetch.py --out-dir "$OUT/epigenomics" ${NO_FETCH:+--no-fetch}

echo "=== Stage 1: select_examples ==="
python select_examples.py --out-dir "$OUT"

echo "=== Stage 1.5: recompute_distances (TSS-aware) ==="
python recompute_distances.py \
    --selected     "$OUT/selected_examples.tsv" \
    --genome-yaml  "$GENOME_YAML" \
    --rbh          "$RBH" \
    --out-dir      "$OUT" || echo "WARN: recompute_distances failed; downstream falls back to mean_in"

echo "=== Stage 1.6: shared_tfs (intersection of TFs at both anchor TSSes) ==="
python shared_tfs.py \
    --distances   "$OUT/selected_examples_distances.tsv" \
    --epigenomics "$OUT/epigenomics" \
    --out-dir     "$OUT" || echo "WARN: shared_tfs failed; Stage 3 will not annotate shared TFs"

echo "=== Stage 2: plot_distribution ==="
python plot_distribution.py \
    --candidates "$OUT/candidate_pairs.tsv" \
    --selected   "$OUT/selected_examples.tsv" \
    --recomputed "$OUT/selected_examples_recomputed.tsv" \
    --focal      Vertebrata Arthropoda \
    --out        "$OUT/stable_pair_distance_distribution.pdf"

echo "=== Stage 3: plot_examples ==="
python plot_examples.py \
    --selected     "$OUT/selected_examples.tsv" \
    --recomputed   "$OUT/selected_examples_recomputed.tsv" \
    --distances    "$OUT/selected_examples_distances.tsv" \
    --top-per-bin  1 \
    --genome-yaml  "$GENOME_YAML" \
    --rbh          "$RBH" \
    --epigenomics  "$OUT/epigenomics" \
    --tree         "$REPO_ROOT/pipeline/step2_divergence_times/BCnS_2026.calibrated_tree.nwk" \
    --focal        Vertebrata Arthropoda \
    --out          "$OUT/stable_pair_example_loci.pdf"

echo "done. outputs:"
ls -lh "$OUT"
