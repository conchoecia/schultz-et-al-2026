#!/bin/bash
#SBATCH --job-name=stable_pair_render_only
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=0-00:20:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Fast iteration: re-render Stage 3 only, reusing all upstream cached outputs.
# Use this when iterating on the tree subplot or other plot-only code paths.

cd /lisc/data/scratch/molevo/dts/git/schultz-et-al-2026/post_analyses/stable_pair_examples
export GENOMES_DIR=/lisc/data/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes
PY=/lisc/data/scratch/molevo/dts/conda_envs/egt/bin/python
GENOME_YAML=/lisc/data/scratch/molevo/dts/git/schultz-et-al-2026/genome_database/genome_list.yaml
RBH=/lisc/data/scratch/molevo/dts/git/schultz-et-al-2026/reference_data/BCnSSimakov2022.rbh
OUT=out

echo "=== Stage 3 only: plot_examples ==="
$PY plot_examples.py \
    --selected     "$OUT/selected_examples.tsv" \
    --recomputed   "$OUT/selected_examples_recomputed.tsv" \
    --distances    "$OUT/selected_examples_distances.tsv" \
    --top-per-bin  1 \
    --genome-yaml  "$GENOME_YAML" \
    --rbh          "$RBH" \
    --epigenomics  "$OUT/epigenomics" \
    --tree         /lisc/data/scratch/molevo/dts/git/schultz-et-al-2026/pipeline/step2_divergence_times/BCnS_2026.calibrated_tree.nwk \
    --focal        Vertebrata Arthropoda \
    --scale        both \
    --tf-enrichment-pairs "$OUT/tf_enrichment_per_pair.tsv" \
    --tf-enrichment-top   "$OUT/tf_enrichment_top_tfs.tsv" \
    --out          "$OUT/stable_pair_example_loci.pdf"

echo done
