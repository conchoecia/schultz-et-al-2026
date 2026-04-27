#!/bin/bash
#SBATCH --job-name=stable_pair_v5
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=0-02:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

cd /lisc/data/scratch/molevo/dts/git/schultz-et-al-2026/post_analyses/stable_pair_examples
export GENOMES_DIR=/lisc/data/scratch/molevo/dts/ODP_genomes/GenDB_annotated_chr/odp_ncbi_genome_db/output/source_data/annotated_genomes
PY=/lisc/data/scratch/molevo/dts/conda_envs/egt/bin/python

GENOME_YAML=/lisc/data/scratch/molevo/dts/git/schultz-et-al-2026/genome_database/genome_list.yaml
RBH=/lisc/data/scratch/molevo/dts/git/schultz-et-al-2026/reference_data/BCnSSimakov2022.rbh
OUT=out

echo "=== Stage 1: select_examples (Vertebrata + Arthropoda percentile-spread) ==="
$PY select_examples.py --out-dir "$OUT" --occupancy-min 0.5 \
    --require-anchor Arthropoda:Drosophilamelanogaster-7227-GCF000001215.4

echo "=== Stage 1.5: recompute_distances across 10+10 anchor species ==="
$PY recompute_distances.py --selected "$OUT/selected_examples.tsv" --out-dir "$OUT"

echo "=== Stage 1.6: shared_tfs (intersection across focal-species) ==="
$PY shared_tfs.py --distances "$OUT/selected_examples_distances.tsv" \
    --epigenomics "$OUT/epigenomics" --out-dir "$OUT"

echo "=== Stage 1.7: tf_enrichment (per-pair z/p/q + top-N TFs) ==="
$PY tf_enrichment.py --distances "$OUT/selected_examples_distances.tsv" \
    --epigenomics "$OUT/epigenomics" --out-dir "$OUT" --top-n 10

echo "=== Stage 2: plot_distribution ==="
$PY plot_distribution.py \
    --candidates "$OUT/candidate_pairs.tsv" \
    --selected   "$OUT/selected_examples.tsv" \
    --recomputed "$OUT/selected_examples_recomputed.tsv" \
    --focal      Vertebrata Arthropoda \
    --out        "$OUT/stable_pair_distance_distribution.pdf"

echo "=== Stage 3: plot_examples (10+10 species locus tracks, normalized + absolute) ==="
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

echo "=== per-pair coverage report ==="
awk -F'\t' 'NR==1 || ($1=="Vertebrata" || $1=="Arthropoda")' \
    "$OUT/selected_examples_recomputed.tsv" | \
    awk -F'\t' 'NR==1 {print "clade\tpair\tn_attempted\tn_same_scaffold\tTSS\torient\tbin"} NR>1 {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $3, $4, $5, $6, $8, $9}' | column -t -s$'\t'

echo done
