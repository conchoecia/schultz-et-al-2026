#!/bin/bash
#SBATCH --job-name=divergence_compute
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=0-06:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Upstream helper for divergence-vs-dispersal. Emits one consensus
# sequence per BCnS family from the HMM library (via hmmemit -c),
# builds a DIAMOND database, then diamond-blastps each species' proteome
# against the consensus DB and records the median per-hit percent
# identity as a per-species divergence metric.
#
# Output: divergence.tsv (species TAB median_pident TAB n_hits)
#
# Keys species by the full RBH-style key (e.g. Branchiostomafloridae-
# 7739-GCF000003815.2) so it joins directly with per_species_ALG_-
# presence_fusions.tsv.
#
# ~10s per species × 1,143 species ~ 3 h single-thread; minutes with
# SLURM_CPUS_PER_TASK.

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
SCRIPT_NAME=divergence_vs_dispersal
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$SCRIPT_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"
# hmmemit + diamond come from LiSC lmod modules; harmless no-op elsewhere.
module load HMMER DIAMOND 2>/dev/null || true

OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT/hits"

: "${HMM_LIB:?set HMM_LIB to the BCnSSimakov2022.hmm file}"
: "${PROTEIN_DIR:?set PROTEIN_DIR to the annotated_genomes root (GCA_/GCF_ subdirs)}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --- 1. One-time: emit BCnS ALG consensus FASTA ---
CONS_FA="$OUT/bcns_consensus.fa"
if [ ! -s "$CONS_FA" ]; then
    echo "Emitting consensus FASTA from $HMM_LIB ..."
    hmmemit -c "$HMM_LIB" > "$CONS_FA"
fi

# --- 2. One-time: build DIAMOND DB ---
CONS_DB="$OUT/bcns_consensus"
if [ ! -s "${CONS_DB}.dmnd" ]; then
    echo "Building DIAMOND DB ..."
    diamond makedb --in "$CONS_FA" -d "$CONS_DB"
fi

# --- 3. Per species: diamond blastp, median %id ---
DIV_TSV="$OUT/divergence.tsv"
echo -e "species\tmedian_pident\tn_hits" > "$DIV_TSV"

# Map each accession-style RBH filename to its (species_key, pep_path).
# RBH filenames look like:
#   BCnSSimakov2022_<species_key>_xy_reciprocal_best_hits.plotted.rbh
# where <species_key> is the full key used in per_species_ALG_presence_fusions.tsv
# (e.g. Branchiostomafloridae-7739-GCF000003815.2). The accession is the last
# hyphen-separated token.
for RBH in "$RBH_DIR"/BCnSSimakov2022_*.rbh; do
    [ -f "$RBH" ] || continue
    BASE=$(basename "$RBH" _xy_reciprocal_best_hits.plotted.rbh)
    SPECIES_KEY="${BASE#BCnSSimakov2022_}"
    ACC=$(echo "$SPECIES_KEY" | awk -F'-' '{print $NF}')

    # Proteome may be staged as either GCF964019385.1/ or GCF_964019385.1/
    # depending on which chrombase extraction was used.
    ACC_UND=$(echo "$ACC" | sed -E 's/^(GCA|GCF)([0-9])/\1_\2/')
    PEP=""
    for TRY in \
        "$PROTEIN_DIR/$ACC/$ACC.chrFilt.pep.gz" \
        "$PROTEIN_DIR/$ACC_UND/$ACC_UND.chrFilt.pep.gz" \
        "$PROTEIN_DIR/$ACC.chrFilt.pep.gz"; do
        if [ -s "$TRY" ]; then PEP="$TRY"; break; fi
    done
    if [ -z "$PEP" ]; then
        # Unannotated genome (no proteome) — skip; downstream NaN.
        continue
    fi

    HITS="$OUT/hits/${SPECIES_KEY}.tsv"
    if [ ! -s "$HITS" ]; then
        diamond blastp \
            --query "$PEP" \
            --db "$CONS_DB" \
            --threads "$THREADS" \
            -k 1 -e 1e-5 \
            -f 6 qseqid sseqid pident bitscore \
            --quiet > "$HITS"
    fi

    # Median percent identity across all hits.
    read -r MEDIAN NHITS <<<"$(awk -F'\t' '{a[NR]=$3} END{
        n=NR;
        if (n==0) { print "NA 0"; exit }
        asort(a);
        if (n%2==1) printf "%.3f %d", a[int(n/2)+1], n;
        else        printf "%.3f %d", (a[n/2]+a[n/2+1])/2, n;
    }' "$HITS")"

    echo -e "${SPECIES_KEY}\t${MEDIAN}\t${NHITS}" >> "$DIV_TSV"
done

echo "Wrote $DIV_TSV with $(( $(wc -l < "$DIV_TSV") - 1 )) species rows."
