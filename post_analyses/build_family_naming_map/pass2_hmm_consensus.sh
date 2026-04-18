#!/bin/bash
#SBATCH --job-name=family_naming_pass2
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Optional pass-2 enrichment for the BCnS family → human gene map:
# for every family without a direct human hit from pass 1, emit the
# HMM consensus and DIAMOND-search it against Swiss-Prot human.
#
# Requires:
#   - HMMER (hmmemit, hmmstat)
#   - diamond
#   - Swiss-Prot human FASTA at $UNIPROT_HUMAN_FASTA
#   - path to the BCnSSimakov2022 HMM library (single concatenated file)

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
SCRIPT_NAME=build_family_naming_map
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$SCRIPT_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT"

: "${HMM_LIB:?set HMM_LIB to the BCnSSimakov2022 HMM file, e.g. /lisc/.../odp/BCnSSimakov2022/BCnSSimakov2022.hmm}"
: "${UNIPROT_HUMAN_FASTA:?set UNIPROT_HUMAN_FASTA to a FASTA of Swiss-Prot human proteins}"

PASS1_TSV="$OUT/bcns_family_to_human_gene.tsv"
if [ ! -f "$PASS1_TSV" ]; then
    echo "ERROR: run pass 1 (run.sh) first to produce $PASS1_TSV" >&2
    exit 1
fi

# 1. Build the DIAMOND DB once.
if [ ! -f "$OUT/uniprot_human.dmnd" ]; then
    diamond makedb --in "$UNIPROT_HUMAN_FASTA" -d "$OUT/uniprot_human"
fi

# 2. For each family with no human_gene in pass 1, emit consensus + DIAMOND.
mkdir -p "$OUT/consensus"
PASS2_TSV="$OUT/pass2_hmm_to_human.tsv"
echo -e "family_id\thuman_gene" > "$PASS2_TSV"

awk -F'\t' 'NR>1 && ($3=="" || $3=="NA") {print $1}' "$PASS1_TSV" | \
while read FAMILY; do
    CONS="$OUT/consensus/${FAMILY}.fa"
    # hmmemit -c emits a consensus sequence for the named HMM within the library
    hmmemit -c -o "$CONS" -n "$FAMILY" "$HMM_LIB" || {
        echo "  no HMM for $FAMILY, skipping" >&2
        continue
    }
    HIT=$(diamond blastp \
            --query "$CONS" \
            --db "$OUT/uniprot_human" \
            -k 1 -f 6 sseqid pident evalue \
            --quiet 2>/dev/null | head -1 | awk '{print $1}')
    if [ -n "$HIT" ]; then
        echo -e "${FAMILY}\t${HIT}" >> "$PASS2_TSV"
    fi
done

echo "Pass 2 mapped $(( $(wc -l < "$PASS2_TSV") - 1 )) additional families."
echo "Now re-run pass 1 with --hmm-map $PASS2_TSV to produce the final map."
