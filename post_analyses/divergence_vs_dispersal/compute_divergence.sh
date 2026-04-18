#!/bin/bash
#SBATCH --job-name=divergence_compute
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Upstream helper for divergence-vs-dispersal. Runs hmmscan of each
# species proteome against the BCnSSimakov2022 HMM library and extracts
# the median per-family bitscore as a proxy for protein divergence.
#
# Output: divergence.tsv (species TAB median_bitscore)

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
# hmmscan is provided by the HMMER lmod module on LiSC; harmless no-op elsewhere.
module load HMMER 2>/dev/null || true

OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT/tbl"

: "${HMM_LIB:?set HMM_LIB to the BCnSSimakov2022.hmm file}"
: "${PROTEIN_DIR:?set PROTEIN_DIR to a directory with {sample}.chrFilt.pep.gz per species}"

DIV_TSV="$OUT/divergence.tsv"
echo -e "species\tmedian_bitscore" > "$DIV_TSV"

for PEP in "$PROTEIN_DIR"/*.chrFilt.pep.gz; do
    SAMPLE=$(basename "$PEP" .chrFilt.pep.gz)
    TBL="$OUT/tbl/${SAMPLE}.tbl"
    if [ ! -s "$TBL" ]; then
        zcat "$PEP" | hmmscan \
            --cpu "${SLURM_CPUS_PER_TASK:-8}" \
            --tblout "$TBL" \
            -E 1e-5 "$HMM_LIB" - > /dev/null
    fi
    MEDIAN=$(awk '!/^#/ && NF>0 {print $6}' "$TBL" | sort -n | \
             awk '{a[NR]=$1} END{ if (NR==0) print "NA"; else if (NR%2==1) print a[int(NR/2)+1]; else printf "%.3f\n", (a[NR/2]+a[NR/2+1])/2 }')
    echo -e "${SAMPLE}\t${MEDIAN}" >> "$DIV_TSV"
done

echo "Wrote $DIV_TSV with $(wc -l < "$DIV_TSV") rows."
