#!/bin/bash
#SBATCH --job-name=go_enrichment_sweep
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Under sbatch, BASH_SOURCE points at /var/spool/slurm/..., so we resolve
# REPO_ROOT by walking up from SLURM_SUBMIT_DIR (or this dir for `bash run.sh`).
STEP_NAME=go_enrichment_sweep
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"
OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT"

# Dryad SupplementaryTable_16 ships inside the workflow bundle; on LiSC
# it's pre-staged at the author's submission directory. Override via env.
DRYAD_ROOT="${DRYAD_ROOT:-/lisc/data/scratch/molevo/dts/manifold/submission_dryad/dryad_repo}"
SUPP_TABLE="${SUPP_TABLE:-$DRYAD_ROOT/SupplementaryTable_16.xlsx}"
FAMILY_MAP="${FAMILY_MAP:-$REPO_ROOT/post_analyses/build_family_naming_map/out/bcns_family_to_human_gene.tsv}"
GENE2ACCESSION="${GENE2ACCESSION:-$SCRIPT_DIR/ncbi_ref/human_gene2accession.tsv.gz}"
GENE2GO="${GENE2GO:-$SCRIPT_DIR/ncbi_ref/human_gene2go.tsv.gz}"

# Auto-download NCBI files if missing.
if [ ! -s "$GENE2ACCESSION" ] || [ ! -s "$GENE2GO" ]; then
    echo "NCBI reference files not present; running download_ncbi.sh ..."
    bash "$SCRIPT_DIR/download_ncbi.sh"
fi

python "$SCRIPT_DIR/sweep.py" \
    --supp-table     "$SUPP_TABLE" \
    --family-map     "$FAMILY_MAP" \
    --gene2accession "$GENE2ACCESSION" \
    --gene2go        "$GENE2GO" \
    --out-dir        "$OUT"
