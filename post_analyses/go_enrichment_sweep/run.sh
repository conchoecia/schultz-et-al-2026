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

# Source table: prefer the freshly-rebuilt unique_pairs.tsv.gz from
# post_analyses/defining_features_pairs (run against the 5,821-species
# 202509 dataset with the 28-clade list) when present; fall back to the
# Dryad SupplementaryTable_16.xlsx otherwise.
REBUILT="$REPO_ROOT/post_analyses/defining_features_pairs/out/unique_pairs.tsv.gz"
DRYAD_ROOT="${DRYAD_ROOT:-/lisc/data/scratch/molevo/dts/manifold/submission_dryad/dryad_repo}"
if [ -s "$REBUILT" ]; then
    DEFAULT_SUPP="$REBUILT"
else
    DEFAULT_SUPP="$DRYAD_ROOT/SupplementaryTable_16.xlsx"
fi
SUPP_TABLE="${SUPP_TABLE:-$DEFAULT_SUPP}"
FAMILY_MAP="${FAMILY_MAP:-$REPO_ROOT/post_analyses/build_family_naming_map/out/bcns_family_to_human_gene.tsv}"
GENE2ACCESSION="${GENE2ACCESSION:-$SCRIPT_DIR/ncbi_ref/human_gene2accession.tsv.gz}"
GENE2GO="${GENE2GO:-$SCRIPT_DIR/ncbi_ref/human_gene2go.tsv.gz}"

# Auto-download NCBI files if missing.
if [ ! -s "$GENE2ACCESSION" ] || [ ! -s "$GENE2GO" ]; then
    echo "NCBI reference files not present; running download_ncbi.sh ..."
    bash "$SCRIPT_DIR/download_ncbi.sh"
fi

# All GO-sweep machinery now lives under `egt go <sub>`. The pre-refactor
# Python scripts in this directory (sweep.py, enrich_plots.py,
# plot_volcano.py, plot_pair_distance.py, plot_pair_coenrich.py,
# pair_coenrichment.py, benchmark_*.py) are retained for historical
# reproducibility but are no longer invoked by run.sh. The `go-sweep-v1`
# git tag in the egt repository pins their pre-refactor state.

egt go sweep \
    --supp-table     "$SUPP_TABLE" \
    --family-map     "$FAMILY_MAP" \
    --gene2accession "$GENE2ACCESSION" \
    --gene2go        "$GENE2GO" \
    --out-dir        "$OUT"

# Downstream plot subcommands consume sweep outputs.
OBO="${OBO:-$REPO_ROOT/post_analyses/entanglement_go_enrich/out/go-basic.obo}"

egt go plot-volcano \
    --significant-terms "$OUT/significant_terms.tsv" \
    --out "$OUT/volcanos.pdf"

egt go plot-volcano \
    --significant-terms "$OUT/significant_terms.tsv" \
    --out "$OUT/volcanos_fold3plus.pdf" \
    --min-fold 3

egt go plot \
    --significant-terms "$OUT/significant_terms.tsv" \
    --obo               "$OBO" \
    --out-dir           "$OUT" \
    --term-gene-lists   "$OUT/term_gene_lists.tsv.gz" \
    --gene-symbols      "$OUT/gene_symbols.tsv"

egt go plot-pair-distance \
    --supp-table "$SUPP_TABLE" \
    --summary    "$OUT/summary.tsv" \
    --out        "$OUT/pair_distance.pdf"
