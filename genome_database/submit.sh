#!/bin/bash
# Build the chromosome-scale genome database from the shipped genome list.
#
# Required env vars (or set in genome_database/config.yaml, then `source` it):
#   CHROMBASE_DIR      — path to a checkout of github.com/conchoecia/chrombase
#   GENOMES_DIR        — where chrombase will stage annotated genomes
#   CHROMBASE_CONFIG   — path to a filled-in copy of chrombase.config.yaml.template
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

: "${CHROMBASE_DIR:?CHROMBASE_DIR must point at a conchoecia/chrombase checkout}"
: "${GENOMES_DIR:?GENOMES_DIR must point at where chrombase will stage genomes}"
: "${CHROMBASE_CONFIG:="$SCRIPT_DIR/chrombase.config.yaml"}"

if [ ! -f "$CHROMBASE_CONFIG" ]; then
    echo "ERROR: $CHROMBASE_CONFIG missing. Copy chrombase.config.yaml.template and edit." >&2
    exit 1
fi

module load Conda/Miniforge3 2>/dev/null || true
source activate "${EGT_CONDA_ENV:-/lisc/data/scratch/molevo/dts/conda_envs/egt}"

# Stage 1: scrape genomes from NCBI driven by chrombase's snakefile +
#          the species list in genome_list.yaml.
snakemake \
    --snakefile "$CHROMBASE_DIR/chrombase_scrape_genomes_NCBI.snakefile" \
    --configfile "$CHROMBASE_CONFIG" \
    --config genome_list="$SCRIPT_DIR/genome_list.yaml" GENOMES_DIR="$GENOMES_DIR" \
    --jobs 100 \
    --rerun-incomplete

# Stage 2: build the annotated-chromosome database from the scraped set.
snakemake \
    --snakefile "$CHROMBASE_DIR/chrombase_build_db_annotated_chr.snakefile" \
    --configfile "$CHROMBASE_CONFIG" \
    --config genome_list="$SCRIPT_DIR/genome_list.yaml" GENOMES_DIR="$GENOMES_DIR" \
    --jobs 100 \
    --rerun-incomplete

echo "Database build complete. Output under $GENOMES_DIR."
