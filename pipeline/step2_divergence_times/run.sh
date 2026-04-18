#!/bin/bash
#SBATCH --job-name=tree_step2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

cd "$SCRIPT_DIR"

# Calibrate the NCBI topology against TimeTree divergence times and emit
# per-node and per-edge TSVs consumed by steps 4, 6, 7, 8 and 9.
#
# TIME_NEWICK defaults to the copy shipped in the Dryad workflow tarball
# (bin/download_data.sh places it at
# dryad_data/newick_and_timetree_20251118/step2_download_newick_from_timetree/newick_timetree.nwk).
# Override in config.yaml if running with a different TimeTree export.
egt newick-to-common-ancestors \
    --topology_newick ../step1_generate_newick/ncbi_tree.nwk \
    --time_newick "${REPO_ROOT}/${TIME_NEWICK}" \
    -c "$GENOME_CONFIG_YAML" \
    -C "$CHROM_COUNTS" \
    -p "$TREE_PREFIX"
