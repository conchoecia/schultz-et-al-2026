#!/bin/bash
#SBATCH --job-name=tree_step1
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

# Build NCBI topology tree from the species list in the genome config.
# --custom_phylogeny: use the Ctenophora-sister arrangement (Myriazoa=-67)
#                    instead of NCBI's default topology.
egt taxids-to-newick -c "$GENOME_CONFIG_YAML" --custom_phylogeny -o ncbi_tree.nwk
