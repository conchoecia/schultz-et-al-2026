#!/bin/bash
#SBATCH --job-name=palette_preview
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

OUT="$SCRIPT_DIR/out"
mkdir -p "$OUT"

# Uses the calibrated tree from step 2. Emits PDF + 3 Newick/NEXUS
# variants so you can pick whichever lays out best in FigTree/iTOL.
TREE="../../pipeline/step2_divergence_times/${TREE_PREFIX}.calibrated_tree.nwk"
cd "$SCRIPT_DIR"

egt palette-preview \
    --tree                   "$TREE" \
    --out                    "$OUT/palette_preview.pdf" \
    --align-tips \
    --emit-colored-newick    "$OUT/palette.colored.nwk" \
    --emit-figtree-nexus     "$OUT/palette.figtree.nex" \
    --emit-collapsed-nexus   "$OUT/palette.collapsed.nex" \
    --collapse-dominance     0.9
