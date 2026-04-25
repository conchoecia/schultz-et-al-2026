#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${HERE}/out"
mkdir -p "${OUTDIR}"

EGT_REPO="${EGT_REPO:-/lisc/data/scratch/molevo/dts/git/egt}"
PYTHON="${PYTHON:-${EGT_REPO}/.venv312/bin/python}"

INPUT_DF="${INPUT_DF:-${HERE}/../figure_umap_inputs/subsample_allsamples.neighbors_250.mind_1.0.missing_large.paper_palette.df.gz}"
TREE_NEWICK="${TREE_NEWICK:-${HERE}/../figure_umap_inputs/20251130Tree.calibrated_tree.nwk}"
PALETTE="${PALETTE:-${EGT_REPO}/src/egt/data/paper_palette_simple.yaml}"
PREFIX="${OUTDIR}/allsamples_250_1.0_paperpalette_simple_featuresplot"

test -s "${TREE_NEWICK}"

cd "${EGT_REPO}"
MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-cache}" "${PYTHON}" -m egt.phylotreeumap_plotdfs \
  -f "${INPUT_DF}" \
  -p "${PREFIX}" \
  --plot-phyla \
  --phyla-clean-output \
  --color-source palette \
  --palette "${PALETTE}" \
  --phyla-manuscript-layout \
  --phyla-main-panel-output "${PREFIX}.phyla_main_panel.pdf" \
  --main-panel-side-buffer 0.10
