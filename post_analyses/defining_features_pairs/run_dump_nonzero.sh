#!/bin/bash
#SBATCH --job-name=dump_nonzero_pairs_annelida
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=0-00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
set -euo pipefail

# Sanity-check NON-zero stored cells in the COO for Annelida: pick a
# random 200 cells with stored value > 0, look up both orthologs in
# the species' RBH, and compare the COO's stored distance vs the
# distance computed from chromosome positions. If these match across
# the sample, the COO is structurally correct and the only weird thing
# about it is the 0.0018% of cells stored as zero (upstream
# placeholders) — our sparse-native fix already handles those.

STEP_NAME=defining_features_pairs
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
while [ "$REPO_ROOT" != "/" ] && [ ! -f "$REPO_ROOT/config.template.yaml" ]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
SCRIPT_DIR="$REPO_ROOT/post_analyses/$STEP_NAME"
source "$REPO_ROOT/config.yaml"

module load Conda/Miniforge3 2>/dev/null || true
source activate "$EGT_CONDA_ENV"

COO_ROOT="${COO_ROOT:-/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP}"

OUT_TSV="$SCRIPT_DIR/nonzero_pairs_chrom_entries_annelida.tsv"
cd "$SCRIPT_DIR"

python dump_zero_pairs.py \
    --coo              "$COO_ROOT/allsamples.coo.npz" \
    --sampledf         "$COO_ROOT/sampledf.tsv" \
    --combo            "$COO_ROOT/combo_to_index.txt" \
    --rbh-dir          "$RBH_DIR" \
    --out-tsv          "$OUT_TSV" \
    --clade-taxid      6340 \
    --nonzero-sample   200

ls -la "$OUT_TSV"
echo
echo "=== quick match/mismatch summary ==="
python - <<PY
import pandas as pd
df = pd.read_csv("$OUT_TSV", sep="\t")
print(f"rows: {len(df)}")
if "matches_coo" in df.columns:
    matches = df["matches_coo"].dropna()
    print(f"checkable (same scaf, not missing): {len(matches)}")
    print(f"  matches: {int(matches.sum())}")
    print(f"  mismatches: {int((~matches).sum())}")
    bad = df[df["matches_coo"] == False].head(5)
    if len(bad):
        print()
        print("first 5 mismatches (abs(fam1_pos - fam2_pos) != coo_value):")
        print(bad[["sample","fam1","fam2","coo_value","fam1_pos","fam2_pos","computed_distance_bp"]].to_string(index=False))
if "same_scaf" in df.columns:
    print()
    same = df["same_scaf"].dropna()
    print(f"same_scaf: True={int(same.sum())}  False={int((~same).sum())}  "
          f"missing/other={len(df)-len(same)}")
PY
