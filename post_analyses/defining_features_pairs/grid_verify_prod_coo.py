#!/usr/bin/env python3
"""Run egt._testing.grid_verify_coo against the production 202509 COO.

Quantifies the row-scramble bug's impact on the current production COO
by sampling an N x N grid of cells and comparing each against the
per-species gb.gz. Writes a full per-cell detail TSV plus a summary
printout.
"""
import json
from pathlib import Path

import pandas as pd
from scipy.sparse import load_npz

from egt._testing import grid_verify_coo

COO_ROOT   = Path("/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP")
COO_FILE   = COO_ROOT / "allsamples.coo.npz"
SAMPLEDF   = COO_ROOT / "sampledf.tsv"
COMBO_FILE = COO_ROOT / "combo_to_index.txt"
GB_DIR     = COO_ROOT / "distance_matrices"

N_ROWS = 20
N_COLS = 20

OUT_DIR = Path(__file__).resolve().parent / "grid_verify_prod_out"
OUT_DIR.mkdir(exist_ok=True)

print(f"[load] sampledf from {SAMPLEDF}")
sdf = pd.read_csv(SAMPLEDF, sep="\t", index_col=0)
# The positional order of sdf must match the COO's row order. The COO
# build writes rows by iloc position after reset_index, so we normalize
# here too.
sdf = sdf.reset_index(drop=True)
print(f"  n_species = {len(sdf)}")

print(f"[load] combo_to_index from {COMBO_FILE}")
pair_to_ix = {}
with open(COMBO_FILE) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        fams = eval(parts[0])
        pair_to_ix[(fams[0], fams[1])] = int(parts[1])
print(f"  n_pairs = {len(pair_to_ix)}")

print(f"[load] COO from {COO_FILE} (5.9 GB, takes a while)...")
coo = load_npz(COO_FILE)
print(f"  shape = {coo.shape}  nnz = {coo.nnz}")

print(f"[grid_verify] sampling {N_ROWS} x {N_COLS} = {N_ROWS*N_COLS} cells, strict=False")
summary = grid_verify_coo(
    coo, sdf, pair_to_ix, GB_DIR,
    n_rows=N_ROWS, n_cols=N_COLS,
    sample_col="sample",
    strict=False,
)

# Summary
print("\n=== SUMMARY ===")
print(f"  n_cells : {summary['n_cells']}")
print(f"  passed  : {summary['passed']}")
print(f"  failed  : {summary['failed']}")
print(f"  skipped : {summary['skipped']}")
if summary["n_cells"]:
    pct_failed = 100.0 * summary["failed"] / summary["n_cells"]
    print(f"  pct_failed : {pct_failed:.1f}%")

# Write full details
details_path = OUT_DIR / "grid_verify_details.tsv"
pd.DataFrame(summary["details"]).to_csv(details_path, sep="\t", index=False)
print(f"\n[write] details -> {details_path}")

# Write summary JSON without the per-cell details
summary_path = OUT_DIR / "grid_verify_summary.json"
slim = {k: v for k, v in summary.items() if k != "details"}
with open(summary_path, "w") as fh:
    json.dump(slim, fh, indent=2)
print(f"[write] summary -> {summary_path}")

# Print a few failing examples for immediate inspection
fails = [d for d in summary["details"] if d.get("outcome") is False]
if fails:
    print(f"\n=== FIRST 10 FAILURES (of {len(fails)}) ===")
    for d in fails[:10]:
        print(f"  row={d['row_idx']:>5} col={d['col_idx']:>8} "
              f"sample={d['sample']} pair={d.get('pair')} "
              f"coo={d['coo_value']} gbgz={d.get('gbgz_value')} "
              f"status={d['status']}")
else:
    print("\nNo failures observed in grid sample.")
