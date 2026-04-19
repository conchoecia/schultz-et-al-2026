#!/usr/bin/env python3
"""Verify the freshly-rebuilt 202509 COO against the per-species gb.gz.

Runs `egt._testing.grid_verify_coo` on a 50x50 grid (2500 cells).
Exits 0 iff 0 cells failed; non-zero otherwise.

Verify-only; no renames. See `promote_rebuild.sh` for the (separate,
user-run) step that promotes rebuild->canonical on pass.
"""
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.sparse import load_npz

from egt._testing import grid_verify_coo

COO_ROOT   = Path("/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP")
REBUILD    = COO_ROOT / "allsamples.coo.REBUILD_20260419.npz"
SAMPLEDF   = COO_ROOT / "sampledf.tsv"
COMBO_FILE = COO_ROOT / "combo_to_index.txt"
GB_DIR     = COO_ROOT / "distance_matrices"
OUT_DIR    = Path(__file__).resolve().parent / "grid_verify_new_out"

N_ROWS = 50
N_COLS = 50

OUT_DIR.mkdir(exist_ok=True)

print(f"[verify] input: {REBUILD}")
if not REBUILD.exists():
    sys.exit(f"FATAL: rebuild output missing at {REBUILD}")

print(f"[verify] loading sampledf ({SAMPLEDF}) ...")
sdf = pd.read_csv(SAMPLEDF, sep="\t", index_col=0).reset_index(drop=True)
print(f"  n_species = {len(sdf)}")

print(f"[verify] loading combo_to_index ({COMBO_FILE}) ...")
pair_to_ix = {}
with open(COMBO_FILE) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        fams = eval(parts[0])
        pair_to_ix[(fams[0], fams[1])] = int(parts[1])
print(f"  n_pairs = {len(pair_to_ix)}")

print(f"[verify] loading COO ({REBUILD}) ...")
coo = load_npz(REBUILD)
print(f"  shape = {coo.shape}  nnz = {coo.nnz}")

print(f"[verify] grid sampling {N_ROWS}x{N_COLS} = {N_ROWS*N_COLS} cells, strict=False")
summary = grid_verify_coo(
    coo, sdf, pair_to_ix, GB_DIR,
    n_rows=N_ROWS, n_cols=N_COLS,
    sample_col="sample",
    strict=False,
)

print("\n=== SUMMARY ===")
print(f"  n_cells : {summary['n_cells']}")
print(f"  passed  : {summary['passed']}")
print(f"  failed  : {summary['failed']}")
print(f"  skipped : {summary['skipped']}")

details_path = OUT_DIR / "grid_verify_details.tsv"
pd.DataFrame(summary["details"]).to_csv(details_path, sep="\t", index=False)
summary_path = OUT_DIR / "grid_verify_summary.json"
slim = {k: v for k, v in summary.items() if k != "details"}
with open(summary_path, "w") as fh:
    json.dump(slim, fh, indent=2)
print(f"[write] details -> {details_path}")
print(f"[write] summary -> {summary_path}")

if summary["failed"] > 0:
    print(f"\n[verify] FAIL: {summary['failed']}/{summary['n_cells']} cells mismatched.")
    fails = [d for d in summary["details"] if d.get("outcome") is False]
    for d in fails[:10]:
        print(f"  row={d['row_idx']:>5} col={d['col_idx']:>8} sample={d['sample']} "
              f"coo={d['coo_value']} gbgz={d.get('gbgz_value')} status={d['status']}")
    print("\n[verify] NOT promoting rebuild; canonical allsamples.coo.npz left untouched.")
    sys.exit(1)

print(f"\n[verify] PASS — 0/{summary['n_cells']} cells mismatched.")
print(f"         Rebuild file remains at {REBUILD} (no rename).")
print(f"         Run promote_rebuild.sh separately when you're ready to make it canonical.")
