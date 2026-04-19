#!/usr/bin/env python3
"""Direct trace of a single COO cell.

For a hardcoded (sample, fam1, fam2), resolve the row_idx and col_idx,
print the COO value, and compare with what the gb.gz for that sample
stores. Then scan ALL gb.gz files in the distance_matrices directory
and report which ones store the COO's value for this pair — if any —
to confirm or rule out a row-scramble bug.
"""
import gzip
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import load_npz

COO_ROOT   = Path("/lisc/data/scratch/molevo/dts/manifold/UMAP_snakemake_202509/GTUMAP")
COO_FILE   = COO_ROOT / "allsamples.coo.npz"
SAMPLEDF   = COO_ROOT / "sampledf.tsv"
COMBO_FILE = COO_ROOT / "combo_to_index.txt"
GB_DIR     = COO_ROOT / "distance_matrices"

TARGET_SAMPLE = "Hirudonipponia-42736-GCA040113095.1"
FAM1 = "Simakov2022BCnS_genefamily_4381"
FAM2 = "Simakov2022BCnS_genefamily_7471"

# --- Resolve indices -------------------------------------------------
sdf = pd.read_csv(SAMPLEDF, sep="\t", index_col=0)
row_idx_candidates = sdf.index[sdf["sample"] == TARGET_SAMPLE].tolist()
assert len(row_idx_candidates) == 1, f"found {len(row_idx_candidates)} rows for {TARGET_SAMPLE}"
row_idx_label = row_idx_candidates[0]
row_idx_pos = sdf.index.get_loc(row_idx_label)
print(f"[sampledf] {TARGET_SAMPLE} -> index_label={row_idx_label}  positional={row_idx_pos}")

# Pair index — key is a Python tuple as written in combo_to_index.txt
pair_to_ix = {}
with open(COMBO_FILE) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 2:
            continue
        fams = eval(parts[0])
        pair_to_ix[(fams[0], fams[1])] = int(parts[1])

keys = [(FAM1, FAM2), (FAM2, FAM1)]
col_idx = None
for k in keys:
    if k in pair_to_ix:
        col_idx = pair_to_ix[k]
        print(f"[combo] pair {k} -> col_idx={col_idx}")
        break
if col_idx is None:
    raise SystemExit("pair not in combo_to_index")

# --- Query COO at (row_idx_pos, col_idx) -----------------------------
print("[load] COO (this takes a bit) ...")
coo = load_npz(COO_FILE).tocsr()
v_coo = coo[row_idx_pos, col_idx]
print(f"[coo] value at ({row_idx_pos}, {col_idx}) = {v_coo}")

# --- Compare with gb.gz for this sample ------------------------------
gbgz = GB_DIR / f"{TARGET_SAMPLE}.gb.gz"
if not gbgz.exists():
    raise SystemExit(f"no gb.gz for {TARGET_SAMPLE}")
df_gb = pd.read_csv(gbgz, sep="\t")
hit = df_gb[(df_gb["rbh1"].isin([FAM1, FAM2])) & (df_gb["rbh2"].isin([FAM1, FAM2]))]
print(f"[gb.gz] distance for ({FAM1}, {FAM2}) in {TARGET_SAMPLE}: "
      f"{hit['distance'].tolist()}")

# --- Scan all gb.gz files: does anyone store COO's value for this pair? -
print(f"\n[scan] which species' gb.gz has distance == {v_coo} for pair ({FAM1},{FAM2})?")
target_val = float(v_coo)
tolerance = 0.5
paths = sorted(GB_DIR.glob("*.gb.gz"))
print(f"  scanning {len(paths)} gb.gz files ...")
matches = []
for p in paths:
    try:
        df = pd.read_csv(p, sep="\t")
    except Exception:
        continue
    rows = df[(df["rbh1"].isin([FAM1, FAM2])) & (df["rbh2"].isin([FAM1, FAM2]))]
    for _, r in rows.iterrows():
        if abs(float(r["distance"]) - target_val) < tolerance:
            matches.append((p.name, float(r["distance"])))
            break
    if len(matches) >= 5:
        break
print(f"  matches (first {len(matches)}):")
for n, d in matches:
    print(f"    {n}  distance={d}")
