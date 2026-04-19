#!/usr/bin/env python3
"""Locate stored-zero entries in the COO and cross-check against the
per-species RBH files.

For each stored cell with value 0 in the pair-distance COO, the
corresponding (ortholog1, ortholog2) should resolve to two BCnS family
genes that are annotated at the exact same genomic position (same
_scaf, same _pos) in that species' RBH file — because the distance was
computed as `abs(pos_x - pos_y)` (see phylotreeumap.py:782).

If that cross-check passes we know our sparse-native path is reading
the "observed zero" cells correctly; if it fails, the COO itself has a
problem upstream.
"""
import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import load_npz


def load_combo_index(path):
    """Return {col_index: (fam1, fam2)}."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            fams_str, idx = parts
            fams = eval(fams_str)  # file stores Python tuple literals
            mapping[int(idx)] = (fams[0], fams[1])
    return mapping


def find_rbh_for_sample(rbh_dir, sample):
    """Locate the RBH file for a given sample key."""
    candidate = Path(rbh_dir) / f"BCnSSimakov2022_{sample}_xy_reciprocal_best_hits.plotted.rbh"
    if candidate.exists():
        return candidate
    # Fallback: glob for any file that contains the sample key.
    hits = list(Path(rbh_dir).glob(f"BCnSSimakov2022_{sample}*.rbh"))
    return hits[0] if hits else None


def lookup_positions(rbh_path, fam1, fam2, sample):
    """Return ALL rows in the sample's RBH table that hit fam1 or fam2.

    Unlike the previous `.iloc[0]` version, this enumerates paralogs so
    the user can see every occurrence. Also computes the distance for
    every (fam1_row, fam2_row) pair on the same scaffold, and reports
    the minimum distance observed — that's the likely candidate for
    what the COO actually stored.
    """
    df = pd.read_csv(rbh_path, sep="\t")
    scaf_col = f"{sample}_scaf"
    pos_col = f"{sample}_pos"
    gene_col = f"{sample}_gene"
    fam_col = "BCnSSimakov2022_gene"
    need = [fam_col, scaf_col, pos_col]
    missing = [c for c in need if c not in df.columns]
    if missing:
        return {"error": f"missing columns {missing}",
                "cols_present": list(df.columns)[:10]}

    def rows_for(fam):
        r = df[df[fam_col] == fam][[fam_col, gene_col, scaf_col, pos_col]]
        return r.reset_index(drop=True) if gene_col in df.columns else (
            df[df[fam_col] == fam][[fam_col, scaf_col, pos_col]]
              .assign(**{gene_col: "(no gene column)"}).reset_index(drop=True)
        )

    rows1 = rows_for(fam1)
    rows2 = rows_for(fam2)

    # Cross-product distance over same-scaffold paralog pairs.
    combos = []
    for _, r1 in rows1.iterrows():
        for _, r2 in rows2.iterrows():
            same = (r1[scaf_col] == r2[scaf_col])
            dist = abs(int(r1[pos_col]) - int(r2[pos_col])) if same else None
            combos.append({
                fam1: {"gene": r1.get(gene_col, "?"),
                       "scaf": r1[scaf_col],
                       "pos": int(r1[pos_col])},
                fam2: {"gene": r2.get(gene_col, "?"),
                       "scaf": r2[scaf_col],
                       "pos": int(r2[pos_col])},
                "same_scaf": same,
                "distance_bp": dist,
            })

    same_scaf_dists = [c["distance_bp"] for c in combos
                        if c["same_scaf"]]
    return {
        "rows1": rows1,
        "rows2": rows2,
        "combos": combos,
        "n_combos_same_scaf": len(same_scaf_dists),
        "min_distance_bp": min(same_scaf_dists) if same_scaf_dists else None,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coo", required=True)
    ap.add_argument("--sampledf", required=True)
    ap.add_argument("--combo", required=True)
    ap.add_argument("--rbh-dir", required=True)
    ap.add_argument("--max-report", type=int, default=15,
                    help="inspect this many stored-zero entries in detail")
    args = ap.parse_args()

    print(f"[load] sampledf: {args.sampledf}")
    sdf = pd.read_csv(args.sampledf, sep="\t", index_col=0)
    print(f"  rows={len(sdf)}")

    print(f"[load] combo_to_index: {args.combo}")
    combo = load_combo_index(args.combo)
    print(f"  n_pairs={len(combo)}")

    print(f"[load] COO: {args.coo}")
    coo = load_npz(args.coo).tocoo()
    print(f"  shape={coo.shape}  nnz={coo.nnz}")

    # Find all stored-zero cells (biologically-real exact-neighbor pairs).
    zero_mask = coo.data == 0
    n_zero = int(zero_mask.sum())
    print(f"[scan] stored-zero cells: {n_zero} ({n_zero / coo.nnz * 100:.4f}% of nnz)")
    if n_zero == 0:
        print("no stored zeros in the COO — nothing to verify")
        return

    zero_rows = coo.row[zero_mask]
    zero_cols = coo.col[zero_mask]

    print()
    print("=== per-species stored-zero counts (top 10) ===")
    row_counts = pd.Series(zero_rows).value_counts().head(10)
    for row_idx, n in row_counts.items():
        sample = sdf.iloc[row_idx]["sample"]
        print(f"  row={row_idx:5d}  n_zeros={n:6d}  sample={sample}")

    print()
    print(f"=== detailed cross-check on first {args.max_report} stored-zero cells ===")
    for i in range(min(args.max_report, n_zero)):
        row_idx = int(zero_rows[i])
        col_idx = int(zero_cols[i])
        sample = sdf.iloc[row_idx]["sample"]
        if col_idx not in combo:
            print(f"[{i}] row={row_idx} col={col_idx} — pair not in combo index")
            continue
        fam1, fam2 = combo[col_idx]
        rbh = find_rbh_for_sample(args.rbh_dir, sample)
        if rbh is None:
            print(f"\n[{i}] sample={sample}  pair=({fam1}, {fam2})  — no RBH file found")
            continue
        hits = lookup_positions(rbh, fam1, fam2, sample)
        print(f"\n[{i}] sample={sample}")
        print(f"    pair=({fam1}, {fam2})  rbh={rbh.name}")
        if "error" in hits:
            print(f"    ERROR  {hits['error']}")
            continue
        rows1 = hits["rows1"]
        rows2 = hits["rows2"]
        print(f"    --- {fam1} rows in RBH: {len(rows1)} ---")
        if rows1.empty:
            print("      (none)")
        else:
            print(rows1.to_string(index=False))
        print(f"    --- {fam2} rows in RBH: {len(rows2)} ---")
        if rows2.empty:
            print("      (none)")
        else:
            print(rows2.to_string(index=False))
        if hits["min_distance_bp"] is None:
            print(f"    → no same-scaf combo across {len(rows1)}x{len(rows2)} "
                  f"paralog pairs; COO stored 0 but RBH gives none  (MISMATCH)")
        else:
            verdict = "MATCH" if hits["min_distance_bp"] == 0 else "MISMATCH (COO says 0)"
            print(f"    → {hits['n_combos_same_scaf']} same-scaf paralog pair(s), "
                  f"min distance = {hits['min_distance_bp']} bp  →  {verdict}")


if __name__ == "__main__":
    main()
