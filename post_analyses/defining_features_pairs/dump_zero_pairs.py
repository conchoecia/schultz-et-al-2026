#!/usr/bin/env python3
"""Dump a TSV with one row per stored-zero cell in the COO.

Columns:
  sample
  fam1, fam2                          — BCnS family IDs for the pair
  fam1_gene, fam1_scaf, fam1_pos      — from that species' RBH file
  fam2_gene, fam2_scaf, fam2_pos
  n_fam1_paralogs, n_fam2_paralogs    — how many rows for each family in the RBH
  same_scaf                           — True if fam1_scaf == fam2_scaf
  computed_distance_bp                — abs(fam1_pos - fam2_pos) if same_scaf, else empty
  note                                — "RBH missing fam1" / "RBH missing fam2" / "ok"

For families with multiple paralogs in the RBH, we pick the (fam1_row,
fam2_row) combination that minimizes computed distance (and prefer
same-scaffold combinations over different-scaffold ones). That's the
"best-case" interpretation — what the COO would have stored if it
chose the closest pair.
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import load_npz


def load_combo_index(path):
    mapping = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2:
                continue
            fams_str, idx = parts
            fams = eval(fams_str)
            mapping[int(idx)] = (fams[0], fams[1])
    return mapping


def find_rbh_for_sample(rbh_dir, sample):
    c = Path(rbh_dir) / f"BCnSSimakov2022_{sample}_xy_reciprocal_best_hits.plotted.rbh"
    if c.exists():
        return c
    hits = list(Path(rbh_dir).glob(f"BCnSSimakov2022_{sample}*.rbh"))
    return hits[0] if hits else None


def best_pair(rows1, rows2, scaf_col, pos_col, gene_col):
    """Choose the (r1, r2) combination with smallest same-scaf distance.
    Falls back to (first, first) if no same-scaf pair exists."""
    if rows1.empty or rows2.empty:
        return None
    best = None
    best_dist = None
    for _, r1 in rows1.iterrows():
        for _, r2 in rows2.iterrows():
            same = (r1[scaf_col] == r2[scaf_col])
            dist = abs(int(r1[pos_col]) - int(r2[pos_col])) if same else None
            cand = (same, dist, r1, r2)
            if best is None:
                best = cand
                best_dist = dist
                continue
            # Prefer same-scaf + smaller distance.
            cur_same, cur_dist, _, _ = best
            if (same and not cur_same) or (
                same and cur_same and dist < cur_dist):
                best = cand
                best_dist = dist
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coo", required=True)
    ap.add_argument("--sampledf", required=True)
    ap.add_argument("--combo", required=True)
    ap.add_argument("--rbh-dir", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--clade-taxid", type=int, default=None,
                    help="if set, restrict to species with this taxid in "
                         "taxid_list (e.g. 6340 for Annelida)")
    ap.add_argument("--nonzero-sample", type=int, default=0,
                    help="if >0, instead of dumping stored-zero cells, "
                         "sample this many NON-zero stored cells at random "
                         "(within the clade filter) and compare the COO's "
                         "stored value against the RBH-computed distance")
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

    if args.nonzero_sample > 0:
        mask = coo.data > 0
        print(f"[mode] sampling {args.nonzero_sample} NON-zero cells "
              f"for COO-vs-RBH distance check")
    else:
        mask = coo.data == 0
    if args.clade_taxid is not None:
        taxid_sets = sdf["taxid_list"].apply(eval).apply(set).to_numpy()
        clade_rows = np.array([args.clade_taxid in s for s in taxid_sets], dtype=bool)
        clade_row_indices = np.where(clade_rows)[0]
        print(f"[filter] clade_taxid={args.clade_taxid}: "
              f"{len(clade_row_indices)} species in clade")
        row_mask = np.isin(coo.row, clade_row_indices)
        mask = mask & row_mask
    n_zero = int(mask.sum())
    if args.nonzero_sample > 0 and n_zero > args.nonzero_sample:
        # Random subsample of non-zero cells within the filter.
        idx_all = np.where(mask)[0]
        rng = np.random.default_rng(42)
        pick = rng.choice(idx_all, size=args.nonzero_sample, replace=False)
        mask = np.zeros_like(mask)
        mask[pick] = True
        n_zero = args.nonzero_sample
        print(f"  subsampled to {n_zero} cells")
    print(f"[scan] stored-zero cells: {n_zero}")
    if n_zero == 0:
        print("no stored zeros — nothing to dump")
        return

    zero_rows = coo.row[mask]
    zero_cols = coo.col[mask]
    stored_values = coo.data[mask].astype(float)

    # Group cells by species so we open each RBH file at most once.
    # Carry the stored value alongside the column index for each cell.
    by_row = {}
    for i in range(n_zero):
        r = int(zero_rows[i]); c = int(zero_cols[i])
        v = float(stored_values[i])
        by_row.setdefault(r, []).append((c, v))

    out_rows = []
    total_species = len(by_row)
    for species_counter, (row_idx, cell_list) in enumerate(by_row.items(), 1):
        sample = sdf.iloc[row_idx]["sample"]
        rbh_path = find_rbh_for_sample(args.rbh_dir, sample)
        if rbh_path is None:
            for c, v in cell_list:
                fam1, fam2 = combo.get(c, (None, None))
                out_rows.append({
                    "sample": sample, "fam1": fam1, "fam2": fam2,
                    "coo_value": v,
                    "note": "RBH file not found",
                })
            continue

        df = pd.read_csv(rbh_path, sep="\t")
        scaf_col = f"{sample}_scaf"
        pos_col  = f"{sample}_pos"
        gene_col = f"{sample}_gene"
        fam_col  = "BCnSSimakov2022_gene"
        need = [fam_col, scaf_col, pos_col]
        missing = [c for c in need if c not in df.columns]
        if missing:
            for c, v in cell_list:
                fam1, fam2 = combo.get(c, (None, None))
                out_rows.append({
                    "sample": sample, "fam1": fam1, "fam2": fam2,
                    "coo_value": v,
                    "note": f"missing columns {missing}",
                })
            continue
        has_gene = gene_col in df.columns

        # Pre-index by family for speed (each family typically appears 1x).
        df_by_fam = {k: g.reset_index(drop=True)
                     for k, g in df.groupby(fam_col, sort=False)}

        for c, v in cell_list:
            fam1, fam2 = combo.get(c, (None, None))
            if fam1 is None:
                out_rows.append({"sample": sample, "coo_value": v,
                                  "note": "col not in combo"})
                continue
            rows1 = df_by_fam.get(fam1, df.iloc[0:0])
            rows2 = df_by_fam.get(fam2, df.iloc[0:0])
            if rows1.empty and rows2.empty:
                note = "both families missing from RBH"
            elif rows1.empty:
                note = "fam1 missing from RBH"
            elif rows2.empty:
                note = "fam2 missing from RBH"
            else:
                note = "ok"

            best = best_pair(rows1, rows2, scaf_col, pos_col,
                              gene_col if has_gene else None)
            rec = {
                "sample": sample,
                "fam1": fam1, "fam2": fam2,
                "coo_value": v,
                "n_fam1_paralogs": len(rows1),
                "n_fam2_paralogs": len(rows2),
                "note": note,
            }
            if best is not None:
                same, dist, r1, r2 = best
                rec["fam1_gene"] = r1[gene_col] if has_gene else ""
                rec["fam1_scaf"] = r1[scaf_col]
                rec["fam1_pos"]  = int(r1[pos_col])
                rec["fam2_gene"] = r2[gene_col] if has_gene else ""
                rec["fam2_scaf"] = r2[scaf_col]
                rec["fam2_pos"]  = int(r2[pos_col])
                rec["same_scaf"] = bool(same)
                rec["computed_distance_bp"] = dist if same else ""
                if same and dist is not None:
                    rec["matches_coo"] = (abs(dist - v) < 0.5)
            out_rows.append(rec)

        if species_counter % 100 == 0 or species_counter == total_species:
            print(f"  processed {species_counter}/{total_species} species "
                  f"({len(out_rows)} rows so far)")

    out_df = pd.DataFrame(out_rows)
    col_order = ["sample", "fam1", "fam2",
                 "coo_value",
                 "fam1_gene", "fam1_scaf", "fam1_pos",
                 "fam2_gene", "fam2_scaf", "fam2_pos",
                 "n_fam1_paralogs", "n_fam2_paralogs",
                 "same_scaf", "computed_distance_bp",
                 "matches_coo", "note"]
    col_order = [c for c in col_order if c in out_df.columns]
    out_df = out_df[col_order]
    out_df.to_csv(args.out_tsv, sep="\t", index=False)

    print()
    print(f"[write] {args.out_tsv}")
    print(f"  total rows: {len(out_df)}")
    if "same_scaf" in out_df.columns:
        print(f"  same_scaf True : {int((out_df['same_scaf'] == True).sum())}")
        print(f"  same_scaf False: {int((out_df['same_scaf'] == False).sum())}")
        if "computed_distance_bp" in out_df.columns:
            d = pd.to_numeric(out_df["computed_distance_bp"], errors="coerce")
            zero = int((d == 0).sum())
            print(f"  computed distance == 0 (TRUE zero-distance cells): {zero}")
            print(f"  computed distance > 0 (COO says 0 but they're apart): {int(((d > 0)).sum())}")


if __name__ == "__main__":
    main()
