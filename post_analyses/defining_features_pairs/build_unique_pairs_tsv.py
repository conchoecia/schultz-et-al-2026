#!/usr/bin/env python3
"""Produce unique_pairs.tsv from 28 per-clade TSVs — just the filter step.

Applies the same flag-derivation + union-filter + ortholog-name-join
that `src/egt/legacy/defining_features_plot2.py::main` does, but
without the 5-figure QC PDFs (those live in `summarize_qc.py`). Reuses
the helpers `add_ratio_columns`, `compute_z_scores`, `assign_flags`
from the legacy module so the math is guaranteed-identical to
SupplementaryTable_16's filter.

Default sigma cutoff = 2 (matches the publication default).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from egt.legacy.defining_features_plot2 import (
    add_ratio_columns,
    compute_z_scores,
    assign_flags,
)
from egt.phylotreeumap import algcomboix_file_to_dict
from egt import rbh_tools

OCC_MIN = 0.5
FLAG_COLS = [
    "close_in_clade", "distant_in_clade",
    "stable_in_clade", "unstable_in_clade",
    "unique_to_clade",
]


def parse_clade_from_filename(path: Path) -> tuple[str, int]:
    m = re.match(r"^(?P<name>.+?)_(?P<taxid>\d+)_unique_pair_df\.tsv\.gz$", path.name)
    if not m:
        raise ValueError(f"unexpected filename: {path}")
    return m.group("name"), int(m.group("taxid"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--clade-stats-dir", type=Path, required=True,
                    help="Directory of *_unique_pair_df.tsv.gz")
    ap.add_argument("--rbh-file", type=Path, required=True,
                    help="ALG database RBH file (BCnSSimakov2022.rbh or similar).")
    ap.add_argument("--pair-combination-path", type=Path, required=True,
                    help="combo_to_index.txt path.")
    ap.add_argument("--output", type=Path, required=True,
                    help="Output unique_pairs.tsv(.gz) path.")
    ap.add_argument("--sigma", type=float, default=2.0,
                    help="Sigma cutoff for non-unique flags (default 2, matches SuppTable_16).")
    args = ap.parse_args(argv)

    print(f"[build] sigma = {args.sigma}")
    print(f"[build] occupancy_in gate for non-unique flags >= {OCC_MIN}")

    # Pair col_idx -> (ortholog1, ortholog2)
    combo_to_ix = algcomboix_file_to_dict(str(args.pair_combination_path))
    ix_to_pair = {v: k for k, v in combo_to_ix.items()}

    # RBH: ortholog name -> ALG family
    print(f"[build] loading RBH: {args.rbh_file}")
    rbh_df = rbh_tools.parse_rbh(str(args.rbh_file))
    rbh_to_alg = dict(zip(rbh_df["rbh"], rbh_df["gene_group"]))

    tsvs = sorted(args.clade_stats_dir.glob("*_unique_pair_df.tsv.gz"))
    print(f"[build] found {len(tsvs)} per-clade TSVs")

    kept_parts = []
    totals = []
    for i, tsv in enumerate(tsvs, 1):
        nodename, taxid = parse_clade_from_filename(tsv)
        df = pd.read_csv(tsv, sep="\t")
        n_total = len(df)

        # Full per-clade pipeline: +pseudocount -> log-ratio -> restrict
        # z-score base to occupancy>=0.5 -> flags.
        # plot2.py::main lines 434-437 filter the df to occupancy_in >= 0.5
        # BEFORE computing the sigma columns, so the mean/std used for
        # z-scoring come from the high-occupancy subset only (low-occupancy
        # pairs have noisy out-clade stats and would skew the distribution).
        # compute_z_scores() alone doesn't apply that filter, so we do it
        # explicitly here to match the original SuppTable_16 filter math.
        df = add_ratio_columns(df)
        df = df[df["occupancy_in"] >= OCC_MIN].copy()
        df = compute_z_scores(df)
        df = assign_flags(df, sd_number=args.sigma)

        kept = df[df[FLAG_COLS].sum(axis=1) > 0].copy()
        kept["nodename"] = nodename
        kept["taxid"] = taxid
        # Ortholog/RBH name joins via pair col_idx
        kept["ortholog1"] = kept["pair"].astype(int).map(lambda ix: ix_to_pair[ix][0])
        kept["ortholog2"] = kept["pair"].astype(int).map(lambda ix: ix_to_pair[ix][1])
        kept["rbh1"] = kept["ortholog1"].map(rbh_to_alg)
        kept["rbh2"] = kept["ortholog2"].map(rbh_to_alg)

        kept_parts.append(kept)
        totals.append((nodename, taxid, n_total, len(kept)))
        print(f"[{i:>2}/{len(tsvs)}] {nodename:<20} n={n_total:>10,}  kept={len(kept):>8,}  "
              f"({100*len(kept)/n_total:.2f}%)")

    result = pd.concat(kept_parts, ignore_index=True)

    # Column order matching SuppTable_16
    priority = ["nodename", "taxid", "ortholog1", "rbh1", "ortholog2", "rbh2"] + FLAG_COLS
    result = result[priority + [c for c in result.columns if c not in priority]]

    # Write
    args.output.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if str(args.output).endswith(".gz") else None
    result.to_csv(args.output, sep="\t", index=False, compression=compression)

    print()
    print(f"[build] wrote {len(result):,} rows -> {args.output}")
    print()
    print("=== per-clade kept counts ===")
    print(f"{'clade':<20}{'taxid':>12}{'total':>12}{'kept':>12}{'%':>8}")
    for nodename, taxid, n_total, n_kept in totals:
        print(f"{nodename:<20}{taxid:>12}{n_total:>12,}{n_kept:>12,}"
              f"{100*n_kept/n_total:>7.2f}%")
    grand_total = sum(t for _, _, t, _ in totals)
    grand_kept = len(result)
    print(f"{'TOTAL':<20}{'':>12}{grand_total:>12,}{grand_kept:>12,}"
          f"{100*grand_kept/grand_total:>7.2f}%")


if __name__ == "__main__":
    main()
