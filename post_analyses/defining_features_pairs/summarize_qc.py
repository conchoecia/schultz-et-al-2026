#!/usr/bin/env python3
"""One-page selection-overview QC PDF per clade (and one cross-clade summary).

The existing full QC PDFs (``*_unique_pair_df_qc.pdf``) contain the 5-panel
distribution plots that the decay-branch code produced — useful for deep
inspection but unwieldy for a quick "what did we select" read.

This companion writes ``*_unique_pair_df_qc_summary.pdf`` — one page per
clade — showing the flag-selection regions overlaid on the actual density
of all pairs. Flagged pairs are highlighted; selection region is shaded;
a barchart reports the count of each flag.

Also writes ``all_clades_selection_summary.pdf`` — one page with per-clade
flag counts, side by side, for cross-clade comparison.

Source of truth for cutoffs: egt/src/egt/legacy/defining_features_plot2.py
(same `sd_number=2` default used in SupplementaryTable_16).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle


SD_NUMBER = 2         # sigma cutoff (matches SuppTable_16)
OCC_MIN   = 0.5       # occupancy_in gate for the 4 non-unique flags


def add_ratio_and_sigma(df, pseudocount=1):
    """Full pipeline: +pseudocount -> log-ratio -> per-clade z-scores on occ>=0.5."""
    df = df.copy()
    df["mean_in_p"]  = df["mean_in"]  + pseudocount
    df["mean_out_p"] = df["mean_out"] + pseudocount
    df["sd_in_p"]    = df["sd_in"]    + pseudocount
    df["sd_out_p"]   = df["sd_out"]   + pseudocount
    df["mean_log"] = np.log10(df["mean_in_p"] / df["mean_out_p"])
    df["sd_log"]   = np.log10(df["sd_in_p"]   / df["sd_out_p"])

    sub = df[df["occupancy_in"] >= OCC_MIN]
    m_mean, s_mean = sub["mean_log"].mean(), sub["mean_log"].std()
    m_sd,   s_sd   = sub["sd_log"].mean(),   sub["sd_log"].std()
    df["mean_sigma"] = (df["mean_log"] - m_mean) / (s_mean if s_mean > 0 else 1)
    df["sd_sigma"]   = (df["sd_log"]   - m_sd)   / (s_sd   if s_sd   > 0 else 1)
    return df


def assign_flags(df):
    close    = (df["mean_sigma"] < -SD_NUMBER) & (df["occupancy_in"] >= OCC_MIN)
    distant  = (df["mean_sigma"] >  SD_NUMBER) & (df["occupancy_in"] >= OCC_MIN)
    stable   = (df["sd_sigma"]   < -SD_NUMBER) & (df["occupancy_in"] >= OCC_MIN)
    unstable = (df["sd_sigma"]   >  SD_NUMBER) & (df["occupancy_in"] >= OCC_MIN)
    unique   = df["notna_out"] == 0
    return {
        "close_in_clade":    close,
        "distant_in_clade":  distant,
        "stable_in_clade":   stable,
        "unstable_in_clade": unstable,
        "unique_to_clade":   unique,
    }


def _hexbin_with_selection(ax, x, y, xlabel, ylabel,
                           neg_select, pos_select, n_neg, n_pos,
                           xlim=None, ylim=None, occ_line=OCC_MIN, sigma=SD_NUMBER):
    """Hexbin of (x, y), vertical lines at ±sigma, horizontal at occ_line,
    red shade over the neg_select region, blue shade over pos_select region.
    Overlay flagged points as dots on top. The y < occ_line region is
    greyed because the σ cutoff is measured only over the y >= occ_line
    subset (matches plot2.py::main line 434).
    """
    ax.hexbin(x, y, gridsize=60, mincnt=1, bins="log", cmap="Greys")
    xmin, xmax = (np.nanmin(x), np.nanmax(x)) if xlim is None else xlim

    # close / stable region (left side) and distant / unstable region (right side)
    ax.add_patch(Rectangle((xmin, occ_line), -sigma - xmin, 1 - occ_line,
                           alpha=0.07, color="red", zorder=0))
    ax.add_patch(Rectangle((sigma, occ_line), xmax - sigma, 1 - occ_line,
                           alpha=0.07, color="blue", zorder=0))

    # Cutoff lines
    ax.axvline(-sigma, color="red",   lw=0.8, ls="--", alpha=0.7, zorder=6)
    ax.axvline( sigma, color="blue",  lw=0.8, ls="--", alpha=0.7, zorder=6)
    ax.axhline(occ_line, color="k",   lw=0.8, ls=":",  alpha=0.7, zorder=6)

    # σ value labels — placed INSIDE the plot at the top so they don't
    # collide with the subplot titles above.
    ax.text(-sigma, 0.97, f"−{sigma}σ", ha="center", va="top",
            color="red",  fontsize=7.5, fontweight="bold",
            transform=ax.get_xaxis_transform(), zorder=7)
    ax.text( sigma, 0.97, f"+{sigma}σ", ha="center", va="top",
            color="blue", fontsize=7.5, fontweight="bold",
            transform=ax.get_xaxis_transform(), zorder=7)

    # Overlay flagged points as scatter
    if n_neg > 0:
        ax.scatter(x[neg_select], y[neg_select], s=2, color="red",
                   alpha=0.4, edgecolors="none", label=f"{n_neg:,} flagged")
    if n_pos > 0:
        ax.scatter(x[pos_select], y[pos_select], s=2, color="blue",
                   alpha=0.4, edgecolors="none", label=f"{n_pos:,} flagged")
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(labelsize=7)


def write_summary_pdf(df_raw, out_pdf, nodename, taxid):
    df = add_ratio_and_sigma(df_raw)
    flags = assign_flags(df)
    total = len(df)
    counts = {k: int(v.sum()) for k, v in flags.items()}
    any_flag = (flags["close_in_clade"] | flags["distant_in_clade"] |
                flags["stable_in_clade"] | flags["unstable_in_clade"] |
                flags["unique_to_clade"])
    counts["any_flag_kept"] = int(any_flag.sum())

    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle(f"{nodename} (taxid {taxid}) — "
                 f"{counts['any_flag_kept']:,} / {total:,} pairs kept "
                 f"({100*counts['any_flag_kept']/total:.2f}%)",
                 fontsize=11, y=0.98)

    # Top-left: mean_sigma vs occupancy_in (close / distant cutoffs)
    ax1 = fig.add_subplot(2, 2, 1)
    _hexbin_with_selection(
        ax1, df["mean_sigma"].to_numpy(), df["occupancy_in"].to_numpy(),
        "log$_{10}$( mean_in / mean_out )  [z-scored]\n"
        "← in-clade distances SHORTER          "
        "in-clade distances LONGER →",
        "occupancy_in\n(fraction of in-clade species with this pair)",
        flags["close_in_clade"],   flags["distant_in_clade"],
        counts["close_in_clade"],  counts["distant_in_clade"],
    )
    # Color-coded title: close (red, left) | distant (blue, right)
    ax1.set_title(f"close = {counts['close_in_clade']:,}",
                  loc="left", color="red", fontsize=10)
    ax1.set_title(f"distant = {counts['distant_in_clade']:,}",
                  loc="right", color="blue", fontsize=10)

    # Top-right: sd_sigma vs occupancy_in (stable / unstable cutoffs)
    ax2 = fig.add_subplot(2, 2, 2)
    _hexbin_with_selection(
        ax2, df["sd_sigma"].to_numpy(), df["occupancy_in"].to_numpy(),
        "log$_{10}$( sd_in / sd_out )  [z-scored]\n"
        "← in-clade distances MORE UNIFORM          "
        "in-clade distances MORE VARIABLE →",
        "occupancy_in\n(fraction of in-clade species with this pair)",
        flags["stable_in_clade"],   flags["unstable_in_clade"],
        counts["stable_in_clade"],  counts["unstable_in_clade"],
    )
    ax2.set_title(f"stable = {counts['stable_in_clade']:,}",
                  loc="left", color="red", fontsize=10)
    ax2.set_title(f"unstable = {counts['unstable_in_clade']:,}",
                  loc="right", color="blue", fontsize=10)

    # Bottom-left: flag count barchart (linear scale; count labels tell the story)
    ax3 = fig.add_subplot(2, 2, 3)
    names  = ["close", "distant", "stable", "unstable", "unique", "any_kept"]
    keys   = ["close_in_clade", "distant_in_clade", "stable_in_clade",
              "unstable_in_clade", "unique_to_clade", "any_flag_kept"]
    vals   = [counts[k] for k in keys]
    colors = ["red", "blue", "red", "blue", "green", "black"]
    bars = ax3.barh(names, vals, color=colors, alpha=0.7)
    max_v = max(vals) if max(vals) > 0 else 1
    # Put count labels INSIDE tall bars, OUTSIDE short ones (so they're always readable)
    for bar, v in zip(bars, vals):
        ha = "right" if v > 0.3 * max_v else "left"
        offset = -0.01 * max_v if ha == "right" else 0.01 * max_v
        text_color = "white" if ha == "right" else "black"
        ax3.text(v + offset, bar.get_y() + bar.get_height()/2,
                 f"{v:,}  ({100*v/total:.2f}%)",
                 va="center", ha=ha, fontsize=7, color=text_color)
    ax3.set_xlabel("pair count (linear)", fontsize=8)
    ax3.set_title("flag counts", fontsize=9)
    ax3.tick_params(labelsize=7)
    ax3.set_xlim(left=0)

    # Bottom-right: text summary
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis("off")
    lines = [
        f"clade:  {nodename}  (taxid {taxid})",
        f"total pairs with >=1 in-clade observation:  {total:,}",
        f"pairs with occupancy_in >= {OCC_MIN}:  {int((df['occupancy_in']>=OCC_MIN).sum()):,}",
        f"pairs with notna_out == 0:  {counts['unique_to_clade']:,}",
        "",
        f"sigma cutoff:  ±{SD_NUMBER}",
        f"occupancy gate (non-unique flags):  >= {OCC_MIN}",
        "",
        "flag counts:",
        f"  close_in_clade     {counts['close_in_clade']:>10,}",
        f"  distant_in_clade   {counts['distant_in_clade']:>10,}",
        f"  stable_in_clade    {counts['stable_in_clade']:>10,}",
        f"  unstable_in_clade  {counts['unstable_in_clade']:>10,}",
        f"  unique_to_clade    {counts['unique_to_clade']:>10,}",
        "",
        f"union (any flag set):  {counts['any_flag_kept']:,}  ({100*counts['any_flag_kept']/total:.3f}%)",
    ]
    ax4.text(0.02, 0.98, "\n".join(lines), va="top", ha="left",
             fontsize=8.5, family="monospace",
             transform=ax4.transAxes)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

    return counts


def write_cross_clade_summary(counts_by_clade, out_pdf):
    """2x3 grid of small barcharts — one panel per flag category, 28 clades
    (alphabetical) as horizontal bars on a linear scale. Each panel is
    readable on its own; no overlapping bars, no symlog.
    """
    clades = sorted(counts_by_clade.keys())
    y = np.arange(len(clades))
    panels = [
        ("close_in_clade",    "red",   "close"),
        ("distant_in_clade",  "blue",  "distant"),
        ("stable_in_clade",   "red",   "stable"),
        ("unstable_in_clade", "blue",  "unstable"),
        ("unique_to_clade",   "green", "unique"),
        ("any_flag_kept",     "black", "any flag kept (union)"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(16, 0.28 * len(clades) + 2),
                             sharey=True)
    fig.suptitle(f"Per-clade flag counts  (σ = ±{SD_NUMBER}, occupancy_in ≥ {OCC_MIN})",
                 fontsize=12, y=0.995)
    for ax, (key, color, label) in zip(axes.ravel(), panels):
        vals = [counts_by_clade[c][key] for c in clades]
        bars = ax.barh(y, vals, color=color, alpha=0.75)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(v, bar.get_y() + bar.get_height()/2,
                        f" {v:,}", va="center", ha="left", fontsize=6.5)
        ax.set_title(label, fontsize=10, color=color)
        ax.set_xlabel("pair count", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.set_xlim(left=0)
        ax.grid(axis="x", alpha=0.25)
    axes[0, 0].set_yticks(y)
    axes[0, 0].set_yticklabels(clades, fontsize=7.5)
    axes[1, 0].set_yticks(y)
    axes[1, 0].set_yticklabels(clades, fontsize=7.5)

    plt.tight_layout(rect=[0, 0, 1, 0.98])
    with PdfPages(out_pdf) as pdf:
        pdf.savefig(fig)
    plt.close(fig)


def parse_clade_from_filename(path):
    m = re.match(r"^(?P<name>.+?)_(?P<taxid>\d+)_unique_pair_df\.tsv\.gz$",
                 Path(path).name)
    if not m:
        raise ValueError(f"unexpected filename: {path}")
    return m.group("name"), int(m.group("taxid"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("in_dir", type=Path,
                    help="Directory of *_unique_pair_df.tsv.gz")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="Output dir (default: in_dir).")
    args = ap.parse_args(argv)
    out_dir = args.out_dir or args.in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    tsvs = sorted(args.in_dir.glob("*_unique_pair_df.tsv.gz"))
    if not tsvs:
        raise SystemExit(f"no *_unique_pair_df.tsv.gz files in {args.in_dir}")
    print(f"[summarize_qc] found {len(tsvs)} per-clade TSVs")

    counts_by_clade = {}
    for i, tsv in enumerate(tsvs, 1):
        nodename, taxid = parse_clade_from_filename(tsv)
        print(f"[{i:>2}/{len(tsvs)}] {nodename} ({taxid}) ... ", end="", flush=True)
        df = pd.read_csv(tsv, sep="\t")
        out_pdf = out_dir / tsv.name.replace(".tsv.gz", "_qc_summary.pdf")
        counts = write_summary_pdf(df, out_pdf, nodename, taxid)
        counts_by_clade[nodename] = counts
        print(f"kept {counts['any_flag_kept']:,}/{len(df):,} -> {out_pdf.name}")

    cross = out_dir / "all_clades_selection_summary.pdf"
    write_cross_clade_summary(counts_by_clade, cross)
    print(f"[summarize_qc] cross-clade summary -> {cross}")

    # Also print a cross-clade table
    print("\n=== per-clade flag counts (cross-clade table) ===")
    cols = ["close_in_clade", "distant_in_clade", "stable_in_clade",
            "unstable_in_clade", "unique_to_clade", "any_flag_kept"]
    print(f"{'clade':<20}" + "".join(f"{c:>18}" for c in cols))
    for c in sorted(counts_by_clade):
        print(f"{c:<20}" + "".join(f"{counts_by_clade[c][k]:>18,}" for k in cols))


if __name__ == "__main__":
    main()
