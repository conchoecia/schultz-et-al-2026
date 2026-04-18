#!/usr/bin/env python3
"""Volcano plots of GO-enrichment results from significant_terms.tsv.

x = log2(fold_enrichment), y = -log10(q). One page per clade, three
subplots per page (stability / closeness / intersection). Each dot is
a (N_threshold, GO-term) config labeled with its N; colored by GO
namespace (BP/MF/CC). Separates genuine enrichment (top-right, high
fold + high significance) from large-N hypergeometric artifacts
(left edge, fold ≈ 1 even at low q).

Uses the conservative sweep_namespace=="all" rows so each term is
plotted at its strictest BH-FDR (the all-family correction).
"""
import argparse
import math
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

NS_COLOR = {"BP": "#1f77b4", "MF": "#ff7f0e", "CC": "#2ca02c", "?": "#999999"}


def log10_safe(q):
    try:
        qf = float(q)
    except (TypeError, ValueError):
        return float("nan")
    if qf <= 0:
        return 300.0
    return -math.log10(qf)


def log2_safe(f):
    try:
        ff = float(f)
    except (TypeError, ValueError):
        return float("nan")
    if not math.isfinite(ff) or ff <= 0:
        return float("nan")
    return math.log2(ff)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--significant-terms", required=True)
    ap.add_argument("--out", required=True, help="output PDF")
    ap.add_argument("--min-fold", type=float, default=None,
                    help="optional: only plot dots with fold >= this "
                         "(use 3.0 to suppress large-N artifacts)")
    args = ap.parse_args()

    df = pd.read_csv(args.significant_terms, sep="\t")
    # Use the conservative "all" sweep_namespace to avoid plotting the
    # same term three times (once per namespace BH family).
    df = df[df["sweep_namespace"] == "all"].copy()
    # Each GO term now appears at multiple N thresholds per (clade, axis).
    # Dedupe to one row per term, keeping the cell with the lowest q
    # (strongest evidence). The N label on that dot then indicates where
    # the term's enrichment peaks.
    df = df.sort_values("q").drop_duplicates(
        subset=["clade", "axis", "go_id"], keep="first")
    df["log2fold"] = df["fold"].map(log2_safe)
    df["mlog10q"] = df["q"].map(log10_safe)
    df = df.dropna(subset=["log2fold", "mlog10q"])
    if args.min_fold is not None:
        df = df[df["fold"] >= args.min_fold]

    clades = sorted(df["clade"].dropna().unique())
    axes_order = ("stability", "closeness", "intersection")

    with PdfPages(args.out) as pdf:
        for clade in clades:
            sub = df[df["clade"] == clade]
            if sub.empty:
                continue
            fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
            for col, axis in enumerate(axes_order):
                ax = axes[col]
                s = sub[sub["axis"] == axis]
                if s.empty:
                    ax.set_title(f"{clade} — {axis}  (no sig terms)")
                    ax.set_xlabel("log2 fold-enrichment")
                    if col == 0:
                        ax.set_ylabel("-log10(q)")
                    continue
                # Plot one color per namespace.
                for ns in ("BP", "MF", "CC", "?"):
                    t = s[s["go_namespace"] == ns]
                    if t.empty:
                        continue
                    ax.scatter(t["log2fold"], t["mlog10q"],
                               c=NS_COLOR[ns], s=28, alpha=0.7,
                               edgecolors="white", linewidths=0.5,
                               label=ns if ns != "?" else None, zorder=3)
                # Annotate each dot with its N_threshold.
                for _, r in s.iterrows():
                    ax.annotate(str(int(r["N_threshold"])),
                                xy=(r["log2fold"], r["mlog10q"]),
                                xytext=(3, 2), textcoords="offset points",
                                fontsize=5.5, color="#333", alpha=0.8,
                                zorder=4)
                ax.axhline(-math.log10(0.05), ls="--", color="red", lw=0.6,
                           label="q=0.05" if col == 0 else None)
                ax.axhline(-math.log10(0.25), ls=":", color="orange", lw=0.6,
                           label="q=0.25" if col == 0 else None)
                ax.axvline(math.log2(3), ls=":", color="#888", lw=0.6,
                           label="fold=3×" if col == 0 else None)
                ax.axvline(0, color="black", lw=0.4, alpha=0.4)
                ax.set_title(f"{clade} — {axis}  (n={len(s)})")
                ax.set_xlabel("log2 fold-enrichment")
                if col == 0:
                    ax.set_ylabel("-log10(q)")
                    ax.legend(fontsize=7, loc="upper left")
            fig.suptitle(clade, fontsize=12, y=1.02)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"[write] {args.out}  clades_plotted={len(clades)}")


if __name__ == "__main__":
    main()
