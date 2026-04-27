"""Stage 2: distribution of stable-pair distances.

Page 1: focal panels (Vertebrata + Diptera by default) — qualifying-pair
distribution drawn on a real-bp log axis with human-readable tick labels
(100 bp / 1 kb / 1 Mb / 10 Mb), with all-clade background overlay and
selected examples annotated.
Page 2: per-clade small multiples on the same axis.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# TrueType fonts (Type 42) so PDF text remains editable in Illustrator/Inkscape.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

BIN_COLORS = {
    "close":        "#1f77b4",
    "intermediate": "#2ca02c",
    "long_range":   "#d62728",
    "ultra_long":   "#9467bd",
}
BASIS_COLORS = {
    "forced":         "#1f77b4",
    "percentile_p5":  "#2ca02c",
    "percentile_p25": "#17becf",
    "percentile_p50": "#9467bd",
    "percentile_p75": "#ff7f0e",
    "percentile_p95": "#d62728",
}

X_LO_BP, X_HI_BP = 1e2, 1e8


def human_bp(x: float, _pos=None) -> str:
    if x < 1e3:
        return f"{int(x)} bp"
    if x < 1e6:
        return f"{int(x/1e3)} kb"
    if x < 1e9:
        return f"{int(x/1e6)} Mb"
    return f"{int(x/1e9)} Gb"


def style_log_bp_axis(ax):
    ax.set_xscale("log")
    ax.set_xlim(X_LO_BP, X_HI_BP)
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(human_bp))
    ax.xaxis.set_minor_locator(
        mticker.LogLocator(base=10, subs=tuple(np.arange(2, 10) * 0.1),
                            numticks=20))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())


def positive_mean_in(df: pd.DataFrame) -> np.ndarray:
    v = df["mean_in"].astype(float).values
    return v[v > 0]


def annotate_examples(ax, sub: pd.DataFrame, recomputed: dict | None,
                       ymax: float):
    """Draw vertical dashed arrows for each example pair, labelled
    pair_label + distance. If TSS-aware recompute is available for that
    pair (recomputed[pair_key]), label shows TSS-TSS in addition."""
    sub = sub.sort_values("mean_in")
    for _, r in sub.iterrows():
        x = float(r["mean_in"]) if r["mean_in"] > 0 else None
        if x is None:
            continue
        basis = r.get("selection_basis", "forced") or "forced"
        c = BASIS_COLORS.get(basis, "0.4")
        ax.axvline(x, color=c, lw=1.0, ls="--", alpha=0.85)
        label = str(r["pair_label"])
        if recomputed and r["pair_key"] in recomputed:
            tss = recomputed[r["pair_key"]].get("median_tss_to_tss_bp")
            if tss is not None and not np.isnan(tss):
                label = f"{label} · TSS {human_bp(tss)} (edge {human_bp(x)})"
            else:
                label = f"{label} · {human_bp(x)}"
        else:
            label = f"{label} · {human_bp(x)}"
        ax.annotate(label, xy=(x, ymax * 0.95),
                    xytext=(3, 0), textcoords="offset points",
                    fontsize=7, rotation=90, va="top", color=c)


def draw_focal_panel(ax, cand: pd.DataFrame, sel: pd.DataFrame, focal: str,
                      recomputed: dict | None):
    bins = np.logspace(np.log10(X_LO_BP), np.log10(X_HI_BP), 60)
    fg = positive_mean_in(cand[cand["nodename"] == focal])
    if len(fg):
        ax.hist(fg, bins=bins, color="#1f77b4", alpha=0.85,
                label=focal, density=True)
    style_log_bp_axis(ax)
    ax.set_xlabel("in-clade mean distance between anchor genes")
    ax.set_ylabel("density")
    ax.set_title(f"{focal} stable-pair distance distribution")
    ax.legend(loc="upper left", fontsize=8)
    sub = sel[sel["nodename"] == focal]
    annotate_examples(ax, sub, recomputed, ax.get_ylim()[1])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--selected", required=True)
    ap.add_argument("--recomputed", default=None,
                    help="optional selected_examples_recomputed.tsv from Stage 1.5")
    ap.add_argument("--focal", nargs="+", default=["Vertebrata", "Diptera"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    cand = pd.read_csv(args.candidates, sep="\t")
    sel = pd.read_csv(args.selected, sep="\t")
    recomputed = None
    if args.recomputed and Path(args.recomputed).exists():
        rdf = pd.read_csv(args.recomputed, sep="\t")
        recomputed = {row["pair_key"]: row.to_dict() for _, row in rdf.iterrows()}

    with PdfPages(args.out) as pdf:
        # Page 1: focal panels stacked.
        fig, axes = plt.subplots(len(args.focal), 1,
                                  figsize=(9.0, 3.4 * len(args.focal)),
                                  squeeze=False)
        for ax, focal in zip(axes[:, 0], args.focal):
            draw_focal_panel(ax, cand, sel, focal, recomputed)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # Page 2: per-clade small multiples, same x-axis style.
        clades = sorted(cand["nodename"].unique())
        ncol = 4
        nrow = int(np.ceil(len(clades) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(11, 1.9 * nrow),
                                  squeeze=False)
        bins = np.logspace(np.log10(X_LO_BP), np.log10(X_HI_BP), 40)
        for i, clade in enumerate(clades):
            ax = axes[i // ncol, i % ncol]
            v = positive_mean_in(cand[cand["nodename"] == clade])
            if len(v):
                ax.hist(v, bins=bins, color="#1f77b4")
            style_log_bp_axis(ax)
            ax.set_title(clade, fontsize=8)
            ax.tick_params(labelsize=6)
            for label in ax.get_xticklabels():
                label.set_rotation(30)
                label.set_horizontalalignment("right")
        for j in range(len(clades), nrow * ncol):
            axes[j // ncol, j % ncol].set_axis_off()
        fig.suptitle("stable-pair distances per clade", fontsize=10)
        fig.tight_layout()
        pdf.savefig(fig); plt.close(fig)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
