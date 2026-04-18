#!/usr/bin/env python3
"""Pair-distance scatter of SupplementaryTable_16 locus pairs per clade.

For each clade, plot every clade-row in (log10 mean_in, log10 mean_out)
space with all pairs as faint background dots. For each ranking axis
(stability / closeness / intersection) separately, highlight the pairs
that fell into the top-N foreground at whichever N produced the most-
significant GO enrichment for that (clade, axis) cell. Points below the
diagonal are closer in-clade than outside; points to the left have
short absolute in-clade distances.

Uses summary.tsv to pick best_N per (clade, axis) and the xlsx to pull
the raw mean_in / mean_out stats.
"""
import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

OCCUPANCY_MIN = 0.5


def pick_best_N(summary, clade, axis, min_fold=3.0):
    """Return (N, top_q, top_term, fold) for the lowest-q cell in (clade, axis)
    under sweep_namespace=='all' whose top term has fold >= min_fold.

    The fold gate suppresses the large-N hypergeometric artifact rebound
    (where broad terms like GO:0005829 "cytosol" flip to q<~1e-3 at
    fold~1.1 because n approaches N/2). Without this, best_N tends
    toward the full row count of the clade, which is not biologically
    meaningful.
    """
    s = summary[(summary["clade"] == clade) & (summary["axis"] == axis)
                & (summary["namespace"] == "all")]
    s = s.dropna(subset=["top_q", "top_term_fold"])
    s = s[s["top_term_fold"] >= min_fold]
    if s.empty:
        return None
    r = s.sort_values("top_q").iloc[0]
    return int(r["N_threshold"]), float(r["top_q"]), str(r["top_term"]), float(r["top_term_fold"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--supp-table", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print("[load] xlsx …")
    supp = pd.read_excel(args.supp_table, engine="openpyxl")
    print("[load] summary.tsv …")
    summary = pd.read_csv(args.summary, sep="\t")

    clades = sorted(supp["nodename"].dropna().unique())
    axes_order = ("stability", "closeness", "intersection")

    with PdfPages(args.out) as pdf:
        for clade in clades:
            sub = supp[supp["nodename"] == clade].copy()
            sub = sub[sub["occupancy_in"].fillna(0) >= OCCUPANCY_MIN]
            sub = sub.dropna(subset=["mean_in", "mean_out",
                                     "sd_in_out_ratio_log_sigma",
                                     "mean_in_out_ratio_log_sigma"])
            sub = sub.reset_index(drop=True)
            if sub.empty:
                continue

            log_mean_in = np.log10(sub["mean_in"].clip(lower=1))
            log_mean_out = np.log10(sub["mean_out"].clip(lower=1))

            # Pick best-N for stability and closeness independently. The
            # intersection panel then uses exactly the overlap of those two
            # foregrounds — NOT the sweep's own "intersection" best_N —
            # so the three panels are mutually consistent and the third
            # reads as "pairs highlighted in both of the first two".
            best_stab = pick_best_N(summary, clade, "stability")
            best_clos = pick_best_N(summary, clade, "closeness")
            if best_stab is None or best_clos is None:
                continue
            N_stab = best_stab[0]
            N_clos = best_clos[0]
            stab_idx = set(sub.sort_values("sd_in_out_ratio_log_sigma").index[:N_stab].tolist())
            clos_idx = set(sub.sort_values("mean_in_out_ratio_log_sigma").index[:N_clos].tolist())
            panel_info = {
                "stability": (N_stab, best_stab[1], best_stab[2], best_stab[3],
                              pd.Index(sorted(stab_idx))),
                "closeness": (N_clos, best_clos[1], best_clos[2], best_clos[3],
                              pd.Index(sorted(clos_idx))),
            }
            # Intersection = overlap of the two red-dot sets above.
            # Pull fold/q from summary only if a matching (clade, axis=
            # "intersection", N=...) row exists for reference; otherwise
            # fall back to "derived from panels 1+2".
            inter_set = stab_idx & clos_idx
            panel_info["intersection"] = (
                None, None, None, None, pd.Index(sorted(inter_set)),
            )

            # Percentile-trim the bulk BUT always include every foreground
            # dot — stable/close pairs often sit in the distribution's
            # tails by design, so clipping them out defeats the figure's
            # purpose. Using the union of stability + closeness + any
            # intersection foregrounds for the "must include" set.
            fg_any = stab_idx | clos_idx | inter_set
            if fg_any:
                fg_in = log_mean_in.iloc[sorted(fg_any)]
                fg_out = log_mean_out.iloc[sorted(fg_any)]
                fg_lo = float(min(fg_in.min(), fg_out.min()))
                fg_hi = float(max(fg_in.max(), fg_out.max()))
            else:
                fg_lo = float("inf")
                fg_hi = float("-inf")
            p1 = min(float(np.percentile(log_mean_in, 1)),
                      float(np.percentile(log_mean_out, 1)),
                      fg_lo)
            p99 = max(float(np.percentile(log_mean_in, 99)),
                       float(np.percentile(log_mean_out, 99)),
                       fg_hi)
            span = p99 - p1
            lim_lo = p1 - 0.05 * span
            lim_hi = p99 + 0.05 * span

            fig, axes = plt.subplots(1, 3, figsize=(15, 5.5),
                                     sharex=True, sharey=True)
            for col, axis in enumerate(axes_order):
                ax = axes[col]
                N_used, best_q, top_term, top_fold, fg_idx = panel_info[axis]
                fg_mask = sub.index.isin(fg_idx)

                # Background pairs (all clade rows not in the foreground).
                ax.scatter(log_mean_in[~fg_mask], log_mean_out[~fg_mask],
                           c="#cccccc", s=5, alpha=0.25,
                           edgecolors="none",
                           label=f"other pairs (n={int((~fg_mask).sum())})",
                           zorder=2)
                # Foreground.
                if axis == "intersection":
                    fg_label = (f"top-{N_stab}-stable ∩ top-{N_clos}-close "
                                f"({fg_mask.sum()} pairs)")
                else:
                    fg_label = f"top-{N_used} ({fg_mask.sum()} pairs)"
                ax.scatter(log_mean_in[fg_mask], log_mean_out[fg_mask],
                           c="#d62728", s=14, alpha=0.8,
                           edgecolors="white", linewidths=0.3,
                           label=fg_label, zorder=3)

                # Diagonal across the equal-limit square.
                ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], ls=":",
                        color="black", lw=0.6, alpha=0.6, label="y = x")
                ax.set_xlim(lim_lo, lim_hi)
                ax.set_ylim(lim_lo, lim_hi)
                ax.set_aspect("equal", adjustable="box")

                # Cosmetics.
                if col == 0:
                    ax.set_ylabel("log10 mean_out  (bp, pair distance outside clade)")
                ax.set_xlabel("log10 mean_in  (bp, pair distance in clade)")
                if axis == "intersection":
                    title = (f"{clade} — intersection\n"
                             f"pairs in both top-{N_stab}-stable and top-{N_clos}-close")
                else:
                    q_label = f"q={best_q:.2e}" if best_q and best_q > 0 else "q=0"
                    title = (f"{clade} — {axis}\nbest_N={N_used}  "
                             f"top={top_term}  fold={top_fold:.1f}×  {q_label}")
                ax.set_title(title, fontsize=9)
                ax.legend(fontsize=7, loc="upper left")
                ax.grid(alpha=0.15)
            fig.suptitle(clade, fontsize=13, y=1.02)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"[write] {args.out}  clades_plotted={len(clades)}")


if __name__ == "__main__":
    main()
