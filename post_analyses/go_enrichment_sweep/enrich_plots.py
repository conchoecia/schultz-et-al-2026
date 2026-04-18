#!/usr/bin/env python3
"""Publication-style GO enrichment plots from the sweep output.

Produces three artifacts conventional for GO-enrichment papers:

1. significant_terms_annotated.tsv -- significant_terms.tsv augmented
   with a human-readable go_name column from go-basic.obo.

2. dotplots.pdf -- one page per clade, clusterProfiler-style dotplot:
   y = GO term name (top-20 unique terms by best q), x = log2 fold,
   dot size = k (gene overlap), dot color = -log10(q). Faceted into
   BP / MF / CC subpanels.

3. heatmap.pdf -- cross-clade heatmap of the most-recurring enriched
   terms (rows) across clades (columns). Cell = -log10(q). Surfaces
   pan-clade vs clade-specific signals. One page per namespace.

For each (clade, go_id) we pick the minimum q across sweep cells -- this
summarizes "best evidence that this term is enriched in this clade."
The q is still the within-cell BH q from the sweep; no extra correction
is applied.
"""
import argparse
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ---------------------------------------------------------------------
# OBO name + namespace parsing
# ---------------------------------------------------------------------
OBO_NS_MAP = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}


def parse_obo(path):
    """Return {go_id: (name, namespace_abbrev)}."""
    term_info = {}
    cur_id = cur_name = cur_ns = None
    in_term = False
    obsolete = False
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line == "[Term]":
                if in_term and cur_id and cur_name and cur_ns and not obsolete:
                    term_info[cur_id] = (cur_name, cur_ns)
                cur_id = cur_name = cur_ns = None
                in_term = True
                obsolete = False
                continue
            if line.startswith("[") and line != "[Term]":
                if in_term and cur_id and cur_name and cur_ns and not obsolete:
                    term_info[cur_id] = (cur_name, cur_ns)
                in_term = False
                continue
            if not in_term or not line:
                continue
            if line.startswith("id: "):
                cur_id = line[4:].strip()
            elif line.startswith("name: "):
                cur_name = line[6:].strip()
            elif line.startswith("namespace: "):
                cur_ns = OBO_NS_MAP.get(line[11:].strip())
            elif line.startswith("is_obsolete: true"):
                obsolete = True
    if in_term and cur_id and cur_name and cur_ns and not obsolete:
        term_info[cur_id] = (cur_name, cur_ns)
    return term_info


# ---------------------------------------------------------------------
# Dotplot (clusterProfiler style)
# ---------------------------------------------------------------------
def draw_dotplot(fig_ax, terms_df, title):
    """Draw one dotplot on an existing matplotlib Axes.

    terms_df: columns ['go_name', 'log2fold', 'mlog10q', 'k'] sorted by
    mlog10q desc (highest-significance first at top of y-axis).
    """
    if terms_df.empty:
        fig_ax.set_title(title + "  (no hits)", fontsize=9)
        fig_ax.set_xticks([])
        fig_ax.set_yticks([])
        return None
    ys = np.arange(len(terms_df))[::-1]  # highest-q at top
    sc = fig_ax.scatter(
        terms_df["log2fold"], ys,
        s=40 + 10 * terms_df["k"].clip(upper=50),
        c=terms_df["mlog10q"],
        cmap="viridis",
        edgecolors="#333", linewidths=0.4,
    )
    fig_ax.set_yticks(ys)
    fig_ax.set_yticklabels(
        [f"{n[:56]}" + ("…" if len(n) > 56 else "")
         for n in terms_df["go_name"]],
        fontsize=7)
    fig_ax.axvline(math.log2(3), ls=":", color="gray", lw=0.5)
    fig_ax.axvline(0, color="black", lw=0.3, alpha=0.3)
    fig_ax.set_xlabel("log2 fold-enrichment", fontsize=8)
    fig_ax.set_title(title, fontsize=9)
    fig_ax.grid(axis="x", alpha=0.15)
    return sc


def make_dotplots(sig_df, out_path, top_n=15, min_fold=3.0):
    """Per clade, three subplots (BP/MF/CC), each showing top-N terms.

    Applies a fold-enrichment gate BEFORE per-term deduplication so the
    large-N hypergeometric artifact (broad terms like GO:0005829 flipping
    to q~1e-5 at fold~1.1 when the foreground approaches N/2) doesn't
    overwrite the real narrow-N signal.
    """
    clades = sorted(sig_df["clade"].dropna().unique())
    # Gate: only keep cells whose fold is high enough to be biologically
    # meaningful. Dedupe afterward.
    gated = sig_df[sig_df["fold"] >= min_fold]
    with PdfPages(out_path) as pdf:
        for clade in clades:
            sub = gated[gated["clade"] == clade].copy()
            if sub.empty:
                continue
            best = (sub.sort_values("q")
                       .drop_duplicates(subset=["go_id"]))
            best["log2fold"] = best["fold"].map(
                lambda f: math.log2(f) if f and f > 0 and math.isfinite(f) else float("nan"))
            best["mlog10q"] = best["q"].map(
                lambda q: 300.0 if q == 0 else (-math.log10(q) if q and q > 0 else float("nan")))
            best = best.dropna(subset=["log2fold", "mlog10q"])

            fig, axes = plt.subplots(1, 3, figsize=(17, 6))
            mappables = []
            per_ns_sizes = []
            for col, ns in enumerate(("BP", "MF", "CC")):
                s = best[best["go_namespace"] == ns]
                s = s.sort_values("q").head(top_n)
                sc = draw_dotplot(axes[col], s, f"{ns}  (top {len(s)})")
                if sc is not None:
                    mappables.append(sc)
                    per_ns_sizes.extend(s["k"].tolist())
            fig.suptitle(f"{clade} — top GO enrichments per namespace  "
                          f"(fold ≥ {min_fold:g}×)",
                          fontsize=11, y=1.01)
            # Right-side colorbar for -log10(q).
            if mappables:
                cbar = fig.colorbar(mappables[-1], ax=axes,
                                     orientation="vertical", fraction=0.03,
                                     pad=0.02)
                cbar.set_label("-log10(q)", fontsize=8)
            # Size legend: k = foreground genes annotated to that term
            # (dot size scales as 40 + 10 * min(k, 50)).
            if per_ns_sizes:
                k_min = int(min(per_ns_sizes))
                k_max = min(50, int(max(per_ns_sizes)))
                k_mid = max(k_min + 1, (k_min + k_max) // 2)
                size_vals = sorted({k_min, k_mid, k_max})
                handles = [plt.scatter([], [], s=40 + 10 * min(k, 50),
                                       c="#777", edgecolors="#333",
                                       linewidths=0.4,
                                       label=f"k = {k}"
                                       + ("+" if k >= 50 else ""))
                           for k in size_vals]
                axes[-1].legend(handles=handles, title="k (foreground\ngenes in term)",
                                 loc="center left",
                                 bbox_to_anchor=(1.25, 0.5),
                                 fontsize=7, title_fontsize=7,
                                 frameon=True)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"[write] {out_path}  clades={len(clades)}  top_n_per_ns={top_n}  "
          f"fold_gate={min_fold:g}")


# ---------------------------------------------------------------------
# Cross-clade heatmap
# ---------------------------------------------------------------------
def make_heatmap(sig_df, out_path, min_clades=2, max_terms_per_ns=40,
                  min_fold=3.0):
    """Pivot (term × clade) with -log10(q). One page per namespace.

    Gates on fold >= min_fold to suppress large-N hypergeometric artifacts.
    """
    gated = sig_df[sig_df["fold"] >= min_fold]
    best = (gated.sort_values("q")
                 .drop_duplicates(subset=["clade", "go_id"]))
    best["mlog10q"] = best["q"].map(
        lambda q: 300.0 if q == 0 else (-math.log10(q) if q and q > 0 else float("nan")))
    best = best.dropna(subset=["mlog10q"])
    clades = sorted(best["clade"].dropna().unique())
    with PdfPages(out_path) as pdf:
        for ns in ("BP", "MF", "CC"):
            nsub = best[best["go_namespace"] == ns]
            if nsub.empty:
                continue
            # Filter to terms hit in at least min_clades clades.
            term_counts = nsub.groupby("go_id")["clade"].nunique()
            terms_kept = term_counts[term_counts >= min_clades].index.tolist()
            if not terms_kept:
                # Fallback: keep the most widely hit terms anyway.
                terms_kept = term_counts.sort_values(ascending=False).head(
                    max_terms_per_ns).index.tolist()
            # Rank remaining terms by count then by median significance.
            if len(terms_kept) > max_terms_per_ns:
                ranking = (nsub[nsub["go_id"].isin(terms_kept)]
                           .groupby("go_id")
                           .agg(n_clades=("clade", "nunique"),
                                med_mlog10q=("mlog10q", "median"))
                           .sort_values(["n_clades", "med_mlog10q"],
                                         ascending=[False, False]))
                terms_kept = ranking.head(max_terms_per_ns).index.tolist()

            mat = nsub[nsub["go_id"].isin(terms_kept)].pivot_table(
                index="go_id", columns="clade",
                values="mlog10q", aggfunc="max")
            mat = mat.reindex(index=terms_kept, columns=clades)
            labels = [f"{gid}  {nsub[nsub.go_id == gid].iloc[0]['go_name'][:44]}"
                      for gid in terms_kept]

            h = max(4, 0.22 * len(terms_kept) + 1.5)
            w = max(7, 0.45 * len(clades) + 4)
            fig, ax = plt.subplots(figsize=(w, h))
            arr = mat.to_numpy(dtype=float)
            vmax = np.nanpercentile(arr, 95) if np.any(np.isfinite(arr)) else 1.0
            im = ax.imshow(arr, aspect="auto", cmap="magma_r",
                            vmin=0, vmax=max(vmax, math.log10(0.05) * -1))
            ax.set_xticks(range(len(clades)))
            ax.set_xticklabels(clades, rotation=45, ha="right", fontsize=8)
            ax.set_yticks(range(len(terms_kept)))
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_title(f"Cross-clade GO enrichment — {ns}  "
                          f"(cell = -log10 best q; terms in \u2265{min_clades} clades)",
                          fontsize=10)
            cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
            cbar.set_label("-log10(q)", fontsize=8)
            # Grid + annotate "significant at q<=0.05" (mlog10q >= 1.301).
            for i in range(arr.shape[0]):
                for j in range(arr.shape[1]):
                    v = arr[i, j]
                    if np.isfinite(v) and v >= 1.301:
                        ax.text(j, i, "*", ha="center", va="center",
                                fontsize=7, color="white")
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"[write] {out_path}")


# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--significant-terms", required=True)
    ap.add_argument("--obo", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--top-n", type=int, default=15,
                    help="top-N terms per (clade, namespace) in dotplots")
    ap.add_argument("--heatmap-min-clades", type=int, default=2,
                    help="terms must be enriched in at least this many clades "
                         "to appear on the heatmap (default: 2)")
    ap.add_argument("--heatmap-max-terms", type=int, default=40,
                    help="cap heatmap rows per namespace (default: 40)")
    ap.add_argument("--min-fold", type=float, default=3.0,
                    help="drop cells whose fold-enrichment is below this "
                         "BEFORE per-term dedupe, so large-N artifacts "
                         "don't mask real signal (default: 3.0)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[load] OBO …")
    term_info = parse_obo(args.obo)
    print(f"  terms_in_obo={len(term_info)}")

    print("[load] significant_terms …")
    sig = pd.read_csv(args.significant_terms, sep="\t")
    print(f"  rows={len(sig)}")

    # Annotate with name + authoritative namespace from OBO.
    def name_of(g):
        v = term_info.get(g)
        return v[0] if v else g
    def ns_of(g):
        v = term_info.get(g)
        return v[1] if v else "?"
    sig["go_name"] = sig["go_id"].map(name_of)
    # Overwrite go_namespace from OBO (authoritative) where available.
    sig["go_namespace"] = sig["go_id"].map(ns_of).where(
        sig["go_id"].map(lambda g: g in term_info),
        sig["go_namespace"])

    annotated_path = out_dir / "significant_terms_annotated.tsv"
    sig.to_csv(annotated_path, sep="\t", index=False)
    print(f"[write] {annotated_path}  rows={len(sig)}")

    make_dotplots(sig, out_dir / "dotplots.pdf", top_n=args.top_n,
                   min_fold=args.min_fold)
    make_heatmap(sig, out_dir / "heatmap.pdf",
                  min_clades=args.heatmap_min_clades,
                  max_terms_per_ns=args.heatmap_max_terms,
                  min_fold=args.min_fold)


if __name__ == "__main__":
    main()
