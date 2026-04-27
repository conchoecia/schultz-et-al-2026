"""Stage 3 v3: per-example multi-species locus tracks (Genomicus-style)
with proper rectangular cladogram subplot, anchor buffer floors, true
absolute-bp scale, and per-pair top-3 enriched TF tracks.

Key v3 changes from v2:
  - Tree subplot uses an external-y-coord rectangular cladogram that
    matches lane centers exactly (no more half-lane offset).
  - Anchors never touch the lane edge (4% buffer floor in normalized mode,
    larger pad floors in absolute mode).
  - Absolute mode uses one figure-wide bp/inch ratio (left-anchored on
    anchor A, with right-context fill out to max_window_bp).
  - ReMap density track removed; replaced with up to 3 per-TF peak tick
    tracks (top-N enriched per `tf_enrichment_top_tfs.tsv` from Stage 1.7).
  - Page header moved to fig.suptitle so it doesn't get clipped.

If GENOMES_DIR is unset write a placeholder PDF and exit 0.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# TrueType fonts (Type 42) so PDF text remains editable in Illustrator/Inkscape.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.path import Path as MplPath

import loci_io

# Figure width = 184 mm = 7.244 in (manuscript single-column wide spec).
FIG_WIDTH_IN = 7.244

ANCHOR_A_COLOR = "#1f77b4"   # ortholog1
ANCHOR_B_COLOR = "#d62728"   # ortholog2
INTERVENING_COLOR = "#bbbbbb"
RIBBON_ALPHA = 0.18
TF_TRACK_COLORS = ["#2ca02c", "#9467bd", "#ff7f0e"]   # top-1, top-2, top-3

FOCAL_SPECIES = {
    "Vertebrata": "Homosapiens-9606-GCF000001405.40",
    "Arthropoda": "Drosophilamelanogaster-7227-GCF000001215.4",
}
FOCAL_BUILD = {
    "Homosapiens-9606-GCF000001405.40": "hg38",
    "Drosophilamelanogaster-7227-GCF000001215.4": "dm6",
}

HG38_REFSEQ_TO_UCSC = {f"NC_{i:06d}.{ver}": f"chr{name}" for (i, name, ver) in [
    (1, 1, 11), (2, 2, 12), (3, 3, 12), (4, 4, 12), (5, 5, 10), (6, 6, 12),
    (7, 7, 14), (8, 8, 11), (9, 9, 12), (10, 10, 11), (11, 11, 10),
    (12, 12, 12), (13, 13, 11), (14, 14, 9), (15, 15, 10), (16, 16, 10),
    (17, 17, 11), (18, 18, 10), (19, 19, 10), (20, 20, 11), (21, 21, 9),
    (22, 22, 11), (23, "X", 11), (24, "Y", 10),
]}
HG38_REFSEQ_TO_UCSC["NC_012920.1"] = "chrM"
DM6_REFSEQ_TO_UCSC = {
    "NT_033779.5": "chr2L", "NT_033778.4": "chr2R",
    "NT_037436.4": "chr3L", "NT_033777.3": "chr3R",
    "NC_004353.4": "chr4",  "NC_004354.4": "chrX",
    "NC_024512.1": "chrY",  "NC_024511.2": "chrM",
}
BUILD_TRANSLATE = {"hg38": HG38_REFSEQ_TO_UCSC, "dm6": DM6_REFSEQ_TO_UCSC}

# Lane geometry — data-unit half-heights. Halved from v6 so anchor blocks
# don't dominate the lane vertically (user feedback: "gene blocks too tall").
GENE_HALFHEIGHT = 0.16           # anchor / intervening gene polygon half-height
INTERVENING_HALFHEIGHT = 0.10
TF_TRACK_HALFHEIGHT = 0.06       # ~1/3 of gene full height (0.32)
TF_TRACK_VERT_SPACING = 0.18
LANE_BOTTOM_PAD = 0.5
LANE_TOP_PAD = 0.5

# ---------------------------------------------------------------------------
# Misc helpers


def write_placeholder(out_path: str, msg: str):
    with PdfPages(out_path) as pdf:
        fig = plt.figure(figsize=(8.5, 4))
        fig.text(0.05, 0.5, msg, fontsize=10, va="center", wrap=True)
        pdf.savefig(fig); plt.close(fig)
    print(f"placeholder -> {out_path}")


def species_key_to_taxid(sk: str) -> str | None:
    parts = sk.split("-")
    return parts[1] if len(parts) >= 2 and parts[1].isdigit() else None


# ---------------------------------------------------------------------------
# .chrom.gz isoform collapse


def collapse_isoforms(slc: pd.DataFrame, max_gene_bp: int = 200_000) -> pd.DataFrame:
    """Collapse same-strand overlapping rows into one block; drop predicted
    readthroughs > max_gene_bp."""
    if slc.empty:
        return slc
    s = slc.sort_values(["scaffold", "strand", "start"]).copy()
    s["start"] = s["start"].astype(int)
    s["end"] = s["end"].astype(int)
    s = s[(s["end"] - s["start"]) <= max_gene_bp]
    if s.empty:
        return s
    out_rows = []
    cur = None
    for _, r in s.iterrows():
        if cur is None:
            cur = r.to_dict(); cur["members"] = [r["gene_id"]]
            continue
        same = (r["scaffold"] == cur["scaffold"] and r["strand"] == cur["strand"])
        if same and r["start"] <= cur["end"] + 500:
            cur["end"] = max(cur["end"], int(r["end"]))
            cur["members"].append(r["gene_id"])
        else:
            out_rows.append(cur); cur = r.to_dict(); cur["members"] = [r["gene_id"]]
    if cur is not None:
        out_rows.append(cur)
    df = pd.DataFrame(out_rows)
    if not df.empty:
        df["gene_id"] = df["members"].apply(lambda m: m[0])
    return df


# ---------------------------------------------------------------------------
# BED loaders


def load_bed(path: Path, keep_chroms: set[str] | None = None,
             with_label: bool = False,
             include_extra_cols: list[int] | None = None) -> pd.DataFrame | None:
    """Read a BED. `with_label` parses col 4 as TF:cell. `include_extra_cols`
    pulls extra named columns by 0-based index (e.g. cCRE category in col 9
    of ENCODE V4 cCRE BED)."""
    if not path.exists() or path.stat().st_size == 0:
        return None
    if with_label:
        cols = [0, 1, 2, 3]
        names = ["chrom", "start", "end", "label"]
    else:
        cols = [0, 1, 2]
        names = ["chrom", "start", "end"]
    dtype = {"chrom": str, "start": "int64", "end": "int64"}
    if with_label:
        dtype["label"] = str
    if include_extra_cols:
        for i in include_extra_cols:
            cols.append(i)
            names.append(f"col{i}")
            dtype[f"col{i}"] = str
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     usecols=cols, names=names, dtype=dtype,
                     compression="infer", low_memory=False,
                     on_bad_lines="skip")
    if keep_chroms is not None:
        df = df[df["chrom"].isin(keep_chroms)].copy()
    if with_label:
        df["tf"] = df["label"].str.split(":").str[0]
    return df


def bed_in_window(bed: pd.DataFrame | None, chrom: str,
                   lo: int, hi: int) -> pd.DataFrame:
    if bed is None:
        return pd.DataFrame()
    return bed[(bed["chrom"] == chrom) & (bed["end"] >= lo) & (bed["start"] <= hi)]


# ---------------------------------------------------------------------------
# Drawing primitives


def draw_gene(ax, x0, x1, y, color, strand, lw=0.4, label=None,
              height=GENE_HALFHEIGHT, notch_frac=0.25, notch_max=0.005):
    if x1 < x0:
        x0, x1 = x1, x0
    width = x1 - x0
    notch = min(width * notch_frac, notch_max)
    if strand == "+":
        verts = [(x0, y - height), (x1 - notch, y - height),
                 (x1, y), (x1 - notch, y + height),
                 (x0, y + height), (x0, y - height)]
    elif strand == "-":
        verts = [(x1, y - height), (x0 + notch, y - height),
                 (x0, y), (x0 + notch, y + height),
                 (x1, y + height), (x1, y - height)]
    else:
        verts = [(x0, y - height), (x1, y - height),
                 (x1, y + height), (x0, y + height), (x0, y - height)]
    poly = mpatches.Polygon(verts, closed=True, facecolor=color,
                             edgecolor="black", linewidth=lw)
    ax.add_patch(poly)
    if label:
        ax.text((x0 + x1) / 2, y + height + 0.06, label, ha="center",
                va="bottom", fontsize=6, color=color)


def draw_ribbon(ax, lane_xy_top, lane_xy_bot, color):
    (x0a, x1a, y_top) = lane_xy_top
    (x0b, x1b, y_bot) = lane_xy_bot
    cy = (y_top + y_bot) / 2
    verts = [(x0a, y_top), (x0a, cy), (x0b, cy), (x0b, y_bot),
             (x1b, y_bot), (x1b, cy), (x1a, cy), (x1a, y_top), (x0a, y_top)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    p = mpatches.PathPatch(MplPath(verts, codes), facecolor=color,
                            edgecolor="none", alpha=RIBBON_ALPHA)
    ax.add_patch(p)


# ---------------------------------------------------------------------------
# Tree subplot — rectangular cladogram, externally fixed tip y-coords
# Mirrors the pattern from EGT/src/egt/newick_to_common_ancestors.py:plot_tree


def load_clade_tree(nwk_path: str | None):
    if not nwk_path or not Path(nwk_path).exists():
        return None
    try:
        from ete4 import Tree
        return Tree(open(nwk_path).read(), parser=1)
    except Exception as exc:  # noqa: BLE001
        print(f"WARN: tree load failed: {exc}")
        return None


def prune_tree_for_species(tree, taxids_in_order: list[str]):
    """Return a pruned copy of `tree` keeping only leaves whose [taxid] is in
    `taxids_in_order`, with children re-ordered so leaves appear in the same
    order as `taxids_in_order` (topmost lane = first in the list).
    """
    if tree is None:
        return None
    target = set(taxids_in_order)
    keep = []
    for L in tree.leaves():
        nm = getattr(L, "name", "") or ""
        if "[" in nm and nm.split("[")[1].rstrip("]") in target:
            keep.append(L)
    if not keep:
        return None
    pruned = tree.copy()
    keep_names = {L.name for L in keep}
    pruned.prune([L for L in pruned.leaves() if L.name in keep_names])

    # Build target rank of each taxid (top = 0).
    tx_rank = {tx: i for i, tx in enumerate(taxids_in_order)}

    def leaf_taxid(L):
        nm = getattr(L, "name", "") or ""
        return nm.split("[")[1].rstrip("]") if "[" in nm else None

    def min_rank(node):
        """Minimum rank (highest priority = lowest int) among descendants."""
        if node.is_leaf:
            tx = leaf_taxid(node)
            return tx_rank.get(tx, 1_000_000)
        return min(min_rank(c) for c in node.children)

    def reorder(node):
        if node.is_leaf:
            return
        node.children = sorted(node.children, key=min_rank)
        for c in node.children:
            reorder(c)

    reorder(pruned)
    return pruned


def draw_taxonomy(ax, pruned_tree, species_keys_in_order: list[str]):
    """Rectangular cladogram on `ax`. Tip y-coords match lane centers
    (`(n - 1 - i) + 0.5` for the i-th species in the top-to-bottom order)."""
    if pruned_tree is None:
        ax.set_axis_off(); return
    n = len(species_keys_in_order)
    sk_to_y = {sk: (n - 1 - i) + 0.5 for i, sk in enumerate(species_keys_in_order)}
    tx_to_sk = {species_key_to_taxid(sk): sk for sk in species_keys_in_order
                if species_key_to_taxid(sk)}

    # Compute depth from root (integer cladogram depth).
    depths: dict = {}
    def walk_depth(node, d=0):
        depths[node] = d
        for c in node.children:
            walk_depth(c, d + 1)
    walk_depth(pruned_tree)
    max_d = max(depths.values()) if depths else 1

    # Assign leaf y from species_keys_in_order; internal y = mean of descendants.
    leaf_y: dict = {}
    for L in pruned_tree.leaves():
        nm = getattr(L, "name", "") or ""
        if "[" not in nm:
            continue
        tx = nm.split("[")[1].rstrip("]")
        sk = tx_to_sk.get(tx)
        if sk and sk in sk_to_y:
            leaf_y[L] = sk_to_y[sk]
    if not leaf_y:
        ax.set_axis_off(); return

    node_y: dict = {}
    def walk_y(node):
        if node in leaf_y:
            node_y[node] = leaf_y[node]
            return node_y[node]
        # Internal node y = midpoint of immediate children's y-positions
        # (NOT the centroid of all descendant leaves — that biases vertical
        # branches toward the larger subtree on unbalanced trees).
        ys = [walk_y(c) for c in node.children]
        ys = [y for y in ys if y is not None]
        node_y[node] = (min(ys) + max(ys)) / 2 if ys else None
        return node_y[node]
    walk_y(pruned_tree)

    # x_node = depth_from_root → root at x=0 (LEFT), tips at x=max_d (RIGHT).
    # Branches coalesce toward x=0 (left). Tip terminal lines + labels live
    # to the right of the tips, between max_d and the gutter.
    def node_x(node):
        return depths[node]

    tip_x = max_d + 0.3   # all leaf terminal lines end at this fixed x
    for n_node in pruned_tree.traverse():
        if n_node.is_leaf:
            # tip terminal: from leaf's own x to the right gutter
            ax.plot([node_x(n_node), tip_x],
                    [node_y[n_node], node_y[n_node]],
                    color="0.35", lw=0.7)
            nm = (n_node.name or "").split("[")[0].rstrip("_")
            parts = nm.split("_")
            # Full species name with abbreviated genus initial (no truncation
            # of species epithet — was cutting "troglodytes" → "troglody").
            short = (f"{parts[0][:1]}. {parts[1]}") if len(parts) >= 2 else nm
            ax.text(tip_x + 0.05, node_y[n_node], short,
                    ha="left", va="center", fontsize=6, color="0.25",
                    style="italic")
        else:
            children = list(n_node.children)
            child_ys = [node_y[c] for c in children if node_y.get(c) is not None]
            if not child_ys:
                continue
            x = node_x(n_node)
            # vertical spine at this internal node
            ax.plot([x, x], [min(child_ys), max(child_ys)], color="0.35", lw=0.7)
            # horizontal branch from this node OUT to each child (x increases)
            for c in children:
                yc = node_y.get(c)
                if yc is None:
                    continue
                ax.plot([x, node_x(c)], [yc, yc], color="0.35", lw=0.7)

    # Branches occupy [0, max_d]; reserve large gutter to right for tip
    # labels so the tree subplot can be physically narrow without label clip.
    ax.set_xlim(-0.3, max_d + max(5.0, max_d * 1.5))
    ax.set_axis_off()


# ---------------------------------------------------------------------------
# Geometry: per-species window (with anchor-buffer floor)


def species_window(sr, scale: str, max_window_bp: int,
                    pad_factor: float, min_pad_bp: int,
                    abs_window_bp: int | None = None) -> dict:
    """Return per-species window parameters.

    For absolute scale, all species share the same `abs_window_bp` width
    (anchor-A is left-anchored at the same physical x; rest is right-context).
    For normalized (or absolute over-budget pairs), use a padded window
    centered on the anchor pair.
    """
    a_lo = min(sr["start1"], sr["start2"])
    a_hi = max(sr["end1"], sr["end2"])
    span = max(1, a_hi - a_lo)
    is_wide = span > max_window_bp

    if scale == "absolute" and abs_window_bp is not None and not is_wide:
        flank_left = max(min_pad_bp, int(0.05 * abs_window_bp))
        win_lo = max(0, a_lo - flank_left)
        win_hi = win_lo + abs_window_bp
        return dict(win_lo=win_lo, win_hi=win_hi, win_w=win_hi - win_lo,
                     wide_pair=False)

    if is_wide:
        # Wide-pair: cap each anchor flank at 2 Mb so we don't try to draw
        # negative coords (Mb-scale spans on real genes — but Mb-scale spans
        # on a single accession are usually NCBI readthrough artifacts).
        anchor_w_each = min(2_000_000,
                             max(2 * min_pad_bp,
                                  int(0.10 * (a_hi - a_lo))))
        win_lo = max(0, a_lo - anchor_w_each)
        win_hi = a_hi + anchor_w_each
        return dict(win_lo=win_lo, win_hi=win_hi, win_w=win_hi - win_lo,
                     wide_pair=True)

    pad = int(max(pad_factor * span, min_pad_bp))
    win_lo = max(0, a_lo - pad)
    win_hi = a_hi + pad
    return dict(win_lo=win_lo, win_hi=win_hi, win_w=win_hi - win_lo,
                 wide_pair=False)


def enforce_anchor_buffer(geom, sr, min_visual_buffer: float = 0.04):
    """Ensure anchor xL/xR (in normalized [0,1] lane coords) are within
    [min_visual_buffer, 1 - min_visual_buffer]. If not, widen the bp window
    symmetrically until they fit. Does nothing for wide_pair (already wide)."""
    if geom["wide_pair"]:
        return geom
    win_lo, win_hi, win_w = geom["win_lo"], geom["win_hi"], geom["win_w"]
    a_lo = min(sr["start1"], sr["start2"])
    a_hi = max(sr["end1"], sr["end2"])
    xL = (a_lo - win_lo) / win_w
    xR = (a_hi - win_lo) / win_w
    if xL >= min_visual_buffer and xR <= 1 - min_visual_buffer:
        return geom
    # Required buffer in bp on each side
    need_left_bp = max(0, (min_visual_buffer * win_w) - (a_lo - win_lo))
    need_right_bp = max(0, (min_visual_buffer * win_w) - (win_hi - a_hi))
    new_lo = win_lo - int(need_left_bp + 1)
    new_hi = win_hi + int(need_right_bp + 1)
    return dict(win_lo=new_lo, win_hi=new_hi, win_w=new_hi - new_lo,
                 wide_pair=False)


# ---------------------------------------------------------------------------
# Per-page renderer


def plot_one_pair(pair_row, dist_rows, recompute_row, genome_cfg,
                   bed_cache: dict[str, pd.DataFrame | None],
                   tree=None, scale: str = "normalized",
                   tf_tracks: dict | None = None,
                   pair_enrichment_row=None,
                   pad_factor: float = 0.7, min_pad_bp: int = 15_000,
                   max_window_bp: int = 600_000):
    """Render one example pair as a multi-species locus figure.

    `tf_tracks` is a per-(pair_key, build) dict of top-N TF rows; each row
    has `tf`, `p_t`, `surprise`, `rank`. Drawn as thin tick lanes under the
    focal-species gene track.
    """
    species_rows = [r for r in dist_rows if r["same_scaffold"]]
    if not species_rows:
        fig = plt.figure(figsize=(11, 3))
        fig.text(0.05, 0.5, f"{pair_row['nodename']} | {pair_row['pair_label']}: "
                            "no same-scaffold species available for plotting.",
                  fontsize=10)
        return fig

    # Re-order species_rows to follow the tree's pruned-leaf order so
    # phylogenetically-adjacent species sit in adjacent lanes (e.g. primate
    # sisters Hsap+Cjac next to each other, not separated by rodents).
    if tree is not None:
        present_taxids = [species_key_to_taxid(r["species_key"])
                           for r in species_rows]
        present_taxids = [t for t in present_taxids if t]
        pruned_for_order = prune_tree_for_species(tree, present_taxids)
        if pruned_for_order is not None:
            tree_taxid_order = []
            for L in pruned_for_order.leaves():
                nm = getattr(L, "name", "") or ""
                if "[" in nm:
                    tree_taxid_order.append(nm.split("[")[1].rstrip("]"))
            order_idx = {tx: i for i, tx in enumerate(tree_taxid_order)}
            species_rows = sorted(species_rows,
                                    key=lambda r: order_idx.get(
                                        species_key_to_taxid(r["species_key"]),
                                        len(order_idx)))

    n_species = len(species_rows)
    pair_key = pair_row.get("pair_key") if isinstance(pair_row, dict) \
        else pair_row["pair_key"]

    # Number of TF lanes per focal species (impacts ylim)
    n_tf_per_focal = 0
    if tf_tracks:
        for (pk, _build), rows in tf_tracks.items():
            if pk == pair_key:
                n_tf_per_focal = max(n_tf_per_focal, min(3, len(rows)))

    # Compute per-species windows up front (needed for absolute-scale max_w).
    geoms = {}
    for sr in species_rows:
        g = species_window(sr, scale="normalized",
                            max_window_bp=max_window_bp,
                            pad_factor=pad_factor, min_pad_bp=min_pad_bp)
        g = enforce_anchor_buffer(g, sr)
        geoms[sr["species_key"]] = g
    if scale == "absolute":
        # Reference window = max non-wide window across species (or fallback
        # to 4 × min_pad_bp if all species are wide_pair).
        non_wide = [g["win_w"] for g in geoms.values() if not g["wide_pair"]]
        abs_window_bp = max(non_wide) if non_wide else 4 * min_pad_bp
        geoms = {}
        for sr in species_rows:
            g = species_window(sr, scale="absolute",
                                max_window_bp=max_window_bp,
                                pad_factor=pad_factor, min_pad_bp=min_pad_bp,
                                abs_window_bp=abs_window_bp)
            g = enforce_anchor_buffer(g, sr)
            geoms[sr["species_key"]] = g

    # Figure: 184 mm wide. Per-lane height stays compressed; figure
    # grows when more TF/regulatory tracks are added (each TF lane gets
    # ~0.20 in of extra height).
    per_lane_in = 0.32
    fig_h = 0.80 + per_lane_in * n_species + 0.20 * n_tf_per_focal
    fig = plt.figure(figsize=(FIG_WIDTH_IN, fig_h))
    if tree is not None:
        # Tree subplot tighter (1.5 / 11.5 → ~12% width). Inside the tree
        # axis, branches now occupy only the LEFT portion (xlim extends
        # well past max_d), leaving the rest as a label gutter that fits
        # full species names without crowding the locus axis.
        gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 11.5], wspace=0.02,
                               left=0.02, right=0.99,
                               top=1 - 0.40 / fig_h, bottom=0.30 / fig_h)
        ax_tree = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1], sharey=ax_tree)
    else:
        ax_tree = None
        ax = fig.add_subplot(1, 1, 1,
                              position=(0.06, 0.30 / fig_h, 0.92,
                                        1 - 0.70 / fig_h))
    ax.set_axis_off()

    lane_y = {sr["species_key"]: n_species - i - 1
               for i, sr in enumerate(species_rows)}

    anchor_footprints_A: dict[str, tuple] = {}
    anchor_footprints_B: dict[str, tuple] = {}

    # x-norm helper. If `flip` is True, the lane is mirrored horizontally
    # so anchor A always lands on the left visually — prevents anchor
    # ribbons from crossing when the genome is flipped/inverted in this
    # species. Mirrored lanes get a "(−)" prefix on the scaffold label.
    def make_normalizer(geom, flip=False):
        win_lo, win_w = geom["win_lo"], geom["win_w"]
        if flip:
            return lambda v: 1.0 - (v - win_lo) / win_w
        return lambda v: (v - win_lo) / win_w

    def flip_strand_char(s):
        return "+" if s == "-" else ("-" if s == "+" else s)

    for sr in species_rows:
        sk = sr["species_key"]
        y = lane_y[sk] + 0.5
        geom = geoms[sk]
        win_lo, win_hi, win_w = geom["win_lo"], geom["win_hi"], geom["win_w"]
        wide_pair = geom["wide_pair"]

        # Flip lane if anchor A is naturally right of anchor B — keeps
        # anchor A on the left across all lanes so ribbons never cross.
        needs_flip = sr["start1"] > sr["start2"]
        norm = make_normalizer(geom, flip=needs_flip)

        ent = genome_cfg["species"].get(sk)
        chrom_path = ent["chrom"] if ent else None
        if not chrom_path or not Path(chrom_path).exists():
            continue
        cdf = loci_io.read_chrom(chrom_path)
        if wide_pair:
            half = (win_hi - win_lo) // 2
            slc_left = cdf[(cdf["scaffold"] == sr["scaffold1"]) &
                            (cdf["end"] >= win_lo) &
                            (cdf["start"] <= win_lo + half)]
            slc_right = cdf[(cdf["scaffold"] == sr["scaffold1"]) &
                             (cdf["end"] >= win_hi - half) &
                             (cdf["start"] <= win_hi)]
            slc = pd.concat([slc_left, slc_right]).copy()
            ax.plot([0.495, 0.505], [y - 0.05, y + 0.05], color="0.4", lw=0.8)
            ax.plot([0.5, 0.51], [y - 0.05, y + 0.05], color="0.4", lw=0.8)
            # Per-species TSS-to-TSS distance (matches the metric in the
            # page title's median TSS-TSS so the numbers are comparable).
            tss_dist_bp = abs(int(sr["tss1"]) - int(sr["tss2"]))
            gap_label = (f"// {tss_dist_bp / 1e6:.1f} Mb TSS-TSS"
                         if tss_dist_bp >= 1e6
                         else f"// {tss_dist_bp / 1e3:.0f} kb TSS-TSS"
                         if tss_dist_bp >= 1e3
                         else f"// {tss_dist_bp:,} bp TSS-TSS")
            ax.text(0.5, y - 0.5, gap_label, ha="center", va="top",
                    fontsize=6, color="0.4", style="italic")
        else:
            slc = cdf[(cdf["scaffold"] == sr["scaffold1"]) &
                      (cdf["end"] >= win_lo) &
                      (cdf["start"] <= win_hi)].copy()
        slc = collapse_isoforms(slc)

        def has_member(row, gid):
            m = row.get("members") if isinstance(row, dict) else row["members"]
            return gid in (m or [])

        # When flipped, draw_gene's xL/xR will be reversed (xL > xR) — its
        # built-in swap handles that; we just have to flip the strand char
        # so the arrow direction stays consistent with the gene's biological
        # 5'→3' polarity in the mirrored frame.
        strand_xform = flip_strand_char if needs_flip else (lambda s: s)
        for _, g in slc.iterrows():
            if has_member(g, sr["gene1"]) or has_member(g, sr["gene2"]):
                continue
            xL = norm(max(g["start"], win_lo))
            xR = norm(min(g["end"], win_hi))
            draw_gene(ax, xL, xR, y, INTERVENING_COLOR, strand_xform(g["strand"]),
                      lw=0.2, height=INTERVENING_HALFHEIGHT)
        a1 = slc[slc["members"].apply(lambda m: sr["gene1"] in (m or []))]
        a2 = slc[slc["members"].apply(lambda m: sr["gene2"] in (m or []))]
        if not a1.empty:
            r = a1.iloc[0]
            xL, xR = norm(r["start"]), norm(r["end"])
            draw_gene(ax, xL, xR, y, ANCHOR_A_COLOR, strand_xform(r["strand"]),
                      lw=0.6, height=GENE_HALFHEIGHT)
            # Store footprint with min/max so ribbon code sees a left→right pair
            anchor_footprints_A[sk] = (min(xL, xR), max(xL, xR), y)
        if not a2.empty:
            r = a2.iloc[0]
            xL, xR = norm(r["start"]), norm(r["end"])
            draw_gene(ax, xL, xR, y, ANCHOR_B_COLOR, strand_xform(r["strand"]),
                      lw=0.6, height=GENE_HALFHEIGHT)
            anchor_footprints_B[sk] = (min(xL, xR), max(xL, xR), y)

        # Species label removed (redundant with the italic tree tip on the
        # left). Chrom + window label moved closer to the gene track.
        win_kb = win_w / 1e3
        # Genome-browser-pasteable region: chrom:start-end with integer
        # coords. Prefix "(−) " when this lane was mirrored so the user
        # can tell the displayed orientation is flipped vs the genome.
        flip_tag = "(−) " if needs_flip else ""
        ax.text(0, y - 0.22,
                f"{flip_tag}{sr['scaffold1']}:{int(win_lo):,}-{int(win_hi):,}  "
                f"({win_kb:.0f} kb shown)",
                ha="left", va="top", fontsize=6, color="0.4")
        bar_bp = 10_000 if win_w >= 100_000 else max(1_000, win_w // 10)
        bar_w = bar_bp / win_w
        bar_label = (f"{bar_bp//1000} kb" if bar_bp >= 1000
                     else f"{bar_bp} bp")
        ax.plot([1.0 - bar_w, 1.0], [y - 0.18, y - 0.18], color="0.3", lw=1.2)
        ax.text(1.0 - bar_w / 2, y - 0.22, bar_label, ha="center", va="top",
                fontsize=6, color="0.3")

        # cCRE + phastCons + TF tracks render ABOVE the focal species lane
        # (so they don't crowd the species below). Order: phastCons signal
        # → cCRE density → top-3 TF tracks.
        build = FOCAL_BUILD.get(sk)
        if build and bed_cache.get(build):
            ucsc_chrom = BUILD_TRANSLATE.get(build, {}).get(sr["scaffold1"])
            ccre = bed_cache[build].get("cCRE")
            phastcons = bed_cache[build].get("phastCons_bw")  # pyBigWig handle
            offset = +0.55
            # phyloP signed conservation track (UCSC 100 vertebrates / 27 insects).
            # Positive = conserved, negative = accelerated. Drawn as a
            # signed bar chart: positive bars go up (green), negative bars
            # go down (orange).
            phylop = bed_cache[build].get("phyloP_bw")
            if ucsc_chrom and phylop is not None:
                try:
                    chrom_size = phylop.chroms().get(ucsc_chrom)
                    q_lo = max(0, int(win_lo))
                    q_hi = min(chrom_size if chrom_size else int(win_hi),
                                int(win_hi))
                    if q_hi - q_lo >= 200:
                        n_bins = 480
                        means = phylop.stats(ucsc_chrom, q_lo, q_hi,
                                                type="mean", nBins=n_bins)
                        means = [(m if m is not None else 0.0) for m in means]
                        # phyloP scale: typically -3..+8 for 100way. Clip to
                        # ±5 for visualization, map to ±0.10 lane height.
                        # Use norm() so the bins respect lane flip.
                        bx_a, bx_b = norm(q_lo), norm(q_hi)
                        bin_x_lo = min(bx_a, bx_b)
                        bin_x_w = abs(bx_b - bx_a)
                        if needs_flip:
                            means = list(reversed(means))
                        baseline = y + offset
                        for i, m in enumerate(means):
                            v = max(-5, min(5, m))
                            h = (v / 5.0) * 0.12
                            if abs(h) < 1e-3: continue
                            color = "#0c5d0e" if h > 0 else "#d97706"
                            ax.add_patch(mpatches.Rectangle(
                                (bin_x_lo + i * (bin_x_w / n_bins),
                                 baseline if h >= 0 else baseline + h),
                                bin_x_w / n_bins, abs(h),
                                facecolor=color, edgecolor="none", alpha=0.9))
                        # zero baseline line
                        ax.plot([bin_x_lo, bin_x_lo + bin_x_w],
                                 [baseline, baseline],
                                 color="0.5", lw=0.3)
                        ax.text(-0.005, baseline, "phyloP",
                                ha="right", va="center", fontsize=6,
                                color="#0c5d0e")
                        offset += TF_TRACK_VERT_SPACING + 0.04
                except Exception as exc:  # noqa: BLE001
                    print(f"WARN phyloP stats failed for {ucsc_chrom} "
                          f"({win_lo}-{win_hi}): {exc}")
            # phastCons signal track (continuous, 0..1; UCSC multiz alignment).
            if ucsc_chrom and phastcons is not None:
                try:
                    chrom_size = phastcons.chroms().get(ucsc_chrom)
                    q_lo = max(0, int(win_lo))
                    q_hi = min(chrom_size if chrom_size else int(win_hi), int(win_hi))
                    if q_hi - q_lo < 200:
                        raise ValueError("window too narrow / off chrom end")
                    n_bins = 480   # finer resolution per user request
                    means = phastcons.stats(ucsc_chrom, q_lo, q_hi,
                                              type="mean", nBins=n_bins)
                    means = [(m if m is not None else 0.0) for m in means]
                    bar_w = 1.0 / n_bins
                    # Flip-aware bin geometry
                    bx_a, bx_b = norm(q_lo), norm(q_hi)
                    bin_x_lo = min(bx_a, bx_b)
                    bin_x_w = abs(bx_b - bx_a)
                    if needs_flip:
                        means = list(reversed(means))
                    for i, m in enumerate(means):
                        h = max(0.0, min(1.0, m)) * 0.13
                        if h <= 0: continue
                        ax.add_patch(mpatches.Rectangle(
                            (bin_x_lo + i * (bin_x_w / n_bins),
                             y + offset - 0.06),
                            bin_x_w / n_bins, h, facecolor="#3d8b37",
                            edgecolor="none", alpha=0.9))
                    ax.text(-0.005, y + offset, "phastCons",
                            ha="right", va="center", fontsize=6, color="#3d8b37")
                    offset += TF_TRACK_VERT_SPACING
                except Exception as exc:  # noqa: BLE001
                    print(f"WARN phastCons stats failed for {ucsc_chrom} "
                          f"({win_lo}-{win_hi}): {exc}")
            if ucsc_chrom and ccre is not None:
                hits = bed_in_window(ccre, ucsc_chrom, win_lo, win_hi)
                if not hits.empty:
                    n_bins = 120
                    density = np.zeros(n_bins)
                    starts = hits["start"].clip(lower=win_lo).values
                    ends = hits["end"].clip(upper=win_hi).values
                    # bin index is bp-based and stays the same regardless
                    # of flip; we apply the flip when drawing.
                    mids = ((starts + ends) // 2 - win_lo) / win_w
                    bin_ix = np.clip((mids * n_bins).astype(int), 0, n_bins - 1)
                    np.add.at(density, bin_ix, 1)
                    if density.max() > 0:
                        norm_d = density / density.max() * 0.12
                        bw = 1.0 / n_bins
                        if needs_flip:
                            norm_d = norm_d[::-1]
                        for i, h in enumerate(norm_d):
                            if h <= 0: continue
                            ax.add_patch(mpatches.Rectangle(
                                (i * bw, y + offset - h / 2),
                                bw, h, facecolor="#888888", edgecolor="none",
                                alpha=0.85))
                # cCRE = ENCODE Registry of candidate cis-Regulatory
                # Elements (PLS/pELS/dELS/CA-only); functional regulatory
                # annotation in human cell lines, NOT sequence conservation.
                ax.text(-0.005, y + offset, "ENCODE cCRE",
                        ha="right", va="center", fontsize=6, color="0.3")
                offset += TF_TRACK_VERT_SPACING

            if tf_tracks:
                key = (pair_key, build)
                top_rows = list(tf_tracks.get(key, []))[:3]
                remap_full = bed_cache[build].get("ReMap_full")
                for ti, tf_row in enumerate(top_rows):
                    color = TF_TRACK_COLORS[ti % len(TF_TRACK_COLORS)]
                    track_y = y + offset
                    if ucsc_chrom and remap_full is not None:
                        peaks = remap_full[(remap_full["tf"] == tf_row["tf"]) &
                                            (remap_full["chrom"] == ucsc_chrom) &
                                            (remap_full["end"] >= win_lo) &
                                            (remap_full["start"] <= win_hi)]
                        for _, pk_row in peaks.iterrows():
                            mid = (max(pk_row["start"], win_lo) +
                                   min(pk_row["end"], win_hi)) / 2
                            x = norm(mid)
                            in_promoter = False
                            for tss, strand in ((sr["tss1"], sr["strand1"]),
                                                  (sr["tss2"], sr["strand2"])):
                                if strand == "+":
                                    plo, phi = tss - 2000, tss
                                else:
                                    plo, phi = tss, tss + 2000
                                if plo <= mid <= phi:
                                    in_promoter = True; break
                            tick_color = color if in_promoter else "0.55"
                            tick_lw = 1.4 if in_promoter else 0.5
                            ax.plot([x, x],
                                    [track_y - TF_TRACK_HALFHEIGHT,
                                     track_y + TF_TRACK_HALFHEIGHT],
                                    color=tick_color, lw=tick_lw,
                                    solid_capstyle="butt")
                    # Restore short inline label: rank + TF name. The full
                    # detail (% promoters bound) lives in the top-right legend.
                    rank = int(tf_row.get("rank", ti + 1))
                    label = f"#{rank} {tf_row['tf']}"
                    ax.text(-0.005, track_y, label, ha="right", va="center",
                            fontsize=6, color=color)
                    offset += TF_TRACK_VERT_SPACING

    # Cross-species ribbons (anchor A + anchor B).
    sk_order = [sr["species_key"] for sr in species_rows]
    for top_sk, bot_sk in zip(sk_order[:-1], sk_order[1:]):
        for fp_dict, color in [(anchor_footprints_A, ANCHOR_A_COLOR),
                                (anchor_footprints_B, ANCHOR_B_COLOR)]:
            if top_sk in fp_dict and bot_sk in fp_dict:
                t = fp_dict[top_sk]; b = fp_dict[bot_sk]
                draw_ribbon(ax, (t[0], t[1], t[2] - GENE_HALFHEIGHT),
                                 (b[0], b[1], b[2] + GENE_HALFHEIGHT), color)

    # Headroom: TF tracks render ABOVE the focal lane(s), so reserve top
    # space proportional to the TF lane count. Bottom only needs the basic
    # caption + scale-bar pad now.
    top_pad = LANE_TOP_PAD + n_tf_per_focal * TF_TRACK_VERT_SPACING + 0.20
    bottom_pad = LANE_BOTTOM_PAD
    # Negative xlim leaves room for left-margin track labels; reduced from
    # -0.22 since the bold species name is gone.
    ax.set_xlim(-0.18, 1.02)
    ax.set_ylim(-bottom_pad, n_species + top_pad)

    if ax_tree is not None:
        # Force tree axis to share the locus axis y-range so tips align.
        ax_tree.set_ylim(ax.get_ylim())
        species_keys_in_order = [sr["species_key"] for sr in species_rows]
        taxids_in_order = [species_key_to_taxid(sk)
                            for sk in species_keys_in_order]
        taxids_in_order = [t for t in taxids_in_order if t]
        pruned = prune_tree_for_species(tree, taxids_in_order)
        draw_taxonomy(ax_tree, pruned, species_keys_in_order)
        ax_tree.set_ylim(ax.get_ylim())

    # Page header (suptitle so it doesn't get clipped by ax ylim).
    # Pair label is "SYM1/SYM2" — colour-tag each half so the gene-block
    # colour key is built into the title.
    pair_label = pair_row["pair_label"]
    title = f"{pair_row['nodename']} | {pair_label}"
    if recompute_row is not None:
        tss = recompute_row["median_tss_to_tss_bp"]
        edge = recompute_row["median_edge_to_edge_bp"]
        ort = recompute_row["modal_orientation"]
        nss = int(recompute_row["n_species_same_scaffold"])
        title += (f"  ·  TSS-TSS {int(tss):,} bp  ·  edge {int(edge):,} bp  "
                  f"·  {ort}  ·  n_species={nss}")
    if scale == "absolute":
        title += "  ·  absolute-bp scale"
    sub = ""
    if pair_enrichment_row is not None:
        obs = int(pair_enrichment_row.get("n_shared", 0))
        exp = float(pair_enrichment_row.get("expected", 0))
        z = float(pair_enrichment_row.get("z", 0))
        q = float(pair_enrichment_row.get("q_indep", 1))
        sub = (f"TF co-binding at TSSes ±2 kb (focal sp.): obs={obs}, "
               f"expected={exp:.1f}, z={z:.2f}, q={q:.2g}")
    # Two-line header: bigger title on top, smaller TF subtitle below it.
    # Y positioned in figure-coords; spacing tuned so they don't overlap.
    # Place title pieces using matplotlib's Text bbox (renderer-aware) so
    # the colored anchor names tile correctly without manual width guesses.
    head_y = 1 - 0.10 / fig_h
    prefix = f"{pair_row['nodename']} | "
    parts = pair_label.split("/", 1)
    rest = title.split(pair_label, 1)[1] if pair_label in title else ""
    fig.canvas.draw()  # ensure renderer exists for bbox queries
    renderer = fig.canvas.get_renderer()

    def place(prev_text, s, color, bold=False):
        """Add `s` immediately to the right of prev_text. Returns new Text."""
        if prev_text is None:
            x = 0.02
        else:
            bb = prev_text.get_window_extent(renderer=renderer)
            x = fig.transFigure.inverted().transform((bb.x1, bb.y0))[0]
        return fig.text(x, head_y, s, fontsize=8, ha="left", va="top",
                          color=color,
                          fontweight=("bold" if bold else "normal"))

    if len(parts) == 2:
        t = place(None, prefix, "black")
        t = place(t, parts[0], ANCHOR_A_COLOR, bold=True)
        t = place(t, " / ", "black")
        t = place(t, parts[1], ANCHOR_B_COLOR, bold=True)
        t = place(t, rest, "black")
    else:
        fig.text(0.02, head_y, title, fontsize=8, ha="left", va="top")
    if sub:
        fig.text(0.02, 1 - 0.30 / fig_h, sub, fontsize=7, ha="left",
                 va="top", color="#1a5490")

    # Top-right floating legend: per-TF detail (% promoters bound). The
    # rank + TF name is already shown inline beside each track; this
    # legend adds the genome-wide occupancy info that's too long inline.
    if tf_tracks:
        page_focal_builds = {FOCAL_BUILD.get(sk)
                              for sk in (sr["species_key"] for sr in species_rows)
                              if FOCAL_BUILD.get(sk)}
        legend_rows = []
        for build in page_focal_builds:
            for r in tf_tracks.get((pair_key, build), [])[:3]:
                p_t = r.get("p_t", float("nan"))
                legend_rows.append((r["tf"], int(r.get("rank", 1)),
                                     p_t, build))
        if legend_rows:
            x_legend = 0.99
            y_legend = 1 - 0.10 / fig_h
            fig.text(x_legend, y_legend, "Top TFs (focal sp.):",
                      fontsize=6, ha="right", va="top", color="0.25")
            for i, (tf, rank, p_t, build) in enumerate(legend_rows):
                color = TF_TRACK_COLORS[(rank - 1) % len(TF_TRACK_COLORS)]
                pct = (f"{p_t*100:.0f}% promoters"
                        if not (p_t != p_t) else "")  # NaN check
                line = f"#{rank} {tf} ({pct})" if pct else f"#{rank} {tf}"
                fig.text(x_legend, y_legend - (0.18 + i * 0.11) / fig_h,
                          line,
                          fontsize=6, ha="right", va="top", color=color,
                          fontweight="bold")
    return fig


# ---------------------------------------------------------------------------
# main()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected", required=True)
    ap.add_argument("--recomputed", default=None)
    ap.add_argument("--distances", required=True)
    ap.add_argument("--top-per-bin", type=int, default=1)
    ap.add_argument("--genome-yaml", required=True)
    ap.add_argument("--genomes-dir", default=os.environ.get("GENOMES_DIR"))
    ap.add_argument("--rbh", required=True)  # back-compat
    ap.add_argument("--epigenomics", default=None)
    ap.add_argument("--tree", default=None,
                    help="calibrated newick for the taxonomy subplot")
    ap.add_argument("--focal", nargs="+", default=["Vertebrata", "Arthropoda"])
    ap.add_argument("--n-species", type=int, default=10)
    ap.add_argument("--scale", choices=["normalized", "absolute", "both"],
                    default="both")
    ap.add_argument("--tf-enrichment-pairs", default=None,
                    help="path to tf_enrichment_per_pair.tsv (Stage 1.7)")
    ap.add_argument("--tf-enrichment-top", default=None,
                    help="path to tf_enrichment_top_tfs.tsv (Stage 1.7)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    if not args.genomes_dir:
        write_placeholder(args.out,
            "Stage 3 requires GENOMES_DIR. Set in env or sibling symlink.\n"
            "See genome_database/annotated_genomes_link.placeholder.")
        return

    sel = pd.read_csv(args.selected, sep="\t")
    sel = sel[sel["nodename"].isin(args.focal)]
    if sel.empty:
        write_placeholder(args.out, "No selected examples in focal clades.")
        return
    sel = (sel.sort_values("score", ascending=False)
              .drop_duplicates(["nodename", "pair_key"])
              .sort_values(["nodename", "mean_in"]))

    dist = pd.read_csv(args.distances, sep="\t")
    recompute = (pd.read_csv(args.recomputed, sep="\t")
                  if args.recomputed and Path(args.recomputed).exists()
                  else None)

    # Stage 1.7 outputs
    tf_pairs = (pd.read_csv(args.tf_enrichment_pairs, sep="\t")
                 if args.tf_enrichment_pairs and Path(args.tf_enrichment_pairs).exists()
                 else None)
    tf_top = (pd.read_csv(args.tf_enrichment_top, sep="\t")
               if args.tf_enrichment_top and Path(args.tf_enrichment_top).exists()
               else None)
    tf_tracks = {}
    if tf_top is not None:
        for (pkey, build), grp in tf_top.groupby(["pair_key", "build"]):
            tf_tracks[(pkey, build)] = grp.to_dict("records")

    cfg = loci_io.load_genome_list(args.genome_yaml, args.genomes_dir)
    tree_obj = load_clade_tree(args.tree)

    # Pre-load BEDs for focal species, restricted to the chroms we'll touch.
    bed_cache: dict[str, dict[str, pd.DataFrame | None]] = {}
    if args.epigenomics and Path(args.epigenomics).exists():
        focal_skeys = set(FOCAL_SPECIES.values())
        focal_dist = dist[dist["species_key"].isin(focal_skeys) &
                          dist["same_scaffold"]]
        hg38_chroms, dm6_chroms = set(), set()
        for _, r in focal_dist.iterrows():
            build = FOCAL_BUILD.get(r["species_key"])
            ucsc = BUILD_TRANSLATE.get(build, {}).get(r["scaffold1"])
            if ucsc:
                (hg38_chroms if build == "hg38" else dm6_chroms).add(ucsc)
        epi = Path(args.epigenomics)
        print(f"loading BEDs (hg38 chroms={len(hg38_chroms)}, "
              f"dm6 chroms={len(dm6_chroms)})")
        bed_cache["hg38"] = {
            "cCRE": load_bed(epi / "hg38" / "ENCODE_cCRE_hg38.bed",
                              keep_chroms=hg38_chroms or None),
            "ReMap_full": load_bed(epi / "hg38" / "ReMap2022_nr_macs2_hg38.bed.gz",
                                     keep_chroms=hg38_chroms or None,
                                     with_label=True),
        }
        bed_cache["dm6"] = {
            "ReMap_full": load_bed(epi / "dm6" / "ReMap2022_nr_macs2_dm6.bed.gz",
                                     keep_chroms=dm6_chroms or None,
                                     with_label=True),
        }
        # Open phastCons + phyloP bigWigs lazily (per-bin stats queried
        # at draw time). hg38 uses 100-vertebrate alignment; dm6 uses
        # 124-arthropod alignment (broader than 27way, matches the user's
        # "wide-taxa conservation" intent).
        for build, key, fname in [
            ("hg38", "phastCons_bw", "phastCons100way_hg38.bw"),
            ("hg38", "phyloP_bw",    "phyloP100way_hg38.bw"),
            ("dm6",  "phastCons_bw", "phastCons124way_dm6.bw"),
            ("dm6",  "phyloP_bw",    "phyloP124way_dm6.bw"),
        ]:
            bw_path = epi / build / fname
            if bw_path.exists():
                try:
                    import pyBigWig
                    bed_cache.setdefault(build, {})[key] = \
                        pyBigWig.open(str(bw_path))
                    print(f"  {build}/{key}: opened {bw_path.name}")
                except Exception as exc:  # noqa: BLE001
                    print(f"  WARN cannot open {bw_path}: {exc}")
            else:
                print(f"  {build}/{key}: not present (skip)")
        for build, beds in bed_cache.items():
            for name, b in beds.items():
                if b is None:
                    print(f"  {build}/{name}: missing")
                elif isinstance(b, pd.DataFrame):
                    print(f"  {build}/{name}: {len(b)} rows")
                else:
                    print(f"  {build}/{name}: opened (bigWig)")

    # Atomic write: render to a sibling .tmp first, then os.rename onto
    # the final path so a concurrent download can never see a 0-byte
    # intermediate state.
    out_final = args.out
    out_tmp = args.out + ".tmp"
    with PdfPages(out_tmp) as pdf:
        n_pages = 0
        for _, pair in sel.iterrows():
            pair_dist = dist[(dist["nodename"] == pair["nodename"]) &
                              (dist["pair_key"] == pair["pair_key"])]
            if pair_dist.empty:
                continue
            r_recompute = None
            if recompute is not None:
                hit = recompute[(recompute["nodename"] == pair["nodename"]) &
                                 (recompute["pair_key"] == pair["pair_key"])]
                if not hit.empty:
                    r_recompute = hit.iloc[0]
            r_pair_enr = None
            if tf_pairs is not None:
                hit = tf_pairs[(tf_pairs["nodename"] == pair["nodename"]) &
                                (tf_pairs["pair_key"] == pair["pair_key"])]
                if not hit.empty:
                    r_pair_enr = hit.iloc[0]
            scales = (["normalized", "absolute"] if args.scale == "both"
                       else [args.scale])
            for sc in scales:
                fig = plot_one_pair(pair, list(pair_dist.to_dict("records")),
                                      r_recompute, cfg, bed_cache,
                                      tree=tree_obj, scale=sc,
                                      tf_tracks=tf_tracks,
                                      pair_enrichment_row=r_pair_enr)
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
                n_pages += 1
        if n_pages == 0:
            write_placeholder(args.out,
                "No same-scaffold species available for any selected pair.")
            return
    os.replace(out_tmp, out_final)
    print(f"wrote {out_final}  pages={n_pages}")


if __name__ == "__main__":
    main()
