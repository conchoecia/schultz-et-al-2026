"""Stage 1.5: orientation-aware TSS-to-TSS recompute for selected examples.

The COO-derived `mean_in` in defining_features_pairs/out/unique_pairs.tsv.gz
is gene-edge / CDS-extent based, not TSS-to-TSS. Divergent head-to-head pairs
like vertebrate EIF2S1/ATP6V1D mis-bin as ~50 kb (edge) when their TSS-TSS
distance is <1 kb. For each selected example pair, this stage opens the
anchor species' .chrom.gz files in place, looks up both anchor genes via the
per-species RBH, and computes:

    tss_to_tss_bp     |TSS_a - TSS_b|, where TSS = start if + else end
    edge_to_edge_bp   gap between the two CDS extents (0 if overlapping)
    orientation       divergent / convergent / tandem_same_strand
    same_scaffold     True iff both genes hit the same scaffold

Outputs:
    out/selected_examples_distances.tsv     long-form, one row per (pair, species)
    out/selected_examples_recomputed.tsv    one row per pair (median across species)
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

import loci_io

REPO_ROOT = Path(__file__).resolve().parents[2]

# Default RBH dir per top-level config.yaml. Holds per-species
# BCnSSimakov2022_<species_key>_xy_reciprocal_best_hits.plotted.rbh files.
DEFAULT_RBH_DIR = Path(
    "/lisc/data/scratch/molevo/dts/manifold/submission_dryad/dryad_repo/"
    "BCnSSimakov2022_current_rbh_202509")

# Curated anchor species per focal clade — the species we resolve TSS for.
# 8-10 species per clade chosen to span major lineages and to maximize the
# probability that both anchor genes resolve on the same scaffold across
# diverse evolutionary distances.
ANCHORS = {
    "Vertebrata": [
        "Homosapiens-9606-GCF000001405.40",        # primate
        "Musmusculus-10090-GCF000001635.27",       # rodent
        "Rattusnorvegicus-10116-GCF000001895.5",   # rodent
        "Pantroglodytes-9598-GCF000001515.7",      # great-ape primate
        "Gallusgallus-9031-GCF000002315.6",        # bird
        "Anoliscarolinensis-28377-GCF000090745.2", # squamate reptile
        "Xenopustropicalis-8364-GCF000004195.4",   # amphibian
        "Daniorerio-7955-GCF000002035.6",          # teleost fish
        "Lepisosteusoculatus-7918-GCF000242695.1", # holostean (gar) — early-branching ray-finned fish
        "Petromyzonmarinus-7757-GCF010993605.1",   # cyclostome (lamprey)
    ],
    "Arthropoda": [
        "Drosophilamelanogaster-7227-GCF000001215.4",  # Diptera (fly)
        "Anophelesfunestus-62324-GCF943734845.2",      # Diptera (mosquito; in tree)
        "Aedesaegypti-7159-GCF002204515.2",            # Diptera (mosquito)
        "Apismellifera-7460-GCF000002195.4",           # Hymenoptera (bee)
        "Bombusterrestris-30195-GCF000214255.1",       # Hymenoptera (bumblebee)
        "Triboliumcastaneum-7070-GCF000002335.3",      # Coleoptera (beetle)
        "Bombyxmori-7091-GCF014905235.1",              # Lepidoptera (silkmoth)
        "Daphniapulex-6669-GCF021134715.1",            # Crustacea (water flea)
        "Penaeusmonodon-6687-GCF015228065.2",          # Crustacea (shrimp)
        "Parasteatodatepidariorum-114398-GCF043381705.1",  # Chelicerata (spider)
    ],
}


def per_species_rbh_path(rbh_dir: Path, species_key: str) -> Path:
    return rbh_dir / f"BCnSSimakov2022_{species_key}_xy_reciprocal_best_hits.plotted.rbh"


def load_species_family_to_gene(rbh_path: Path,
                                  species_key: str) -> dict[str, str]:
    """Return {family_id: gene_id} for one species' RBH file."""
    if not rbh_path.exists():
        return {}
    gene_col = f"{species_key}_gene"
    df = pd.read_csv(rbh_path, sep="\t", usecols=["rbh", gene_col],
                     dtype={"rbh": str, gene_col: str})
    df = df.dropna(subset=[gene_col])
    return dict(zip(df["rbh"], df[gene_col]))


def lookup_gene(chrom_df: pd.DataFrame, gene_id: str,
                 max_isoform_bp: int = 200_000):
    """Look up a gene by exact gene_id. Rejects entries whose CDS span
    exceeds `max_isoform_bp` — these are NCBI-predicted readthrough
    transcripts (e.g. multi-Mb XP_ entries spanning multiple real genes)."""
    hit = chrom_df[chrom_df["gene_id"] == gene_id]
    if hit.empty:
        return None
    r = hit.iloc[0]
    span = int(r["end"]) - int(r["start"])
    if span > max_isoform_bp:
        return None  # treat as missing → triggers alt-isoform fallback
    return dict(scaffold=r["scaffold"], strand=r["strand"],
                start=int(r["start"]), end=int(r["end"]))


def find_alt_isoform_on_scaffold(chrom_df: pd.DataFrame,
                                   gene_id: str,
                                   target_scaffold: str,
                                   accession_to_symbol: dict[str, str],
                                   partner_pos: int | None = None,
                                   max_isoform_bp: int = 200_000
                                   ) -> dict | None:
    """When the RBH-given gene_id resolves to a different scaffold than the
    partner, look for any other accession with the SAME gene symbol that IS
    on the target scaffold. Useful when BCnSSimakov2022.rbh picked an
    alt-haplotype isoform.

    Filters: rejects accessions whose CDS span exceeds `max_isoform_bp`
    (these are NCBI-predicted readthrough transcripts spanning multiple
    real genes — same artifact `collapse_isoforms` handles in plotting).

    Selection: among the surviving same-symbol accessions on
    `target_scaffold`, picks the one whose midpoint is closest to
    `partner_pos`. If `partner_pos` is None, picks the longest sane CDS.
    """
    sym = accession_to_symbol.get(gene_id)
    if not sym:
        return None
    same_sym_accs = {a for a, s in accession_to_symbol.items() if s == sym}
    if not same_sym_accs:
        return None
    cands = chrom_df[(chrom_df["scaffold"] == target_scaffold) &
                      (chrom_df["gene_id"].isin(same_sym_accs))].copy()
    if cands.empty:
        return None
    cands = cands[(cands["end"] - cands["start"]) <= max_isoform_bp]
    if cands.empty:
        return None
    if partner_pos is not None:
        cands = cands.assign(
            _mid=(cands["start"] + cands["end"]) // 2,
            _dist=lambda d: (d["_mid"] - partner_pos).abs())
        r = cands.sort_values("_dist").iloc[0]
    else:
        cands = cands.assign(_w=cands["end"] - cands["start"])
        r = cands.sort_values("_w", ascending=False).iloc[0]
    return dict(scaffold=r["scaffold"], strand=r["strand"],
                start=int(r["start"]), end=int(r["end"]),
                gene_id=str(r["gene_id"]))


def load_accession_to_symbol(path: Path) -> dict[str, str]:
    """NCBI human_gene2accession.tsv.gz -> {NP_/XP_ accession: HGNC symbol}."""
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t",
                      usecols=["protein_accession.version", "Symbol"],
                      dtype=str)
    df = df[df["protein_accession.version"].fillna("-") != "-"]
    return dict(zip(df["protein_accession.version"], df["Symbol"]))


def orientation(s1, s2, tss1, tss2) -> str:
    if s1 == s2:
        return "tandem_same_strand"
    # + strand TSS at start (left edge), - strand TSS at end (right edge).
    # divergent: TSSes face away from each other (genes pointing outward)
    # convergent: TSSes face each other (genes pointing inward, 3' ends meet)
    if s1 == "+" and s2 == "-":
        return "convergent" if tss1 < tss2 else "divergent"
    if s1 == "-" and s2 == "+":
        return "convergent" if tss2 < tss1 else "divergent"
    return "unknown"


def recompute_pair_species(pair_row, species_key, fam2gene, chrom_df,
                              acc_to_sym: dict[str, str] | None = None) -> dict:
    fam1, fam2 = pair_row["ortholog1"], pair_row["ortholog2"]
    g1, g2 = fam2gene.get(fam1), fam2gene.get(fam2)
    out = {
        "nodename": pair_row["nodename"],
        "pair_label": pair_row["pair_label"],
        "pair_key": pair_row["pair_key"],
        "species_key": species_key,
        "ortholog1": fam1, "ortholog2": fam2,
        "gene1": g1 or "", "gene2": g2 or "",
        "scaffold1": "", "scaffold2": "",
        "strand1": "", "strand2": "",
        "start1": pd.NA, "end1": pd.NA,
        "start2": pd.NA, "end2": pd.NA,
        "tss1": pd.NA, "tss2": pd.NA,
        "tss_to_tss_bp": pd.NA,
        "edge_to_edge_bp": pd.NA,
        "orientation": "",
        "same_scaffold": False,
        "status": "ok",
    }
    if not g1 or not g2:
        out["status"] = "missing_rbh_gene"
        return out
    h1, h2 = lookup_gene(chrom_df, g1), lookup_gene(chrom_df, g2)
    # If either gene was rejected (predicted readthrough) or missing, try
    # the symbol-based alt-isoform fallback BEFORE giving up.
    if (h1 is None or h2 is None) and acc_to_sym:
        # need a partner_pos to anchor the search; take the median chrom
        # position of any same-symbol non-readthrough isoform.
        for which, missing_id in [("h1", g1), ("h2", g2)]:
            if (which == "h1" and h1 is not None) or (which == "h2" and h2 is not None):
                continue
            sym = acc_to_sym.get(missing_id)
            if not sym:
                continue
            same_sym = {a for a, s in acc_to_sym.items() if s == sym}
            cands = chrom_df[chrom_df["gene_id"].isin(same_sym)].copy()
            cands = cands[(cands["end"] - cands["start"]) <= 200_000]
            if cands.empty:
                continue
            # Prefer the most-frequent scaffold (canonical chromosome).
            top_scaf = cands["scaffold"].value_counts().index[0]
            cands = cands[cands["scaffold"] == top_scaf]
            r = cands.assign(_w=cands["end"] - cands["start"]) \
                       .sort_values("_w", ascending=False).iloc[0]
            new_h = dict(scaffold=r["scaffold"], strand=r["strand"],
                          start=int(r["start"]), end=int(r["end"]),
                          gene_id=str(r["gene_id"]))
            if which == "h1":
                h1 = new_h
                out.update(gene1=new_h["gene_id"], status="ok_alt_isoform_g1")
            else:
                h2 = new_h
                out.update(gene2=new_h["gene_id"], status="ok_alt_isoform_g2")
    if not h1 or not h2:
        out["status"] = "missing_chrom_gene"
        return out
    out.update(scaffold1=h1["scaffold"], scaffold2=h2["scaffold"],
               strand1=h1["strand"], strand2=h2["strand"],
               start1=h1["start"], end1=h1["end"],
               start2=h2["start"], end2=h2["end"])
    if h1["scaffold"] != h2["scaffold"]:
        # Try the alt-isoform fallback: maybe one of the genes is on an
        # alt-haplotype contig in the RBH but has a canonical isoform on
        # the partner's scaffold (common for human GCF assembly).
        if acc_to_sym:
            # Pass partner midpoint so the fallback picks the closest isoform
            # on the partner's scaffold (not a far-away alt that happens to
            # share a chromosome).
            partner2_mid = (h2["start"] + h2["end"]) // 2
            alt1 = find_alt_isoform_on_scaffold(chrom_df, g1, h2["scaffold"],
                                                  acc_to_sym, partner_pos=partner2_mid)
            if alt1:
                h1 = alt1
                out.update(gene1=alt1["gene_id"], scaffold1=alt1["scaffold"],
                            strand1=alt1["strand"], start1=alt1["start"],
                            end1=alt1["end"], status="ok_alt_isoform_g1")
            else:
                partner1_mid = (h1["start"] + h1["end"]) // 2
                alt2 = find_alt_isoform_on_scaffold(chrom_df, g2, h1["scaffold"],
                                                      acc_to_sym, partner_pos=partner1_mid)
                if alt2:
                    h2 = alt2
                    out.update(gene2=alt2["gene_id"], scaffold2=alt2["scaffold"],
                                strand2=alt2["strand"], start2=alt2["start"],
                                end2=alt2["end"], status="ok_alt_isoform_g2")
        if h1["scaffold"] != h2["scaffold"]:
            out["status"] = "different_scaffold"
            return out
    tss1 = h1["start"] if h1["strand"] == "+" else h1["end"]
    tss2 = h2["start"] if h2["strand"] == "+" else h2["end"]
    edge = max(0, max(h1["start"], h2["start"]) - min(h1["end"], h2["end"]))
    out.update(same_scaffold=True, tss1=tss1, tss2=tss2,
               tss_to_tss_bp=int(abs(tss1 - tss2)),
               edge_to_edge_bp=int(edge),
               orientation=orientation(h1["strand"], h2["strand"], tss1, tss2))
    return out


def re_bin(tss_bp: float) -> str:
    if pd.isna(tss_bp):
        return "unknown"
    if tss_bp < 1e3: return "close"
    if tss_bp < 1e5: return "intermediate"
    if tss_bp < 1e6: return "long_range"
    return "ultra_long"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--selected", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--genome-yaml",
                    default=str(REPO_ROOT / "genome_database/genome_list.yaml"))
    ap.add_argument("--rbh", default=None,
                    help="(unused; per-species RBH dir is --rbh-dir)")
    ap.add_argument("--rbh-dir", default=str(DEFAULT_RBH_DIR))
    ap.add_argument("--anchors", action="append", default=None,
                    help="extra species_key to add to the anchor set "
                         "(applies to all clades)")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    sel = pd.read_csv(args.selected, sep="\t")
    if sel.empty:
        print("no selected examples; nothing to do")
        return
    # de-dupe pairs in case selection picked the same pair multiple times
    sel = sel.drop_duplicates(["nodename", "pair_key"], keep="first")

    cfg = loci_io.load_genome_list(args.genome_yaml)
    species_cfg = cfg.get("species", {})
    rbh_dir = Path(args.rbh_dir)
    if not rbh_dir.exists():
        print(f"WARN: rbh_dir not found ({rbh_dir}); cannot recompute distances")
        return

    # Expand anchor list per clade with any --anchors additions.
    anchors_by_clade = {k: list(v) for k, v in ANCHORS.items()}
    if args.anchors:
        for clade, lst in anchors_by_clade.items():
            lst.extend(args.anchors)

    # Load NCBI human accession→symbol once for the alt-isoform fallback
    # (only used for human; other species don't have an analogous bridge here).
    ncbi_acc_path = (REPO_ROOT
        / "post_analyses/go_enrichment_sweep/ncbi_ref/human_gene2accession.tsv.gz")
    acc_to_sym = load_accession_to_symbol(ncbi_acc_path)
    if acc_to_sym:
        print(f"loaded {len(acc_to_sym):,} NCBI human accession→symbol bridge "
              f"entries (for RBH alt-isoform fallback)")

    # Cache per-species RBH and chrom.gz reads (one per species, not per pair).
    fam2gene_cache: dict[str, dict[str, str]] = {}
    chrom_cache: dict[str, pd.DataFrame] = {}

    rows = []
    for _, pair in sel.iterrows():
        clade = pair["nodename"]
        anchors = anchors_by_clade.get(clade, [])
        if not anchors:
            continue
        for sk in anchors:
            if sk not in fam2gene_cache:
                fam2gene_cache[sk] = load_species_family_to_gene(
                    per_species_rbh_path(rbh_dir, sk), sk)
            if sk not in chrom_cache:
                ent = species_cfg.get(sk)
                if not ent or not Path(ent.get("chrom", "")).exists():
                    chrom_cache[sk] = pd.DataFrame(columns=loci_io.CHROM_COLUMNS)
                else:
                    chrom_cache[sk] = loci_io.read_chrom(ent["chrom"])
            # Apply alt-isoform fallback only for human (it's the species
            # whose RBH most often picks alt-haplotype contigs).
            sym_bridge = acc_to_sym if sk.startswith("Homosapiens-") else None
            rows.append(recompute_pair_species(
                pair, sk, fam2gene_cache[sk], chrom_cache[sk], sym_bridge))

    if not rows:
        print("no recompute rows produced")
        return

    long_df = pd.DataFrame(rows)
    long_path = out_dir / "selected_examples_distances.tsv"
    long_df.to_csv(long_path, sep="\t", index=False)
    print(f"wrote {long_path}  rows={len(long_df)}")

    # Per-pair median across species_with same_scaffold==True.
    ok = long_df[long_df["same_scaffold"]]
    summary = []
    for (clade, pkey), g in long_df.groupby(["nodename", "pair_key"]):
        ok_g = g[g["same_scaffold"]]
        if len(ok_g):
            tss = float(np.nanmedian(ok_g["tss_to_tss_bp"].astype(float)))
            edge = float(np.nanmedian(ok_g["edge_to_edge_bp"].astype(float)))
            modal = ok_g["orientation"].mode().iloc[0] if len(ok_g) else ""
        else:
            tss = float("nan"); edge = float("nan"); modal = ""
        summary.append(dict(
            nodename=clade, pair_key=pkey,
            pair_label=g["pair_label"].iloc[0],
            n_species_attempted=len(g),
            n_species_same_scaffold=int(g["same_scaffold"].sum()),
            median_tss_to_tss_bp=tss,
            median_edge_to_edge_bp=edge,
            modal_orientation=modal,
            bin_tss=re_bin(tss),
        ))
    summary_df = pd.DataFrame(summary)
    sum_path = out_dir / "selected_examples_recomputed.tsv"
    summary_df.to_csv(sum_path, sep="\t", index=False)
    print(f"wrote {sum_path}  rows={len(summary_df)}")

    # Smoke check: EIF2S1/ATP6V1D@Vertebrata.
    eif = summary_df[(summary_df["nodename"] == "Vertebrata") &
                      summary_df["pair_label"].isin(
                          ["EIF2S1/ATP6V1D", "ATP6V1D/EIF2S1"])]
    if len(eif):
        r = eif.iloc[0]
        print(f"  CHECK EIF2S1/ATP6V1D@Vertebrata: "
              f"TSS-TSS={r['median_tss_to_tss_bp']:.0f} bp  "
              f"edge={r['median_edge_to_edge_bp']:.0f} bp  "
              f"orientation={r['modal_orientation']}  bin_tss={r['bin_tss']}  "
              f"n_same_scaf={r['n_species_same_scaffold']}/{r['n_species_attempted']}")


if __name__ == "__main__":
    main()
