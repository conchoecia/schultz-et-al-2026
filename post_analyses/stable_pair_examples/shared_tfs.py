"""Stage 1.6: shared transcription factors between anchor TSSes.

For each selected example pair × anchor species, identify TFs whose ChIP-seq
peaks fall within ±W bp of EACH anchor TSS, then take the intersection. The
shared set is what justifies a "co-regulation" interpretation (the same TF
binds both promoters in the same species).

Inputs:
    out/selected_examples_distances.tsv   per-(pair, species) anchor coords
    out/epigenomics/{hg38,dm6}/ReMap2022_nr_macs2_*.bed.gz   TF:cell rows
Outputs:
    out/shared_tfs_per_pair_species.tsv   long-form (pair × species)
    out/shared_tfs_summary.tsv            per pair: union over species

ReMap column 4 is "TF:cell_type" — we strip the cell suffix and report
unique TF names.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import plot_examples  # for FOCAL_BUILD + BUILD_TRANSLATE


def load_remap_with_tf(path: Path,
                       keep_chroms: set[str] | None) -> pd.DataFrame | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     usecols=[0, 1, 2, 3],
                     names=["chrom", "start", "end", "label"],
                     dtype={"chrom": str, "start": "int64", "end": "int64",
                             "label": str},
                     compression="infer", low_memory=False,
                     on_bad_lines="skip")
    if keep_chroms is not None:
        df = df[df["chrom"].isin(keep_chroms)].copy()
    df["tf"] = df["label"].str.split(":").str[0]
    return df


def tfs_at_locus(remap: pd.DataFrame, ucsc_chrom: str,
                  pos: int, half_window_bp: int = 2000) -> set[str]:
    if remap is None or not ucsc_chrom:
        return set()
    lo, hi = pos - half_window_bp, pos + half_window_bp
    sub = remap[(remap["chrom"] == ucsc_chrom) &
                (remap["end"] >= lo) & (remap["start"] <= hi)]
    return set(sub["tf"].dropna().unique())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--distances", required=True)
    ap.add_argument("--epigenomics", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--half-window-bp", type=int, default=2000,
                    help="±bp around each TSS to scan for ChIP peaks")
    args = ap.parse_args(argv)

    dist = pd.read_csv(args.distances, sep="\t")
    dist = dist[dist["same_scaffold"] == True].copy()
    if dist.empty:
        print("no same-scaffold rows in distances; nothing to do")
        return

    # Pre-filter ReMap BEDs to the chromosomes we need.
    epi = Path(args.epigenomics)
    remap_cache: dict[str, pd.DataFrame | None] = {}
    needed_chroms: dict[str, set[str]] = {"hg38": set(), "dm6": set()}
    for _, r in dist.iterrows():
        build = plot_examples.FOCAL_BUILD.get(r["species_key"])
        if not build:
            continue
        ucsc = plot_examples.BUILD_TRANSLATE.get(build, {}).get(r["scaffold1"])
        if ucsc:
            needed_chroms[build].add(ucsc)
    for build, chroms in needed_chroms.items():
        if not chroms:
            remap_cache[build] = None
            continue
        path = epi / build / f"ReMap2022_nr_macs2_{build}.bed.gz"
        print(f"loading {build} ReMap with {len(chroms)} chroms ...")
        remap_cache[build] = load_remap_with_tf(path, keep_chroms=chroms)
        n = 0 if remap_cache[build] is None else len(remap_cache[build])
        print(f"  {build}: {n:,} peak rows kept")

    rows = []
    for _, r in dist.iterrows():
        build = plot_examples.FOCAL_BUILD.get(r["species_key"])
        if build is None:
            continue
        remap = remap_cache.get(build)
        ucsc = plot_examples.BUILD_TRANSLATE.get(build, {}).get(r["scaffold1"])
        if remap is None or not ucsc:
            continue
        tss1 = int(r["tss1"]); tss2 = int(r["tss2"])
        tfs1 = tfs_at_locus(remap, ucsc, tss1, args.half_window_bp)
        tfs2 = tfs_at_locus(remap, ucsc, tss2, args.half_window_bp)
        shared = tfs1 & tfs2
        rows.append(dict(
            nodename=r["nodename"], pair_label=r["pair_label"],
            pair_key=r["pair_key"], species_key=r["species_key"],
            build=build, ucsc_chrom=ucsc,
            tss1=tss1, tss2=tss2,
            n_tfs_at_anchor1=len(tfs1),
            n_tfs_at_anchor2=len(tfs2),
            n_shared_tfs=len(shared),
            shared_tfs=";".join(sorted(shared)),
        ))

    out_dir = Path(args.out_dir)
    long_df = pd.DataFrame(rows)
    long_path = out_dir / "shared_tfs_per_pair_species.tsv"
    long_df.to_csv(long_path, sep="\t", index=False)
    print(f"wrote {long_path}  rows={len(long_df)}")

    if long_df.empty:
        return

    # Summary: per pair, union of shared TFs across all species, plus
    # intersection (TFs shared in EVERY focal species).
    summary = []
    for (clade, pkey), g in long_df.groupby(["nodename", "pair_key"]):
        all_shared_sets = [set(s.split(";")) - {""} for s in g["shared_tfs"]]
        if not all_shared_sets:
            continue
        union = set().union(*all_shared_sets)
        intersect = set.intersection(*all_shared_sets) if len(all_shared_sets) > 1 else all_shared_sets[0]
        summary.append(dict(
            nodename=clade, pair_key=pkey,
            pair_label=g["pair_label"].iloc[0],
            n_focal_species=len(g),
            n_shared_tfs_intersection=len(intersect),
            n_shared_tfs_union=len(union),
            intersection_tfs=";".join(sorted(intersect)),
            union_tfs=";".join(sorted(union)),
        ))
    sum_df = pd.DataFrame(summary)
    sum_path = out_dir / "shared_tfs_summary.tsv"
    sum_df.to_csv(sum_path, sep="\t", index=False)
    print(f"wrote {sum_path}  rows={len(sum_df)}")

    # Quick check for the EIF2S1/ATP6V1D + NRF1 hypothesis.
    eif = sum_df[sum_df["pair_label"].isin(["EIF2S1/ATP6V1D", "ATP6V1D/EIF2S1"])]
    if not eif.empty:
        r = eif.iloc[0]
        print(f"  CHECK EIF2S1/ATP6V1D: intersection across "
              f"{r['n_focal_species']} focal species = "
              f"{r['n_shared_tfs_intersection']} TFs")
        if r["intersection_tfs"]:
            print(f"    {r['intersection_tfs']}")


if __name__ == "__main__":
    main()
