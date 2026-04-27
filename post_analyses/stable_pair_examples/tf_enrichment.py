"""Stage 1.7: per-pair TF-binding enrichment test.

The shared-TF intersection produced by `shared_tfs.py` reports the *raw count*
of TFs that bind both anchor TSSes. That isn't a hypothesis test — many TFs
bind many promoters, so some intersection is expected by chance.

This stage computes a per-pair enrichment statistic using each TF's
genome-wide occupancy as the prior. For each pair we report:
  - obs   = TFs binding both anchors' promoter windows
  - exp   = Σ p_t² where p_t is per-TF prior (fraction of all promoters bound)
  - z, p, q (BH-FDR) under independence
  - per-TF surprise (low p_t × co-binding) so the most informative TFs surface

The promoter window is **asymmetric upstream-only**: -2 kb..0 from TSS on +
strand, 0..+2 kb on - strand. This matches the user's preference for
"right upstream of the TSS".

Outputs:
    out/promoter_universe_<build>.tsv   one TSS per gene-locus
    out/tf_enrichment_per_pair.tsv      pair-level z/p/q
    out/tf_enrichment_top_tfs.tsv       per-pair top-N TFs by surprise
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import loci_io
import plot_examples  # for FOCAL_BUILD, BUILD_TRANSLATE, collapse_isoforms

REPO_ROOT = Path(__file__).resolve().parents[2]
PROMOTER_UPSTREAM_BP = 2000
DEFAULT_TOP_N = 10


def build_promoter_universe(species_key: str, genome_cfg: dict,
                              build: str, translate: dict[str, str]
                              ) -> pd.DataFrame:
    """Return one row per gene-locus: chrom (UCSC), tss, strand, gene_id."""
    ent = genome_cfg["species"].get(species_key)
    if not ent or not Path(ent.get("chrom", "")).exists():
        return pd.DataFrame(columns=["chrom", "tss", "strand", "gene_id"])
    cdf = loci_io.read_chrom(ent["chrom"])
    cdf = plot_examples.collapse_isoforms(cdf)
    cdf["tss"] = np.where(cdf["strand"] == "+", cdf["start"], cdf["end"])
    cdf["chrom"] = cdf["scaffold"].map(translate)
    cdf = cdf[cdf["chrom"].notna()].copy()
    return cdf[["chrom", "tss", "strand", "gene_id"]].reset_index(drop=True)


def upstream_window(tss: int, strand: str,
                     half_bp: int = PROMOTER_UPSTREAM_BP) -> tuple[int, int]:
    if strand == "+":
        return tss - half_bp, tss
    return tss, tss + half_bp


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


def per_tf_priors(promoters: pd.DataFrame, remap: pd.DataFrame
                   ) -> pd.DataFrame:
    """For each TF, count promoters with at least one peak in their upstream
    window, then derive p_t = K_t / N_promoters."""
    n = len(promoters)
    promoters = promoters.copy()
    promoters["win_lo"] = np.where(promoters["strand"] == "+",
                                    promoters["tss"] - PROMOTER_UPSTREAM_BP,
                                    promoters["tss"])
    promoters["win_hi"] = np.where(promoters["strand"] == "+",
                                    promoters["tss"],
                                    promoters["tss"] + PROMOTER_UPSTREAM_BP)
    # For speed, group both by chrom and do interval intersection per chrom.
    counts: dict[str, set[int]] = {}  # tf -> set of promoter indices bound
    for chrom, prom_chrom in promoters.groupby("chrom"):
        peaks_chrom = remap[remap["chrom"] == chrom]
        if peaks_chrom.empty:
            continue
        # build a per-TF index list of bound promoters by sweep
        peak_mids = ((peaks_chrom["start"] + peaks_chrom["end"]) // 2).values
        peak_tfs = peaks_chrom["tf"].values
        prom_lo = prom_chrom["win_lo"].values
        prom_hi = prom_chrom["win_hi"].values
        prom_idx = prom_chrom.index.values
        # naive O(N_peaks * log N_promoters) via sort + searchsorted
        order_lo = np.argsort(prom_lo)
        sorted_lo = prom_lo[order_lo]
        sorted_hi = prom_hi[order_lo]
        sorted_idx = prom_idx[order_lo]
        for mid, tf in zip(peak_mids, peak_tfs):
            # find candidate promoters with win_lo <= mid
            i = np.searchsorted(sorted_lo, mid, side="right")
            if i == 0:
                continue
            # check the candidates whose win_hi >= mid; in practice few
            for k in range(i - 1, max(-1, i - 50), -1):
                if sorted_hi[k] >= mid and sorted_lo[k] <= mid:
                    counts.setdefault(tf, set()).add(int(sorted_idx[k]))
                if sorted_lo[k] < mid - PROMOTER_UPSTREAM_BP:
                    break
    rows = [{"tf": tf, "K_t": len(prom_set), "p_t": len(prom_set) / n}
            for tf, prom_set in counts.items()]
    rows.sort(key=lambda r: -r["K_t"])
    return pd.DataFrame(rows)


def tfs_at_window(remap: pd.DataFrame, chrom: str,
                   lo: int, hi: int) -> set[str]:
    if remap is None:
        return set()
    sub = remap[(remap["chrom"] == chrom) &
                (remap["end"] >= lo) & (remap["start"] <= hi)]
    return set(sub["tf"].dropna().unique())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--distances", required=True,
                    help="out/selected_examples_distances.tsv")
    ap.add_argument("--epigenomics", required=True,
                    help="out/epigenomics dir")
    ap.add_argument("--genome-yaml",
                    default=str(REPO_ROOT / "genome_database/genome_list.yaml"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--top-n", type=int, default=DEFAULT_TOP_N)
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    epi = Path(args.epigenomics)

    dist = pd.read_csv(args.distances, sep="\t")
    dist = dist[dist["same_scaffold"] == True].copy()
    if dist.empty:
        print("no same-scaffold rows; nothing to do")
        return

    cfg = loci_io.load_genome_list(args.genome_yaml)

    # Per build: build promoter universe, load ReMap, compute priors. Cache to
    # disk so re-runs are fast.
    build_data: dict[str, dict] = {}
    for build in {plot_examples.FOCAL_BUILD.get(sk)
                  for sk in dist["species_key"].unique()
                  if plot_examples.FOCAL_BUILD.get(sk)}:
        focal_sk = next(sk for sk, b in plot_examples.FOCAL_BUILD.items()
                        if b == build)
        translate = plot_examples.BUILD_TRANSLATE.get(build, {})
        prom_path = out_dir / f"promoter_universe_{build}.tsv"
        if prom_path.exists() and prom_path.stat().st_size > 0:
            print(f"[{build}] using cached {prom_path}")
            promoters = pd.read_csv(prom_path, sep="\t")
        else:
            print(f"[{build}] building promoter universe from {focal_sk} ...")
            promoters = build_promoter_universe(focal_sk, cfg, build, translate)
            promoters.to_csv(prom_path, sep="\t", index=False)
            print(f"  wrote {prom_path}  rows={len(promoters)}")
        if promoters.empty:
            print(f"[{build}] empty promoter universe; skipping")
            continue
        # ReMap restricted to chroms in our promoter universe.
        remap_chroms = set(promoters["chrom"].unique())
        remap_path = epi / build / f"ReMap2022_nr_macs2_{build}.bed.gz"
        print(f"[{build}] loading ReMap ({len(remap_chroms)} chroms) ...")
        remap = load_remap_with_tf(remap_path, keep_chroms=remap_chroms)
        if remap is None or remap.empty:
            print(f"[{build}] ReMap empty; skipping")
            continue
        print(f"[{build}] {len(remap):,} peak rows kept")
        priors_path = out_dir / f"tf_priors_{build}.tsv"
        if priors_path.exists() and priors_path.stat().st_size > 0:
            print(f"[{build}] using cached priors {priors_path}")
            priors = pd.read_csv(priors_path, sep="\t")
        else:
            print(f"[{build}] computing per-TF priors (this is the slow part) ...")
            priors = per_tf_priors(promoters, remap)
            priors.to_csv(priors_path, sep="\t", index=False)
            print(f"  wrote {priors_path}  rows={len(priors)}")
        build_data[build] = {"promoters": promoters, "remap": remap,
                              "priors": priors,
                              "p_t_lookup": dict(zip(priors["tf"], priors["p_t"]))}

    if not build_data:
        print("no usable build data; nothing to test")
        return

    # Per-pair test.
    per_pair_rows = []
    per_tf_rows = []
    for _, r in dist.iterrows():
        sk = r["species_key"]
        build = plot_examples.FOCAL_BUILD.get(sk)
        if build is None or build not in build_data:
            continue
        bd = build_data[build]
        translate = plot_examples.BUILD_TRANSLATE.get(build, {})
        ucsc = translate.get(r["scaffold1"])
        if not ucsc:
            continue
        lo1, hi1 = upstream_window(int(r["tss1"]), r["strand1"])
        lo2, hi2 = upstream_window(int(r["tss2"]), r["strand2"])
        tfs1 = tfs_at_window(bd["remap"], ucsc, lo1, hi1)
        tfs2 = tfs_at_window(bd["remap"], ucsc, lo2, hi2)
        shared = tfs1 & tfs2
        all_tfs = tfs1 | tfs2 | set(bd["p_t_lookup"].keys())  # universe = TFs in priors
        p_lookup = bd["p_t_lookup"]
        # Aggregate independence test across the universe of TFs assayed.
        ps = np.array([p_lookup.get(tf, 0.0) for tf in p_lookup if p_lookup.get(tf, 0.0) > 0])
        # NB: p_lookup keys ARE the universe (TFs that bound at least one promoter).
        exp = float(np.sum(ps ** 2))
        var = float(np.sum(ps ** 2 * (1 - ps ** 2)))
        obs = len(shared)
        z = (obs - exp) / np.sqrt(var) if var > 0 else 0.0
        # one-sided p (upper tail) via standard normal
        from math import erf, sqrt
        p_indep = 0.5 * (1 - erf(z / sqrt(2)))
        per_pair_rows.append(dict(
            nodename=r["nodename"], pair_label=r["pair_label"],
            pair_key=r["pair_key"], species_key=sk, build=build,
            ucsc_chrom=ucsc, n_tfs_assayed=len(p_lookup),
            n_at_anchor1=len(tfs1), n_at_anchor2=len(tfs2),
            n_shared=obs, expected=exp, z=z, p_indep=p_indep,
        ))
        # Per-TF surprise table (only for TFs in the intersection so far).
        for tf in shared:
            p_t = p_lookup.get(tf, 0.0)
            if p_t <= 0:
                continue
            surprise = -2.0 * np.log10(p_t)
            per_tf_rows.append(dict(
                nodename=r["nodename"], pair_label=r["pair_label"],
                pair_key=r["pair_key"], species_key=sk, build=build,
                tf=tf, p_t=p_t, surprise=surprise,
            ))

    if not per_pair_rows:
        print("no enrichment rows produced")
        return

    pp_df = pd.DataFrame(per_pair_rows)
    # BH-FDR over p_indep within the table.
    pp_df = pp_df.sort_values("p_indep")
    m = len(pp_df)
    pp_df["q_indep"] = (pp_df["p_indep"].values * m
                         / np.arange(1, m + 1)).clip(0, 1)
    # enforce monotone q (BH)
    qs = pp_df["q_indep"].values[::-1]
    qs = np.minimum.accumulate(qs)
    pp_df["q_indep"] = qs[::-1]
    pp_path = out_dir / "tf_enrichment_per_pair.tsv"
    pp_df.to_csv(pp_path, sep="\t", index=False)
    print(f"wrote {pp_path}  rows={len(pp_df)}")

    tf_df = pd.DataFrame(per_tf_rows)
    if not tf_df.empty:
        tf_df = (tf_df.sort_values(["pair_key", "surprise"], ascending=[True, False])
                       .groupby("pair_key", group_keys=False)
                       .head(args.top_n)
                       .reset_index(drop=True))
        tf_df["rank"] = tf_df.groupby("pair_key").cumcount() + 1
    tf_path = out_dir / "tf_enrichment_top_tfs.tsv"
    tf_df.to_csv(tf_path, sep="\t", index=False)
    print(f"wrote {tf_path}  rows={len(tf_df)}")

    # Smoke check: EIF2S1/ATP6V1D top TFs
    eif = tf_df[tf_df["pair_label"].isin(["EIF2S1/ATP6V1D", "ATP6V1D/EIF2S1"])]
    if not eif.empty:
        print("  CHECK EIF2S1/ATP6V1D top-3 TFs by surprise:")
        for _, x in eif.head(3).iterrows():
            print(f"    {x['tf']}  p_t={x['p_t']:.4f}  surprise={x['surprise']:.2f}")


if __name__ == "__main__":
    main()
