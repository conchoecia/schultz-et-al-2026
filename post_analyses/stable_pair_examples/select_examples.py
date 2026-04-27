"""Stage 1: select stable-pair examples per clade and distance bin.

Reads the canonical unique_pairs.tsv.gz from defining_features_pairs/out,
joins family→human-symbol, GO co-membership, entanglement fold-enrichment,
and reviewer-fusion overlap. Emits:

    out/candidate_pairs.tsv     all qualifying stable pairs, all clades
    out/selected_examples.tsv   top-N per (clade, distance_bin) by composite score
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DEF_INPUT = REPO_ROOT / "post_analyses/defining_features_pairs/out/unique_pairs.tsv.gz"
DEF_FAMMAP = REPO_ROOT / "post_analyses/build_family_naming_map/out/bcns_family_to_human_gene.tsv"
DEF_GOSIG = REPO_ROOT / "post_analyses/go_enrichment_sweep/out/significant_terms_annotated.tsv"
DEF_ENTANG = REPO_ROOT / "post_analyses/entanglement_browse/out/entangled_pairs_per_clade.tsv"
DEF_FUSION = REPO_ROOT / "post_analyses/reviewer_common_fusions/pair_metrics.tsv"
DEF_NCBI_ACC = REPO_ROOT / "post_analyses/go_enrichment_sweep/ncbi_ref/human_gene2accession.tsv.gz"


def load_np_to_symbol(path: Path) -> dict[str, str]:
    """NP_xxx (and XP_xxx) accession → HGNC symbol from NCBI gene2accession."""
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", usecols=["protein_accession.version", "Symbol"],
                     dtype=str)
    df = df[df["protein_accession.version"].fillna("-") != "-"]
    return dict(zip(df["protein_accession.version"], df["Symbol"]))

DISTANCE_BINS = [
    ("close",        0,        1_000),
    ("gap_short",    1_000,    10_000),
    ("intermediate", 10_000,   100_000),
    ("long_range",   100_000,  1_000_000),
    ("ultra_long",   1_000_000, np.inf),
]


def assign_bin(d: float) -> str:
    for name, lo, hi in DISTANCE_BINS:
        if lo <= d < hi:
            return name
    return "unknown"


def pair_key(a: str, b: str) -> str:
    return "|".join(sorted([str(a), str(b)]))


def add_symbols(df: pd.DataFrame, fammap_path: Path,
                np_to_sym: dict[str, str]) -> pd.DataFrame:
    fam = pd.read_csv(fammap_path, sep="\t")[["family_id", "alg", "human_gene"]]
    fam["symbol"] = fam["human_gene"].map(np_to_sym).fillna(fam["human_gene"])
    fam = fam.rename(columns={"family_id": "ortholog", "alg": "alg_code"})
    fam = fam[["ortholog", "alg_code", "human_gene", "symbol"]]
    df = df.merge(fam.add_suffix("1"), how="left",
                  left_on="ortholog1", right_on="ortholog1")
    df = df.merge(fam.add_suffix("2"), how="left",
                  left_on="ortholog2", right_on="ortholog2")
    df["pair_label"] = (df["symbol1"].fillna(df["ortholog1"]).astype(str)
                        + "/" +
                        df["symbol2"].fillna(df["ortholog2"]).astype(str))
    return df


def add_go_support(df: pd.DataFrame, gosig_path: Path) -> pd.DataFrame:
    if not gosig_path.exists():
        df["go_terms_shared"] = 0
        df["go_hits_pair"] = 0
        return df
    go = pd.read_csv(gosig_path, sep="\t", dtype=str)
    go["bcns_families"] = go["bcns_families"].fillna("")
    # For each (clade, term) row, list of families hit.
    go["fam_set"] = go["bcns_families"].apply(
        lambda s: frozenset(x for x in s.split(";") if x))
    # Map (clade, family) -> count of significant terms it appears in.
    long = go[["clade", "fam_set"]].explode("fam_set").dropna()
    hits_per_fam = long.groupby(["clade", "fam_set"]).size().rename("hits").reset_index()
    hits_per_fam = hits_per_fam.rename(columns={"fam_set": "family"})
    h1 = hits_per_fam.rename(columns={"clade": "nodename", "family": "ortholog1",
                                       "hits": "go_hits1"})
    h2 = hits_per_fam.rename(columns={"clade": "nodename", "family": "ortholog2",
                                       "hits": "go_hits2"})
    df = df.merge(h1, on=["nodename", "ortholog1"], how="left")
    df = df.merge(h2, on=["nodename", "ortholog2"], how="left")
    df[["go_hits1", "go_hits2"]] = df[["go_hits1", "go_hits2"]].fillna(0)
    df["go_hits_pair"] = df["go_hits1"] + df["go_hits2"]
    # Co-membership: count terms whose fam_set contains BOTH ortholog1 and ortholog2.
    df["pair_key"] = [pair_key(a, b) for a, b in zip(df["ortholog1"], df["ortholog2"])]
    co_counts = {}
    for clade, fset in zip(go["clade"], go["fam_set"]):
        if len(fset) < 2:
            continue
        fams = sorted(fset)
        for i in range(len(fams)):
            for j in range(i + 1, len(fams)):
                co_counts[(clade, pair_key(fams[i], fams[j]))] = \
                    co_counts.get((clade, pair_key(fams[i], fams[j])), 0) + 1
    df["go_terms_shared"] = [co_counts.get((c, k), 0)
                              for c, k in zip(df["nodename"], df["pair_key"])]
    return df


def add_entanglement(df: pd.DataFrame, ent_path: Path) -> pd.DataFrame:
    if not ent_path.exists():
        df["entangle_fold"] = np.nan
        return df
    ent = pd.read_csv(ent_path, sep="\t")
    ent = ent[ent["clade"] != "clade"].copy()  # drop accidental header rows
    ent["alg_key"] = [pair_key(a, b) for a, b in zip(ent["alg_a"], ent["alg_b"])]
    ent_slim = (ent[["clade", "alg_key", "fold_enrichment"]]
                .drop_duplicates(["clade", "alg_key"])
                .rename(columns={"clade": "nodename",
                                 "fold_enrichment": "entangle_fold"}))
    df["alg_key"] = [pair_key(a, b) for a, b in zip(df["alg_code1"], df["alg_code2"])]
    df = df.merge(ent_slim, on=["nodename", "alg_key"], how="left")
    return df


def add_fusion(df: pd.DataFrame, fusion_path: Path) -> pd.DataFrame:
    if not fusion_path.exists():
        df["fusion_category"] = ""
        return df
    fus = pd.read_csv(fusion_path, sep="\t")
    fus["alg_key"] = [pair_key(a, b) for a, b in zip(fus["alg1"], fus["alg2"])]
    slim = fus[["alg_key", "category"]].drop_duplicates("alg_key").rename(
        columns={"category": "fusion_category"})
    return df.merge(slim, on="alg_key", how="left")


def zscore(s: pd.Series) -> pd.Series:
    s = s.astype(float)
    sd = s.std(ddof=0)
    if sd == 0 or np.isnan(sd):
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.mean()) / sd


def composite_score(df: pd.DataFrame) -> pd.Series:
    # Group within (clade, distance_bin); within each group, z-score components.
    out = pd.Series(np.zeros(len(df)), index=df.index, dtype=float)
    weights = {
        "occupancy_in":      0.40,
        "abs_sigma":         0.30,
        "go_terms_shared":   0.20,
        "go_hits_pair":      0.05,
        "entangle_fold":     0.05,
    }
    df = df.copy()
    df["abs_sigma"] = df["sd_in_out_ratio_log_sigma"].abs()
    for (_, _), idx in df.groupby(["nodename", "distance_bin"]).groups.items():
        sub = df.loc[idx]
        s = pd.Series(0.0, index=idx)
        for col, w in weights.items():
            v = sub[col].fillna(0)
            s = s + w * zscore(v).fillna(0)
        out.loc[idx] = s
    return out


def percentile_pick(df: pd.DataFrame, clade: str,
                     percentiles: list[float],
                     log_window: float = 0.10,
                     min_anchors_resolved: int = 0,
                     min_subclades_resolved: int = 0) -> pd.DataFrame:
    """For one focal clade, return one row per target percentile of mean_in.

    Within each target's log-space ±log_window window pick the highest-scoring
    pair. Falls back to the nearest qualifying pair to the target if the
    window is empty. log_window=0.10 ≈ ±26% in linear space (log10±0.1).
    """
    sub = df[(df["nodename"] == clade) & (df["mean_in"] > 0)].copy()
    cov_col = f"n_anchors_resolved_{clade}"
    sub_col = f"n_subclades_resolved_{clade}"
    if cov_col in sub.columns and min_anchors_resolved > 0:
        sub = sub[sub[cov_col] >= min_anchors_resolved]
    if sub_col in sub.columns and min_subclades_resolved > 0:
        sub = sub[sub[sub_col] >= min_subclades_resolved]
    if sub.empty:
        return sub
    log_d = np.log10(sub["mean_in"].astype(float).values)
    out_rows = []
    for p in percentiles:
        target_log = float(np.percentile(log_d, p))
        in_window = sub[(log_d >= target_log - log_window) &
                         (log_d <= target_log + log_window)]
        if in_window.empty:
            # fallback: nearest pair by absolute log distance
            nearest_idx = np.argmin(np.abs(log_d - target_log))
            row = sub.iloc[[nearest_idx]].copy()
        else:
            row = in_window.sort_values("score", ascending=False).head(1).copy()
        row["selection_basis"] = f"percentile_p{int(p)}"
        row["target_mean_in_bp"] = int(10 ** target_log)
        out_rows.append(row)
    if not out_rows:
        return sub.head(0)
    return pd.concat(out_rows, ignore_index=True)


def force_keep(df: pd.DataFrame, force_pairs: list[str]) -> pd.DataFrame:
    """force_pairs: list of 'SYMBOL1/SYMBOL2@CLADE' strings."""
    if not force_pairs:
        return df
    df = df.copy()
    df["forced"] = False
    for spec in force_pairs:
        try:
            sym, clade = spec.split("@")
            a, b = sym.split("/")
        except ValueError:
            print(f"WARN: bad --force-pair spec {spec!r}, expected 'A/B@Clade'")
            continue
        mask = ((df["nodename"] == clade) &
                (((df["symbol1"] == a) & (df["symbol2"] == b)) |
                 ((df["symbol1"] == b) & (df["symbol2"] == a))))
        df.loc[mask, "forced"] = True
        df.loc[mask, "score"] = df["score"].max() + 1.0
    return df


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DEF_INPUT))
    ap.add_argument("--fammap", default=str(DEF_FAMMAP))
    ap.add_argument("--gosig", default=str(DEF_GOSIG))
    ap.add_argument("--entang", default=str(DEF_ENTANG))
    ap.add_argument("--fusion", default=str(DEF_FUSION))
    ap.add_argument("--np-acc", default=str(DEF_NCBI_ACC),
                    help="NCBI human_gene2accession.tsv.gz for NP→symbol")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--clade", default=None,
                    help="restrict to a single clade (nodename)")
    ap.add_argument("--occupancy-min", type=float, default=0.5)
    ap.add_argument("--top-n", type=int, default=5,
                    help="rows kept per (clade, distance_bin) in selected_examples.tsv")
    ap.add_argument("--dry-run", action="store_true",
                    help="sample the first 100k rows of the input table")
    ap.add_argument("--force-pair", action="append", default=[
        "EIF2S1/ATP6V1D@Vertebrata",
        "SLC25A26/TMF1@Vertebrata",
    ])
    ap.add_argument("--strategy", choices=["bin", "percentile"],
                    default="percentile",
                    help="how to pick selected_examples.tsv rows for focal clades")
    ap.add_argument("--percentile-clade", action="append",
                    default=["Vertebrata", "Arthropoda"],
                    help="focal clade(s) for percentile-spread strategy")
    ap.add_argument("--percentiles", default="5,15,25,35,45,55,65,75,85,95",
                    help="comma-separated target percentiles for percentile mode")
    ap.add_argument("--require-anchor", action="append", default=[],
                    help="<clade>:<species_key> — require both anchor genes "
                         "to resolve in this species before the pair is "
                         "eligible. Repeatable.")
    ap.add_argument("--min-anchors-resolved", type=int, default=7,
                    help="require ≥N of the focal-clade anchor species to "
                         "have both anchor genes on the same scaffold "
                         "(needs candidate_pairs_anchor_coverage.tsv)")
    ap.add_argument("--min-subclades-resolved", type=int, default=3,
                    help="require ≥N distinct sublineages among resolved anchors")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"reading {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    if args.dry_run:
        df = df.head(100_000)
    if args.clade:
        df = df[df["nodename"] == args.clade].copy()

    df = df[(df["stable_in_clade"] == 1) &
            (df["occupancy_in"] >= args.occupancy_min)].copy()
    df["distance_bin"] = df["mean_in"].apply(assign_bin)
    print(f"  qualifying stable pairs: {len(df)}  clades: {df['nodename'].nunique()}")

    np_to_sym = load_np_to_symbol(Path(args.np_acc))
    df = add_symbols(df, Path(args.fammap), np_to_sym)
    df = add_go_support(df, Path(args.gosig))
    df = add_entanglement(df, Path(args.entang))
    df = add_fusion(df, Path(args.fusion))

    # Anchor-coverage scoring (only for focal clades — heavy disk read).
    if args.strategy == "percentile" and args.percentile_clade:
        try:
            import score_anchor_coverage
            import recompute_distances
            cfg_yaml = score_anchor_coverage.REPO_ROOT / "genome_database/genome_list.yaml"
            anchor_set = sorted({sk for c in args.percentile_clade
                                 for sk in recompute_distances.ANCHORS.get(c, [])})
            print(f"computing anchor coverage for {len(anchor_set)} species ...")
            cfg = __import__('loci_io').load_genome_list(str(cfg_yaml))
            fam2scaf = score_anchor_coverage.build_lookup(
                anchor_set, recompute_distances.DEFAULT_RBH_DIR, cfg)
            # Also emit per-anchor-species presence boolean columns so
            # selection can require specific species (e.g. Drosophila).
            present_by_species: dict[str, list[bool]] = {}
            for clade in args.percentile_clade:
                anchors = recompute_distances.ANCHORS.get(clade, [])
                n_resolved, n_subclades = [], []
                per_sp_lists = {sk: [] for sk in anchors}
                for fam1, fam2 in zip(df["ortholog1"], df["ortholog2"]):
                    cnt = 0; sub = set()
                    for sk in anchors:
                        s = fam2scaf.get(sk, {})
                        a, b = s.get(fam1), s.get(fam2)
                        present = (a is not None and a == b)
                        per_sp_lists[sk].append(present)
                        if present:
                            cnt += 1
                            sc = score_anchor_coverage.SUBCLADE.get(sk)
                            if sc:
                                sub.add(sc)
                    n_resolved.append(cnt); n_subclades.append(len(sub))
                df[f"n_anchors_resolved_{clade}"] = n_resolved
                df[f"n_subclades_resolved_{clade}"] = n_subclades
                for sk, lst in per_sp_lists.items():
                    present_by_species[f"present_in_{sk}"] = lst
            for col, lst in present_by_species.items():
                df[col] = lst
        except Exception as exc:  # noqa: BLE001
            print(f"WARN: anchor-coverage scoring failed: {exc}; "
                  "proceeding without coverage filter")

    df["score"] = composite_score(df)
    df = force_keep(df, args.force_pair)

    keep_cols = [
        "nodename", "taxid", "distance_bin", "pair_label",
        "symbol1", "symbol2", "alg_code1", "alg_code2",
        "ortholog1", "ortholog2", "pair_key",
        "mean_in", "sd_in", "mean_out", "sd_out",
        "occupancy_in", "occupancy_out",
        "sd_in_out_ratio_log_sigma", "stable_in_clade",
        "go_hits_pair", "go_terms_shared",
        "entangle_fold", "fusion_category",
        "score",
    ]
    # Add per-clade anchor coverage columns if computed.
    for c in df.columns:
        if c.startswith("n_anchors_resolved_") or c.startswith("n_subclades_resolved_"):
            keep_cols.append(c)
    keep_cols = [c for c in keep_cols if c in df.columns]
    if "forced" in df.columns:
        keep_cols.append("forced")
    df = df[keep_cols].sort_values(["nodename", "distance_bin", "score"],
                                    ascending=[True, True, False])

    cand_path = out_dir / "candidate_pairs.tsv"
    df.to_csv(cand_path, sep="\t", index=False)
    print(f"wrote {cand_path}  rows={len(df)}")

    # Selected examples — strategy depends on CLI flag.
    headline = {"close", "intermediate", "long_range", "ultra_long"}
    pieces = []
    pcts = [float(x) for x in args.percentiles.split(",") if x.strip()]
    focal = set(args.percentile_clade) if args.strategy == "percentile" else set()
    # Parse --require-anchor specs.
    required_per_clade: dict[str, list[str]] = {}
    for spec in args.require_anchor:
        try:
            cl, sk = spec.split(":", 1)
        except ValueError:
            print(f"WARN bad --require-anchor {spec!r} (expect <clade>:<species>)")
            continue
        required_per_clade.setdefault(cl, []).append(sk)

    for clade in sorted(df["nodename"].unique()):
        if clade in focal:
            df_for_clade = df
            req = required_per_clade.get(clade, [])
            for sk in req:
                col = f"present_in_{sk}"
                if col in df_for_clade.columns:
                    df_for_clade = df_for_clade[df_for_clade[col] == True]
                    print(f"  filter [{clade}] requires {sk}: "
                          f"{len(df_for_clade)} pairs remaining")
                else:
                    print(f"  WARN no {col} column; "
                          "anchor-coverage step did not run for this species")
            picks = percentile_pick(df_for_clade, clade, pcts,
                                     min_anchors_resolved=args.min_anchors_resolved,
                                     min_subclades_resolved=args.min_subclades_resolved)
            pieces.append(picks)
        else:
            sub = (df[(df["nodename"] == clade) &
                       df["distance_bin"].isin(headline)]
                   .groupby("distance_bin", group_keys=False)
                   .head(args.top_n)
                   .copy())
            sub["selection_basis"] = "bin_" + sub["distance_bin"].astype(str)
            sub["target_mean_in_bp"] = pd.NA
            pieces.append(sub)
    sel = pd.concat(pieces, ignore_index=True)
    # Always include forced pairs (e.g. EIF2S1/ATP6V1D@Vertebrata).
    if "forced" in df.columns and df["forced"].any():
        forced = df[df["forced"]].copy()
        forced["selection_basis"] = "forced"
        forced["target_mean_in_bp"] = pd.NA
        # de-dupe vs sel by (nodename, pair_key)
        sel_keys = set(zip(sel["nodename"], sel["pair_key"]))
        new = forced[~forced.apply(
            lambda r: (r["nodename"], r["pair_key"]) in sel_keys, axis=1)]
        if len(new):
            sel = pd.concat([sel, new], ignore_index=True)
    # De-dupe (clade, pair_key); keep the first selection_basis seen.
    sel = (sel.drop_duplicates(["nodename", "pair_key"], keep="first")
              .sort_values(["nodename", "mean_in"]))
    sel_path = out_dir / "selected_examples.tsv"
    sel.to_csv(sel_path, sep="\t", index=False)
    print(f"wrote {sel_path}  rows={len(sel)}")

    # Smoke check: EIF2S1/ATP6V1D in Vertebrata close bin.
    if not args.clade or args.clade == "Vertebrata":
        hit = df[(df["nodename"] == "Vertebrata") &
                 (df["pair_label"].isin(["EIF2S1/ATP6V1D", "ATP6V1D/EIF2S1"]))]
        if len(hit):
            row = hit.iloc[0]
            print(f"  CHECK EIF2S1/ATP6V1D@Vertebrata: bin={row['distance_bin']} "
                  f"mean_in={row['mean_in']:.0f} occupancy_in={row['occupancy_in']:.3f}")
        else:
            print("  CHECK EIF2S1/ATP6V1D@Vertebrata: NOT FOUND in qualifying set")


if __name__ == "__main__":
    main()
